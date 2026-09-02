using DataFrames

function compute_los_spacing(box_length_pc, depth)
    depth > 0 || error("Cannot compute line-of-sight spacing for an empty cube.")
    los_pixel_length_pc = Float64(box_length_pc) / depth
    los_pixel_length_cm = los_pixel_length_pc * PARSEC_TO_CM
    los_distance_array = range(start=0.0, step=los_pixel_length_pc, length=depth)
    return los_pixel_length_pc, los_pixel_length_cm, los_distance_array
end

_stage(message) = @info message

function _validate_processing_cube_shapes(simu, LOS, expected_shape, cubes::Pair...)
    expected_shape === nothing && return nothing
    expected = Tuple(Int.(expected_shape))

    for (name, cube) in cubes
        cube === nothing && continue
        actual = size(cube)
        actual == expected || throw_config_error(
            "Cube shape mismatch for $(name) in $(simu) after LOS=$(LOS): expected $(expected) from `BoxLength_pix`, got $(actual). " *
            "Update `BoxLength_pix` or check the FITS cube dimensions.";
            code=:cube_shape_mismatch,
        )
    end

    return nothing
end

function _write_integrated_quantities(resultspath, B1, B2, BLOS, T, ne, Bperpcube, los_pixel_length_cm; metadata=nothing)
    Btotal = Btot(B1, B2, BLOS)
    intBtotal = intLOS(Btotal, los_pixel_length_cm)
    sigmaBtotal = sigmaLOS(Btotal)
    Ne = intLOS(ne, los_pixel_length_cm)
    sigmane = sigmaLOS(ne)
    sigmaT = sigmaLOS(T)
    intBLOS = intLOS(BLOS, los_pixel_length_cm)
    sigmaBLOS = sigmaLOS(BLOS)
    intBperp = intLOS(Bperpcube, los_pixel_length_cm)

    WriteData2D(resultspath, intBtotal, "intBtotal"; ensure_path=false, metadata=metadata)
    WriteData2D(resultspath, sigmaBtotal, "sigmaBtotal"; ensure_path=false, metadata=metadata)
    WriteData2D(resultspath, Ne, "intne"; ensure_path=false, metadata=metadata)
    WriteData2D(resultspath, sigmane, "sigmane"; ensure_path=false, metadata=metadata)
    WriteData2D(resultspath, sigmaT, "sigmaT"; ensure_path=false, metadata=metadata)
    WriteData2D(resultspath, intBLOS, "intBLOS"; ensure_path=false, metadata=metadata)
    WriteData2D(resultspath, sigmaBLOS, "sigmaBLOS"; ensure_path=false, metadata=metadata)
    WriteData2D(resultspath, intBperp, "intBperp"; ensure_path=false, metadata=metadata)

    return nothing
end

function _max_finite_cube(cube::AbstractArray{<:Real, 3})
    nx, ny = size(cube, 1), size(cube, 2)
    best = fill(-Inf, nx, ny)

    # Slice-major traversal keeps the innermost index on the contiguous
    # dimension instead of striding by nx*ny along each sightline.
    @inbounds for k in axes(cube, 3)
        for j in axes(cube, 2), i in axes(cube, 1)
            value = cube[i, j, k]
            isfinite(value) && (best[i, j] = max(best[i, j], Float64(value)))
        end
    end

    out = Matrix{eltype(cube)}(undef, nx, ny)
    nan = convert(eltype(out), NaN)
    # `best` only ever receives finite values, so a non-finite entry means the
    # sightline had no valid cell at all.
    @inbounds for idx in eachindex(out, best)
        out[idx] = isfinite(best[idx]) ? convert(eltype(out), best[idx]) : nan
    end

    return out
end

function _healpix_unit(DataName::String)
    return haskey(DictHeader, DataName) ? String(DictHeader[DataName]["bunit"]) : ""
end

function _healpix_extname(DataName::String)
    return uppercase(replace(DataName, "_" => ""))
end

function _healpix_map_vector(data)
    if data isa AbstractVector
        return collect(data)
    elseif ndims(data) == 2 && size(data, 2) == 1
        return collect(view(data, :, 1))
    end

    error("Expected a HEALPix map vector or Npix x 1 array, got size $(size(data)).")
end

function _healpix_stack_matrix(data)
    ndims(data) == 3 && size(data, 2) == 1 ||
        error("Expected a HEALPix stack encoded as Npix x 1 x Nslice, got size $(size(data)).")
    return reshape(data, size(data, 1), size(data, 3))
end

_healpix_coordsys(hp_meta) = get(hp_meta, :coordsys, nothing)

function _write_healpix_map_quantity(resultspath, data, DataName::String, hp_meta)
    path = joinpath(resultspath, "$(DataName).fits")
    write_healpix_map(path, _healpix_map_vector(data);
        nside = hp_meta.nside,
        order = hp_meta.order,
        unit = _healpix_unit(DataName),
        extname = _healpix_extname(DataName),
        coordsys = _healpix_coordsys(hp_meta),
    )
    return path
end

function _write_healpix_stack_quantity(resultspath, data, DataName::String, coordinates, hp_meta; basename::Union{Nothing,String}=nothing)
    matrix = _healpix_stack_matrix(data)
    if size(matrix, 2) == 1
        return [_write_healpix_map_quantity(resultspath, view(matrix, :, 1), DataName, hp_meta)]
    end

    return write_healpix_stack(resultspath, matrix, basename === nothing ? DataName : basename, collect(coordinates);
        nside = hp_meta.nside,
        order = hp_meta.order,
        unit = _healpix_unit(DataName),
        extname = _healpix_extname(DataName),
        coordsys = _healpix_coordsys(hp_meta),
    )
end

function _write_healpix_rmclean_result(resultspath, result, hp_meta)
    coordsys = _healpix_coordsys(hp_meta)
    write_healpix_stack(resultspath, result.cleanFDF, "cleanFDF", result.phi;
        nside = hp_meta.nside,
        order = hp_meta.order,
        unit = _healpix_unit("cleanFDF"),
        extname = _healpix_extname("cleanFDF"),
        coordsys = coordsys,
    )
    write_healpix_stack(resultspath, result.realCleanFDF, "realCleanFDF", result.phi;
        nside = hp_meta.nside,
        order = hp_meta.order,
        unit = _healpix_unit("realCleanFDF"),
        extname = _healpix_extname("realCleanFDF"),
        coordsys = coordsys,
    )
    write_healpix_stack(resultspath, result.imagCleanFDF, "imagCleanFDF", result.phi;
        nside = hp_meta.nside,
        order = hp_meta.order,
        unit = _healpix_unit("imagCleanFDF"),
        extname = _healpix_extname("imagCleanFDF"),
        coordsys = coordsys,
    )
    write_healpix_stack(resultspath, abs.(result.residual), "residualFDF", result.phi;
        nside = hp_meta.nside,
        order = hp_meta.order,
        unit = _healpix_unit("residualFDF"),
        extname = _healpix_extname("residualFDF"),
        coordsys = coordsys,
    )

    return nothing
end

function _write_rmclean_cartesian(resultspath, result; metadata=nothing)
    WriteData3D(resultspath, result.cleanFDF, "cleanFDF", result.phi; ensure_path=false, metadata=metadata)
    WriteData3D(resultspath, result.realCleanFDF, "realCleanFDF", result.phi; ensure_path=false, metadata=metadata)
    WriteData3D(resultspath, result.imagCleanFDF, "imagCleanFDF", result.phi; ensure_path=false, metadata=metadata)
    WriteData3D(resultspath, abs.(result.residual), "residualFDF", result.phi; ensure_path=false, metadata=metadata)
    return nothing
end

function _write_integrated_quantities_healpix(resultspath, B1, B2, BLOS, T, ne, Bperpcube, los_pixel_length_cm, hp_meta)
    Btotal = Btot(B1, B2, BLOS)
    quantities = (
        "intBtotal" => intLOS(Btotal, los_pixel_length_cm),
        "sigmaBtotal" => sigmaLOS(Btotal),
        "intne" => intLOS(ne, los_pixel_length_cm),
        "sigmane" => sigmaLOS(ne),
        "sigmaT" => sigmaLOS(T),
        "intBLOS" => intLOS(BLOS, los_pixel_length_cm),
        "sigmaBLOS" => sigmaLOS(BLOS),
        "intBperp" => intLOS(Bperpcube, los_pixel_length_cm),
    )

    for (name, data) in quantities
        _write_healpix_map_quantity(resultspath, data, name, hp_meta)
    end

    return nothing
end

function _apply_synchrotron_filter!(Qnu, Unu, T_nu, Llarge_filter_pix)
    # The filter works entirely in PIXEL units (Δx = Δy = 1 pixel), matching
    # the documented convention: `kernel_size_synchrotron` is the largest
    # retained spatial scale in pixels. No small-scale cut is applied beyond
    # the Nyquist limit (Lcut_small = 2 pixels ⇔ f = fNy = 0.5 cycle/pixel).
    n, m = size(Qnu, 1), size(Qnu, 2)
    H, _ = instrument_bandpass_L(
        n,
        m;
        Δx = 1.0,
        Δy = 1.0,
        Lcut_small = 2.0,
        Llarge = Float64(Llarge_filter_pix),
        fNy = 0.5,
    )

    Qnu .= apply_to_array_xy(Qnu, H; n = n, m = m)
    Unu .= apply_to_array_xy(Unu, H; n = n, m = m)
    T_nu .= apply_to_array_xy(T_nu, H; n = n, m = m)

    return nothing
end

"""
    _add_noise!(Qnu, Unu, SNR_nu, rng)

Add gaussian noise to the Q and U cubes so that the per-channel polarized
signal-to-noise ratio equals `SNR_nu`. The reference signal in each frequency
channel is the rms of the polarized intensity, P_rms = sqrt(<Q²> + <U²>), and
the same standard deviation σ = P_rms / SNR_nu is applied to Q and U
(σ in the same unit as the cubes, i.e. Kelvin).
"""
function _add_noise!(Qnu, Unu, SNR_nu, rng)
    isfinite(SNR_nu) && SNR_nu > 0 || error("SNR_nu must be finite and > 0, got $SNR_nu")
    noiseQ = similar(Qnu[:, :, 1])
    noiseU = similar(Unu[:, :, 1])

    @views for i in axes(Qnu, 3)
        Qch = Qnu[:, :, i]
        Uch = Unu[:, :, i]
        sum_p2 = 0.0
        nfinite = 0
        @inbounds for idx in eachindex(Qch, Uch)
            q = Qch[idx]
            u = Uch[idx]
            if isfinite(q) && isfinite(u)
                sum_p2 += abs2(q) + abs2(u)
                nfinite += 1
            end
        end
        # An entirely masked channel contains no signal from which a noise
        # level can be estimated. Leave it masked instead of contaminating
        # the other channels (or attempting randn! with a NaN scale).
        nfinite == 0 && continue
        P_rms = sqrt(sum_p2 / nfinite)
        sigma = P_rms / SNR_nu
        isfinite(sigma) && sigma > 0 || continue
        randn!(rng, noiseQ)
        randn!(rng, noiseU)
        @inbounds for idx in eachindex(Qch, Uch, noiseQ, noiseU)
            if isfinite(Qch[idx]) && isfinite(Uch[idx])
                Qch[idx] += sigma * noiseQ[idx]
                Uch[idx] += sigma * noiseU[idx]
            end
        end
    end

    return nothing
end

"""Return the channels whose centre frequency falls inside an RFI range (MHz)."""
function _rfi_channel_mask(nuArray, rfi_ranges)
    isempty(rfi_ranges) && return falses(length(nuArray))
    return BitVector(any(lo <= nu <= hi for (lo, hi) in rfi_ranges) for nu in nuArray)
end

function _apply_rfi_mask!(mask, arrays...)
    any(mask) || return arrays
    for array in arrays
        array === nothing && continue
        @views array[:, :, mask] .= eltype(array)(NaN)
    end
    return arrays
end

# Convert an array to the requested working precision. `Float64` (the default)
# returns the input unchanged, preserving the historical behaviour exactly.
_to_precision(::Type{T}, x::Nothing) where {T} = nothing
_to_precision(::Type{T}, x::AbstractArray) where {T} = eltype(x) === T ? x : T.(x)

function _physical_mask_threshold(mask_config, key::AbstractString)
    mask_config === nothing && return nothing
    mask_config isa AbstractDict || return nothing
    return get(mask_config, key, nothing)
end

function _has_physical_mask(mask_config)
    mask_config isa AbstractDict || return false
    return any(_physical_mask_threshold(mask_config, key) !== nothing for key in ("T_min", "T_max", "n_min", "n_max"))
end

function _physical_cell_mask(T, n, mask_config)
    mask = trues(size(T))
    Tmin = _physical_mask_threshold(mask_config, "T_min")
    Tmax = _physical_mask_threshold(mask_config, "T_max")
    nmin = _physical_mask_threshold(mask_config, "n_min")
    nmax = _physical_mask_threshold(mask_config, "n_max")

    Tmin === nothing || (mask .&= T .>= Tmin)
    Tmax === nothing || (mask .&= T .<= Tmax)
    nmin === nothing || (mask .&= n .>= nmin)
    nmax === nothing || (mask .&= n .<= nmax)

    return mask
end

function _apply_physical_mask!(B1, B2, BLOS, T, n, nHp, mask_config)
    _has_physical_mask(mask_config) || return nothing

    mask = _physical_cell_mask(T, n, mask_config)
    rejected = .!mask
    B1[rejected] .= zero(eltype(B1))
    B2[rejected] .= zero(eltype(B2))
    BLOS[rejected] .= zero(eltype(BLOS))
    T[rejected] .= zero(eltype(T))
    n[rejected] .= zero(eltype(n))
    nHp === nothing || (nHp[rejected] .= zero(eltype(nHp)))
    return mask
end

function _density_to_number_density!(n, density_kind::AbstractString, mean_molecular_weight::Real, hydrogen_mass_g::Real)
    kind = lowercase(String(density_kind))
    kind == "number_density" && return n
    kind == "mass_density" || error("Unknown density_kind $(density_kind). Use \"number_density\" or \"mass_density\".")
    n ./= (Float64(mean_molecular_weight) * Float64(hydrogen_mass_g))
    return n
end

function _process_synchrotron_common(
    simu::AbstractString,
    LOS,
    FaradayRotation::AbstractString,
    responseSynchrotron::AbstractString,
    df::DataFrame,
    add_noise,
    SNR_nu,
    kernel_size_synchrotron,
    nuArray::AbstractArray,
    PhiArray,
    BoxLength_pc,
    conversionn,
    conversionT,
    conversionB,
    electron_density_builder;
    write_ne::Bool = true,
    log_progress::Bool = false,
    rng = Random.default_rng(),
    expected_shape = nothing,
    metadata = nothing,
    rm_clean_enabled::Bool = false,
    rm_clean_gain::Real = 0.1,
    rm_clean_niter::Integer = 1000,
    rm_clean_threshold::Real = 0.0,
    float_type::Type{<:AbstractFloat} = Float64,
    tile_rows::Union{Nothing, Integer} = nothing,
    field_sources=nothing,
    physical_mask=nothing,
    density_kind::AbstractString="number_density",
    mean_molecular_weight::Real=1.0,
    hydrogen_mass_g::Real=M_p,
    outputs=Set(["integrated", "stokes", "rm", "fdf", "spectral_index", "diagnostics"]),
    checkpoint_signature=nothing,
    rfi_ranges=Tuple{Float64, Float64}[],
)
    selected_outputs = Set(String.(outputs))
    want(name) = name in selected_outputs
    need_qu = want("stokes") || want("fdf") || want("diagnostics")
    need_t = want("stokes") || want("spectral_index") || want("diagnostics")
    resultspath = joinpath(simu, LOS, "Synchrotron")
    mkpath(resultspath)
    fits_metadata = metadata === nothing ? Dict{String, Any}() : copy(metadata)
    fits_metadata["LOS"] = String(LOS)
    fits_metadata["BLENPC"] = Float64(BoxLength_pc)
    grid_kind = simulation_grid_kind(simu, field_sources)
    grid_kind == :healpix && _validate_healpix_los(LOS)
    hp_meta = grid_kind == :healpix ? healpix_simulation_metadata(simu, field_sources) : nothing
    fits_metadata["GRID"] = grid_kind == :healpix ? "HEALPIX" : grid_kind == :amr ? "AMR" : "IMAGE"

    if float_type !== Float64
        hp_meta === nothing || throw_config_error(
            "`precision = \"float32\"` is not supported with HEALPix inputs yet. Use the default precision.";
            code=:unsupported_grid_operation)
        fits_metadata["PRECIS"] = lowercase(string(float_type))
    end

    if tile_rows !== nothing
        grid_kind == :amr && throw_config_error(
            "`tile_size` is not supported with AMR inputs because rasterization must validate the complete leaf-cell coverage.";
            code=:unsupported_grid_operation)
        if hp_meta !== nothing
            return _process_synchrotron_tiled_healpix(
                simu, LOS, FaradayRotation, df, nuArray, PhiArray, BoxLength_pc,
                conversionn, conversionT, conversionB, electron_density_builder;
                write_ne = write_ne,
                log_progress = log_progress,
                tile_rows = Int(tile_rows),
                resultspath = resultspath,
                hp_meta = hp_meta,
                field_sources = field_sources,
                physical_mask = physical_mask,
                density_kind = density_kind,
                mean_molecular_weight = mean_molecular_weight,
                hydrogen_mass_g = hydrogen_mass_g,
                outputs = selected_outputs,
                checkpoint_signature = checkpoint_signature,
                rfi_ranges = rfi_ranges,
            )
        end
        return _process_synchrotron_tiled(
            simu, LOS, FaradayRotation, df, nuArray, PhiArray, BoxLength_pc,
            conversionn, conversionT, conversionB, electron_density_builder;
            write_ne = write_ne,
            log_progress = log_progress,
            expected_shape = expected_shape,
            fits_metadata = fits_metadata,
            float_type = float_type,
            tile_rows = Int(tile_rows),
            resultspath = resultspath,
            field_sources = field_sources,
            physical_mask = physical_mask,
            density_kind = density_kind,
            mean_molecular_weight = mean_molecular_weight,
            hydrogen_mass_g = hydrogen_mass_g,
            outputs = selected_outputs,
            checkpoint_signature = checkpoint_signature,
            rfi_ranges = rfi_ranges,
        )
    end

    SimuParameters = ReadSimulation(simu, LOS, conversionn, conversionT, conversionB; field_sources = field_sources)
    B1, B2, BLOS = SimuParameters[1], SimuParameters[2], SimuParameters[3]
    T, n = SimuParameters[4], SimuParameters[5]
    nHp = SimuParameters[7]

    B1 = _to_precision(float_type, B1)
    B2 = _to_precision(float_type, B2)
    BLOS = _to_precision(float_type, BLOS)
    T = _to_precision(float_type, T)
    n = _to_precision(float_type, n)
    nHp = _to_precision(float_type, nHp)
    _density_to_number_density!(n, density_kind, mean_molecular_weight, hydrogen_mass_g)
    _validate_processing_cube_shapes(simu, LOS, expected_shape,
        "B1" => B1,
        "B2" => B2,
        "BLOS" => BLOS,
        "temperature" => T,
        "density" => n,
        "densityHp" => nHp,
    )
    _apply_physical_mask!(B1, B2, BLOS, T, n, nHp, physical_mask)

    los_pixel_length_pc, los_pixel_length_cm, los_distance_array = compute_los_spacing(BoxLength_pc, size(B1, 3))

    Bperpcube = Bperp(B1, B2)
    psi_src = IntrinsicAngle(B1, B2)

    faraday_enabled = uppercase(FaradayRotation) == "Y"
    filtering_enabled = uppercase(responseSynchrotron) == "Y"
    noise_enabled = uppercase(add_noise) == "Y"

    _stage("Computing electron density")
    ne = _to_precision(float_type, electron_density_builder(T, n, nHp))
    if write_ne && want("integrated")
        if hp_meta === nothing
            WriteData3D(resultspath, ne, "ne", los_distance_array; ensure_path=false, metadata=fits_metadata)
        else
            _write_healpix_stack_quantity(resultspath, ne, "ne", los_distance_array, hp_meta)
        end
    end

    _stage("Computing integrated quantities")
    if want("integrated")
        if hp_meta === nothing
            _write_integrated_quantities(resultspath, B1, B2, BLOS, T, ne, Bperpcube, los_pixel_length_cm; metadata=fits_metadata)
        else
            _write_integrated_quantities_healpix(resultspath, B1, B2, BLOS, T, ne, Bperpcube, los_pixel_length_cm, hp_meta)
        end
    end
    B1 = nothing
    B2 = nothing
    T = nothing
    n = nothing

    RMcube = nothing
    if faraday_enabled && (want("rm") || need_qu)
        resultspath = joinpath(resultspath, "WithFaraday")
        mkpath(resultspath)
        _stage("Computing RM")
        dRM = deltaRM(BLOS, ne, los_pixel_length_pc)
        RMcube = RM(_to_precision(float_type, dRM))
        RMmap = RMcube[:, :, end]
        if want("rm")
            if hp_meta === nothing
                WriteData2D(resultspath, RMmap, "RMmap"; ensure_path=false, metadata=fits_metadata)
            else
                _write_healpix_map_quantity(resultspath, RMmap, "RMmap", hp_meta)
            end
        end
    else
        resultspath = joinpath(resultspath, "noFaraday")
        mkpath(resultspath)
        @info "No Faraday rotation included"
    end

    if need_qu && faraday_enabled
        _stage("Computing Qnu and Unu with Faraday rotation")
        Qnu, Unu = QUnu3D(Bperpcube, psi_src, RMcube, nuArray, df, los_pixel_length_cm; log_progress = log_progress)
    elseif need_qu
        _stage("Computing Qnu and Unu without Faraday rotation")
        Qnu, Unu = QUnuNoFaraday3D(Bperpcube, psi_src, nuArray, df, los_pixel_length_cm; log_progress = log_progress)
    end

    if !need_qu
        Qnu = nothing
        Unu = nothing
    end
    T_nu = need_t ? Tnu3D(Bperpcube, nuArray, df, los_pixel_length_cm) : nothing
    Bperpcube = nothing

    if filtering_enabled && (need_qu || need_t)
        hp_meta === nothing ||
            throw_config_error("HEALPix inputs do not support the current cartesian Fourier filtering (`responseSynchrotron=Y`). Disable filtering or add a spherical-harmonic HEALPix filter.";
                code=:unsupported_grid_operation)
        _stage("Applying interferometric Fourier mask")
        need_qu && need_t || throw_config_error("Filtering selective outputs requires both Q/U and T; include `stokes` or `diagnostics`."; code=:invalid_outputs)
        _apply_synchrotron_filter!(Qnu, Unu, T_nu, kernel_size_synchrotron)
        resultspath = joinpath(resultspath, "filtered")
        mkpath(resultspath)
    else
        @info "No filtering performed"
    end

    if noise_enabled && need_qu
        _stage("Adding gaussian noise to Q and U")
        _add_noise!(Qnu, Unu, SNR_nu, rng)
    end


    rfi_mask = _rfi_channel_mask(nuArray, rfi_ranges)
    if any(rfi_mask)
        _stage("Flagging $(count(rfi_mask)) RFI frequency channel(s)")
        _apply_rfi_mask!(rfi_mask, Qnu, Unu, T_nu)
    end
    valid_channels = .!rfi_mask
    valid_nu = nuArray[valid_channels]

    # The internal frequency axis is in MHz; FITS headers declare CUNIT3 = "Hz",
    # so the spectral axis is converted to Hz once, here, at write time.
    nuArray_Hz = collect(Float64, nuArray) .* 1e6

    if want("stokes")
        if hp_meta === nothing
            WriteQUnu3D(resultspath, Qnu, Unu, nuArray_Hz; ensure_path=false, metadata=fits_metadata)
            WriteData3D(resultspath, T_nu, "T_nu", nuArray_Hz; ensure_path=false, metadata=fits_metadata, filename="Tnu.fits")
        else
            _write_healpix_stack_quantity(resultspath, Qnu, "Qnu", nuArray_Hz, hp_meta)
            _write_healpix_stack_quantity(resultspath, Unu, "Unu", nuArray_Hz, hp_meta)
            _write_healpix_stack_quantity(resultspath, T_nu, "T_nu", nuArray_Hz, hp_meta; basename="Tnu")
        end
    end

    Pnumax = nothing
    if want("stokes") || want("diagnostics")
        _stage("Computing polarized intensity")
        Pnucube = Pnu(Qnu, Unu)
        Pnumax = _max_finite_cube(Pnucube)
        if want("stokes")
            polfrac = PolarizationFraction(Pnucube, T_nu)
            polfracmax = _max_finite_cube(polfrac)
            if hp_meta === nothing
                WriteData3D(resultspath, Pnucube, "Pnu", nuArray_Hz; ensure_path=false, metadata=fits_metadata)
                WriteData2D(resultspath, Pnumax, "Pnumax"; ensure_path=false, metadata=fits_metadata)
                WriteData3D(resultspath, polfrac, "polfrac", nuArray_Hz; ensure_path=false, metadata=fits_metadata)
                WriteData2D(resultspath, polfracmax, "polfracmax"; ensure_path=false, metadata=fits_metadata)
            else
                _write_healpix_stack_quantity(resultspath, Pnucube, "Pnu", nuArray_Hz, hp_meta)
                _write_healpix_map_quantity(resultspath, Pnumax, "Pnumax", hp_meta)
                _write_healpix_stack_quantity(resultspath, polfrac, "polfrac", nuArray_Hz, hp_meta)
                _write_healpix_map_quantity(resultspath, polfracmax, "polfracmax", hp_meta)
            end
        end
    end
    if want("spectral_index")
      _stage("Computing spectral index map")
      try
        if length(nuArray) >= 2
            # T_nu is a brightness temperature, so the fitted log-log slope is
            # the temperature index beta (T ∝ ν^β). The map is written in the
            # flux-density convention alpha = beta + 2 (S_ν ∝ ν^α); the slope
            # error is identical for both conventions. min_channels = 2 keeps
            # short frequency axes usable (the error map is then NaN).
            beta, alpha_err = spectral_index_map(T_nu, nuArray; min_channels = 2)
            alpha = beta .+ 2.0
            alpha_metadata = copy(fits_metadata)
            alpha_metadata["ALPHADEF"] = "S_nu ~ nu^alpha; alpha = beta_Tb + 2"
            if hp_meta === nothing
                WriteData2D(resultspath, alpha, "alpha"; ensure_path=false, metadata=alpha_metadata)
                WriteData2D(resultspath, alpha_err, "alpha_err"; ensure_path=false, metadata=alpha_metadata)
            else
                _write_healpix_map_quantity(resultspath, alpha, "alpha", hp_meta)
                _write_healpix_map_quantity(resultspath, alpha_err, "alpha_err", hp_meta)
            end
        else
            @info "Skipping spectral index map: at least two frequency channels are required"
        end
    catch err
        @warn "Failed to compute or write the spectral index map" path = resultspath exception = (err, catch_backtrace())
      end
    end

    if want("diagnostics")
      try
        write_polarization_diagnostic_plots(resultspath, Qnu, Unu, T_nu, nuArray_Hz; Pnumax = Pnumax)
    catch err
        @warn "Failed to write polarization diagnostic plots" path = resultspath exception = (err, catch_backtrace())
      end
    end

    if faraday_enabled && want("fdf")
        _stage("Performing RM synthesis")
        rmsf = nothing
        if hp_meta === nothing
            FDF, realFDF, imagFDF = RMSynthesis(Qnu[:, :, valid_channels], Unu[:, :, valid_channels], valid_nu * 1e6, PhiArray; log_progress = log_progress)
            Pmax = maxCube(FDF)
            WriteData3D(resultspath, FDF, "FDF", PhiArray; ensure_path=false, metadata=fits_metadata)
            WriteData3D(resultspath, realFDF, "realFDF", PhiArray; ensure_path=false, metadata=fits_metadata)
            WriteData3D(resultspath, imagFDF, "imagFDF", PhiArray; ensure_path=false, metadata=fits_metadata)
            WriteData2D(resultspath, Pmax, "Pmax"; ensure_path=false, metadata=fits_metadata)
        else
            q_stack = HealpixStack(_healpix_stack_matrix(Qnu)[:, valid_channels]; nside=hp_meta.nside, order=hp_meta.order)
            u_stack = HealpixStack(_healpix_stack_matrix(Unu)[:, valid_channels]; nside=hp_meta.nside, order=hp_meta.order)
            hp_result = RMSynthesisHealpix(q_stack, u_stack, valid_nu * 1e6, PhiArray; log_progress = log_progress)
            write_healpix_rm_result(resultspath, hp_result; overwrite=true)
            Pmax = maximum(hp_result.fdf, dims=2)
            _write_healpix_map_quantity(resultspath, Pmax, "Pmax", hp_meta)
        end

        # RMSF diagnostics: write the spread function and its resolution metrics
        # alongside the FDF products. If RM-CLEAN is requested, diagnostics are
        # required and failures should abort the run; otherwise they only warn.
        if rm_clean_enabled
            rmsf = rmsf_diagnostics(valid_nu * 1e6, PhiArray)
            @info "RMSF diagnostics" fwhm = rmsf.fwhm delta_phi_theory = rmsf.fwhm_theoretical phi_max = rmsf.phi_max max_scale = rmsf.max_scale
            write_rmsf(resultspath, rmsf; ensure_path=false, metadata=fits_metadata)
        else
            try
                rmsf = rmsf_diagnostics(valid_nu * 1e6, PhiArray)
                @info "RMSF diagnostics" fwhm = rmsf.fwhm delta_phi_theory = rmsf.fwhm_theoretical phi_max = rmsf.phi_max max_scale = rmsf.max_scale
                write_rmsf(resultspath, rmsf; ensure_path=false, metadata=fits_metadata)
            catch err
                @warn "Failed to compute or write RMSF diagnostics" exception = err
            end
        end

        if rm_clean_enabled
            _stage("Performing RM-CLEAN")
            if hp_meta === nothing
                clean_result = rmclean(realFDF, imagFDF, PhiArray, rmsf;
                    gain = rm_clean_gain,
                    threshold = rm_clean_threshold,
                    niter = rm_clean_niter,
                    log_progress = log_progress,
                )
                _write_rmclean_cartesian(resultspath, clean_result; metadata=fits_metadata)
            else
                clean_result = RMClean(q_stack.pixels, u_stack.pixels, valid_nu * 1e6, PhiArray;
                    gain = rm_clean_gain,
                    threshold = rm_clean_threshold,
                    niter = rm_clean_niter,
                    diagnostics = rmsf,
                    log_progress = log_progress,
                )
                _write_healpix_rmclean_result(resultspath, clean_result, hp_meta)
            end
        end
    else
        @info "No Faraday tomography performed"
    end

    return nothing
end

function ProcessSynchrotron(simu::AbstractString, LOS, FaradayRotation::AbstractString, responseSynchrotron::AbstractString,
                       df::DataFrame, add_noise, SNR_nu, kernel_size_synchrotron, zeta::Float64, Geff::Float64,
                       omegaPAH::Float64, XC::Float64, nuArray::AbstractArray, PhiArray,
                       PixelLength_pc, PixelLength_cm, BoxLength_pc,
                       DistanceArray, conversionn, conversionT, conversionB; log_progress::Bool = false, rng = Random.default_rng(), expected_shape = nothing, metadata = nothing,
                       rm_clean_enabled::Bool = false, rm_clean_gain::Real = 0.1, rm_clean_niter::Integer = 1000, rm_clean_threshold::Real = 0.0,
                       float_type::Type{<:AbstractFloat} = Float64, tile_rows::Union{Nothing, Integer} = nothing, field_sources=nothing,
                       physical_mask=nothing, density_kind::AbstractString="number_density", mean_molecular_weight::Real=1.0,
                       hydrogen_mass_g::Real=M_p, outputs=Set(["integrated", "stokes", "rm", "fdf", "spectral_index", "diagnostics"]),
                       checkpoint_signature=nothing, rfi_ranges=Tuple{Float64, Float64}[])
    return _process_synchrotron_common(
        simu,
        LOS,
        FaradayRotation,
        responseSynchrotron,
        df,
        add_noise,
        SNR_nu,
        kernel_size_synchrotron,
        nuArray,
        PhiArray,
        BoxLength_pc,
        conversionn,
        conversionT,
        conversionB,
        (T, n, nHp) -> Wolfire_ne(zeta, Geff, omegaPAH, XC, T, n);
        write_ne = true,
        log_progress = log_progress,
        rng = rng,
        expected_shape = expected_shape,
        metadata = metadata,
        rm_clean_enabled = rm_clean_enabled,
        rm_clean_gain = rm_clean_gain,
        rm_clean_niter = rm_clean_niter,
        rm_clean_threshold = rm_clean_threshold,
        float_type = float_type,
        tile_rows = tile_rows,
        field_sources = field_sources,
        physical_mask = physical_mask,
        density_kind = density_kind,
        mean_molecular_weight = mean_molecular_weight,
        hydrogen_mass_g = hydrogen_mass_g,
        outputs = outputs,
        checkpoint_signature = checkpoint_signature,
        rfi_ranges = rfi_ranges,
    )
end

function ProcessSynchrotron(simu::String, LOS, FaradayRotation::String, responseSynchrotron::String,
                       df::DataFrame, add_noise, SNR_nu, kernel_size_synchrotron, IonizationFraction::Float64,
                       nuArray::AbstractArray, PhiArray, PixelLength_pc, PixelLength_cm,
                       BoxLength_pc, DistanceArray, conversionn, conversionT, conversionB; log_progress::Bool = false, rng = Random.default_rng(), expected_shape = nothing, metadata = nothing,
                       rm_clean_enabled::Bool = false, rm_clean_gain::Real = 0.1, rm_clean_niter::Integer = 1000, rm_clean_threshold::Real = 0.0,
                       float_type::Type{<:AbstractFloat} = Float64, tile_rows::Union{Nothing, Integer} = nothing, field_sources=nothing,
                       physical_mask=nothing, density_kind::AbstractString="number_density", mean_molecular_weight::Real=1.0,
                       hydrogen_mass_g::Real=M_p, outputs=Set(["integrated", "stokes", "rm", "fdf", "spectral_index", "diagnostics"]),
                       checkpoint_signature=nothing, rfi_ranges=Tuple{Float64, Float64}[])
    return _process_synchrotron_common(
        simu,
        LOS,
        FaradayRotation,
        responseSynchrotron,
        df,
        add_noise,
        SNR_nu,
        kernel_size_synchrotron,
        nuArray,
        PhiArray,
        BoxLength_pc,
        conversionn,
        conversionT,
        conversionB,
        (T, n, nHp) -> ne_propto_nH(n, IonizationFraction);
        write_ne = true,
        log_progress = log_progress,
        rng = rng,
        expected_shape = expected_shape,
        metadata = metadata,
        rm_clean_enabled = rm_clean_enabled,
        rm_clean_gain = rm_clean_gain,
        rm_clean_niter = rm_clean_niter,
        rm_clean_threshold = rm_clean_threshold,
        float_type = float_type,
        tile_rows = tile_rows,
        field_sources = field_sources,
        physical_mask = physical_mask,
        density_kind = density_kind,
        mean_molecular_weight = mean_molecular_weight,
        hydrogen_mass_g = hydrogen_mass_g,
        outputs = outputs,
        checkpoint_signature = checkpoint_signature,
        rfi_ranges = rfi_ranges,
    )
end

function ProcessSynchrotron(simu::String, LOS, FaradayRotation::String, responseSynchrotron::String,
                       df::DataFrame,  add_noise, SNR_nu, kernel_size_synchrotron, nuArray::AbstractArray, PhiArray,
                       PixelLength_pc, PixelLength_cm, BoxLength_pc,
                       DistanceArray, conversionn, conversionT, conversionB; log_progress::Bool = false, rng = Random.default_rng(), expected_shape = nothing, metadata = nothing,
                       rm_clean_enabled::Bool = false, rm_clean_gain::Real = 0.1, rm_clean_niter::Integer = 1000, rm_clean_threshold::Real = 0.0,
                       float_type::Type{<:AbstractFloat} = Float64, tile_rows::Union{Nothing, Integer} = nothing, field_sources=nothing,
                       physical_mask=nothing, density_kind::AbstractString="number_density", mean_molecular_weight::Real=1.0,
                       hydrogen_mass_g::Real=M_p, outputs=Set(["integrated", "stokes", "rm", "fdf", "spectral_index", "diagnostics"]),
                       checkpoint_signature=nothing, rfi_ranges=Tuple{Float64, Float64}[])
    return _process_synchrotron_common(
        simu,
        LOS,
        FaradayRotation,
        responseSynchrotron,
        df,
        add_noise,
        SNR_nu,
        kernel_size_synchrotron,
        nuArray,
        PhiArray,
        BoxLength_pc,
        conversionn,
        conversionT,
        conversionB,
        (T, n, nHp) -> nHp;
        write_ne = false,
        log_progress = log_progress,
        rng = rng,
        expected_shape = expected_shape,
        metadata = metadata,
        rm_clean_enabled = rm_clean_enabled,
        rm_clean_gain = rm_clean_gain,
        rm_clean_niter = rm_clean_niter,
        rm_clean_threshold = rm_clean_threshold,
        float_type = float_type,
        tile_rows = tile_rows,
        field_sources = field_sources,
        physical_mask = physical_mask,
        density_kind = density_kind,
        mean_molecular_weight = mean_molecular_weight,
        hydrogen_mass_g = hydrogen_mass_g,
        outputs = outputs,
        checkpoint_signature = checkpoint_signature,
        rfi_ranges = rfi_ranges,
    )
end
