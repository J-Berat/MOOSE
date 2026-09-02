module MooseFromConfig

using CSV
using Crayons
using DataFrames
using Dates
using JSON

using ..Moose: PARSEC_TO_CM, RunConfig, ValidationResult, ensure_directory_access,
               build_density_parameters,
               ensure_readable_file, normalize_los_float_values, normalize_los_int_values,
               normalize_rng_seed, normalize_field_sources, normalize_physical_mask,
               run_moose_processing, throw_config_error,
               validate_los_float_values, validate_los_int_values,
               validation_failure, validation_success

function normalize_base_dir(cfg, config_path)
    raw_dir = get(cfg, "base_dir") do
        return validation_failure(String, "Config file $(config_path) must define `base_dir`.")
    end

    raw_dir isa ValidationResult && return raw_dir

    normalized_dir = abspath(expanduser(String(raw_dir)))
    validation_error = ensure_directory_access(normalized_dir)
    validation_error === nothing || return validation_failure(String, validation_error)

    resolved_dir = try
        realpath(normalized_dir)
    catch
        normalized_dir
    end

    return validation_success(resolved_dir)
end

function collect_simulations(cfg, base_dir, config_path)
    sims = haskey(cfg, "simulations") ? cfg["simulations"] : get(cfg, "chosen_simu") do
        return validation_failure(Vector{String}, "Config file $(config_path) must define either `simulations` or `chosen_simu` (array of simulations).")
    end
    sims isa ValidationResult && return sims
    sims isa AbstractVector || return validation_failure(Vector{String}, "`simulations` must be an array of paths or folder names.")

    paths = String[]
    errors = String[]

    for simu in sims
        candidate = isabspath(simu) ? String(simu) : joinpath(base_dir, String(simu))
        expanded = abspath(expanduser(candidate))

        validation_error = ensure_directory_access(expanded)
        if validation_error !== nothing
            push!(errors, validation_error)
            continue
        end

        resolved = try
            realpath(expanded)
        catch
            expanded
        end

        push!(paths, resolved)
    end

    if !isempty(errors)
        return validation_failure(Vector{String}, join(errors, "\n"))
    end
    isempty(paths) && return validation_failure(Vector{String}, "Config file $(config_path) must reference at least one simulation.")

    return validation_success(paths)
end

function collect_los(cfg)
    los_values = get(cfg, "chosen_LOS", ["x", "y", "z"])
    los_values isa AbstractVector ||
        throw_config_error("`chosen_LOS` must be an array of LOS identifiers (x, y, z)."; code=:invalid_los)

    normalized = unique(lowercase.(String.(los_values)))
    invalid = filter(los -> !(los in ("x", "y", "z")), normalized)
    !isempty(invalid) &&
        throw_config_error("Invalid line(s) of sight: $(join(invalid, ", ")). Allowed values are x, y, or z."; code=:invalid_los)

    return normalized
end

function normalize_yes_no_flag(value, field_name)
    if value isa Bool
        return value ? "Y" : "N"
    end

    normalized = uppercase(strip(String(value)))
    if normalized in ("Y", "YES", "TRUE", "1")
        return "Y"
    elseif normalized in ("N", "NO", "FALSE", "0")
        return "N"
    end

    throw_config_error("`$(field_name)` must be one of Y/N (or true/false). Got: $(value)"; code=:invalid_flag)
end

function parse_config_float(value, field_name)
    numeric_value = try
        value isa AbstractString ? tryparse(Float64, strip(value)) : Float64(value)
    catch
        nothing
    end
    numeric_value === nothing && throw_config_error("`$(field_name)` must be numeric. Got: $(value)"; code=:invalid_numeric)
    return numeric_value
end

function validate_positive_finite(value, field_name)
    numeric_value = parse_config_float(value, field_name)
    isfinite(numeric_value) || throw_config_error("`$(field_name)` must be finite. Got: $(value)"; code=:invalid_numeric)
    numeric_value > 0 || throw_config_error("`$(field_name)` must be > 0. Got: $(value)"; code=:invalid_numeric)
    return numeric_value
end

function validate_nonnegative_finite(value, field_name)
    numeric_value = parse_config_float(value, field_name)
    isfinite(numeric_value) || throw_config_error("`$(field_name)` must be finite. Got: $(value)"; code=:invalid_numeric)
    numeric_value >= 0 || throw_config_error("`$(field_name)` must be >= 0. Got: $(value)"; code=:invalid_numeric)
    return numeric_value
end

function build_frequency_array(cfg)
    freq_cfg = get(cfg, "freq", nothing)
    if freq_cfg isa AbstractDict
        start_val = get(freq_cfg, "start", 120.0)
        end_val = get(freq_cfg, "end", 167.0)
        step_val = get(freq_cfg, "step", 0.098)
    else
        start_val = get(cfg, "nustart", 120.0)
        end_val = get(cfg, "nuend", 167.0)
        step_val = get(cfg, "dnu", 0.098)
    end

    start_numeric = parse_config_float(start_val, "freq.start")
    end_numeric = parse_config_float(end_val, "freq.end")
    step_numeric = parse_config_float(step_val, "freq.step")

    step_numeric <= 0 && throw_config_error("Frequency step must be positive (received $(step_val))."; code=:invalid_frequency)
    end_numeric <= start_numeric && throw_config_error(
        "The end frequency ($(end_val)) must be greater than the start frequency ($(start_val)).";
        code=:invalid_frequency,
    )

    return start_numeric, end_numeric, step_numeric
end

function normalize_rfi_ranges(cfg)
    raw = get(cfg, "rfi", get(cfg, "rfi_ranges_mhz", nothing))
    raw === nothing && return Tuple{Float64, Float64}[]
    if raw isa AbstractDict
        enabled = normalize_yes_no_flag(get(raw, "enabled", true), "rfi.enabled")
        enabled == "N" && return Tuple{Float64, Float64}[]
        raw = get(raw, "ranges_mhz", get(raw, "ranges", nothing))
        raw === nothing && throw_config_error(
            "`rfi.ranges_mhz` is required when RFI flagging is enabled."; code=:invalid_rfi)
    end
    raw isa AbstractVector || throw_config_error(
        "`rfi.ranges_mhz` must be an array of [min, max] frequency ranges."; code=:invalid_rfi)

    ranges = Tuple{Float64, Float64}[]
    for (index, item) in enumerate(raw)
        if item isa AbstractVector && length(item) == 2
            lo = parse_config_float(item[1], "rfi.ranges_mhz[$index][1]")
            hi = parse_config_float(item[2], "rfi.ranges_mhz[$index][2]")
        elseif item isa AbstractDict
            lo = parse_config_float(get(item, "min", get(item, "start", nothing)), "rfi.ranges_mhz[$index].min")
            hi = parse_config_float(get(item, "max", get(item, "end", nothing)), "rfi.ranges_mhz[$index].max")
        else
            throw_config_error(
                "RFI range $index must be [min, max] or an object with min/max (MHz)."; code=:invalid_rfi)
        end
        isfinite(lo) && isfinite(hi) && lo > 0 && hi > 0 || throw_config_error(
            "RFI range $index must contain positive finite frequencies."; code=:invalid_rfi)
        lo <= hi || throw_config_error(
            "RFI range $index is reversed: $lo MHz > $hi MHz."; code=:invalid_rfi)
        push!(ranges, (lo, hi))
    end
    sort!(ranges; by=first)
    return ranges
end

function build_faraday(cfg)
    faraday_cfg = get(cfg, "faraday", nothing)
    if faraday_cfg isa AbstractDict
        enabled = normalize_yes_no_flag(get(faraday_cfg, "enabled", true), "faraday.enabled")
        phimin = get(faraday_cfg, "phimin", -20.0)
        phimax = get(faraday_cfg, "phimax", 20.0)
        dphi = get(faraday_cfg, "dphi", 0.1)
    else
        enabled = normalize_yes_no_flag(get(cfg, "FaradayRotation", "Y"), "FaradayRotation")
        phimin = get(cfg, "phimin", -20.0)
        phimax = get(cfg, "phimax", 20.0)
        dphi = get(cfg, "dphi", 0.1)
    end

    phimin = parse_config_float(phimin, "phimin")
    phimax = parse_config_float(phimax, "phimax")
    dphi = parse_config_float(dphi, "dphi")
    isfinite(phimin) || throw_config_error("`phimin` must be finite. Got: $(phimin)"; code=:invalid_faraday_range)
    isfinite(phimax) || throw_config_error("`phimax` must be finite. Got: $(phimax)"; code=:invalid_faraday_range)
    isfinite(dphi) || throw_config_error("`dphi` must be finite. Got: $(dphi)"; code=:invalid_faraday_range)
    dphi <= 0 && throw_config_error("The Faraday step dphi must be positive (received $(dphi))."; code=:invalid_faraday_range)
    phimax <= phimin && throw_config_error(
        "Faraday rotation range is invalid: phimax ($(phimax)) must be greater than phimin ($(phimin)).";
        code=:invalid_faraday_range,
    )
    return enabled, phimin, phimax, dphi
end

function validate_nonnegative_int(value, field_name)
    int_value = try
        value isa AbstractString ? tryparse(Int, strip(value)) : Int(value)
    catch
        nothing
    end
    int_value === nothing && throw_config_error("`$(field_name)` must be an integer. Got: $(value)"; code=:invalid_numeric)
    int_value >= 0 || throw_config_error("`$(field_name)` must be >= 0. Got: $(value)"; code=:invalid_numeric)
    return int_value
end

function build_rm_clean(cfg)
    clean_cfg = get(cfg, "rm_clean", nothing)
    if clean_cfg isa AbstractDict
        enabled = normalize_yes_no_flag(get(clean_cfg, "enabled", false), "rm_clean.enabled")
        gain = get(clean_cfg, "gain", 0.1)
        niter = get(clean_cfg, "niter", 1000)
        threshold = get(clean_cfg, "threshold", 0.0)
    else
        enabled = normalize_yes_no_flag(get(cfg, "do_rm_clean", get(cfg, "RMClean", false)), "do_rm_clean")
        gain = get(cfg, "rm_clean_gain", 0.1)
        niter = get(cfg, "rm_clean_niter", 1000)
        threshold = get(cfg, "rm_clean_threshold", 0.0)
    end

    gain = validate_positive_finite(gain, "rm_clean.gain")
    gain <= 1.0 || throw_config_error("`rm_clean.gain` must be <= 1.0. Got: $(gain)"; code=:invalid_numeric)
    niter = validate_nonnegative_int(niter, "rm_clean.niter")
    threshold = validate_nonnegative_finite(threshold, "rm_clean.threshold")

    return enabled == "Y", gain, niter, threshold
end

function normalize_box_lengths(box_length)
    if box_length isa AbstractVector
        length(box_length) == 3 || throw_config_error("Box length array must have three elements (x, y, z).")
        return (; x = Float64(box_length[1]), y = Float64(box_length[2]), z = Float64(box_length[3]))
    elseif box_length isa AbstractDict
        return (; x = Float64(get(box_length, "x", get(box_length, "X", get(box_length, "size_pc", 50.0)))),
                 y = Float64(get(box_length, "y", get(box_length, "Y", get(box_length, "size_pc", 50.0)))),
                 z = Float64(get(box_length, "z", get(box_length, "Z", get(box_length, "size_pc", 50.0)))))
    else
        return (; x = Float64(box_length), y = Float64(box_length), z = Float64(box_length))
    end
end

function normalize_box_pixels(box_length_pix)
    if box_length_pix isa AbstractVector
        length(box_length_pix) == 3 || throw_config_error("Box pixel array must have three elements (x, y, z).")
        return (; x = Int(box_length_pix[1]), y = Int(box_length_pix[2]), z = Int(box_length_pix[3]))
    elseif box_length_pix isa AbstractDict
        return (; x = Int(get(box_length_pix, "x", get(box_length_pix, "X", get(box_length_pix, "npix", 256)))),
                 y = Int(get(box_length_pix, "y", get(box_length_pix, "Y", get(box_length_pix, "npix", 256)))),
                 z = Int(get(box_length_pix, "z", get(box_length_pix, "Z", get(box_length_pix, "npix", 256)))))
    else
        return (; x = Int(box_length_pix), y = Int(box_length_pix), z = Int(box_length_pix))
    end
end

function build_distance_parameters(cfg)
    box_cfg = get(cfg, "box", nothing)
    raw_box_length_pc = box_cfg isa AbstractDict ? get(cfg, "BoxLength_pc", box_cfg) : get(cfg, "BoxLength_pc", 50.0)
    raw_box_length_pix = if haskey(cfg, "BoxLength_pix")
        cfg["BoxLength_pix"]
    elseif box_cfg isa AbstractDict
        get(box_cfg, "npix", get(box_cfg, "pixels", 256))
    else
        256
    end

    box_length_pc = normalize_los_float_values(raw_box_length_pc; default = 50.0, fallback_keys = (:size_pc,))
    box_length_pix = normalize_los_int_values(raw_box_length_pix; default = 256, fallback_keys = (:npix,))

    return box_length_pc, box_length_pix
end

"""
    MOOSE_from_config(config_path::AbstractString; quiet::Bool = false)

Run the MOOSE pipeline non-interactively using parameters defined in a JSON configuration
file. The configuration keys mirror those written by the interactive `MOOSE` function
(e.g. `base_dir`, `chosen_simu`, `chosen_LOS`, `conversionB`, `FaradayRotation`, etc.),
and also supports a streamlined schema with nested keys as produced by the Streamlit
frontend example (`freq`, `box`, `faraday`, `ne`, `emissivity`).

# Arguments
- `config_path`: path to the JSON configuration file.
- `quiet`: when `true`, skips the startup logo from `MOOSE`.

The function will raise an error if required keys are missing from the configuration.
"""
function build_config(cfg, config_path)
    base_dir_result = normalize_base_dir(cfg, config_path)
    base_dir_result.error === nothing || throw_config_error(base_dir_result.error; code=:missing_base_dir)
    base_dir = base_dir_result.value

    simu_paths_result = collect_simulations(cfg, base_dir, config_path)
    simu_paths_result.error === nothing || throw_config_error(simu_paths_result.error; code=:missing_simulation)
    simu_paths = simu_paths_result.value

    chosen_LOS = map(x -> lowercase(String(x)), collect_los(cfg))
    isempty(chosen_LOS) && throw_config_error("`chosen_LOS` must contain at least one value among x, y, z."; code=:invalid_los)
    invalid_los = [los for los in chosen_LOS if !(los in ("x", "y", "z"))]
    isempty(invalid_los) || throw_config_error("`chosen_LOS` contains invalid values: $(join(invalid_los, ", ")). Allowed values: x, y, z."; code=:invalid_los)

    conversionB = validate_positive_finite(get(cfg, "conversionB", 1.0), "conversionB")
    conversionn = validate_positive_finite(get(cfg, "conversionn", 1.0), "conversionn")
    conversionT = validate_positive_finite(get(cfg, "conversionT", 1.0), "conversionT")
    log_progress = get(cfg, "log_progress", true)
    field_sources = normalize_field_sources(get(cfg, "field_sources", get(cfg, "input_fields", nothing)))
    physical_mask = normalize_physical_mask(get(cfg, "physical_mask", get(cfg, "mask", nothing)))
    density_kind, mean_molecular_weight, hydrogen_mass_g = build_density_parameters(cfg)
    responseSynchrotron = normalize_yes_no_flag(get(cfg, "responseSynchrotron", "N"), "responseSynchrotron")
    kernel_size_synchrotron = get(cfg, "kernel_size_synchrotron", nothing)
    add_noise = normalize_yes_no_flag(get(cfg, "add_noise", "N"), "add_noise")
    SNR_nu = get(cfg, "SNR_nu", nothing)

    interpolation_file_path = get(cfg, "interpolation_file_path") do
        emiss_cfg = get(cfg, "emissivity", nothing)
        emiss_cfg isa AbstractDict ? get(emiss_cfg, "path") : nothing
    end
    interpolation_file_path === nothing && throw_config_error(
        "Config file $(config_path) must define `interpolation_file_path` or `emissivity.path`.";
        code=:missing_interpolation_path,
    )
    interpolation_file_path = isabspath(interpolation_file_path) ? interpolation_file_path : joinpath(base_dir, interpolation_file_path)

    ne_option = string(get(cfg, "ne_option", get(get(cfg, "ne", Dict()), "mode", 1)))
    IonizationFraction = get(cfg, "IonizationFraction", get(get(cfg, "ne", Dict()), "ion_fraction", 0.01))
    wolfire_constants = collect_wolfire_constants(cfg)

    BoxLength_pc, BoxLength_pix = build_distance_parameters(cfg)
    BoxLength_pc = validate_los_float_values(BoxLength_pc, "BoxLength_pc")
    BoxLength_pix = validate_los_int_values(BoxLength_pix, "BoxLength_pix")

    nustart, nuend, dnu = build_frequency_array(cfg)
    nustart = validate_positive_finite(nustart, "nustart")
    nuend = validate_positive_finite(nuend, "nuend")
    dnu = validate_positive_finite(dnu, "dnu")
    nuend > nustart || throw_config_error("`nuend` must be strictly greater than `nustart`."; code=:invalid_frequency)
    rfi_ranges = normalize_rfi_ranges(cfg)
    frequencies = collect(range(start=nustart, stop=nuend, step=dnu))
    all(nu -> any(lo <= nu <= hi for (lo, hi) in rfi_ranges), frequencies) && throw_config_error(
        "The configured RFI ranges flag every frequency channel."; code=:invalid_rfi)

    FaradayRotation, phimin, phimax, dphi = build_faraday(cfg)
    FaradayRotation = normalize_yes_no_flag(FaradayRotation, "FaradayRotation")
    dphi = validate_positive_finite(dphi, "dphi")
    phimax > phimin || throw_config_error("`phimax` must be strictly greater than `phimin`."; code=:invalid_faraday_range)
    rm_clean_enabled, rm_clean_gain, rm_clean_niter, rm_clean_threshold = build_rm_clean(cfg)
    if rm_clean_enabled && FaradayRotation != "Y"
        throw_config_error("`rm_clean.enabled` requires Faraday rotation to be enabled."; code=:invalid_rm_clean)
    end

    ne_option in ("1", "2", "3") || throw_config_error("`ne_option` must be one of \"1\", \"2\", or \"3\"."; code=:invalid_ne_option)
    if ne_option == "2"
        IonizationFraction = validate_nonnegative_finite(IonizationFraction, "IonizationFraction")
        IonizationFraction <= 1.0 || throw_config_error("`IonizationFraction` must be <= 1.0."; code=:invalid_ne_option)
    end

    if responseSynchrotron == "Y"
        kernel_size_synchrotron === nothing && throw_config_error("`kernel_size_synchrotron` is required when `responseSynchrotron` is enabled."; code=:missing_filter_kernel)
        kernel_size_synchrotron = validate_positive_finite(kernel_size_synchrotron, "kernel_size_synchrotron")
    end

    if add_noise == "Y"
        SNR_nu === nothing && throw_config_error("`SNR_nu` is required when `add_noise` is enabled."; code=:missing_noise_snr)
        SNR_nu = validate_positive_finite(SNR_nu, "SNR_nu")
    elseif SNR_nu !== nothing
        SNR_nu = validate_positive_finite(SNR_nu, "SNR_nu")
    end

    interpolation_validation = ensure_readable_file(interpolation_file_path)
    interpolation_validation === nothing || throw_config_error(interpolation_validation; code=:missing_interpolation_file)

    rng_seed = normalize_rng_seed(get(cfg, "rng_seed", nothing))

    precision = lowercase(string(get(cfg, "precision", "float64")))
    precision in ("float64", "float32") ||
        throw_config_error("`precision` must be \"float64\" or \"float32\". Got: $(precision)"; code=:invalid_precision)

    tile_size = get(cfg, "tile_size", nothing)
    if tile_size !== nothing
        tile_size = validate_nonnegative_int(tile_size, "tile_size")
        tile_size > 0 || throw_config_error("`tile_size` must be a positive integer. Got: $(tile_size)"; code=:invalid_tile_size)
        responseSynchrotron == "Y" && throw_config_error(
            "`tile_size` is incompatible with `responseSynchrotron = Y`: the interferometric Fourier mask needs the full sky plane. Disable filtering or remove `tile_size`.";
            code=:invalid_tile_size)
        add_noise == "Y" && throw_config_error(
            "`tile_size` is incompatible with `add_noise = Y`: the per-channel noise level is derived from the full-map rms. Disable noise or remove `tile_size`.";
            code=:invalid_tile_size)
        rm_clean_enabled && throw_config_error(
            "`tile_size` does not support RM-CLEAN yet. Disable `rm_clean` or remove `tile_size`.";
            code=:invalid_tile_size)
    end

    resume = lowercase(strip(String(get(cfg, "resume", "off"))))
    resume in ("off", "safe") || throw_config_error(
        "`resume` must be \"off\" or \"safe\". Got: $(resume)"; code=:invalid_resume)
    resume == "safe" && add_noise == "Y" && throw_config_error(
        "`resume = safe` is not compatible with noise injection because skipped LOS runs would change the shared random-number sequence."; code=:invalid_resume)

    raw_outputs = get(cfg, "outputs", ["all"])
    raw_outputs isa AbstractVector || throw_config_error("`outputs` must be an array of output groups."; code=:invalid_outputs)
    outputs = unique(lowercase.(String.(raw_outputs)))
    allowed_outputs = Set(["all", "integrated", "stokes", "rm", "fdf", "spectral_index", "diagnostics"])
    invalid_outputs = [name for name in outputs if !(name in allowed_outputs)]
    isempty(invalid_outputs) || throw_config_error("Unknown output groups: $(join(invalid_outputs, ", "))."; code=:invalid_outputs)
    isempty(outputs) && throw_config_error("`outputs` must select at least one output group."; code=:invalid_outputs)
    all_outputs = "all" in outputs
    selected_outputs = all_outputs ? Set(["integrated", "stokes", "rm", "fdf", "spectral_index", "diagnostics"]) : Set(outputs)
    !all_outputs && !isempty(intersect(selected_outputs, Set(["rm", "fdf"]))) && FaradayRotation != "Y" && throw_config_error(
        "The `rm` and `fdf` output groups require Faraday rotation."; code=:invalid_outputs)
    rm_clean_enabled && !("fdf" in selected_outputs) && throw_config_error(
        "RM-CLEAN requires the `fdf` output group."; code=:invalid_outputs)
    tile_size !== nothing && selected_outputs == Set(["diagnostics"]) && throw_config_error(
        "The `diagnostics` output group cannot be used alone with `tile_size`: tiled processing does not retain full Q/U/T cubes for plots."; code=:invalid_outputs)

    return RunConfig(
        base_dir,
        simu_paths,
        chosen_LOS,
        conversionB,
        conversionn,
        conversionT,
        FaradayRotation,
        phimin,
        phimax,
        dphi,
        responseSynchrotron,
        kernel_size_synchrotron,
        add_noise,
        SNR_nu,
        interpolation_file_path,
        ne_option,
        IonizationFraction,
        wolfire_constants,
        nustart,
        nuend,
        dnu,
        BoxLength_pc,
        BoxLength_pix,
        config_path,
        log_progress,
        rng_seed;
        rm_clean_enabled = rm_clean_enabled,
        rm_clean_gain = rm_clean_gain,
        rm_clean_niter = rm_clean_niter,
        rm_clean_threshold = rm_clean_threshold,
        precision = precision,
        tile_size = tile_size,
        field_sources = field_sources,
        physical_mask = physical_mask,
        density_kind = density_kind,
        mean_molecular_weight = mean_molecular_weight,
        hydrogen_mass_g = hydrogen_mass_g,
        resume = resume,
        outputs = outputs,
        rfi_ranges = rfi_ranges,
    ), simu_paths
end

function collect_wolfire_constants(cfg)
    ne_cfg = get(cfg, "ne", nothing)
    zeta = ne_cfg isa AbstractDict ? get(ne_cfg, "zeta", get(cfg, "zeta", nothing)) : get(cfg, "zeta", nothing)
    Geff = ne_cfg isa AbstractDict ? get(ne_cfg, "Geff", get(cfg, "Geff", nothing)) : get(cfg, "Geff", nothing)
    phiPAH = ne_cfg isa AbstractDict ? get(ne_cfg, "phiPAH", get(cfg, "phiPAH", nothing)) : get(cfg, "phiPAH", nothing)
    XC = ne_cfg isa AbstractDict ? get(ne_cfg, "XC", get(cfg, "XC", nothing)) : get(cfg, "XC", nothing)

    values = (zeta, Geff, phiPAH, XC)
    all_missing = all(value === nothing for value in values)
    any_missing = any(value === nothing for value in values)

    if all_missing
        return nothing
    elseif any_missing
        throw_config_error("Wolfire constants must include zeta, Geff, phiPAH, and XC when provided."; code=:invalid_wolfire_constants)
    end

    return (
        parse_config_float(zeta, "zeta"),
        parse_config_float(Geff, "Geff"),
        parse_config_float(phiPAH, "phiPAH"),
        parse_config_float(XC, "XC"),
    )
end

function MOOSE_from_config_dict(
    cfg::AbstractDict;
    config_path::AbstractString = "<in-memory>",
    quiet::Bool = false,
    source_config_path = config_path,
    saved_config_path = config_path,
    write_config_file::Bool = true,
)
    run_config, _ = build_config(cfg, config_path)
    run_moose_processing(
        run_config;
        quiet = quiet,
        persisted_config = cfg,
        source_config_path = source_config_path,
        saved_config_path = saved_config_path,
        write_config_file = write_config_file,
    )

    return nothing
end

function MOOSE_from_config(config_path::AbstractString; quiet::Bool = false)
    cfg = JSON.parsefile(config_path)
    MOOSE_from_config_dict(
        cfg;
        config_path = config_path,
        quiet = quiet,
        source_config_path = config_path,
        saved_config_path = config_path,
        write_config_file = true,
    )

    return nothing
end

end # module
