"""
    QUnu(Bperp::AbstractArray, psi_src::AbstractArray, RM::AbstractArray, nuArray::AbstractArray, df::DataFrame, PixelLength_cm::Float64) -> Tuple{AbstractArray, AbstractArray}

Calculate the Stokes parameters Q and U as functions of frequency.

# Arguments
- `Bperp::AbstractArray`: Array representing the perpendicular component of the magnetic field.
- `psi_src::AbstractArray`: Array representing the source polarization angles.
- `RM::AbstractArray`: Array representing the Rotation Measure.
- `nuArray::AbstractArray`: Array of frequencies at which to compute the Stokes parameters.
- `df::DataFrame`: DataFrame containing columns `B` (magnetic field values), `nu` (frequency values), `e_perp` (perpendicular emissivity), and `e_para` (parallel emissivity).
- `PixelLength_cm::Float64`: The pixel length in centimeters.

# Returns
- `Tuple{AbstractArray, AbstractArray}`: Two arrays, Qnu and Unu, representing the Stokes parameters Q and U as functions of frequency.

# Example
```julia
using DataFrames

## The dataframe "df" should be the dataframe computed by the emissivity interpolation code

# Example input arrays
Bperp = [1.0, 2.0]
psi_src = [0.1, 0.2]
RM = [0.001, 0.002]
nuArray = [1e9, 1.1e9]
PixelLength_cm = 1.0

# Function call
Qnu, Unu = QUnu(Bperp, psi_src, RM, nuArray, df, PixelLength_cm)

"""

"""
    emissivity_grid(df, values) -> (B, nu, eps)

Reorder a flat emissivity table into the matrix layout expected by
`Spline2D(x, y, z)`, i.e. `eps[i, j] = f(B[i], nu[j])` with `B` and `nu`
sorted in increasing order. The table must contain exactly one row per
(B, ν) pair of the full cartesian product — any missing, duplicated or
extra row raises an explicit error instead of silently permuting
emissivities.
"""
function emissivity_grid(df::DataFrame, values::AbstractVector)
    B = sort(unique(df.B))
    nu = sort(unique(df.nu))

    nrow(df) == length(B) * length(nu) || error(
        "Emissivity table is not a complete (B, nu) grid: $(nrow(df)) rows ≠ " *
        "$(length(B)) B values × $(length(nu)) nu values.")
    allunique(zip(df.B, df.nu)) || error(
        "Emissivity table contains duplicated (B, nu) pairs.")

    # Sort rows with nu as the outer (slow) key and B as the inner (fast) key,
    # so the reshape below fills the (B, nu) matrix in column-major order.
    perm = sortperm(collect(zip(df.nu, df.B)))
    eps = reshape(Float64.(values[perm]), (length(B), length(nu)))

    return Float64.(B), Float64.(nu), eps
end

struct EmissivityInterpolator
    B::Vector{Float64}
    eps_interp::Spline2D
end

function EmissivityInterpolator(df::DataFrame)
    B, nu, eps = emissivity_grid(df, df.e_perp .- df.e_para)
    eps_interp = Spline2D(B, nu, eps)
    return EmissivityInterpolator(B, eps_interp)
end

function build_emissivity_frequency_cache(interpolator::EmissivityInterpolator, nuArray)
    Nfreq = length(nuArray)
    B = interpolator.B
    cache = Matrix{Float64}(undef, length(B), Nfreq)
    @inbounds for i in 1:Nfreq
        nui = nuArray[i]
        for j in eachindex(B)
            cache[j, i] = interpolator.eps_interp(B[j], nui)
        end
    end
    return cache
end

@inline function linear_interp_extrapolated(xgrid::Vector{Float64}, ygrid::AbstractVector{Float64}, x::Float64)
    n = length(xgrid)
    @inbounds if x <= xgrid[1]
        x1, x2 = xgrid[1], xgrid[2]
        y1, y2 = ygrid[1], ygrid[2]
        return y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    elseif x >= xgrid[n]
        x1, x2 = xgrid[n - 1], xgrid[n]
        y1, y2 = ygrid[n - 1], ygrid[n]
        return y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    else
        idx = searchsortedlast(xgrid, x)
        x1, x2 = xgrid[idx], xgrid[idx + 1]
        y1, y2 = ygrid[idx], ygrid[idx + 1]
        return y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    end
end

@inline function linear_interp_extrapolated(xgrid::Vector{Float64}, ygrid::AbstractMatrix{Float64}, col::Int, x::Float64)
    n = length(xgrid)
    @inbounds if x <= xgrid[1]
        x1, x2 = xgrid[1], xgrid[2]
        y1, y2 = ygrid[1, col], ygrid[2, col]
        return y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    elseif x >= xgrid[n]
        x1, x2 = xgrid[n - 1], xgrid[n]
        y1, y2 = ygrid[n - 1, col], ygrid[n, col]
        return y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    else
        idx = searchsortedlast(xgrid, x)
        x1, x2 = xgrid[idx], xgrid[idx + 1]
        y1, y2 = ygrid[idx, col], ygrid[idx + 1, col]
        return y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    end
end

@inline function _linear_interp_at_index(xgrid::Vector{Float64}, ygrid::AbstractMatrix{Float64},
                                         col::Int, x::Float64, idx::Int)
    @inbounds begin
        x1, x2 = xgrid[idx], xgrid[idx + 1]
        y1, y2 = ygrid[idx, col], ygrid[idx + 1, col]
        return y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    end
end

# `Threads.@threads` may migrate a task between threads mid-iteration, so
# `threadid()` is not a stable index into per-thread scratch space: two tasks can
# end up sharing one buffer and silently corrupt each other's sightline. Partition
# the pixel grid instead and let every task own its buffers for a whole chunk.
function _threaded_pixel_chunks(npixels::Integer)
    npixels <= 0 && return UnitRange{Int}[]
    n = Int(npixels)
    nchunks = max(1, min(Threads.nthreads(), n))
    return collect(Iterators.partition(1:n, cld(n, nchunks)))
end

function _emissivity_brackets!(indices, xgrid::Vector{Float64}, values)
    last = length(xgrid) - 1
    @inbounds for k in eachindex(indices, values)
        x = Float64(values[k])
        indices[k] = x <= xgrid[1] ? 1 : x >= xgrid[end] ? last : searchsortedlast(xgrid, x)
    end
    return indices
end

function emissivity_at_frequency!(buffer, B::Vector{Float64}, eps_interp::Spline2D, Bperp::AbstractArray, nui;
    eps_cache_col=nothing, eps_line_buffer=nothing)
    eps_i = eps_cache_col
    if eps_i === nothing
        eps_i = eps_line_buffer
        @inbounds for j in eachindex(B)
            eps_i[j] = eps_interp(B[j], nui)
        end
    end
    @inbounds for idx in eachindex(Bperp, buffer)
        buffer[idx] = linear_interp_extrapolated(B, eps_i, Float64(Bperp[idx]))
    end
    return buffer
end

function _QUnu!(Qnu, Unu, Bperp, psi_src, RM, nuArray, PixelLength_cm, interpolator;
                 emissivity_cache=nothing, interp_indices=nothing)
    Nfreq = length(nuArray)
    eps_line_buffer = emissivity_cache === nothing ? similar(interpolator.B, Float64) : nothing
    indices = interp_indices === nothing ? similar(Bperp, Int) : interp_indices
    emissivity_cache === nothing || _emissivity_brackets!(indices, interpolator.B, Bperp)

    for i in 1:Nfreq
        nui = nuArray[i]
        if emissivity_cache === nothing
            @inbounds for j in eachindex(interpolator.B)
                eps_line_buffer[j] = interpolator.eps_interp(interpolator.B[j], nui)
            end
        end

        faraday_factor = (C_m / (nui * 1e6))^2
        sum_u = 0.0
        sum_q = 0.0
        @inbounds for idx in eachindex(Bperp, psi_src, RM)
            arg = 2.0 * (psi_src[idx] + RM[idx] * faraday_factor)
            eps_val = emissivity_cache === nothing ?
                linear_interp_extrapolated(interpolator.B, eps_line_buffer, Float64(Bperp[idx])) :
                _linear_interp_at_index(interpolator.B, emissivity_cache, i, Float64(Bperp[idx]), indices[idx])
            # One argument reduction for both quadratures instead of two.
            sin_arg, cos_arg = sincos(arg)
            sum_u += eps_val * sin_arg
            sum_q += eps_val * cos_arg
        end

        Unu[i] = BrightnessTemperature(nui, sum_u * PixelLength_cm)
        Qnu[i] = BrightnessTemperature(nui, sum_q * PixelLength_cm)
    end

    return Qnu, Unu
end

function QUnu(Bperp, psi_src, RM, nuArray, df, PixelLength_cm; precomputed_interp=nothing, emissivity_cache=nothing)
    interpolator = precomputed_interp === nothing ? EmissivityInterpolator(df) : precomputed_interp
    Nfreq = length(nuArray)
    Qnu = zeros(Nfreq)
    Unu = zeros(Nfreq)
    return _QUnu!(Qnu, Unu, Bperp, psi_src, RM, nuArray, PixelLength_cm, interpolator; emissivity_cache=emissivity_cache)
end

"""
    QUnu3D(Bperpcube::AbstractArray, psi_src::AbstractArray, RM::AbstractArray, nuArray::AbstractArray, df::DataFrame, PixelLength_cm::Float64) -> Tuple{AbstractArray, AbstractArray}

Calculate the Stokes parameters Q and U for a 3D cube as functions of frequency.

# Arguments
- `Bperpcube::AbstractArray`: 3D array representing the perpendicular component of the magnetic field for each pixel in the cube.
- `psi_src::AbstractArray`: 3D array representing the source polarization angles for each pixel in the cube.
- `RM::AbstractArray`: 3D array representing the Rotation Measure for each pixel in the cube.
- `nuArray::AbstractArray`: Array of frequencies at which to compute the Stokes parameters.
- `df::DataFrame`: DataFrame containing columns `B` (magnetic field values), `nu` (frequency values), `e_perp` (perpendicular emissivity), and `e_para` (parallel emissivity).
- `PixelLength_cm::Float64`: The pixel length in centimeters.

# Returns
- `Tuple{AbstractArray, AbstractArray}`: Two 3D arrays, Qnu and Unu, representing the Stokes parameters Q and U as functions of frequency for each pixel in the cube.

# Example
```julia
using DataFrames

## The dataframe "df" should be the dataframe computed by the emissivity interpolation code

# Example input arrays
Bperpcube = rand(10, 10, 2)  # 10x10 pixels, 2 depth slices
psi_src = rand(10, 10, 2)
RM = rand(10, 10, 2)
nuArray = [1e9, 1.1e9]
PixelLength_cm = 1.0

# Function call
Qnu, Unu = QUnu3D(Bperpcube, psi_src, RM, nuArray, df, PixelLength_cm)
"""
function QUnu3D(Bperpcube, psi_src, RM, nuArray, df, PixelLength_cm; log_progress::Bool=false)
    Nfreq = length(nuArray)
    nx, ny = size(Bperpcube, 1), size(Bperpcube, 2)
    # Output cubes follow the working precision of the input cube (the
    # per-pixel accumulation below still runs in Float64 scalars).
    T = float(eltype(Bperpcube))
    Qnu = zeros(T, nx, ny, Nfreq)
    Unu = zeros(T, nx, ny, Nfreq)
    interpolator = EmissivityInterpolator(df)
    emissivity_cache = build_emissivity_frequency_cache(interpolator, nuArray)
    depth = size(Bperpcube, 3)
    pixels = CartesianIndices((1:nx, 1:ny))
    total_pixels = length(pixels)
    progress_counter = Threads.Atomic{Int}(0)
    progress_step = max(floor(Int, total_pixels / 100), 1)
    progress_lock = ReentrantLock()

    @sync for part in _threaded_pixel_chunks(total_pixels)
        Threads.@spawn begin
            indices = Vector{Int}(undef, depth)
            for p in part
                i, j = Tuple(pixels[p])
                @views Bperp_vec = Bperpcube[i, j, :]
                @views RM_vec = RM[i, j, :]
                @views psi_src_vec = psi_src[i, j, :]
                @views qdest = Qnu[i, j, :]
                @views udest = Unu[i, j, :]
                _QUnu!(qdest, udest, Bperp_vec, psi_src_vec, RM_vec, nuArray, PixelLength_cm,
                    interpolator; emissivity_cache=emissivity_cache, interp_indices=indices)

                if log_progress
                    done = Threads.atomic_add!(progress_counter, 1) + 1
                    if done % progress_step == 0
                        lock(progress_lock) do
                            print_progress(done, total_pixels; label="Computing Q/U (Faraday)")
                        end
                    end
                end
            end
        end
    end

    if log_progress && total_pixels > 0
        print_progress(total_pixels, total_pixels; label="Computing Q/U (Faraday)")
    end

    return Qnu, Unu
end

"""
    QUnuNoFaraday(Bperp::AbstractArray, psi_src::AbstractArray, nuArray::AbstractArray, df::DataFrame, PixelLength_cm::Float64) -> Tuple{AbstractArray, AbstractArray}

Calculate the Stokes parameters Q and U as functions of frequency without considering Faraday rotation.

# Arguments
- `Bperp::AbstractArray`: Array representing the perpendicular component of the magnetic field.
- `psi_src::AbstractArray`: Array representing the source polarization angles.
- `nuArray::AbstractArray`: Array of frequencies at which to compute the Stokes parameters.
- `df::DataFrame`: DataFrame containing columns `B` (magnetic field values), `nu` (frequency values), `e_perp` (perpendicular electric field), and `e_para` (parallel electric field).
- `PixelLength_cm::Float64`: The pixel length in centimeters.

# Returns
- `Tuple{AbstractArray, AbstractArray}`: Two arrays, Qnu and Unu, representing the Stokes parameters Q and U as functions of frequency.

# Example
```julia
using DataFrames

## The dataframe "df" should be the dataframe computed by the emissivity interpolation code

# Example input arrays
Bperp = [1.0, 2.0]
psi_src = [0.1, 0.2]
nuArray = [1e9, 1.1e9]
PixelLength_cm = 1.0

# Function call
Qnu, Unu = QUnuNoFaraday(Bperp, psi_src, nuArray, df, PixelLength_cm)
"""
function _QUnuNoFaraday!(Qnu, Unu, Bperp, psi_src, nuArray, PixelLength_cm, interpolator;
                          emissivity_cache=nothing, interp_indices=nothing,
                          sin_angles=nothing, cos_angles=nothing)
    Nfreq = length(nuArray)
    eps_line_buffer = emissivity_cache === nothing ? similar(interpolator.B, Float64) : nothing
    indices = interp_indices === nothing ? similar(Bperp, Int) : interp_indices
    sins = sin_angles === nothing ? similar(Bperp, Float64) : sin_angles
    coss = cos_angles === nothing ? similar(Bperp, Float64) : cos_angles
    if emissivity_cache !== nothing
        _emissivity_brackets!(indices, interpolator.B, Bperp)
    end
    @inbounds for idx in eachindex(psi_src, sins, coss)
        sins[idx], coss[idx] = sincos(2.0 * psi_src[idx])
    end

    for i in 1:Nfreq
        nui = nuArray[i]
        if emissivity_cache === nothing
            @inbounds for j in eachindex(interpolator.B)
                eps_line_buffer[j] = interpolator.eps_interp(interpolator.B[j], nui)
            end
        end

        sum_u = 0.0
        sum_q = 0.0
        @inbounds for idx in eachindex(Bperp, psi_src)
            eps_val = emissivity_cache === nothing ?
                linear_interp_extrapolated(interpolator.B, eps_line_buffer, Float64(Bperp[idx])) :
                _linear_interp_at_index(interpolator.B, emissivity_cache, i, Float64(Bperp[idx]), indices[idx])
            sum_u += eps_val * sins[idx]
            sum_q += eps_val * coss[idx]
        end

        Unu[i] = BrightnessTemperature(nui, sum_u * PixelLength_cm)
        Qnu[i] = BrightnessTemperature(nui, sum_q * PixelLength_cm)
    end

    return Qnu, Unu
end

function QUnuNoFaraday(Bperp, psi_src, nuArray, df, PixelLength_cm; precomputed_interp=nothing, emissivity_cache=nothing)

    interpolator = precomputed_interp === nothing ? EmissivityInterpolator(df) : precomputed_interp

    Nfreq = length(nuArray)
    Qnu = zeros(Nfreq)
    Unu = zeros(Nfreq)
    return _QUnuNoFaraday!(Qnu, Unu, Bperp, psi_src, nuArray, PixelLength_cm, interpolator; emissivity_cache=emissivity_cache)

end

"""
    QUnuNoFaraday3D(Bperpcube::AbstractArray, psi_src::AbstractArray, nuArray::AbstractArray, df::DataFrame, PixelLength_cm::Float64) -> Tuple{AbstractArray, AbstractArray}

Calculate the Stokes parameters Q and U for a 3D cube as functions of frequency without considering Faraday rotation.

# Arguments
- `Bperpcube::AbstractArray`: 3D array representing the perpendicular component of the magnetic field for each pixel in the cube.
- `psi_src::AbstractArray`: 3D array representing the source polarization angles for each pixel in the cube.
- `nuArray::AbstractArray`: Array of frequencies at which to compute the Stokes parameters.
- `df::DataFrame`: DataFrame containing columns `B` (magnetic field values), `nu` (frequency values), `e_perp` (perpendicular emissivity), and `e_para` (parallel emissivity).
- `PixelLength_cm::Float64`: The pixel length in centimeters.

# Returns
- `Tuple{AbstractArray, AbstractArray}`: Two 3D arrays, Qnu and Unu, representing the Stokes parameters Q and U as functions of frequency for each pixel in the cube.

# Example
```julia
using DataFrames

## The dataframe "df" should be the dataframe computed by the emissivity interpolation code

# Example input arrays
Bperpcube = rand(10, 10, 2)  # 10x10 pixels, 2 depth slices
psi_src = rand(10, 10, 2)
nuArray = [1e9, 1.1e9]
PixelLength_cm = 1.0

# Function call
Qnu, Unu = QUnuNoFaraday3D(Bperpcube, psi_src, nuArray, df, PixelLength_cm)
"""
function QUnuNoFaraday3D(Bperpcube, psi_src, nuArray, df, PixelLength_cm; log_progress::Bool=false)

    Nfreq = length(nuArray)
    T = float(eltype(Bperpcube))
    Qnu = zeros(T, size(Bperpcube,1), size(Bperpcube,2), Nfreq)
    Unu = zeros(T, size(Bperpcube,1), size(Bperpcube,2), Nfreq)
    interpolator = EmissivityInterpolator(df)
    emissivity_cache = build_emissivity_frequency_cache(interpolator, nuArray)
    depth = size(Bperpcube, 3)
    pixels = CartesianIndices((1:size(Bperpcube,1), 1:size(Bperpcube,2)))
    total_pixels = length(pixels)
    progress_counter = Threads.Atomic{Int}(0)
    progress_step = max(floor(Int, total_pixels / 100), 1)
    progress_lock = ReentrantLock()

    @sync for part in _threaded_pixel_chunks(total_pixels)
        Threads.@spawn begin
            indices = Vector{Int}(undef, depth)
            sins = Vector{Float64}(undef, depth)
            coss = Vector{Float64}(undef, depth)
            for p in part
                i, j = Tuple(pixels[p])
                @views Bperp_vec = Bperpcube[i,j,:]
                @views psi_src_vec = psi_src[i,j,:]
                @views qdest = Qnu[i, j, :]
                @views udest = Unu[i, j, :]
                _QUnuNoFaraday!(qdest, udest, Bperp_vec, psi_src_vec, nuArray, PixelLength_cm,
                    interpolator; emissivity_cache=emissivity_cache,
                    interp_indices=indices, sin_angles=sins, cos_angles=coss)

                if log_progress
                    done = Threads.atomic_add!(progress_counter, 1) + 1
                    if done % progress_step == 0
                        lock(progress_lock) do
                            print_progress(done, total_pixels; label="Computing Q/U")
                        end
                    end
                end
            end
        end
    end

    if log_progress && total_pixels > 0
        print_progress(total_pixels, total_pixels; label="Computing Q/U")
    end

    return Qnu, Unu

end
