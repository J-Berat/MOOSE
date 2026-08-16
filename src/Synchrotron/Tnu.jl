"""
    Tnu(Bperp::AbstractArray, nuArray::AbstractArray, df::DataFrame, PixelLength_cm::Float64) -> AbstractArray

Calculate the brightness temperature T as a function of frequency.

# Arguments
- `Bperp::AbstractArray`: Array representing the perpendicular component of the magnetic field.
- `nuArray::AbstractArray`: Array of frequencies at which to compute the brightness temperature.
- `df::DataFrame`: DataFrame containing columns `B` (magnetic field values), `nu` (frequency values), `e_perp` (perpendicular emissivity), and `e_para` (parallel emissivity).
- `PixelLength_cm::Float64`: The pixel length in centimeters.

# Returns
- `AbstractArray`: An array representing the brightness temperature T as a function of frequency.

# Example
```julia
using DataFrames

## The dataframe "df" should be the dataframe computed by the emissivity interpolation code

# Example input arrays
Bperp = [1.0, 2.0]
nuArray = [1e9, 1.1e9]
PixelLength_cm = 1.0

# Function call
T_nu = Tnu(Bperp, nuArray, df, PixelLength_cm)
"""

struct TemperatureInterpolator
    B::Vector{Float64}
    eps_interp::Spline2D
end

function TemperatureInterpolator(df::DataFrame)
    B, nu, eps = emissivity_grid(df, df.e_para .+ df.e_perp)
    eps_interp = Spline2D(B, nu, eps)
    return TemperatureInterpolator(B, eps_interp)
end

function build_emissivity_frequency_cache(interpolator::TemperatureInterpolator, nuArray)
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

function emissivity_total_at_frequency!(buffer, B::Vector{Float64}, eps_interp::Spline2D, Bperp::AbstractArray, nui;
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

function _Tnu!(T_nu, Bperp, nuArray, PixelLength_cm, interpolator;
              emissivity_cache=nothing, interp_indices=nothing)
    B = interpolator.B
    eps_line_buffer = emissivity_cache === nothing ? similar(B, Float64) : nothing
    # The bracketing interval of each cell in the B grid is frequency independent,
    # so the binary search runs once per sightline rather than once per channel.
    indices = interp_indices === nothing ? similar(Bperp, Int) : interp_indices
    emissivity_cache === nothing || _emissivity_brackets!(indices, B, Bperp)

    for i in eachindex(nuArray)
        nui = nuArray[i]
        if emissivity_cache === nothing
            @inbounds for j in eachindex(B)
                eps_line_buffer[j] = interpolator.eps_interp(B[j], nui)
            end
        end

        # Accumulate the line-of-sight emissivity directly: the intermediate
        # per-cell buffer was only ever consumed by `sum`.
        Inui = 0.0
        @inbounds for idx in eachindex(Bperp)
            Inui += emissivity_cache === nothing ?
                linear_interp_extrapolated(B, eps_line_buffer, Float64(Bperp[idx])) :
                _linear_interp_at_index(B, emissivity_cache, i, Float64(Bperp[idx]), indices[idx])
        end

        T_nu[i] = BrightnessTemperature(nui, Inui * PixelLength_cm)
    end

    return T_nu
end

function Tnu(Bperp, nuArray, df, PixelLength_cm; precomputed_interp = nothing, emissivity_cache=nothing)
    interpolator = precomputed_interp === nothing ? TemperatureInterpolator(df) : precomputed_interp
    Nfreq = length(nuArray)
    T_nu = zeros(Nfreq)
    return _Tnu!(T_nu, Bperp, nuArray, PixelLength_cm, interpolator; emissivity_cache=emissivity_cache)

end

"""
    Tnu3D(Bperpcube::AbstractArray, nuArray::AbstractArray, df::DataFrame, PixelLength_cm::Float64) -> AbstractArray

Calculate the brightness temperature T for a 3D cube as a function of frequency.

# Arguments
- `Bperpcube::AbstractArray`: 3D array representing the perpendicular component of the magnetic field for each pixel in the cube.
- `nuArray::AbstractArray`: Array of frequencies at which to compute the brightness temperature.
- `df::DataFrame`: DataFrame containing columns `B` (magnetic field values), `nu` (frequency values), `e_perp` (perpendicular electric field), and `e_para` (parallel electric field).
- `PixelLength_cm::Float64`: The pixel length in centimeters.

# Returns
- `AbstractArray`: A 3D array representing the brightness temperature T as a function of frequency for each pixel in the cube.

# Example
```julia
using DataFrames

## The dataframe "df" should be the dataframe computed by the emissivity interpolation code

# Example input arrays
Bperpcube = rand(10, 10, 2)  # 10x10 pixels, 2 depth slices
nuArray = [1e9, 1.1e9]
PixelLength_cm = 1.0

# Function call
T_nu = Tnu3D(Bperpcube, nuArray, df, PixelLength_cm)
"""

function Tnu3D(Bperpcube, nuArray, df, PixelLength_cm)
    nx, ny = size(Bperpcube, 1), size(Bperpcube, 2)
    Nfreq = length(nuArray)
    # Output cube follows the working precision of the input cube (the
    # per-pixel accumulation still runs in Float64 scalars).
    T_nu = zeros(float(eltype(Bperpcube)), nx, ny, Nfreq)
    interpolator = TemperatureInterpolator(df)
    emissivity_cache = build_emissivity_frequency_cache(interpolator, nuArray)

    depth = size(Bperpcube, 3)
    pixels = CartesianIndices((1:nx, 1:ny))

    @sync for part in _threaded_pixel_chunks(length(pixels))
        Threads.@spawn begin
            indices = Vector{Int}(undef, depth)
            for p in part
                i, j = Tuple(pixels[p])
                @views Bperp_vec = Bperpcube[i, j, :]
                @views tdest = T_nu[i, j, :]
                _Tnu!(tdest, Bperp_vec, nuArray, PixelLength_cm, interpolator;
                    emissivity_cache=emissivity_cache, interp_indices=indices)
            end
        end
    end

    return T_nu
end
