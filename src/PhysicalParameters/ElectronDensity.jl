"""
    constant_ne(constant::Float64, SizeCube::Tuple{Int, Int, Int}) -> Array{Float64, 3}

Create a 3D array with constant electron density.

# Arguments
- `constant::Float64`: The constant value for the electron density.
- `SizeCube::Tuple{Int, Int, Int}`: The size of the 3D array.

# Returns
- `Array{Float64, 3}`: A 3D array with all elements set to the specified constant value.

# Description
This function creates a 3D array with the specified size, where all elements are set to the given constant value. This can be used to represent a homogeneous electron density distribution.

# Example
```julia
# Example usage
constant_value = 1.0e-3  # Example constant electron density
size_cube = (100, 100, 100)  # Size of the 3D array
electron_density = constant_ne(constant_value, size_cube)
println(electron_density)
"""
constant_ne(constant::Float64, SizeCube::Tuple{Int, Int, Int}) = fill(constant, SizeCube)

"""
    ne_propto_nH(n::AbstractArray, IonizationFraction::Float64) -> AbstractArray

Calculate the electron density proportional to the hydrogen density and ionization fraction.

# Arguments
- `n::AbstractArray`: The array representing hydrogen density.
- `IonizationFraction::Float64`: The constant ionization fraction.

# Returns
- `AbstractArray`: An array representing the electron density.

# Description
This function calculates the electron density (ne) as a product of the hydrogen density (n) and a constant ionization fraction. The formula used is:
n_e = IonizationFraction x n 
where:
- n_e is the electron density.
- n is the hydrogen density.
- IonizationFraction is the constant ionization fraction.

# Example
```julia
# Example usage
hydrogen_density = rand(100, 100, 100)  # Example hydrogen density data
ionization_fraction = 0.1  # Example ionization fraction
electron_density = ne_propto_nH(hydrogen_density, ionization_fraction)
println(electron_density)
"""
ne_propto_nH(n::AbstractArray, IonizationFraction::Float64) = IonizationFraction .* n
"""
    Wolfire_ne(zeta::Float64, Geff::Float64, omegaPAH::Float64, XC::Float64, T::AbstractArray, n::AbstractArray) -> AbstractArray

Calculate the electron density using the Wolfire et al. 2003 model.

# Arguments
- `zeta::Float64`: Cosmic-ray ionization rate.
- `Geff::Float64`: Effective FUV field strength.
- `omegaPAH::Float64`: PAH abundance factor.
- `XC::Float64`: Conversion factor from hydrogen density to electron density.
- `T::AbstractArray`: Array representing the gas temperature.
- `n::AbstractArray`: Array representing the hydrogen density.

# Returns
- `AbstractArray`: An array representing the electron density calculated using the Wolfire et al. model.

# Description
This function calculates the electron density (ne) using the Wolfire et al. model. The formula used is:
 n_e = 2.4e-3*sqrt(zeta/1e-16)*(T/100)^(0.25)*sqrt(Geff)/omegaPAH+n*XC
where :
- n_e is the electron density.
- zeta is the cosmic-ray ionization rate.
- Geff is the effective FUV field strength.
- omegaPAH is the PAH abundance factor.
- XC is the conversion factor from hydrogen density to Carbon.
- T is the gas temperature.
- n is the hydrogen density.

# Example
```julia
# Example usage
zeta = 1.0e-16  # Example cosmic-ray ionization rate
Geff = 1.0  # Example effective FUV field strength
omegaPAH = 1.0  # Example PAH abundance factor
XC = 0.1  # Example conversion factor
T = rand(100, 100, 100)  # Example temperature data
n = rand(100, 100, 100)  # Example hydrogen density data

electron_density = Wolfire_ne(zeta, Geff, omegaPAH, XC, T, n)
println(electron_density)
"""
Wolfire_ne(zeta::Float64, Geff::Float64, omegaPAH::Float64, XC::Float64, T::AbstractArray, n::AbstractArray) = @. 2.4e-3*sqrt(zeta/1e-16)*(T/100)^(0.25)*sqrt(Geff)/omegaPAH+n*XC

"""
    DM(ne::AbstractArray{T, 1}, PixelLength_pc::Float64) -> Float64
    DM(ne::AbstractArray{T, 3}, PixelLength_pc::Float64) -> AbstractArray{Float64, 2}

Calculate the dispersion measure (DM) for electron density along the line of sight.

# Arguments
- `ne::AbstractArray{T, 1}`: A 1D array representing the electron density along the line of sight.
- OR `ne::AbstractArray{T, 3}`: A 3D array representing the electron density in a volume.
- `PixelLength_pc::Float64`: The length of a pixel in parsecs.

# Returns
- For 1D input: `Float64` - The total dispersion measure along the line of sight.
- For 3D input: `AbstractArray{Float64, 2}` - A 2D array representing the dispersion measure along the third dimension.

# Description
This function calculates the dispersion measure (DM), which is the integral of the electron density (ne) along the line of sight, multiplied by the pixel length in parsecs. The function is overloaded to handle both 1D and 3D input arrays.

For a 1D array, the dispersion measure is calculated as:
DM = sum(ne .* PixelLength_pc)

For a 3D array, the dispersion measure is calculated along the third dimension and the resulting 2D array is returned:
DM = sum(ne .* PixelLength_pc, dims=3)

# Example
```julia
# Example usage for 1D array
ne_1d = rand(100)  # Example electron density data
PixelLength_pc = 0.1  # Example pixel length in parsecs
dm_1d = DM(ne_1d, PixelLength_pc)
println("1D DM: ", dm_1d)

# Example usage for 3D array
ne_3d = rand(100, 100, 100)  # Example 3D electron density data
PixelLength_pc = 0.1  # Example pixel length in parsecs
dm_3d = DM(ne_3d, PixelLength_pc)
println("3D DM: ", dm_3d)
"""
function _integrate_density_vector(ne::AbstractVector, PixelLength_pc::Real, transform;
                                   nan_policy::Symbol)
    _check_nan_policy(nan_policy)
    isfinite(PixelLength_pc) || throw(ArgumentError("PixelLength_pc must be finite."))
    if nan_policy === :error
        all(isfinite, ne) || throw(ArgumentError("Density contains non-finite values."))
    end
    nan_policy !== :omit && return sum(x -> transform(x) * PixelLength_pc, ne)

    values = (transform(x) * PixelLength_pc for x in ne if isfinite(x))
    result = sum(values; init=0.0)
    any(isfinite, ne) ? result : NaN
end

DM(ne::AbstractVector, PixelLength_pc; nan_policy::Symbol=:propagate) =
    _integrate_density_vector(ne, PixelLength_pc, identity; nan_policy)

DM(ne::AbstractArray{<:Real, 3}, PixelLength_pc; nan_policy::Symbol=:propagate) =
    intLOS(ne, PixelLength_pc; nan_policy)

"""
    EM(ne::Vector{T}, PixelLength_pc::Float64) -> Float64
    EM(ne::Array{T, 3}, PixelLength_pc::Float64) -> AbstractArray{Float64, 2}

Calculate the emission measure (EM) for electron density along the line of sight.

# Arguments
- `ne::Vector{T}`: A 1D array representing the electron density along the line of sight.
- OR `ne::Array{T, 3}`: A 3D array representing the electron density in a volume.
- `PixelLength_pc::Float64`: The length of a pixel in parsecs.

# Returns
- For 1D input: `Float64` - The total emission measure along the line of sight.
- For 3D input: `AbstractArray{Float64, 2}` - A 2D array representing the emission measure along the third dimension.

# Description
This function calculates the emission measure (EM), which is the integral of the square of the electron density (ne) along the line of sight, multiplied by the pixel length in parsecs. The function is overloaded to handle both 1D and 3D input arrays.

For a 1D array, the emission measure is calculated as:
EM = sum(ne .^ 2 .* PixelLength_pc)

For a 3D array, the emission measure is calculated along the third dimension and the resulting 2D array is returned:
EM = sum(ne .^ 2 .* PixelLength_pc, dims=3)

# Example
```julia
# Example usage for 1D array
ne_1d = rand(100)  # Example electron density data
PixelLength_pc = 0.1  # Example pixel length in parsecs
em_1d = EM(ne_1d, PixelLength_pc)
println("1D EM: ", em_1d)

# Example usage for 3D array
ne_3d = rand(100, 100, 100)  # Example 3D electron density data
em_3d = EM(ne_3d, PixelLength_pc)
println("3D EM: ", em_3d)
"""
EM(ne::AbstractVector, PixelLength_pc; nan_policy::Symbol=:propagate) =
    _integrate_density_vector(ne, PixelLength_pc, abs2; nan_policy)

function EM(ne::AbstractArray{<:Real, 3}, PixelLength_pc;
            nan_policy::Symbol=:propagate)
    _check_nan_policy(nan_policy)
    isfinite(PixelLength_pc) || throw(ArgumentError("PixelLength_pc must be finite."))
    if nan_policy === :error
        all(isfinite, ne) || throw(ArgumentError("Density contains non-finite values."))
    end
    nan_policy !== :omit &&
        return dropdims(sum(x -> abs2(x) * PixelLength_pc, ne; dims=3), dims=3)

    T = float(promote_type(eltype(ne), typeof(PixelLength_pc)))
    out = fill(T(NaN), size(ne, 1), size(ne, 2))
    @inbounds for j in axes(ne, 2), i in axes(ne, 1)
        acc = zero(T)
        found = false
        for k in axes(ne, 3)
            value = ne[i, j, k]
            if isfinite(value)
                acc += abs2(value) * PixelLength_pc
                found = true
            end
        end
        found && (out[i, j] = acc)
    end
    return out
end

"""
    WolfireConstants() -> Tuple{Float64, Float64, Float64, Float64}

Prompt the user to enter values for the Wolfire model constants and return them.

# Returns
- `Tuple{Float64, Float64, Float64, Float64}`: A tuple containing:
  - `zeta`: Ionization rate by Cosmic Rays.
  - `Geff`: Effective radiation field.
  - `phiPAH`: PAH grain alignment efficiency.
  - `XC`: Conversion factor of H into C.

# Description
This function prompts the user to enter values for the constants used in the Wolfire et al. model, specifically the ionization rate by cosmic rays (zeta), the effective radiation field (Geff), collision rate parameter for PAH (phiPAH), and the conversion factor of H into C (XC). Default values are provided for each constant, which the user can override.

# Example
```julia
# Example usage
zeta, Geff, omegaPAH, XC = WolfireConstants()
println("zeta: ", zeta)
println("Geff: ", Geff)
println("phiPAH: ", phiPAH)
println("XC: ", XC)
"""
function WolfireConstants()
    
    println("Please enter the values for the constants:")
    zeta = ask_user("zeta (ionization rate by Cosmic Rays)", 2.5e-16)
    Geff = ask_user("Geff (effective radiation field)", 1.0)
    phiPAH = ask_user("phiPAH (collision rate parameter for PAH)", 0.5)
    XC = ask_user("XC (Conversion factor of H into C)", 1.4e-4)
    
    return zeta, Geff, phiPAH, XC
end
