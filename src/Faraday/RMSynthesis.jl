"""
    RMSynthesis(Q::AbstractArray, U::AbstractArray, nuArray::AbstractArray, PhiArray::AbstractArray) -> Tuple{AbstractArray, AbstractArray, AbstractArray}

Perform Rotation Measure (RM) synthesis on Stokes Q and U parameters.

# Arguments
- `Q::AbstractArray`: Array representing the Stokes Q parameter. Can be 1D, 2D, or 3D.
- `U::AbstractArray`: Array representing the Stokes U parameter. Can be 1D, 2D, or 3D.
- `nuArray::AbstractArray`: Array of frequency values in Hz.
- `PhiArray::AbstractArray`: Array of Faraday depths in rad/m².

# Returns
- `Tuple{AbstractArray, AbstractArray, AbstractArray}`: A tuple containing:
  - `absF::AbstractArray`: Absolute value of the Faraday dispersion function.
  - `realF::AbstractArray`: Real part of the Faraday dispersion function.
  - `imagF::AbstractArray`: Imaginary part of the Faraday dispersion function.

# Description
The `RMSynthesis` function is a technique used in radio astronomy to study the Faraday rotation effect. It combines the Stokes Q and U parameters into the complex polarization P, and computes the Faraday dispersion function F for each value in `PhiArray`. The result is provided in absolute, real, and imaginary components.

# Example
```julia
# Example input arrays
Q = [0.1, 0.2, 0.3]
U = [0.4, 0.5, 0.6]
nuArray = [1e9, 1.1e9, 1.2e9]
PhiArray = [-100, 0, 100]

# Function call
absF, realF, imagF = RMSynthesis(Q, U, nuArray, PhiArray)

# Output
absF = [...]
realF = [...]
imagF = [...]
```
"""

const _RM_SYNTHESIS_BLOCK_SIZE = 64

function _validate_faraday_axes(nuArray::AbstractArray, PhiArray::AbstractArray;
                                min_frequencies::Int = 1,
                                require_frequency_span::Bool = false)
    length(nuArray) >= min_frequencies || throw(ArgumentError(
        "At least $min_frequencies frequency channel(s) are required."))
    isempty(PhiArray) && throw(ArgumentError("PhiArray must not be empty."))
    all(nu -> isfinite(nu) && nu > 0, nuArray) || throw(ArgumentError(
        "nuArray must contain only positive, finite frequencies."))
    all(isfinite, PhiArray) || throw(ArgumentError(
        "PhiArray must contain only finite Faraday depths."))
    if require_frequency_span
        extrema(nuArray)[1] < extrema(nuArray)[2] || throw(ArgumentError(
            "At least two distinct frequencies are required."))
    end
    return nothing
end

function _validate_rmsynthesis_inputs(Q::AbstractArray, U::AbstractArray,
                                      nuArray::AbstractArray, PhiArray::AbstractArray)
    size(Q) == size(U) || throw(ArgumentError(
        "Q and U must have the same size (got $(size(Q)) and $(size(U)))."))
    1 <= ndims(Q) <= 3 || throw(ArgumentError(
        "Q and U must be 1D, 2D, or 3D arrays."))
    _validate_faraday_axes(nuArray, PhiArray)
    size(Q, ndims(Q)) == length(nuArray) || throw(ArgumentError(
        "The last Q/U dimension ($(size(Q, ndims(Q)))) must match nuArray ($(length(nuArray)))."))
    return nothing
end

function _rmsynthesis_mul!(F::AbstractMatrix{CT}, P::AbstractMatrix{CT}, a,
                           PhiArray, scale::T; log_progress::Bool = false) where {T<:Real, CT<:Complex{T}}
    nlambda = size(P, 2)
    nphi = length(PhiArray)
    block_size = min(_RM_SYNTHESIS_BLOCK_SIZE, max(nphi, 1))
    phases = Matrix{CT}(undef, nlambda, block_size)

    # Process several Faraday depths at once. This turns the former sequence
    # of BLAS level-2 matrix-vector products into level-3 matrix-matrix
    # products while bounding temporary storage to nlambda * block_size.
    for first_phi in 1:block_size:nphi
        last_phi = min(first_phi + block_size - 1, nphi)
        width = last_phi - first_phi + 1
        @inbounds for column in 1:width
            phi = Float64(PhiArray[first_phi + column - 1])
            for l in eachindex(a)
                phases[l, column] = CT(cis(-2.0 * phi * a[l]))
            end
        end

        LinearAlgebra.mul!(view(F, :, first_phi:last_phi), P,
                           view(phases, :, 1:width), scale, zero(T))
        if log_progress
            print_progress(last_phi, nphi; label="RM synthesis")
            @debug "RM synthesis accumulation" idx = last_phi total = nphi
        end
    end

    return F
end

function RMSynthesis(Q::AbstractArray, U::AbstractArray, nuArray::AbstractArray, PhiArray::AbstractArray; log_progress::Bool = false)

    _validate_rmsynthesis_inputs(Q, U, nuArray, PhiArray)
    log_progress && @info "Starting RM synthesis" n_phi = length(PhiArray) n_lambda = length(nuArray)

    LambdaSqArray = @. (C_m/nuArray)^2

    nPhi = length(PhiArray)
    nLambda = length(LambdaSqArray)
    nDims = length(size(Q))

    # The FDF is accumulated in the working precision of the Q/U cubes
    # (Float64 inputs reproduce the historical behaviour bit for bit;
    # Float32 inputs halve the memory footprint of the FDF products).
    T = float(real(promote_type(eltype(Q), eltype(U))))
    CT = Complex{T}

    K = 1.0 / nLambda

    if nDims == 1
        Q = reshape(Q, (1,1,size(Q,1)))
        U = reshape(U, (1,1,size(U,1)))
    elseif nDims == 2
        Q = reshape(Q, (1,size(Q,1), size(Q,2)))
        U = reshape(U, (1,size(U,1), size(U,2)))
    end

    P = @. (Q + 1im * U)

    nx, ny = size(Q,1), size(Q,2)

    Lambda0Sq = sum(LambdaSqArray) * K
    a = (LambdaSqArray .- Lambda0Sq)

    Pmat = reshape(P, nx * ny, nLambda)
    # Masked pixels (NaN from HEALPix UNSEEN or partial-sky maps) are skipped:
    # their FDF is NaN by definition, so there is no point running the
    # matrix product over them.
    npix_total = nx * ny
    valid = Vector{Bool}(undef, npix_total)
    @inbounds for r in 1:npix_total
        ok = true
        for l in 1:nLambda
            z = Pmat[r, l]
            if !(isfinite(real(z)) && isfinite(imag(z)))
                ok = false
                break
            end
        end
        valid[r] = ok
    end
    nvalid = count(valid)

    if nvalid == npix_total
        F = Matrix{CT}(undef, npix_total, nPhi)
        _rmsynthesis_mul!(F, Pmat, a, PhiArray, T(K); log_progress)
    else
        log_progress && @info "Skipping masked pixels in RM synthesis" masked = npix_total - nvalid total = npix_total
        rows = findall(valid)
        Pvalid = Pmat[rows, :]
        Fvalid = Matrix{CT}(undef, nvalid, nPhi)
        _rmsynthesis_mul!(Fvalid, Pvalid, a, PhiArray, T(K); log_progress)
        F = fill(CT(T(NaN), T(NaN)), npix_total, nPhi)
        F[rows, :] = Fvalid
    end

    F = reshape(F, nx, ny, nPhi)
    
    if nDims == 1
        F = dropdims(dropdims(F,dims=1),dims=1)
    elseif nDims == 2
        F = dropdims(F,dims=1)
    end
     
    log_progress && @info "RM synthesis complete" output_size = size(F)

    return(abs.(F),real.(F),imag.(F))
end

"""
    getRMSF(nuArray::AbstractArray, PhiArray::AbstractArray) -> Tuple{AbstractArray, Float64}

Calculate the Rotation Measure Spread Function (RMSF) and its full width at half maximum (FWHM).

# Arguments
- `nuArray::AbstractArray`: Array of frequency values in Hz.
- `PhiArray::AbstractArray`: Array of Faraday depths in rad/m².

# Returns
- `Tuple{AbstractArray, Float64}`: A tuple containing:
  - `absRMSF::AbstractArray`: Absolute value of the Rotation Measure Spread Function.
  - `fwhmRMSF::Float64`: Full width at half maximum (FWHM) of the RMSF.

# Description
The `getRMSF` function calculates the Rotation Measure Spread Function (RMSF), which describes the response of an RM synthesis to a single Faraday depth component. It uses the input frequency values (`nuArray`) and Faraday depths (`PhiArray`) to compute the RMSF and its FWHM.

# Example
```julia
# Example input arrays
nuArray = [1e9, 1.1e9, 1.2e9]
PhiArray = [-100, 0, 100]

# Function call
absRMSF, fwhmRMSF = getRMSF(nuArray, PhiArray)

# Output
absRMSF = [...]
fwhmRMSF = ...
```
"""
function getRMSF(nuArray::AbstractArray, PhiArray::AbstractArray; log_progress::Bool = false)

    _validate_faraday_axes(nuArray, PhiArray; min_frequencies=2, require_frequency_span=true)
    log_progress && @info "Starting RMSF computation" n_phi = length(PhiArray)
    
    LambdaSqArray = @. (C_m/nuArray)^2
    
    nPhi = length(PhiArray)
    nLambda = length(LambdaSqArray)

    K = 1.0 / nLambda

    Lambda0Sq = sum(LambdaSqArray) * K
    a = (LambdaSqArray .- Lambda0Sq)
    
    fwhmRMSF = 3.8 / (maximum(LambdaSqArray) - minimum(LambdaSqArray))

    RMSF = Vector{ComplexF64}(undef, nPhi)
    phase = Vector{ComplexF64}(undef, nLambda)
    for i in 1:nPhi
        phi = Float64(PhiArray[i])
        @inbounds for l in eachindex(a)
            phase[l] = cis(-2.0 * phi * a[l])
        end
        RMSF[i] = K * sum(phase)
        if log_progress
            print_progress(i, nPhi; label="RMSF computation")
            @debug "RMSF accumulation" idx = i total = nPhi
        end
    end

    log_progress && @info "RMSF computation complete" output_size = length(RMSF) fwhm = fwhmRMSF

    return(abs.(RMSF),fwhmRMSF)
end

"""
    RMSFDiagnostics

Container for the Rotation Measure Spread Function (RMSF) and the resolution
metrics that characterise an RM-synthesis experiment.

# Fields
- `phi::Vector{Float64}`: Faraday-depth lags at which the RMSF is sampled
  (rad/m²). The grid is symmetric about `0` and shares the spacing of the
  Faraday-depth array passed to [`rmsf_diagnostics`](@ref).
- `rmsf::Vector{ComplexF64}`: complex RMSF `R(φ)`, normalised so that
  `R(0) ≈ 1`.
- `fwhm::Float64`: full width at half maximum of `|R(φ)|`, measured directly
  from the sampled main lobe (rad/m²). This is the effective Faraday
  resolution and the natural restoring-beam width for RM-CLEAN.
- `fwhm_theoretical::Float64`: analytic resolution `2√3 / Δλ²` (rad/m²).
- `phi_max::Float64`: largest recoverable `|φ|` set by the channel width in
  `λ²` (rad/m²).
- `max_scale::Float64`: largest Faraday-thick structure that stays sensitive,
  `π / λ²_min` (rad/m²).
- `lambda2_min`, `lambda2_max`, `lambda0_2::Float64`: minimum, maximum and
  weighted-mean `λ²` of the observation (m²).
"""
struct RMSFDiagnostics
    phi::Vector{Float64}
    rmsf::Vector{ComplexF64}
    fwhm::Float64
    fwhm_theoretical::Float64
    phi_max::Float64
    max_scale::Float64
    lambda2_min::Float64
    lambda2_max::Float64
    lambda0_2::Float64
end

# Uniform spacing of a Faraday-depth sampling array.
function _uniform_step(values::AbstractArray; name::AbstractString = "sample grid")
    n = length(values)
    n >= 2 || error("At least two samples are required to infer a step size.")
    all(isfinite, values) || throw(ArgumentError("$(name) must contain only finite values."))

    step = Float64(values[2]) - Float64(values[1])
    !iszero(step) || throw(ArgumentError("$(name) must have a non-zero step."))
    tol = max(1e-10, 1e-8 * max(abs(step), maximum(abs, Float64.(values))))
    for i in 3:n
        local_step = Float64(values[i]) - Float64(values[i - 1])
        isapprox(local_step, step; atol = tol, rtol = 1e-8) ||
            error("$(name) must be uniformly spaced; step $(i - 1) is $(local_step), expected $(step).")
    end

    return step
end

# FWHM of a sampled, single-peaked curve via linear interpolation of the
# half-maximum crossings on either side of the global maximum.
function _measure_fwhm(x::AbstractVector, y::AbstractVector)
    peak, ipk = findmax(y)
    (isfinite(peak) && peak > 0) || return NaN
    half = peak / 2

    cross(i, j) = begin
        yi, yj = y[i], y[j]
        yi == yj ? x[i] : x[i] + (half - yi) * (x[j] - x[i]) / (yj - yi)
    end

    l = ipk
    while l > 1 && y[l] > half
        l -= 1
    end
    left = l == ipk ? x[ipk] : (y[l] > half ? x[l] : cross(l, l + 1))

    r = ipk
    while r < length(y) && y[r] > half
        r += 1
    end
    right = r == ipk ? x[ipk] : (y[r] > half ? x[r] : cross(r, r - 1))

    return right - left
end

"""
    rmsf_diagnostics(nuArray, PhiArray) -> RMSFDiagnostics

Compute the complex Rotation Measure Spread Function and the standard
RM-synthesis resolution metrics for a given frequency coverage.

# Arguments
- `nuArray::AbstractArray`: observed frequencies in **Hz**.
- `PhiArray::AbstractArray`: the (uniformly spaced) Faraday-depth array used for
  RM synthesis, in rad/m². Its spacing sets the RMSF sampling.

# Returns
- [`RMSFDiagnostics`](@ref) with the sampled complex RMSF and the resolution
  metrics (`fwhm`, `fwhm_theoretical`, `phi_max`, `max_scale`).

The RMSF is evaluated on a symmetric lag grid spanning `±(N-1)·δφ`, where `N`
is `length(PhiArray)` and `δφ` its spacing, so the same object can act as the
convolution kernel for [`RMClean`](@ref).

# Example
```julia
nu = collect(range(1.0e9, 1.5e9, length = 64))
phi = collect(range(-100.0, 100.0, length = 201))
diag = rmsf_diagnostics(nu, phi)
diag.fwhm          # measured Faraday resolution (rad/m²)
diag.phi_max       # maximum recoverable |φ| (rad/m²)
```
"""
function rmsf_diagnostics(nuArray::AbstractArray, PhiArray::AbstractArray)
    length(PhiArray) >= 2 || error("rmsf_diagnostics requires at least two Faraday depths.")
    _validate_faraday_axes(nuArray, PhiArray; min_frequencies=2, require_frequency_span=true)

    LambdaSqArray = @. (C_m / nuArray)^2
    nLambda = length(LambdaSqArray)
    K = 1.0 / nLambda
    Lambda0Sq = sum(LambdaSqArray) * K
    a = LambdaSqArray .- Lambda0Sq

    dphi = _uniform_step(PhiArray; name = "PhiArray")
    n = length(PhiArray) - 1
    phigrid = collect(-n:n) .* dphi

    rmsf = Vector{ComplexF64}(undef, length(phigrid))
    @inbounds for (idx, phi) in enumerate(phigrid)
        acc = zero(ComplexF64)
        for l in eachindex(a)
            acc += cis(-2.0 * phi * a[l])
        end
        rmsf[idx] = K * acc
    end

    fwhm_measured = _measure_fwhm(phigrid, abs.(rmsf))
    delta_phi, max_scale, phi_max = rmsynthesis_parameters(nuArray)
    fwhm = (isfinite(fwhm_measured) && fwhm_measured > 0) ? fwhm_measured : delta_phi

    return RMSFDiagnostics(
        phigrid,
        rmsf,
        fwhm,
        delta_phi,
        phi_max,
        max_scale,
        minimum(LambdaSqArray),
        maximum(LambdaSqArray),
        Lambda0Sq,
    )
end
