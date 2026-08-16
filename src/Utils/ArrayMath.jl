"""
Collection of array-level helpers used throughout the pipeline.
"""

maxCube(cube::AbstractArray) = dropdims(maximum(cube, dims = 3), dims = 3)

function _check_nan_policy(nan_policy::Symbol)
    nan_policy in (:propagate, :omit, :error) || throw(ArgumentError(
        "nan_policy must be :propagate, :omit, or :error (got :$nan_policy)."))
    return nan_policy
end

"""
    intLOS(cube, pixel_length; nan_policy=:propagate)

Integrate a 3D cube along its third dimension. `nan_policy=:propagate` keeps
masked sightlines masked, `:omit` sums only finite cells (all-invalid
sightlines return `NaN`), and `:error` rejects non-finite input.
"""
function intLOS(cube::AbstractArray{<:Real, 3}, pixel_length::Real;
                nan_policy::Symbol = :propagate)
    _check_nan_policy(nan_policy)
    isfinite(pixel_length) || throw(ArgumentError("pixel_length must be finite."))
    if nan_policy === :error
        all(isfinite, cube) || throw(ArgumentError("intLOS input contains non-finite values."))
    end
    nan_policy !== :omit &&
        return dropdims(sum(x -> x * pixel_length, cube; dims=3), dims=3)

    T = float(promote_type(eltype(cube), typeof(pixel_length)))
    out = zeros(T, size(cube, 1), size(cube, 2))
    seen = falses(size(cube, 1), size(cube, 2))
    # Slice-major traversal: the innermost index runs along the contiguous
    # dimension, whereas accumulating one sightline at a time strides through
    # the cube by nx*ny elements on every step. The cells of a given sightline
    # are still visited in increasing k, so the sums are unchanged.
    @inbounds for k in axes(cube, 3)
        for j in axes(cube, 2), i in axes(cube, 1)
            value = cube[i, j, k]
            if isfinite(value)
                out[i, j] += value * pixel_length
                seen[i, j] = true
            end
        end
    end
    @inbounds for idx in eachindex(out, seen)
        seen[idx] || (out[idx] = T(NaN))
    end
    return out
end

"""
    sigmaLOS(cube; nan_policy=:propagate)

Compute the sample standard deviation along the third dimension. With
`:omit`, only finite cells are used and fewer than two valid cells yield
`NaN`.
"""
function sigmaLOS(cube::AbstractArray{<:Real, 3}; nan_policy::Symbol = :propagate)
    _check_nan_policy(nan_policy)
    if nan_policy === :error
        all(isfinite, cube) || throw(ArgumentError("sigmaLOS input contains non-finite values."))
    end
    nan_policy !== :omit && return dropdims(std(cube, dims=3), dims=3)

    T = float(eltype(cube))
    nx, ny = size(cube, 1), size(cube, 2)
    counts = zeros(Int, nx, ny)
    means = zeros(T, nx, ny)
    m2 = zeros(T, nx, ny)
    # Slice-major traversal (see `intLOS`). Each sightline still runs its Welford
    # recurrence over increasing k, so the accumulated variance is unchanged.
    @inbounds for k in axes(cube, 3)
        for j in axes(cube, 2), i in axes(cube, 1)
            value = cube[i, j, k]
            if isfinite(value)
                count = counts[i, j] + 1
                counts[i, j] = count
                delta = value - means[i, j]
                means[i, j] += delta / count
                m2[i, j] += delta * (value - means[i, j])
            end
        end
    end

    out = fill(T(NaN), nx, ny)
    @inbounds for idx in eachindex(out, counts, m2)
        count = counts[idx]
        count >= 2 && (out[idx] = sqrt(max(m2[idx], zero(T)) / (count - 1)))
    end
    return out
end
