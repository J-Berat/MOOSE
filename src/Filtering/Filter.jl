"""
    apply_instrument_2d(img, H)

Apply a 2D Fourier-domain instrumental mask. `H` is a binary transfer
function with ones for sampled spatial frequencies and zeros for rejected
frequencies, as in the instrumental filtering used by the Depolarization
pipeline.
"""
function apply_instrument_2d(img::AbstractMatrix, H::AbstractMatrix)
    size(img) == size(H) || error("Filter shape mismatch: image size=$(size(img)) filter size=$(size(H))")
    all(isfinite, img) || throw(ArgumentError(
        "apply_instrument_2d requires finite image values; fill, crop, or inpaint masked pixels before the FFT."))
    all(isfinite, H) || throw(ArgumentError("The Fourier filter H contains non-finite values."))
    return real.(ifft(fft(img) .* H))
end

"""
    apply_to_array_xy(data, H; n=size(H, 1), m=size(H, 2))

Apply a Fourier-domain instrumental mask to a 2D image or to every sky-plane
slice of a Stokes cube. Supported cube layouts are `(n, m, nν)` and
`(nν, n, m)`.
"""
function apply_to_array_xy(data, H; n::Int=size(H, 1), m::Int=size(H, 2))
    size(H) == (n, m) || error("Filter mask H must have size ($n,$m), got $(size(H))")
    nd = ndims(data)

    if nd == 2
        size(data) == (n, m) || error("2D input must have size ($n,$m), got $(size(data))")
        return apply_instrument_2d(data, H)
    elseif nd == 3
        sz = size(data)
        Tout = float(eltype(data))
        out = similar(data, Tout, sz)

        if sz[1] == n && sz[2] == m
            @views for k in axes(data, 3)
                out[:, :, k] = apply_instrument_2d(data[:, :, k], H)
            end
            return out
        elseif sz[2] == n && sz[3] == m
            @views for k in axes(data, 1)
                out[k, :, :] = apply_instrument_2d(data[k, :, :], H)
            end
            return out
        else
            error("Unsupported 3D shape $(sz). Expected (n,m,nν) or (nν,n,m) with n=$n m=$m.")
        end
    else
        error("Unsupported ndims(data)=$nd")
    end
end

"""
    instrument_bandpass_L(n, m; Δx, Δy=Δx, Lcut_small, Llarge, fNy)

Build a hard 0/1 spatial-frequency band-pass mask. The mask removes scales
larger than `Llarge` and smaller than `Lcut_small`, capped at the Nyquist
frequency `fNy`.
"""
function instrument_bandpass_L(n::Int, m::Int;
                               Δx::Real, Δy::Real=Δx,
                               Lcut_small::Real,
                               Llarge::Real,
                               fNy::Real)
    n > 0 || error("n must be positive, got $n")
    m > 0 || error("m must be positive, got $m")
    isfinite(Δx) && Δx > 0 || error("Δx must be positive and finite, got $Δx")
    isfinite(Δy) && Δy > 0 || error("Δy must be positive and finite, got $Δy")
    isfinite(Lcut_small) && Lcut_small > 0 || error("Lcut_small must be positive and finite, got $Lcut_small")
    isfinite(Llarge) && Llarge > 0 || error("Llarge must be positive and finite, got $Llarge")
    isfinite(fNy) && fNy > 0 || error("fNy must be positive and finite, got $fNy")

    # fftfreq(n, fs) expects the *sampling frequency* fs = 1/Δx, not the step Δx.
    # Spatial frequencies are then in cycles per unit of Δx (same unit as
    # Lcut_small/Llarge, which must be expressed in that same length unit).
    fx = FFTW.fftfreq(n, 1 / Δx)
    fy = FFTW.fftfreq(m, 1 / Δy)

    flo = 1 / Llarge
    fhi_raw = 1 / Lcut_small
    fhi = min(fhi_raw, fNy)

    @debug "Band-pass filter frequencies" Lcut_small=Lcut_small Llarge=Llarge flo=flo fhi=fhi fhi_raw=fhi_raw

    flo2 = flo^2
    fhi2 = fhi^2
    H = Matrix{Float32}(undef, n, m)
    @inbounds for j in 1:m
        fy2 = float(fy[j])^2
        for i in 1:n
            f2 = float(fx[i])^2 + fy2
            H[i, j] = (f2 >= flo2 && f2 <= fhi2) ? 1f0 : 0f0
        end
    end

    return H, fftshift(H)
end
