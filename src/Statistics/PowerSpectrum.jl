using Statistics

const _cairomakie_pkgid = Base.PkgId(Base.UUID("13f3f980-e62b-5c42-98c6-ff1f3baf88f0"), "CairoMakie")

function _ensure_cairomakie()
    pkg_path = Base.find_package(_cairomakie_pkgid.name)
    pkg_path === nothing && error("CairoMakie is required for plotting. Install it with `Pkg.add(\"CairoMakie\")` and ensure it is available in the active environment.")

    return Base.require(_cairomakie_pkgid)
end

"""
    power_spectrum_2d(field; pixel_size = 1.0, center = true, detrend_mean = true, normalize = true)

Compute the 2D power spectral density (PSD) of a 2D field.

# Arguments
- `field::AbstractMatrix`: Real- or complex-valued image/map.
- `pixel_size`: Physical size of a pixel. Frequencies are returned in cycles per `pixel_size`.
- `center`: When `true`, apply `fftshift` so the zero-frequency component is centered.
- `detrend_mean`: Subtract the global mean before transforming to remove the DC peak.
- `normalize`: Divide by the number of pixels so Parseval-like scaling is preserved.

# Returns
`(kx, ky, psd)`, where `kx` and `ky` are frequency axes and `psd` is the 2D PSD.
"""
function power_spectrum_2d(field::AbstractMatrix; pixel_size::Real = 1.0, center::Bool = true,
    detrend_mean::Bool = true, normalize::Bool = true, log_progress::Bool = false)

    isfinite(pixel_size) && pixel_size > 0 ||
        throw(ArgumentError("pixel_size must be positive and finite (got $pixel_size)."))
    isempty(field) && throw(ArgumentError("field must not be empty."))
    all(isfinite, field) || throw(ArgumentError(
        "power_spectrum_2d requires a fully finite rectangular field; masked NaN/Inf pixels cannot be passed to an FFT without an explicit inpainting or windowing strategy."))

    log_progress && @info "Computing 2D power spectrum" size = size(field) pixel_size = pixel_size detrend_mean = detrend_mean normalize = normalize

    data = detrend_mean ? field .- mean(field) : field
    nx, ny = size(data)

    ft = fft(data)
    psd = abs.(ft).^2
    normalize && (psd ./= nx * ny)

    kx = FFTW.fftfreq(nx, 1 / pixel_size)
    ky = FFTW.fftfreq(ny, 1 / pixel_size)

    if center
        result = (fftshift(kx), fftshift(ky), fftshift(psd))
    else
        result = (kx, ky, psd)
    end

    log_progress && @info "2D power spectrum ready" centered = center result_sizes = map(size, result)
    return result
end

"""
    radial_psd(field; pixel_size = 1.0, nbins = nothing, detrend_mean = true, normalize = true)

Compute the radially averaged 1D PSD of a 2D field.

# Arguments
- `field::AbstractMatrix`: Real- or complex-valued image/map.
- `pixel_size`: Physical size of a pixel. Frequencies are returned in cycles per `pixel_size`.
- `nbins`: Number of radial bins. Defaults to half the smallest field dimension.
- `detrend_mean`: Subtract the global mean before transforming to remove the DC peak.
- `normalize`: Divide by the number of pixels so Parseval-like scaling is preserved.

# Returns
`(k, psd_1d)`, where `k` are bin centers and `psd_1d` is the mean PSD per radial bin.
"""
function radial_psd(field::AbstractMatrix; pixel_size::Real = 1.0, nbins::Union{Int, Nothing} = nothing,
    detrend_mean::Bool = true, normalize::Bool = true, log_progress::Bool = false)

    kx, ky, psd2d = power_spectrum_2d(field; pixel_size = pixel_size, center = true,
        detrend_mean = detrend_mean, normalize = normalize, log_progress = log_progress)

    nx, ny = size(psd2d)

    nbins = isnothing(nbins) ? max(floor(Int, min(nx, ny) / 2), 1) : nbins
    max_radius = hypot(maximum(abs, kx), maximum(abs, ky))
    edges = collect(range(0, max_radius, length = nbins + 1))
    radius_type = real(promote_type(eltype(kx), eltype(ky)))
    bin_width = max(edges[2] - edges[1], eps(radius_type))
    inv_bin_width = inv(bin_width)

    bin_sums = zeros(Float64, nbins)
    bin_counts = zeros(Int, nbins)

    total_points = nx * ny
    progress_step = max(floor(Int, total_points / 10), 1)
    processed = 0
    @inbounds for j in eachindex(ky)
        kyj = ky[j]
        for i in eachindex(kx)
            radius = hypot(kx[i], kyj)
            bin = clamp(floor(Int, radius * inv_bin_width) + 1, 1, nbins)
            bin_sums[bin] += psd2d[i, j]
            bin_counts[bin] += 1
            processed += 1
            if log_progress && processed % progress_step == 0
                print_progress(processed, total_points; label="Power-spectrum radial binning")
                @debug "Radial PSD binning" processed = processed total = total_points
            end
        end
    end

    bin_centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    psd1d = bin_sums ./ max.(bin_counts, 1)

    log_progress && @info "Radial PSD complete" nbins = nbins

    return bin_centers, psd1d
end

"""
    plot_power_spectrum_figure(field; pixel_size = 1.0, nbins = nothing, detrend_mean = true,
        normalize = true, log2d = true, log1d = true, slope = nothing,
        slope_k_range = nothing, slope_color = :red, colormap = :viridis, fig_kwargs...)

Create a CairoMakie figure showing the 2D power spectrum and its radially averaged 1D PSD
side by side.

# Arguments
- `field::AbstractMatrix`: Real- or complex-valued image/map.
- `pixel_size`: Physical size of a pixel. Frequencies are returned in cycles per `pixel_size`.
- `nbins`: Number of radial bins for the 1D PSD. Defaults to half the smallest field dimension.
- `detrend_mean`: Subtract the global mean before transforming to remove the DC peak.
- `normalize`: Divide by the number of pixels so Parseval-like scaling is preserved.
- `log2d`: When `true`, apply a base-10 logarithm to the 2D PSD heatmap values.
- `log1d`: When `true`, set the radial PSD axes to log-scale and safeguard the values.
- `slope`: Optional power-law slope to overlay on the 1D PSD. When `nothing`, the
  slope is estimated via linear regression in log-log space using the available bins.
- `slope_k_range`: Optional `(kmin, kmax)` tuple limiting the frequencies used for slope
  estimation and the overlay line. When `nothing`, all positive bins are used.
- `slope_color`: Color to use for the slope overlay line.
- `colormap`: Colormap to use for the 2D PSD heatmap.
- `fig_kwargs`: Additional keyword arguments forwarded to `CairoMakie.Figure`.

# Returns
`CairoMakie.Figure` containing the 2D PSD heatmap (with colorbar) and the 1D PSD line plot.
`CairoMakie` must be available in the active environment.
"""
function plot_power_spectrum_figure(field::AbstractMatrix; pixel_size::Real = 1.0,
    nbins::Union{Int, Nothing} = nothing, detrend_mean::Bool = true,
    normalize::Bool = true, log2d::Bool = true, log1d::Bool = true,
    slope::Union{Nothing, Real} = nothing,
    slope_k_range::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
    slope_color = :red, colormap = :viridis, fig_kwargs...)

    CairoMakie = _ensure_cairomakie()

    kx, ky, psd2d = power_spectrum_2d(field; pixel_size = pixel_size, center = true,
        detrend_mean = detrend_mean, normalize = normalize)
    k, psd1d = radial_psd(field; pixel_size = pixel_size, nbins = nbins,
        detrend_mean = detrend_mean, normalize = normalize)

    eps_val = eps(real(eltype(psd2d)))
    heatmap_values = log2d ? log10.(psd2d .+ eps_val) : psd2d
    k_safe = log1d ? max.(k, eps_val) : k
    psd1d_safe = log1d ? max.(psd1d, eps_val) : psd1d

    fig = CairoMakie.Figure(; fig_kwargs...)
    left = fig[1, 1] = CairoMakie.GridLayout()

    ax2d = CairoMakie.Axis(left[1, 1]; title = "2D power spectrum",
        xlabel = "kₓ (cycles/pixel)", ylabel = "k_y (cycles/pixel)")
    hm = CairoMakie.heatmap!(ax2d, kx, ky, heatmap_values; colormap = colormap)
    CairoMakie.Colorbar(left[1, 2], hm; label = log2d ? "log₁₀ PSD" : "PSD")
    CairoMakie.tightlimits!(ax2d)

    ax1d = CairoMakie.Axis(fig[1, 2]; title = "Radially averaged PSD",
        xlabel = "|k| (cycles/pixel)", ylabel = "PSD")
    if log1d
        ax1d.xscale = CairoMakie.log10
        ax1d.yscale = CairoMakie.log10
    end
    CairoMakie.lines!(ax1d, k_safe, psd1d_safe; color = :black, label = "PSD")

    # Overlay an estimated or user-provided slope on the 1D PSD
    slope_mask = (k_safe .> eps_val) .& (psd1d_safe .> eps_val)
    if !isnothing(slope_k_range)
        kmin, kmax = slope_k_range
        slope_mask .&= (k_safe .>= kmin) .& (k_safe .<= kmax)
    end

    if any(slope_mask)
        xs = log10.(k_safe[slope_mask])
        ys = log10.(psd1d_safe[slope_mask])
        fitted_slope = isnothing(slope) ? (sum((xs .- mean(xs)) .* (ys .- mean(ys))) /
            max(sum((xs .- mean(xs)) .^ 2), eps_val)) : slope

        # Anchor the slope line to the PSD at a representative frequency
        k_line_min = isnothing(slope_k_range) ? minimum(k_safe[slope_mask]) : first(slope_k_range)
        k_line_max = isnothing(slope_k_range) ? maximum(k_safe[slope_mask]) : last(slope_k_range)
        k_line = range(k_line_min, k_line_max; length = 100)
        anchor_k = median(k_safe[slope_mask])
        anchor_idx = argmin(abs.(k_safe .- anchor_k))
        anchor_psd = psd1d_safe[anchor_idx]
        slope_line = anchor_psd .* (k_line ./ anchor_k) .^ fitted_slope

        CairoMakie.lines!(ax1d, k_line, slope_line; color = slope_color, linestyle = :dash,
            label = "pente = $(round(fitted_slope, digits = 2))")
        CairoMakie.axislegend(ax1d; position = :rb)
    end
    CairoMakie.tightlimits!(ax1d)

    return fig
end
