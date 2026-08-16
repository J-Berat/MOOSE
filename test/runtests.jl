using Test
using Moose
using FITSIO
using HDF5
using JSON
using CSV
include(joinpath(@__DIR__, "..", "src", "MOOSE_cli.jl"))
using Moose.MooseFromConfig: build_config

function write_test_fits(path, data)
    FITS(path, "w") do fits
        write(fits, data)
    end
end

function write_test_hdf5(path, data; dataset=splitext(basename(path))[1])
    h5open(path, "w") do h5
        h5[dataset] = data
    end
end

function write_test_healpix_cube(path, pixels; nside=1, order=:ring)
    # Store one row per HEALPix pixel and one vector element per LOS shell.
    data = permutedims(Matrix(pixels))
    ordering = order == :nested ? "NESTED" : "RING"
    header = FITSHeader(
        ["PIXTYPE", "ORDERING", "NSIDE", "FIRSTPIX", "LASTPIX", "INDXSCHM"],
        ["HEALPIX", ordering, Int(nside), 1, size(pixels, 1), "IMPLICIT"],
        ["", "", "", "", "", ""],
    )
    FITS(path, "w") do fits
        write(fits, ["PIXELS"], [data]; header=header, name="MAP")
    end
end

function write_test_emissivity(path)
    open(path, "w") do io
        write(io, "B\tnu\te_para\te_perp\n")
        for nu in (99.0, 100.0, 101.0, 102.0)
            for B in (0.0, 1.0, 2.0, 3.0)
                write(io, "$(B)\t$(nu)\t$(1.0e-40 * (1 + B))\t$(2.0e-40 * (1 + B + nu / 1000))\n")
            end
        end
    end
end

@testset "RMS" begin
    x = [1.0, 2.0, 3.0]
    @test isapprox(Moose.RMS(x), sqrt(2 / 3))

    x2 = [1.0, 2.0]
    y2 = [3.0, 4.0]
    @test isapprox(Moose.RMS(x2, y2), sqrt(0.5))

    z2 = [5.0, 6.0]
    @test isapprox(Moose.RMS(x2, y2, z2), sqrt(0.75))

    @test Moose.RMS([1.0, NaN, 3.0]) == 1.0
    @test Moose.RMS([1.0, NaN, 3.0], [2.0, 50.0, 4.0]) == sqrt(2.0)
    @test isnan(Moose.RMS([NaN, Inf]))
end

@testset "Pnu" begin
    q = [1.0, 2.0]
    u = [3.0, 4.0]
    expected = sqrt.(q .^ 2 .+ u .^ 2)
    @test Moose.Pnu(q, u) == expected
end

@testset "Polarization fraction" begin
    p = [2.0, 1.0, NaN, 4.0]
    i = [4.0, 0.0, 3.0, -1.0]
    frac = Moose.PolarizationFraction(p, i)
    @test frac[1] == 0.5
    @test all(isnan, frac[2:end])
end

@testset "Conversion Jy/beam ↔ K" begin
    intensity = [0.0, 1.5]
    nu = 1.0e9
    theta = 10.0
    expected = 1.222e3 .* intensity ./ (nu^2 * theta^2)
    @test all(isapprox.(Moose.ConversionJyK(intensity, nothing, nu, theta), expected; rtol = 1e-12))
    @test iszero(Moose.ConversionJyK(0.0, nothing, nu, theta))
end

@testset "Rotation Measure" begin
    BLOS = [1.0, -2.0, 0.5]
    ne = [0.5, 1.0, 0.0]
    pixel_length = 2.0
    expected_delta = 0.81 .* ne .* BLOS .* pixel_length
    @test Moose.deltaRM(BLOS, ne, pixel_length) ≈ expected_delta atol = 0 rtol = 1e-12

    @test Moose.RM([1.0, 2.0, 3.0]) == [1.0, 3.0, 6.0]

    cube = reshape(1.0:8.0, 2, 2, 2)
    rm_cube = Moose.RM(cube)
    @test rm_cube[:, :, 1] == cube[:, :, 1]
    @test rm_cube[:, :, 2] == cube[:, :, 1] .+ cube[:, :, 2]
end

@testset "Line-of-sight reductions" begin
    cube = reshape(Float32.(1:12), 2, 2, 3)
    pixel_length = 0.5

    @test Moose.constant_ne(1.25, (2, 3, 4)) == fill(1.25, 2, 3, 4)
    @test Moose.DM(@view(cube[1, 1, :]), pixel_length) ≈ sum(cube[1, 1, :] .* pixel_length)
    @test Moose.EM(@view(cube[1, 1, :]), pixel_length) ≈ sum(cube[1, 1, :] .^ 2 .* pixel_length)
    @test Moose.DM(cube, pixel_length) ≈ dropdims(sum(cube .* pixel_length, dims = 3), dims = 3)
    @test Moose.EM(cube, pixel_length) ≈ dropdims(sum(cube .^ 2 .* pixel_length, dims = 3), dims = 3)
    @test Moose.intLOS(cube, pixel_length) ≈ dropdims(sum(cube .* pixel_length, dims = 3), dims = 3)

    masked = copy(cube)
    masked[1, 1, 2] = NaN
    masked[2, 2, :] .= NaN
    @test isnan(Moose.intLOS(masked, pixel_length)[1, 1])
    omitted_int = Moose.intLOS(masked, pixel_length; nan_policy=:omit)
    @test omitted_int[1, 1] ≈ (masked[1, 1, 1] + masked[1, 1, 3]) * pixel_length
    @test isnan(omitted_int[2, 2])
    @test_throws ArgumentError Moose.intLOS(masked, pixel_length; nan_policy=:error)

    omitted_sigma = Moose.sigmaLOS(masked; nan_policy=:omit)
    @test omitted_sigma[1, 1] ≈ Moose.std([masked[1, 1, 1], masked[1, 1, 3]])
    @test isnan(omitted_sigma[2, 2])
    @test_throws ArgumentError Moose.sigmaLOS(masked; nan_policy=:invalid)

    profile = [1.0, NaN, 3.0]
    @test isnan(Moose.DM(profile, pixel_length))
    @test Moose.DM(profile, pixel_length; nan_policy=:omit) == 2.0
    @test Moose.EM(profile, pixel_length; nan_policy=:omit) == 5.0
    @test_throws ArgumentError Moose.EM(profile, pixel_length; nan_policy=:error)
end

@testset "Physical cell mask" begin
    B1 = reshape(Float64.(1:8), 2, 2, 2)
    B2 = copy(B1) .+ 10
    BLOS = copy(B1) .+ 20
    T = reshape([50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0], 2, 2, 2)
    n = reshape([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], 2, 2, 2)
    nHp = copy(n) .+ 1.0

    mask_cfg = Dict("T_min" => 150.0, "T_max" => 350.0, "n_min" => 0.3, "n_max" => 0.7)
    expected = (T .>= 150.0) .& (T .<= 350.0) .& (n .>= 0.3) .& (n .<= 0.7)
    mask = Moose._apply_physical_mask!(B1, B2, BLOS, T, n, nHp, mask_cfg)

    @test mask == expected
    @test all(iszero, B1[.!expected])
    @test all(iszero, B2[.!expected])
    @test all(iszero, BLOS[.!expected])
    @test all(iszero, T[.!expected])
    @test all(iszero, n[.!expected])
    @test all(iszero, nHp[.!expected])
    @test all(!iszero, B1[expected])
end

@testset "Density kind conversion" begin
    mu = 1.4
    mH = Moose.M_p
    nH_expected = reshape([1.0, 2.0, 3.0, 4.0], 2, 2, 1)
    rho = nH_expected .* (mu * mH)

    out = Moose._density_to_number_density!(copy(rho), "mass_density", mu, mH)
    @test out ≈ nH_expected

    same = Moose._density_to_number_density!(copy(nH_expected), "number_density", mu, mH)
    @test same == nH_expected
    @test_throws ErrorException Moose._density_to_number_density!(copy(rho), "rho", mu, mH)
end

@testset "Moments" begin
    y = [1.0, 2.0, 3.0]
    x = [1.0, 2.0, 3.0]
    m0, m1, m2 = Moose.moments(y, x = x)
    @test isapprox(m0, 6.0; atol = 1e-12)
    @test isapprox(m1, 14 / 6; atol = 1e-12)
    expected_m2 = sqrt(sum(y .* (x .- m1) .^ 2) / m0)
    @test isapprox(m2, expected_m2; atol = 1e-12)

    empty_moments = Moose.moments([0.0, 0.0]; x = [1.0, 2.0], threshold = 0.5)
    @test all(isnan, empty_moments)

    finite_moments = Moose.moments([1.0, NaN, 3.0]; x = [1.0, 2.0, 3.0])
    @test finite_moments[1] == 4.0
    @test finite_moments[2] == 2.5
    @test isfinite(finite_moments[3])
    @test_throws ArgumentError Moose.moments([1.0, 2.0]; x = [1.0])
end

@testset "Spectral index map" begin
    nu = [100.0, 120.0, 150.0, 200.0, 300.0]  # MHz
    nx, ny = 3, 2
    beta_true = [-2.7 -0.5; -3.1 0.0; -2.0 -1.5]
    amp = [1.0e-3 2.0; 5.0 0.7; 3.0e2 1.0]
    cube = [amp[i, j] * nu[k]^beta_true[i, j] for i in 1:nx, j in 1:ny, k in eachindex(nu)]

    alpha, alpha_err = spectral_index_map(cube, nu)
    @test size(alpha) == (nx, ny)
    @test size(alpha_err) == (nx, ny)
    @test all(isapprox.(alpha, beta_true; atol = 1e-8))
    # Exact power laws leave only round-off residuals; the slope error floor
    # in double precision is ~1e-8 for |slope| of a few.
    @test all(err -> isfinite(err) && err < 1e-6, alpha_err)

    # The log-log slope is invariant under a rescaling of the frequency axis.
    alpha_hz, _ = spectral_index_map(cube, nu .* 1e6)
    @test all(isapprox.(alpha_hz, alpha; atol = 1e-8))

    # Non-positive / non-finite channels are excluded pixel by pixel.
    dirty = copy(cube)
    dirty[1, 1, 2] = -1.0
    dirty[1, 1, 4] = NaN
    alpha_dirty, err_dirty = spectral_index_map(dirty, nu)
    @test isapprox(alpha_dirty[1, 1], beta_true[1, 1]; atol = 1e-8)
    @test isfinite(err_dirty[1, 1])

    # Pixels with fewer than min_channels valid samples become NaN.
    dirty[2, 2, :] .= 0.0
    dirty[2, 2, 1] = 1.0
    dirty[2, 2, 5] = 1.0
    alpha_dirty, err_dirty = spectral_index_map(dirty, nu; min_channels = 3)
    @test isnan(alpha_dirty[2, 2])
    @test isnan(err_dirty[2, 2])

    # Exactly two points fit a slope but leave no residual dof for the error.
    two_pt, two_err = spectral_index_map(dirty, nu; min_channels = 2)
    @test isapprox(two_pt[2, 2], 0.0; atol = 1e-8)
    @test isnan(two_err[2, 2])

    # The slope error tracks injected log-space scatter.
    rng = Moose.Random.MersenneTwister(11)
    noisy = [10.0 * nuk^-2.7 * 10.0^(0.01 * randn(rng)) for _ in 1:1, _ in 1:1, nuk in nu]
    alpha_noisy, err_noisy = spectral_index_map(noisy, nu)
    @test isapprox(alpha_noisy[1, 1], -2.7; atol = 0.5)
    @test 0.0 < err_noisy[1, 1] < 0.5

    # Input validation.
    @test_throws Exception spectral_index_map(cube, nu[1:3])
    @test_throws Exception spectral_index_map(cube, replace(nu, 100.0 => -100.0))
    @test_throws Exception spectral_index_map(cube, nu; min_channels = 1)
    @test_throws Exception spectral_index_map(cube[:, :, 1:2], nu[1:2]; min_channels = 3)

    # HEALPix-style stacks (Npix x 1 x Nchan) are supported as-is.
    hp_cube = reshape(cube[:, 1, :], nx, 1, length(nu))
    hp_alpha, _ = spectral_index_map(hp_cube, nu)
    @test size(hp_alpha) == (nx, 1)
    @test all(isapprox.(hp_alpha[:, 1], beta_true[:, 1]; atol = 1e-8))
end

@testset "Power spectrum" begin
    delta_field = [1.0 0.0; 0.0 0.0]

    kx, ky, psd2d = Moose.power_spectrum_2d(delta_field; detrend_mean = false, normalize = true)
    @test kx == [-0.5, 0.0]
    @test ky == [-0.5, 0.0]
    @test all(isapprox.(psd2d, 0.25; atol = 1e-12))

    k, psd1d = Moose.radial_psd(delta_field; detrend_mean = false, normalize = true, nbins = 1)
    @test length(k) == 1
    @test isapprox(psd1d[1], 0.25; atol = 1e-12)
    @test_throws ArgumentError Moose.power_spectrum_2d([1.0 NaN; 2.0 3.0])
end

@testset "Interferometric filtering" begin
    H, Hshift = Moose.instrument_bandpass_L(8, 8; Δx = 1.0, Δy = 1.0, Lcut_small = 1.0, Llarge = 4.0, fNy = 0.5)
    @test size(H) == (8, 8)
    @test size(Hshift) == (8, 8)
    @test Set(vec(H)) ⊆ Set(Float32[0, 1])
    @test H[1, 1] == 0f0
    @test any(==(1f0), H)

    img = reshape(Float64.(1:64), 8, 8)
    out = Moose.apply_instrument_2d(img, H)
    @test size(out) == size(img)
    @test eltype(out) <: Real

    cube = repeat(img, 1, 1, 3)
    filtered = Moose.apply_to_array_xy(cube, H; n = 8, m = 8)
    @test size(filtered) == size(cube)
    @test filtered[:, :, 1] ≈ out
    @test_throws ArgumentError Moose.apply_instrument_2d([1.0 NaN; 2.0 3.0], ones(2, 2))
    @test_throws ArgumentError Moose.apply_instrument_2d(ones(2, 2), [1.0 Inf; 1.0 1.0])
    @test_throws ErrorException Moose.instrument_bandpass_L(2, 2; Δx=NaN,
        Lcut_small=1.0, Llarge=2.0, fNy=0.5)
end

@testset "Regression — LOS basis is cyclic for all three LOS (BUG-1)" begin
    @test Moose.los_basis(:Ax, :Ay, :Az, "z") == (:Ax, :Ay, :Az)
    @test Moose.los_basis(:Ax, :Ay, :Az, "x") == (:Ay, :Az, :Ax)
    @test Moose.los_basis(:Ax, :Ay, :Az, "y") == (:Az, :Ax, :Ay)
    @test_throws ErrorException Moose.los_basis(1, 2, 3, "w")
end

@testset "Regression — noise σ derives from SNR (BUG-2)" begin
    rng = Moose.Random.MersenneTwister(1234)
    Q = fill(3.0, 32, 32, 2)
    U = fill(4.0, 32, 32, 2)
    Q0, U0 = copy(Q), copy(U)
    snr = 10.0
    Moose._add_noise!(Q, U, snr, rng)
    # P_rms = sqrt(3² + 4²) = 5 per channel → expected σ = 5 / 10 = 0.5
    sigmaQ = sqrt(sum(abs2, Q .- Q0) / length(Q))
    sigmaU = sqrt(sum(abs2, U .- U0) / length(U))
    @test isapprox(sigmaQ, 0.5; rtol = 0.15)
    @test isapprox(sigmaU, 0.5; rtol = 0.15)
    @test_throws ErrorException Moose._add_noise!(Q, U, -1.0, rng)

    # Masked pixels neither poison the channel-wide noise estimate nor get
    # spuriously unmasked by the random draw.
    Qmasked = fill(3.0, 32, 32, 1)
    Umasked = fill(4.0, 32, 32, 1)
    Qmasked[1, 1, 1] = NaN
    Umasked[1, 1, 1] = NaN
    Moose._add_noise!(Qmasked, Umasked, snr, Moose.Random.MersenneTwister(1234))
    @test isnan(Qmasked[1, 1, 1]) && isnan(Umasked[1, 1, 1])
    @test all(isfinite, @view Qmasked[2:end, :, :])
    @test all(isfinite, @view Umasked[2:end, :, :])
end

@testset "NaN-aware descriptive statistics" begin
    data = [1.0 NaN; 4.0 -Inf]
    stats = Moose.CalculateStatistics(data)
    @test stats[1] == 4.0
    @test stats[2] == CartesianIndex(2, 1)
    @test stats[3] == 1.0
    @test stats[4] == CartesianIndex(1, 1)
    @test stats[5] == 2.5
    @test_throws ArgumentError Moose.CalculateStatistics(fill(NaN, 2, 2))

    @test Moose.EffectiveWidth([1.0, NaN, 2.0], [0.0, 1.0, 2.0]) == 1.5
    @test isnan(Moose.EffectiveWidth(fill(NaN, 3), 1:3))

    summary_input = reshape([1.0, 2.0, 3.0, 10.0], 2, 2)
    summary = Moose.SummarizeStats("z", summary_input, summary_input,
        summary_input, summary_input, summary_input,
        summary_input, summary_input, summary_input)
    expected_stats = Moose.CalculateStatistics(summary_input)
    @test summary.Quantity == ["Max", "IndMax", "Min", "IndMin", "Mean", "Median", "Std", "Skew", "Kurt"]
    @test summary.n[6] == expected_stats[6]
    @test summary.n[7] == expected_stats[7]
    @test summary.n[9] == expected_stats[9]
end

@testset "Regression — band-pass filter selects requested scales (BUG-3)" begin
    n = 64
    # Pure sine of wavelength 8 pixels along axis 1 (exact FFT bin).
    img = [sin(2π * i / 8) for i in 1:n, j in 1:n]
    Hkeep, _ = Moose.instrument_bandpass_L(n, n; Δx = 1.0, Δy = 1.0,
                                           Lcut_small = 2.0, Llarge = 16.0, fNy = 0.5)
    Hkill, _ = Moose.instrument_bandpass_L(n, n; Δx = 1.0, Δy = 1.0,
                                           Lcut_small = 2.0, Llarge = 4.0, fNy = 0.5)
    out_keep = Moose.apply_instrument_2d(img, Hkeep)
    out_kill = Moose.apply_instrument_2d(img, Hkill)
    @test out_keep ≈ img atol = 1e-8       # 2 ≤ 8 ≤ 16 → survives
    @test maximum(abs, out_kill) < 1e-8    # 8 > 4 → annihilated
end

@testset "Regression — emissivity grid handles rectangular/shuffled tables (BUG-4)" begin
    Bvals = [0.0, 1.0, 2.0, 3.0, 4.0]
    nuvals = [10.0, 20.0, 30.0, 40.0]
    rows = [(b, ν) for ν in nuvals for b in Bvals]
    reverse!(rows)  # scramble row order: result must not depend on it
    df = Moose.DataFrame(
        B = [r[1] for r in rows],
        nu = [r[2] for r in rows],
        e_perp = [r[1] + 1000 * r[2] for r in rows],
        e_para = zeros(length(rows)),
    )
    B, nu, eps = Moose.emissivity_grid(df, df.e_perp .- df.e_para)
    @test size(eps) == (length(Bvals), length(nuvals))
    @test all(eps[i, j] == B[i] + 1000 * nu[j] for i in eachindex(B), j in eachindex(nu))

    df_incomplete = df[1:end-1, :]
    @test_throws ErrorException Moose.emissivity_grid(df_incomplete, df_incomplete.e_perp)
end

@testset "Regression — header cache is thread-safe (BUG-5)" begin
    empty!(Moose._HEADER_PARAMS_CACHE)
    tasks = [Threads.@spawn Moose._header_params_cached("Qnu") for _ in 1:8]
    results = fetch.(tasks)
    @test all(r -> r["bunit"] == "K" && r["cunit3"] == "Hz", results)
end

@testset "Regression — FITS header units (BUG-6)" begin
    @test Moose.DictHeader["Pnumax"]["bunit"] == "K"
    @test Moose.DictHeader["polfrac"]["bunit"] == ""
    @test Moose.DictHeader["polfracmax"]["bunit"] == ""
    @test Moose.DictHeader["intBLOS"]["bunit"] == "muG cm"
    @test Moose.DictHeader["intBtotal"]["bunit"] == "muG cm"
    @test Moose.DictHeader["intBperp"]["bunit"] == "muG cm"
    @test Moose.DictHeader["ne"]["ctype3"] == "DIST"
    # Single-channel spectral axis must not throw (CDELT3 falls back to 0.0).
    h = Moose.buildHeader3D(3, (2, 2, 1), "", "", "FREQ", "", "", "Hz", "K", [1.0e8])
    @test h["CDELT3"] == 0.0
    @test h["CRVAL3"] == 1.0e8

    h_meta = Moose.buildHeader3D(
        3,
        (2, 2, 1),
        "",
        "",
        "FREQ",
        "",
        "",
        "Hz",
        "K",
        [1.0e8];
        metadata = Dict("MOOSEV" => "1.0.0", "CFGHASH" => repeat("a", 64), "SKIPME" => nothing),
    )
    @test h_meta["MOOSEV"] == "1.0.0"
    @test h_meta["CFGHASH"] == repeat("a", 64)
    @test_throws KeyError h_meta["SKIPME"]
end

@testset "Unit constants" begin
    @test Moose.C_m == 2.99792458e8
    @test Moose.C == 2.99792458e10
end

@testset "Regression — FITS cube reading fails early with actionable validation" begin
    mktempdir() do dir
        two_d_path = joinpath(dir, "map.fits")
        write_test_fits(two_d_path, ones(2, 2))
        err = try
            Moose.read_file(two_d_path, 1.0; expected_ndims = 3)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("expected 3D", sprint(showerror, err))

        nan_path = joinpath(dir, "nan_cube.fits")
        cube = ones(2, 2, 2)
        cube[1, 1, 1] = NaN
        write_test_fits(nan_path, cube)
        err = try
            Moose.read_file(nan_path, 1.0; expected_ndims = 3)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("non-finite", sprint(showerror, err))
    end
end

@testset "HDF5 simulation cube reading" begin
    cube = reshape(Float64.(1:8), 2, 2, 2)

    mktempdir() do sim_dir
        bx = cube
        by = cube .+ 10
        bz = cube .+ 20
        density = cube .+ 30
        temperature = cube .+ 40
        density_hp = cube .+ 50

        write_test_hdf5(joinpath(sim_dir, "Bx.h5"), bx)
        write_test_hdf5(joinpath(sim_dir, "By.hdf5"), by)
        write_test_hdf5(joinpath(sim_dir, "Bz.h5"), bz)
        write_test_hdf5(joinpath(sim_dir, "density.h5"), density)
        write_test_hdf5(joinpath(sim_dir, "temperature.h5"), temperature)
        write_test_hdf5(joinpath(sim_dir, "densityHp.h5"), density_hp)

        @test Moose.simulation_grid_kind(sim_dir) == :image
        @test Moose.read_file(joinpath(sim_dir, "Bx.h5"), 2.0; expected_ndims=3) == 2.0 .* bx

        B1, B2, BLOS, T, n, nH2, nHp = Moose.ReadSimulation(sim_dir, "y", 2.0, 3.0, 4.0)
        @test B1 == permutedims(4.0 .* bz, [3, 1, 2])
        @test B2 == permutedims(4.0 .* bx, [3, 1, 2])
        @test BLOS == permutedims(4.0 .* by, [3, 1, 2])
        @test T == permutedims(3.0 .* temperature, [3, 1, 2])
        @test n == permutedims(2.0 .* density, [3, 1, 2])
        @test nH2 === nothing
        @test nHp == permutedims(2.0 .* density_hp, [3, 1, 2])
    end

    mktempdir() do sim_dir
        shared_path = joinpath(sim_dir, "simulation.hdf5")
        h5open(shared_path, "w") do h5
            h5["fields/Bx"] = cube
            h5["fields/By"] = cube .+ 1
            h5["fields/Bz"] = cube .+ 2
            h5["gas/density"] = cube .+ 3
            h5["gas/temperature"] = cube .+ 4
        end

        @test Moose.source_label(Moose.simulation_field_source(sim_dir, "Bx")) == "$(shared_path):fields/Bx"
        B1, B2, BLOS, T, n, nH2, nHp = Moose.ReadSimulation(sim_dir, "z", 1.0, 1.0, 1.0)
        @test B1 == cube
        @test B2 == cube .+ 1
        @test BLOS == cube .+ 2
        @test T == cube .+ 4
        @test n == cube .+ 3
        @test nH2 === nothing
        @test nHp === nothing
    end
end

@testset "Configurable simulation field sources" begin
    cube = reshape(Float64.(1:8), 2, 2, 2)

    mktempdir() do sim_dir
        write_test_fits(joinpath(sim_dir, "magnetic_x.fits"), cube)
        write_test_fits(joinpath(sim_dir, "magnetic_y.fits"), cube .+ 10)
        write_test_fits(joinpath(sim_dir, "magnetic_z.fits"), cube .+ 20)
        write_test_fits(joinpath(sim_dir, "nH.fits"), cube .+ 30)
        write_test_fits(joinpath(sim_dir, "temp.fits"), cube .+ 40)
        write_test_fits(joinpath(sim_dir, "electron.fits"), cube .+ 50)

        field_sources = Dict(
            "Bx" => "magnetic_x.fits",
            "By" => "magnetic_y.fits",
            "Bz" => "magnetic_z.fits",
            "density" => "nH.fits",
            "temperature" => "temp.fits",
            "densityHp" => "electron.fits",
        )

        @test Moose.source_label(Moose.simulation_field_source(sim_dir, "Bx", field_sources)) == joinpath(sim_dir, "magnetic_x.fits")
        @test Moose.simulation_grid_kind(sim_dir, field_sources) == :image

        B1, B2, BLOS, T, n, nH2, nHp = Moose.ReadSimulation(sim_dir, "z", 2.0, 3.0, 4.0; field_sources = field_sources)
        @test B1 == 4.0 .* cube
        @test B2 == 4.0 .* (cube .+ 10)
        @test BLOS == 4.0 .* (cube .+ 20)
        @test T == 3.0 .* (cube .+ 40)
        @test n == 2.0 .* (cube .+ 30)
        @test nH2 === nothing
        @test nHp == 2.0 .* (cube .+ 50)
    end

    mktempdir() do sim_dir
        shared_path = joinpath(sim_dir, "simulation.h5")
        h5open(shared_path, "w") do h5
            h5["mhd/mag_x"] = cube
            h5["mhd/mag_y"] = cube .+ 1
            h5["mhd/mag_z"] = cube .+ 2
            h5["gas/nH"] = cube .+ 3
            h5["gas/temp"] = cube .+ 4
        end

        field_sources = Dict(
            "Bx" => Dict("file" => "simulation.h5", "dataset" => "mhd/mag_x"),
            "By" => Dict("file" => "simulation.h5", "dataset" => "mhd/mag_y"),
            "Bz" => Dict("file" => "simulation.h5", "dataset" => "mhd/mag_z"),
            "density" => Dict("file" => "simulation.h5", "dataset" => "gas/nH"),
            "temperature" => Dict("file" => "simulation.h5", "dataset" => "gas/temp"),
        )

        @test Moose.source_label(Moose.simulation_field_source(sim_dir, "Bx", field_sources)) == "$(shared_path):mhd/mag_x"
        B1, B2, BLOS, T, n, nH2, nHp = Moose.ReadSimulation(sim_dir, "z", 1.0, 1.0, 1.0; field_sources = field_sources)
        @test B1 == cube
        @test B2 == cube .+ 1
        @test BLOS == cube .+ 2
        @test T == cube .+ 4
        @test n == cube .+ 3
        @test nH2 === nothing
        @test nHp === nothing
    end
end

@testset "AMR leaf-cell HDF5 inputs" begin
    mktempdir() do sim_dir
        path = joinpath(sim_dir, "amr.h5")
        centers = collect(Iterators.product((0.25, 0.75), (0.25, 0.75), (0.25, 0.75)))
        x = [point[1] for point in centers]
        y = [point[2] for point in centers]
        z = [point[3] for point in centers]
        base = Float64.(1:length(centers))
        h5open(path, "w") do h5
            h5["cells/x"] = x
            h5["cells/y"] = y
            h5["cells/z"] = z
            h5["cells/level"] = fill(1, length(centers))
            h5["cells/dx"] = fill(0.5, length(centers))
            h5["cells/Bx"] = base
            h5["cells/By"] = base .+ 10
            h5["cells/Bz"] = base .+ 20
            h5["cells/density"] = base .+ 30
            h5["cells/temperature"] = base .+ 40
        end

        field_sources = Dict{String, Any}(
            field => Dict("path" => "amr.h5", "dataset" => "cells/$(dataset)")
            for (field, dataset) in (
                "Bx" => "Bx", "By" => "By", "Bz" => "Bz",
                "density" => "density", "temperature" => "temperature",
            )
        )
        field_sources["amr"] = Dict(
            "file" => "amr.h5",
            "x" => "cells/x", "y" => "cells/y", "z" => "cells/z",
            "level" => "cells/level",
            "bounds" => [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]],
            "shape" => [2, 2, 2],
        )

        @test Moose.simulation_grid_kind(sim_dir, field_sources) == :amr
        B1, B2, BLOS, T, n, nH2, nHp = Moose.ReadSimulation(
            sim_dir, "z", 2.0, 3.0, 4.0; field_sources = field_sources)
        expected = reshape(base, 2, 2, 2)
        @test B1 == 4.0 .* expected
        @test B2 == 4.0 .* (expected .+ 10)
        @test BLOS == 4.0 .* (expected .+ 20)
        @test T == 3.0 .* (expected .+ 40)
        @test n == 2.0 .* (expected .+ 30)
        @test nH2 === nothing
        @test nHp === nothing

        # The same AMR volume can be observed along any cartesian LOS.
        B1x, B2x, BLOSx, _, _, _, _ = Moose.ReadSimulation(
            sim_dir, "x", 1.0, 1.0, 1.0; field_sources = field_sources)
        @test B1x == permutedims(expected .+ 10, [2, 3, 1])
        @test B2x == permutedims(expected .+ 20, [2, 3, 1])
        @test BLOSx == permutedims(expected, [2, 3, 1])

        # Alternate geometry file, explicit sizes, level offsets, and textual
        # booleans all use the same exact cell-to-voxel assignment.
        geometry_path = joinpath(sim_dir, "geometry.h5")
        h5open(geometry_path, "w") do h5
            h5["cells/x"] = x
            h5["cells/y"] = y
            h5["cells/z"] = z
            h5["cells/dx"] = fill(0.5, length(centers))
        end
        size_sources = deepcopy(field_sources)
        size_sources["amr"] = Dict(
            "file" => "geometry.h5",
            "x" => "cells/x", "y" => "cells/y", "z" => "cells/z",
            "size" => "cells/dx", "shape" => [2, 2, 2], "strict" => "false",
        )
        size_plan = Moose.load_amr_raster_plan(sim_dir, size_sources["amr"], path)
        @test !size_plan.geometry.strict
        @test Moose.rasterize_amr_field(base, size_plan) == expected

        offset_sources = deepcopy(field_sources)
        offset_sources["amr"]["level_offset"] = 1
        h5open(path, "r+") do h5
            h5["cells/level"][:] = fill(2, length(centers))
        end
        offset_plan = Moose.load_amr_raster_plan(sim_dir, offset_sources["amr"], path)
        @test Moose.rasterize_amr_field(base, offset_plan) == expected

        # Structured AMR errors survive the generic grid-reader boundary.
        misaligned_sources = deepcopy(field_sources)
        misaligned_sources["amr"]["level_offset"] = 1
        misaligned_sources["amr"]["shape"] = [3, 3, 3]
        err = try
            Moose.ReadSimulation(sim_dir, "z", 1.0, 1.0, 1.0; field_sources = misaligned_sources)
            nothing
        catch caught
            caught
        end
        @test err isa Moose.MooseError
        @test err.code == :amr_grid_misaligned

        # Full AMR pipeline through science-map FITS output and preflight.
        interpolation_path = joinpath(sim_dir, "emissivity.csv")
        write_test_emissivity(interpolation_path)
        config_path = joinpath(sim_dir, "amr_config.json")
        cfg = Dict(
            "base_dir" => sim_dir,
            "simulations" => [sim_dir],
            "chosen_LOS" => ["z"],
            "field_sources" => offset_sources,
            "interpolation_file_path" => interpolation_path,
            "faraday" => Dict("enabled" => false, "phimin" => -10.0, "phimax" => 10.0, "dphi" => 1.0),
            "ne_option" => "2", "IonizationFraction" => 0.1,
            "freq" => Dict("start" => 100.0, "end" => 101.0, "step" => 1.0),
            "BoxLength_pc" => 2.0, "BoxLength_pix" => 2,
            "log_progress" => false,
        )
        write(config_path, JSON.json(cfg))
        run_cfg, _ = build_config(cfg, config_path)
        plan_io = IOBuffer()
        preflight = Moose.preflight_plan(run_cfg; io=plan_io)
        @test only(preflight.entries).grid == :amr
        @test occursin("AMR geometry/lookup", String(take!(plan_io)))
        Moose.MOOSE_from_config(config_path; quiet=true)
        q_path = joinpath(sim_dir, "z", "Synchrotron", "noFaraday", "Qnu.fits")
        @test isfile(q_path)
        @test FITS(q_path) do fits
            read_header(fits[1])["GRID"] == "AMR"
        end
    end

    # Strict coverage catches a non-leaf hierarchy (coarse and fine cells
    # occupying the same volume) instead of silently double-counting it.
    geometry = Moose.AMRGeometry(
        [0.5 0.5 0.5; 0.25 0.25 0.25],
        [1.0 1.0 1.0; 0.5 0.5 0.5],
        ((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
        (2, 2, 2), true, 1e-8,
    )
    @test_throws Moose.MooseError Moose.rasterize_amr_field([1.0, 2.0], geometry)

    coarse = Moose.AMRGeometry(
        reshape([0.5, 0.5, 0.5], 1, 3),
        reshape([1.0, 1.0, 1.0], 1, 3),
        ((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
        (4, 3, 2), true, 1e-8,
    )
    @test Moose.rasterize_amr_field([7.5], coarse) ≈ fill(7.5, 4, 3, 2)

    # A genuine mixed leaf hierarchy: seven level-1 octants and eight level-2
    # leaves replacing the upper corner.
    coarse_centers = [collect(point) for point in Iterators.product((0.25, 0.75), (0.25, 0.75), (0.25, 0.75)) if point != (0.75, 0.75, 0.75)]
    fine_centers = [collect(point) for point in Iterators.product((0.625, 0.875), (0.625, 0.875), (0.625, 0.875))]
    mixed_centers = permutedims(hcat(coarse_centers..., fine_centers...))
    mixed_widths = vcat(fill(0.5, length(coarse_centers), 3), fill(0.25, length(fine_centers), 3))
    mixed = Moose.AMRGeometry(mixed_centers, mixed_widths,
        ((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)), (4, 4, 4), true, 1e-8)
    mixed_plan = Moose.build_amr_raster_plan(mixed)
    @test all(>(0), mixed_plan.cell_index)
    mixed_map = Moose.rasterize_amr_field(collect(1.0:15.0), mixed_plan)
    @test mixed_map[1:2, 1:2, 1:2] == fill(1.0, 2, 2, 2)
    @test sort(unique(vec(mixed_map[3:4, 3:4, 3:4]))) == collect(8.0:15.0)

    partial = Moose.AMRGeometry(
        reshape([0.25, 0.5, 0.5], 1, 3), reshape([0.5, 1.0, 1.0], 1, 3),
        ((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)), (2, 2, 2), false, 1e-8)
    partial_map = Moose.rasterize_amr_field([3.0], partial)
    @test all(==(3.0), partial_map[1, :, :])
    @test all(isnan, partial_map[2, :, :])

    fine_grid = Moose.AMRGeometry(
        permutedims(hcat(([point...] for point in Iterators.product((0.125, 0.375, 0.625, 0.875), (0.125, 0.375, 0.625, 0.875), (0.125, 0.375, 0.625, 0.875)))...)),
        fill(0.25, 64, 3), ((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
        (2, 2, 2), true, 1e-8)
    err = try
        Moose.build_amr_raster_plan(fine_grid)
        nothing
    catch caught
        caught
    end
    @test err isa Moose.MooseError
    @test err.code == :amr_resolution_too_coarse
end

@testset "Regression — Tnu3D cached path preserves column values" begin
    mktempdir() do dir
        emissivity_path = joinpath(dir, "emissivity.csv")
        write_test_emissivity(emissivity_path)
        df = CSV.File(emissivity_path) |> Moose.DataFrame
        Bperpcube = reshape(collect(range(0.1, 2.5, length = 12)), 2, 3, 2)
        nuArray = [99.0, 100.0, 101.0]
        pixel_length_cm = 2.0

        interpolator = Moose.TemperatureInterpolator(df)
        emissivity_cache = Moose.build_emissivity_frequency_cache(interpolator, nuArray)
        t3d = Moose.Tnu3D(Bperpcube, nuArray, df, pixel_length_cm)

        for i in axes(Bperpcube, 1), j in axes(Bperpcube, 2)
            expected = Moose.Tnu(@view(Bperpcube[i, j, :]), nuArray, df, pixel_length_cm;
                precomputed_interp = interpolator,
                emissivity_cache = emissivity_cache,
            )
            @test t3d[i, j, :] ≈ expected
        end
    end
end

@testset "Regression — QUnu3D cached path preserves column values" begin
    mktempdir() do dir
        emissivity_path = joinpath(dir, "emissivity.csv")
        write_test_emissivity(emissivity_path)
        df = CSV.File(emissivity_path) |> Moose.DataFrame
        Bperpcube = reshape(collect(range(0.1, 2.5, length = 12)), 2, 3, 2)
        psi_src = reshape(collect(range(0.2, 0.8, length = 12)), 2, 3, 2)
        RM = reshape(collect(range(-0.2, 0.3, length = 12)), 2, 3, 2)
        nuArray = [99.0, 100.0, 101.0]
        pixel_length_cm = 2.0

        interpolator = Moose.EmissivityInterpolator(df)
        emissivity_cache = Moose.build_emissivity_frequency_cache(interpolator, nuArray)
        q3d, u3d = Moose.QUnu3D(Bperpcube, psi_src, RM, nuArray, df, pixel_length_cm)
        qnf3d, unf3d = Moose.QUnuNoFaraday3D(Bperpcube, psi_src, nuArray, df, pixel_length_cm)

        function old_buffered_qu!(Qnu, Unu, Bperp, psi, rm)
            eps_buffer = similar(Bperp, Float64)
            for k in eachindex(nuArray)
                nui = nuArray[k]
                cache_col = @view emissivity_cache[:, k]
                Moose.emissivity_at_frequency!(
                    eps_buffer,
                    interpolator.B,
                    interpolator.eps_interp,
                    Bperp,
                    nui;
                    eps_cache_col = cache_col,
                )

                faraday_factor = (Moose.C_m / (nui * 1e6))^2
                sum_u = 0.0
                sum_q = 0.0
                @inbounds for idx in eachindex(eps_buffer, psi, rm)
                    arg = 2.0 * (psi[idx] + rm[idx] * faraday_factor)
                    eps_val = eps_buffer[idx]
                    sum_u += eps_val * sin(arg)
                    sum_q += eps_val * cos(arg)
                end

                Unu[k] = Moose.BrightnessTemperature(nui, sum_u * pixel_length_cm)
                Qnu[k] = Moose.BrightnessTemperature(nui, sum_q * pixel_length_cm)
            end
            return Qnu, Unu
        end

        function old_buffered_qu_no_faraday!(Qnu, Unu, Bperp, psi)
            eps_buffer = similar(Bperp, Float64)
            for k in eachindex(nuArray)
                nui = nuArray[k]
                cache_col = @view emissivity_cache[:, k]
                Moose.emissivity_at_frequency!(
                    eps_buffer,
                    interpolator.B,
                    interpolator.eps_interp,
                    Bperp,
                    nui;
                    eps_cache_col = cache_col,
                )

                sum_u = 0.0
                sum_q = 0.0
                @inbounds for idx in eachindex(eps_buffer, psi)
                    arg = 2.0 * psi[idx]
                    eps_val = eps_buffer[idx]
                    sum_u += eps_val * sin(arg)
                    sum_q += eps_val * cos(arg)
                end

                Unu[k] = Moose.BrightnessTemperature(nui, sum_u * pixel_length_cm)
                Qnu[k] = Moose.BrightnessTemperature(nui, sum_q * pixel_length_cm)
            end
            return Qnu, Unu
        end

        old_q3d = zeros(size(q3d))
        old_u3d = zeros(size(u3d))
        old_qnf3d = zeros(size(qnf3d))
        old_unf3d = zeros(size(unf3d))
        for i in axes(Bperpcube, 1), j in axes(Bperpcube, 2)
            old_buffered_qu!(
                @view(old_q3d[i, j, :]),
                @view(old_u3d[i, j, :]),
                @view(Bperpcube[i, j, :]),
                @view(psi_src[i, j, :]),
                @view(RM[i, j, :]),
            )
            old_buffered_qu_no_faraday!(
                @view(old_qnf3d[i, j, :]),
                @view(old_unf3d[i, j, :]),
                @view(Bperpcube[i, j, :]),
                @view(psi_src[i, j, :]),
            )
        end
        @test q3d == old_q3d
        @test u3d == old_u3d
        @test qnf3d == old_qnf3d
        @test unf3d == old_unf3d

        for i in axes(Bperpcube, 1), j in axes(Bperpcube, 2)
            q, u = Moose.QUnu(
                @view(Bperpcube[i, j, :]),
                @view(psi_src[i, j, :]),
                @view(RM[i, j, :]),
                nuArray,
                df,
                pixel_length_cm;
                precomputed_interp = interpolator,
                emissivity_cache = emissivity_cache,
            )
            @test q3d[i, j, :] ≈ q
            @test u3d[i, j, :] ≈ u

            qnf, unf = Moose.QUnuNoFaraday(
                @view(Bperpcube[i, j, :]),
                @view(psi_src[i, j, :]),
                nuArray,
                df,
                pixel_length_cm;
                precomputed_interp = interpolator,
                emissivity_cache = emissivity_cache,
            )
            @test qnf3d[i, j, :] ≈ qnf
            @test unf3d[i, j, :] ≈ unf
        end
    end
end

@testset "CLI argument validation" begin
    err = try
        parse_cli_args(["--faraday", "maybe"])
        nothing
    catch e
        e
    end

    @test err isa Moose.MooseError
    @test err.code == :cli_invalid_argument
    @test occursin("expects Y or N", err.message)

    _, _, _, overrides = parse_cli_args(["--faraday", "true", "--noise", "0"])
    @test overrides["FaradayRotation"] == "Y"
    @test overrides["add_noise"] == "N"

    _, _, _, overrides = parse_cli_args(["--density-kind", "mass_density", "--mean-molecular-weight", "1.4", "--hydrogen-mass-g", "1.6726231e-24"])
    @test overrides["density_kind"] == "mass_density"
    @test overrides["mean_molecular_weight"] == 1.4
    @test overrides["hydrogen_mass_g"] == 1.6726231e-24

    _, _, _, overrides = parse_cli_args(["--resume", "safe"])
    @test overrides["resume"] == "safe"

    _, _, _, overrides = parse_cli_args(["--plan"])
    @test overrides["__plan_only"] === true

    _, _, _, overrides = parse_cli_args(["--outputs", "rm,fdf"])
    @test overrides["outputs"] == ["rm", "fdf"]

    err = try
        parse_cli_args(["--write-back"])
        nothing
    catch e
        e
    end

    @test err isa Moose.MooseError
    @test err.code == :cli_invalid_argument

    missing_config = joinpath(tempdir(), "moose-missing-config-$(rand(UInt)).json")
    err = try
        parse_cli_args([missing_config])
        nothing
    catch e
        e
    end

    @test err isa Moose.MooseError
    @test err.code == :cli_invalid_argument
    @test occursin("Config file not found", err.message)
end

@testset "Config validation" begin
    @test first(Moose.MooseFromConfig.build_faraday(Dict{String,Any}())) == "Y"
    @test first(Moose.MooseFromConfig.build_faraday(Dict("faraday" => Dict{String,Any}()))) == "Y"

    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "simu1")
        mkdir(sim_dir)

        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "interpolation_file_path" => joinpath(base_dir, "emissivity.csv"),
            "chosen_LOS" => ["x", "invalid"],
        )

        err = try
            build_config(cfg, "config.json")
            nothing
        catch e
            e
        end

        @test err isa Moose.MooseError
        @test occursin("Invalid line(s) of sight", err.message)
    end

    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "simu2")
        mkdir(sim_dir)

        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "interpolation_file_path" => joinpath(base_dir, "emissivity.csv"),
            "freq" => Dict("start" => 150.0, "end" => 140.0, "step" => -1.0),
        )

        err = try
            build_config(cfg, "config.json")
            nothing
        catch e
            e
        end

        @test err isa Moose.MooseError
        @test err.code == :invalid_frequency
        @test occursin("Frequency step", err.message) || occursin("greater than the start frequency", err.message)
    end

    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "simu3")
        mkdir(sim_dir)
        interpolation_path = joinpath(base_dir, "emissivity.csv")
        write(interpolation_path, "B\tnu\te_para\te_perp\n1.0\t120.0\t1.0\t2.0\n")

        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "interpolation_file_path" => interpolation_path,
            "box" => Dict(
                "x" => 30.0,
                "y" => 60.0,
                "z" => 90.0,
                "npix" => 128,
            ),
            "rng_seed" => 1234,
        )

        run_cfg, _ = build_config(cfg, "config.json")
        @test run_cfg.BoxLength_pc.x == 30.0
        @test run_cfg.BoxLength_pc.y == 60.0
        @test run_cfg.BoxLength_pc.z == 90.0
        @test run_cfg.BoxLength_pix == (; x = 128, y = 128, z = 128)
        @test run_cfg.rng_seed == 1234
    end

    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "simu_fields")
        mkdir(sim_dir)
        interpolation_path = joinpath(base_dir, "emissivity.csv")
        write(interpolation_path, "B\tnu\te_para\te_perp\n1.0\t120.0\t1.0\t2.0\n")

        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "interpolation_file_path" => interpolation_path,
            "field_sources" => Dict(
                "Bx" => "mag_x.fits",
                "By" => "mag_y.fits",
                "Bz" => "mag_z.fits",
                "density" => "nH.fits",
                "temperature" => Dict("file" => "gas.h5", "dataset" => "thermal/T"),
            ),
        )

        run_cfg, _ = build_config(cfg, "config.json")
        @test run_cfg.field_sources["Bx"] == "mag_x.fits"
        @test run_cfg.field_sources["temperature"]["dataset"] == "thermal/T"
    end

    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "simu_mask")
        mkdir(sim_dir)
        interpolation_path = joinpath(base_dir, "emissivity.csv")
        write(interpolation_path, "B\tnu\te_para\te_perp\n1.0\t120.0\t1.0\t2.0\n")

        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "interpolation_file_path" => interpolation_path,
            "physical_mask" => Dict(
                "temperature_min" => 100.0,
                "temperature_max" => 1.0e6,
                "nH_min" => 1.0e-3,
                "nH_max" => 10.0,
            ),
        )

        run_cfg, _ = build_config(cfg, "config.json")
        @test run_cfg.physical_mask["T_min"] == 100.0
        @test run_cfg.physical_mask["T_max"] == 1.0e6
        @test run_cfg.physical_mask["n_min"] == 1.0e-3
        @test run_cfg.physical_mask["n_max"] == 10.0

        cfg["physical_mask"]["temperature_min"] = 2.0e6
        err = try
            build_config(cfg, "config.json")
            nothing
        catch e
            e
        end
        @test err isa Moose.MooseError
        @test err.code == :invalid_physical_mask
    end

    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "simu_density_kind")
        mkdir(sim_dir)
        interpolation_path = joinpath(base_dir, "emissivity.csv")
        write(interpolation_path, "B\tnu\te_para\te_perp\n1.0\t120.0\t1.0\t2.0\n")

        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "interpolation_file_path" => interpolation_path,
            "density_kind" => "mass_density",
            "mean_molecular_weight" => 1.4,
            "hydrogen_mass_g" => Moose.M_p,
        )

        run_cfg, _ = build_config(cfg, "config.json")
        @test run_cfg.density_kind == "mass_density"
        @test run_cfg.mean_molecular_weight == 1.4
        @test run_cfg.hydrogen_mass_g == Moose.M_p

        cfg["density_kind"] = "bad"
        err = try
            build_config(cfg, "config.json")
            nothing
        catch e
            e
        end
        @test err isa Moose.MooseError
        @test err.code == :invalid_density_kind
    end

    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "simu4")
        mkdir(sim_dir)
        interpolation_path = joinpath(base_dir, "emissivity.csv")
        write(interpolation_path, "B\tnu\te_para\te_perp\n1.0\t120.0\t1.0\t2.0\n")

        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "interpolation_file_path" => interpolation_path,
            "faraday" => Dict("enabled" => "true", "phimin" => 20.0, "phimax" => -20.0, "dphi" => 0.1),
        )

        err = try
            build_config(cfg, "config.json")
            nothing
        catch e
            e
        end

        @test err isa Moose.MooseError
        @test err.code == :invalid_faraday_range
    end

    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "simu5")
        mkdir(sim_dir)
        interpolation_path = joinpath(base_dir, "emissivity.csv")
        write(interpolation_path, "B\tnu\te_para\te_perp\n1.0\t120.0\t1.0\t2.0\n")

        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "interpolation_file_path" => interpolation_path,
            "faraday" => Dict("enabled" => "N", "phimin" => -20.0, "phimax" => 20.0, "dphi" => 1.0),
            "rm_clean" => Dict("enabled" => true),
        )

        err = try
            build_config(cfg, "config.json")
            nothing
        catch e
            e
        end

        @test err isa Moose.MooseError
        @test err.code == :invalid_rm_clean
    end
end

@testset "Regression — BoxLength_pix must match FITS cube dimensions" begin
    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "tiny_simu")
        mkdir(sim_dir)

        cube = fill(1.0, 2, 2, 2)
        write_test_fits(joinpath(sim_dir, "Bx.fits"), cube)
        write_test_fits(joinpath(sim_dir, "By.fits"), fill(0.5, 2, 2, 2))
        write_test_fits(joinpath(sim_dir, "Bz.fits"), fill(0.25, 2, 2, 2))
        write_test_fits(joinpath(sim_dir, "density.fits"), cube)
        write_test_fits(joinpath(sim_dir, "temperature.fits"), fill(100.0, 2, 2, 2))

        interpolation_path = joinpath(base_dir, "emissivity.csv")
        write_test_emissivity(interpolation_path)

        config_path = joinpath(base_dir, "moose_config.json")
        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "chosen_LOS" => ["z"],
            "interpolation_file_path" => interpolation_path,
            "faraday" => Dict("enabled" => "N", "phimin" => -10.0, "phimax" => 10.0, "dphi" => 1.0),
            "responseSynchrotron" => "N",
            "add_noise" => "N",
            "ne_option" => "2",
            "IonizationFraction" => 0.1,
            "freq" => Dict("start" => 100.0, "end" => 101.0, "step" => 1.0),
            "BoxLength_pc" => 2.0,
            "BoxLength_pix" => Dict("x" => 3, "y" => 2, "z" => 2),
            "log_progress" => false,
        )
        write(config_path, JSON.json(cfg))

        err = try
            Moose.MOOSE_from_config(config_path; quiet = true)
            nothing
        catch e
            e
        end

        @test err isa Moose.MooseError
        @test err.code == :cube_shape_mismatch
        @test occursin("BoxLength_pix", err.message)
    end
end

@testset "Minimal config end-to-end FITS run" begin
    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "tiny_simu")
        mkdir(sim_dir)

        cube = fill(1.0, 2, 2, 2)
        write_test_fits(joinpath(sim_dir, "Bx.fits"), cube)
        write_test_fits(joinpath(sim_dir, "By.fits"), fill(0.5, 2, 2, 2))
        write_test_fits(joinpath(sim_dir, "Bz.fits"), fill(0.25, 2, 2, 2))
        write_test_fits(joinpath(sim_dir, "density.fits"), cube)
        write_test_fits(joinpath(sim_dir, "temperature.fits"), fill(100.0, 2, 2, 2))

        interpolation_path = joinpath(base_dir, "emissivity.csv")
        write_test_emissivity(interpolation_path)

        config_path = joinpath(base_dir, "moose_config.json")
        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "chosen_LOS" => ["z"],
            "interpolation_file_path" => interpolation_path,
            "faraday" => Dict("enabled" => "N", "phimin" => -10.0, "phimax" => 10.0, "dphi" => 1.0),
            "responseSynchrotron" => "false",
            "add_noise" => "N",
            "ne_option" => "2",
            "IonizationFraction" => 0.1,
            "freq" => Dict("start" => 100.0, "end" => 101.0, "step" => 1.0),
            "BoxLength_pc" => 2.0,
            "BoxLength_pix" => 2,
            "log_progress" => false,
            "rng_seed" => 7,
            "resume" => "safe",
        )
        write(config_path, JSON.json(cfg))

        run_cfg, _ = build_config(cfg, config_path)
        plan_io = IOBuffer()
        plan = Moose.preflight_plan(run_cfg; io=plan_io)
        @test plan.frequency_channels == 2
        @test plan.faraday_channels == 0
        @test length(plan.entries) == 1
        @test only(plan.entries).shape == (2, 2, 2)
        @test plan.disk_bytes > 0
        @test occursin("peak RAM", String(take!(plan_io)))

        Moose.MOOSE_from_config(config_path; quiet = true)

        result_dir = joinpath(sim_dir, "z", "Synchrotron", "noFaraday")
        @test isdir(result_dir)
        @test isfile(joinpath(result_dir, "Qnu.fits"))
        @test isfile(joinpath(result_dir, "Unu.fits"))
        @test isfile(joinpath(result_dir, "Pnu.fits"))
        @test isfile(joinpath(result_dir, "Tnu.fits"))
        @test isfile(joinpath(result_dir, "Pnumax.fits"))
        @test isfile(joinpath(result_dir, "polfrac.fits"))
        @test isfile(joinpath(result_dir, "polfracmax.fits"))
        @test isfile(joinpath(result_dir, "alpha.fits"))
        @test isfile(joinpath(result_dir, "alpha_err.fits"))
        @test isfile(joinpath(result_dir, "polarization_angle_vs_lambda2.png"))
        @test isfile(joinpath(result_dir, "fractional_polarization_vs_lambda2.png"))
        @test isfile(joinpath(result_dir, "stokes_qu_diagram.png"))
        @test isfile(joinpath(result_dir, "polarization_diagnostics.png"))
        @test isfile(joinpath(result_dir, "polarization_diagnostics.pdf"))
        @test isfile(joinpath(sim_dir, "z", "Synchrotron", "ne.fits"))
        @test isfile(joinpath(base_dir, "MOOSE_summary.log"))

        manifest_path = joinpath(sim_dir, "z", "Synchrotron", Moose.COMPLETION_MANIFEST)
        @test isfile(manifest_path)
        manifest_before = read(manifest_path, String)
        manifest = JSON.parse(manifest_before)
        @test manifest["status"] == "complete"
        @test manifest["los"] == "z"
        @test !isempty(manifest["outputs"])

        # An identical safe-resume run must skip the completed LOS. The
        # completion timestamp therefore remains byte-for-byte unchanged.
        Moose.MOOSE_from_config(config_path; quiet = true)
        @test read(manifest_path, String) == manifest_before

        qnu = read(FITS(joinpath(result_dir, "Qnu.fits"))[1])
        @test size(qnu) == (2, 2, 2)
        @test all(isfinite, qnu)
        unu = read(FITS(joinpath(result_dir, "Unu.fits"))[1])
        tnu = read(FITS(joinpath(result_dir, "Tnu.fits"))[1])
        polfrac = read(FITS(joinpath(result_dir, "polfrac.fits"))[1])
        polfracmax = read(FITS(joinpath(result_dir, "polfracmax.fits"))[1])
        expected_polfrac = sqrt.(qnu .^ 2 .+ unu .^ 2) ./ tnu
        @test all(isapprox.(polfrac, expected_polfrac; rtol = 1e-10))
        @test polfracmax ≈ dropdims(maximum(polfrac, dims = 3), dims = 3)

        alpha = read(FITS(joinpath(result_dir, "alpha.fits"))[1])
        alpha_err = read(FITS(joinpath(result_dir, "alpha_err.fits"))[1])
        @test size(alpha) == (2, 2)
        @test all(isfinite, alpha)
        # Two frequency channels leave no residual dof for the slope error.
        @test all(isnan, alpha_err)
        alpha_header = FITS(joinpath(result_dir, "alpha.fits")) do fits
            read_header(fits[1])
        end
        @test alpha_header["ALPHADEF"] == "S_nu ~ nu^alpha; alpha = beta_Tb + 2"

        q_header = FITS(joinpath(result_dir, "Qnu.fits")) do fits
            read_header(fits[1])
        end
        @test q_header["MOOSEV"] == Moose.moose_version()
        @test q_header["CFGHASH"] == Moose.moose_config_hash(cfg)
        @test q_header["LOS"] == "z"
        @test q_header["FARADAY"] == "N"
        @test q_header["FILTER"] == "N"
        @test q_header["NOISE"] == "N"
        @test q_header["NUNIT"] == "MHz input; Hz FITS"

        summary = read(joinpath(base_dir, "MOOSE_summary.log"), String)
        @test occursin("Config read: $(config_path)", summary)
        @test occursin("Config effective: $(config_path)", summary)
        @test occursin("Config saved: $(config_path)", summary)
        @test occursin("Config hash: $(Moose.moose_config_hash(cfg))", summary)

        config_before_cli = read(config_path, String)
        run_with_config(config_path, true, false, Dict{String, Any}("chosen_LOS" => ["z"]))
        @test read(config_path, String) == config_before_cli

        summary_after_cli = read(joinpath(base_dir, "MOOSE_summary.log"), String)
        @test occursin("Config effective: $(config_path) + CLI overrides", summary_after_cli)
        @test occursin("Config saved: <not written>", summary_after_cli)
    end
end

@testset "Selective RM-only outputs" begin
    mktempdir() do dir
        demo = make_demo_data(joinpath(dir, "selective"); npix=2)
        cfg = JSON.parsefile(demo.config_path)
        cfg["outputs"] = ["rm"]
        cfg["rm_clean"] = Dict("enabled" => false)
        write(demo.config_path, JSON.json(cfg))

        Moose.MOOSE_from_config(demo.config_path; quiet=true)
        root = joinpath(demo.simulation_dir, "z", "Synchrotron")
        result = joinpath(root, "WithFaraday")
        @test isfile(joinpath(result, "RMmap.fits"))
        @test !isfile(joinpath(result, "Qnu.fits"))
        @test !isfile(joinpath(result, "FDF.fits"))
        @test !isfile(joinpath(result, "alpha.fits"))
        @test !isfile(joinpath(root, "intBtotal.fits"))
    end
end

@testset "Config-driven RM-CLEAN FITS outputs" begin
    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "tiny_clean")
        mkdir(sim_dir)

        cube = fill(1.0, 2, 2, 2)
        write_test_fits(joinpath(sim_dir, "Bx.fits"), cube)
        write_test_fits(joinpath(sim_dir, "By.fits"), fill(0.5, 2, 2, 2))
        write_test_fits(joinpath(sim_dir, "Bz.fits"), fill(0.25, 2, 2, 2))
        write_test_fits(joinpath(sim_dir, "density.fits"), cube)
        write_test_fits(joinpath(sim_dir, "temperature.fits"), fill(100.0, 2, 2, 2))

        interpolation_path = joinpath(base_dir, "emissivity.csv")
        write_test_emissivity(interpolation_path)

        config_path = joinpath(base_dir, "moose_config.json")
        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "chosen_LOS" => ["z"],
            "interpolation_file_path" => interpolation_path,
            "faraday" => Dict("enabled" => "Y", "phimin" => -2.0, "phimax" => 2.0, "dphi" => 2.0),
            "rm_clean" => Dict("enabled" => true, "gain" => 0.1, "niter" => 5, "threshold" => 0.0),
            "responseSynchrotron" => "N",
            "add_noise" => "N",
            "ne_option" => "2",
            "IonizationFraction" => 0.1,
            "freq" => Dict("start" => 100.0, "end" => 101.0, "step" => 1.0),
            "BoxLength_pc" => 2.0,
            "BoxLength_pix" => 2,
            "log_progress" => false,
        )
        write(config_path, JSON.json(cfg))

        Moose.MOOSE_from_config(config_path; quiet = true)

        result_dir = joinpath(sim_dir, "z", "Synchrotron", "WithFaraday")
        @test isfile(joinpath(result_dir, "FDF.fits"))
        @test isfile(joinpath(result_dir, "RMSF.fits"))
        @test isfile(joinpath(result_dir, "cleanFDF.fits"))
        @test isfile(joinpath(result_dir, "realCleanFDF.fits"))
        @test isfile(joinpath(result_dir, "imagCleanFDF.fits"))
        @test isfile(joinpath(result_dir, "residualFDF.fits"))

        clean = read(FITS(joinpath(result_dir, "cleanFDF.fits"))[1])
        residual = read(FITS(joinpath(result_dir, "residualFDF.fits"))[1])
        @test size(clean) == (2, 2, 3)
        @test size(residual) == (2, 2, 3)
        @test all(isfinite, clean)
        @test all(isfinite, residual)

        header = FITS(joinpath(result_dir, "cleanFDF.fits")) do fits
            read_header(fits[1])
        end
        @test header["RMCLEAN"] == true
        @test header["RMCNITER"] == 5
    end
end

@testset "Minimal config end-to-end HEALPix run" begin
    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "tiny_hp")
        mkdir(sim_dir)

        npix = 12
        write_test_healpix_cube(joinpath(sim_dir, "Bx.fits"), hcat(fill(1.0, npix), fill(1.2, npix)))
        write_test_healpix_cube(joinpath(sim_dir, "By.fits"), hcat(fill(0.5, npix), fill(0.4, npix)))
        write_test_healpix_cube(joinpath(sim_dir, "Bz.fits"), hcat(fill(0.25, npix), fill(0.3, npix)))
        write_test_healpix_cube(joinpath(sim_dir, "density.fits"), hcat(fill(1.0, npix), fill(1.1, npix)))
        write_test_healpix_cube(joinpath(sim_dir, "temperature.fits"), hcat(fill(100.0, npix), fill(110.0, npix)))

        interpolation_path = joinpath(base_dir, "emissivity.csv")
        write_test_emissivity(interpolation_path)

        config_path = joinpath(base_dir, "moose_config.json")
        cfg = Dict(
            "base_dir" => base_dir,
            "simulations" => [sim_dir],
            "chosen_LOS" => ["z"],
            "interpolation_file_path" => interpolation_path,
            "faraday" => Dict("enabled" => "Y", "phimin" => -2.0, "phimax" => 2.0, "dphi" => 2.0),
            "rm_clean" => Dict("enabled" => true, "gain" => 0.1, "niter" => 5, "threshold" => 0.0),
            "responseSynchrotron" => "N",
            "add_noise" => "N",
            "ne_option" => "2",
            "IonizationFraction" => 0.1,
            "freq" => Dict("start" => 100.0, "end" => 101.0, "step" => 1.0),
            "BoxLength_pc" => 2.0,
            "BoxLength_pix" => 2,
            "log_progress" => false,
        )
        write(config_path, JSON.json(cfg))

        Moose.MOOSE_from_config(config_path; quiet = true)

        root_dir = joinpath(sim_dir, "z", "Synchrotron")
        result_dir = joinpath(root_dir, "WithFaraday")
        stack = Moose.read_fits_grid_stack(joinpath(sim_dir, "Bx.fits"); column=:all)
        @test stack isa Moose.HealpixStack
        @test size(stack) == (npix, 2)
        @test Moose.detect_fits_grid(joinpath(root_dir, "intBtotal.fits")) == :healpix
        @test isfile(joinpath(root_dir, "ne_0001_p0p0.fits"))
        @test isfile(joinpath(root_dir, "ne_0002_p1p0.fits"))
        @test isfile(joinpath(result_dir, "Qnu_0001_p1p0e8.fits"))
        @test isfile(joinpath(result_dir, "Unu_0002_p1p01e8.fits"))
        @test isfile(joinpath(result_dir, "Pnumax.fits"))
        @test isfile(joinpath(result_dir, "polfrac_0001_p1p0e8.fits"))
        @test isfile(joinpath(result_dir, "polfracmax.fits"))
        @test isfile(joinpath(result_dir, "alpha.fits"))
        @test isfile(joinpath(result_dir, "alpha_err.fits"))
        @test isfile(joinpath(result_dir, "FDF_0002_p0p0.fits"))
        @test isfile(joinpath(result_dir, "realFDF_0001_m2p0.fits"))
        @test isfile(joinpath(result_dir, "imagFDF_0003_p2p0.fits"))
        @test isfile(joinpath(result_dir, "cleanFDF_0002_p0p0.fits"))
        @test isfile(joinpath(result_dir, "realCleanFDF_0001_m2p0.fits"))
        @test isfile(joinpath(result_dir, "imagCleanFDF_0003_p2p0.fits"))
        @test isfile(joinpath(result_dir, "residualFDF_0002_p0p0.fits"))

        q_map = Moose.read_healpix_map(joinpath(result_dir, "Qnu_0001_p1p0e8.fits"))
        @test length(q_map) == npix
        @test all(isfinite, collect(q_map))
    end
end

@testset "Demo data generator with known results" begin
    # Faraday-enabled demo (Faraday screen + uniform emitter): RM map,
    # spectral index, Tnu, Q/U per channel and integrated maps are exact;
    # the FDF peak position is quantized on the dphi grid.
    mktempdir() do dir
        demo = make_demo_data(dir; npix = 4)
        @test isfile(demo.config_path)
        @test isfile(demo.emissivity_path)
        @test isfile(joinpath(demo.base_dir, "expected_results.json"))

        Moose.MOOSE_from_config(demo.config_path; quiet = true)

        result_dir = joinpath(demo.simulation_dir, "z", "Synchrotron", "WithFaraday")
        @test isdir(result_dir)

        rmmap = read(FITS(joinpath(result_dir, "RMmap.fits"))[1])
        @test all(isapprox.(rmmap, demo.expected.rm; rtol = 1e-10))

        alpha = read(FITS(joinpath(result_dir, "alpha.fits"))[1])
        @test all(isapprox.(alpha, demo.expected.alpha; atol = 1e-8))

        tnu = read(FITS(joinpath(result_dir, "Tnu.fits"))[1])
        qnu = read(FITS(joinpath(result_dir, "Qnu.fits"))[1])
        unu = read(FITS(joinpath(result_dir, "Unu.fits"))[1])
        polfrac = read(FITS(joinpath(result_dir, "polfrac.fits"))[1])
        @test size(tnu, 3) == length(demo.expected.nu_MHz)
        for k in axes(tnu, 3)
            @test all(isapprox.(tnu[:, :, k], demo.expected.Tnu[k]; rtol = 1e-8))
            @test all(isapprox.(qnu[:, :, k] ./ tnu[:, :, k], demo.expected.qnu_over_tnu[k]; atol = 1e-8))
            @test all(isapprox.(unu[:, :, k] ./ tnu[:, :, k], demo.expected.unu_over_tnu[k]; atol = 1e-8))
        end

        # Faraday-thin emitter: no depolarization at any frequency.
        pol_fraction = sqrt.(qnu .^ 2 .+ unu .^ 2) ./ tnu
        @test all(isapprox.(pol_fraction, demo.expected.pol_fraction; rtol = 1e-8))
        @test all(isapprox.(polfrac, demo.expected.pol_fraction; rtol = 1e-8))

        root_dir = joinpath(demo.simulation_dir, "z", "Synchrotron")
        intne = read(FITS(joinpath(root_dir, "intne.fits"))[1])
        @test all(isapprox.(intne, demo.expected.intne; rtol = 1e-10))
        intblos = read(FITS(joinpath(root_dir, "intBLOS.fits"))[1])
        @test all(isapprox.(intblos, demo.expected.intBLOS; rtol = 1e-10))

        # The |FDF| peak of the Faraday-thin emitter sits at RM_screen,
        # within one Faraday-depth grid step.
        fdf = read(FITS(joinpath(result_dir, "FDF.fits"))[1])
        phi = collect(-5.0:0.25:5.0)
        @test size(fdf, 3) == length(phi)
        peak_phi = phi[argmax(fdf[1, 1, :])]
        @test abs(peak_phi - demo.expected.fdf_peak_phi) <= 0.25
    end

    # Faraday-disabled demo: the intrinsic polarization observables are exact.
    mktempdir() do dir
        demo = make_demo_data(dir; npix = 4, faraday = false)
        @test demo.expected.fdf_peak_phi === nothing

        Moose.MOOSE_from_config(demo.config_path; quiet = true)

        result_dir = joinpath(demo.simulation_dir, "z", "Synchrotron", "noFaraday")
        qnu = read(FITS(joinpath(result_dir, "Qnu.fits"))[1])
        unu = read(FITS(joinpath(result_dir, "Unu.fits"))[1])
        tnu = read(FITS(joinpath(result_dir, "Tnu.fits"))[1])
        polfrac = read(FITS(joinpath(result_dir, "polfrac.fits"))[1])

        # psi_src = π modulo π ⇒ Q/T = +pol_fraction, U ≈ 0.
        @test all(isapprox.(qnu ./ tnu, demo.expected.pol_fraction; atol = 1e-8))
        @test all(abs.(unu) .<= 1e-10 .* tnu)

        pol_fraction = sqrt.(qnu .^ 2 .+ unu .^ 2) ./ tnu
        @test all(isapprox.(pol_fraction, demo.expected.pol_fraction; rtol = 1e-8))
        @test all(isapprox.(polfrac, demo.expected.pol_fraction; rtol = 1e-8))

        pol_angle = 0.5 .* atan.(unu, qnu)
        @test all(isapprox.(abs.(pol_angle), demo.expected.intrinsic_pol_angle; atol = 1e-8))
    end

    # Input validation.
    mktempdir() do dir
        @test_throws Exception make_demo_data(dir; npix = 1)
        @test_throws Exception make_demo_data(dir; B_perp_uG = 2.5)
        @test_throws Exception make_demo_data(dir; nu_MHz = [100.0, 125.0, 150.0])
        @test_throws Exception make_demo_data(dir; pol_fraction = 1.5)
        @test_throws Exception make_demo_data(dir; B_nodes_uG = 1.0:1.0:6.0)
    end
end

@testset "Precision option (float32)" begin
    mktempdir() do dir
        demo = make_demo_data(dir; npix = 4)
        cfg = JSON.parsefile(demo.config_path)
        cfg["precision"] = "float32"
        write(demo.config_path, JSON.json(cfg))

        Moose.MOOSE_from_config(demo.config_path; quiet = true)

        result_dir = joinpath(demo.simulation_dir, "z", "Synchrotron", "WithFaraday")
        qnu = read(FITS(joinpath(result_dir, "Qnu.fits"))[1])
        tnu = read(FITS(joinpath(result_dir, "Tnu.fits"))[1])
        fdf = read(FITS(joinpath(result_dir, "FDF.fits"))[1])
        rmmap = read(FITS(joinpath(result_dir, "RMmap.fits"))[1])

        # The cubes are stored in single precision...
        @test eltype(qnu) == Float32
        @test eltype(tnu) == Float32
        @test eltype(fdf) == Float32
        @test eltype(rmmap) == Float32

        # ... and still match the analytic expectations at float32 accuracy.
        @test all(isapprox.(rmmap, demo.expected.rm; rtol = 1e-5))
        for k in axes(tnu, 3)
            @test all(isapprox.(tnu[:, :, k], demo.expected.Tnu[k]; rtol = 1e-5))
            @test all(isapprox.(qnu[:, :, k] ./ tnu[:, :, k], demo.expected.qnu_over_tnu[k]; atol = 1e-4))
        end
        alpha = read(FITS(joinpath(result_dir, "alpha.fits"))[1])
        @test all(isapprox.(alpha, demo.expected.alpha; atol = 1e-3))

        header = FITS(joinpath(result_dir, "Qnu.fits")) do fits
            read_header(fits[1])
        end
        @test header["PRECIS"] == "float32"
    end

    mktempdir() do dir
        demo = make_demo_data(dir; npix = 4)
        cfg = JSON.parsefile(demo.config_path)
        cfg["precision"] = "float16"
        write(demo.config_path, JSON.json(cfg))
        @test_throws Exception Moose.MOOSE_from_config(demo.config_path; quiet = true)
    end
end

@testset "Tiled processing option" begin
    # A tiled run must reproduce the plain run (same per-pixel math), and
    # therefore also the analytic demo expectations. npix = 5 with
    # tile_size = 2 exercises an uneven final band.
    mktempdir() do dir
        ref = make_demo_data(joinpath(dir, "ref"); npix = 5)
        Moose.MOOSE_from_config(ref.config_path; quiet = true)

        tiled = make_demo_data(joinpath(dir, "tiled"); npix = 5)
        cfg = JSON.parsefile(tiled.config_path)
        cfg["tile_size"] = 2
        write(tiled.config_path, JSON.json(cfg))
        Moose.MOOSE_from_config(tiled.config_path; quiet = true)

        products = [
            joinpath("WithFaraday", name) for name in
            ("Qnu.fits", "Unu.fits", "Tnu.fits", "Pnu.fits", "polfrac.fits", "FDF.fits",
             "realFDF.fits", "imagFDF.fits", "RMmap.fits", "Pnumax.fits", "polfracmax.fits",
             "Pmax.fits", "alpha.fits", "alpha_err.fits")
        ]
        append!(products, ["ne.fits", "intBtotal.fits", "sigmaBtotal.fits", "intne.fits",
                           "sigmane.fits", "sigmaT.fits", "intBLOS.fits", "sigmaBLOS.fits", "intBperp.fits"])

        for rel in products
            a = read(FITS(joinpath(ref.simulation_dir, "z", "Synchrotron", rel))[1])
            b = read(FITS(joinpath(tiled.simulation_dir, "z", "Synchrotron", rel))[1])
            @test size(a) == size(b)
            scale = maximum(abs, a[isfinite.(a)]; init = 0.0)
            @test all(@. (isnan(a) & isnan(b)) | (abs(a - b) <= 1e-10 * scale + 1e-10 * abs(a)))
        end

        # The tiled run also matches the analytic expectations directly.
        result_dir = joinpath(tiled.simulation_dir, "z", "Synchrotron", "WithFaraday")
        rmmap = read(FITS(joinpath(result_dir, "RMmap.fits"))[1])
        @test all(isapprox.(rmmap, tiled.expected.rm; rtol = 1e-10))
        alpha = read(FITS(joinpath(result_dir, "alpha.fits"))[1])
        @test all(isapprox.(alpha, tiled.expected.alpha; atol = 1e-8))

        header = FITS(joinpath(result_dir, "Qnu.fits")) do fits
            read_header(fits[1])
        end
        @test header["TILESIZE"] == 2
        @test header["CFGHASH"] isa String

        # No leftover streaming scratch files.
        @test isempty(filter(name -> endswith(name, ".part"), readdir(result_dir)))
    end

    # A safe tiled run resumes from the last durable band checkpoint. Build a
    # realistic interrupted state from a completed small run, then verify the
    # checkpoint is consumed and the final products remain identical.
    mktempdir() do dir
        demo = make_demo_data(dir; npix = 5)
        cfg = JSON.parsefile(demo.config_path)
        cfg["tile_size"] = 2
        cfg["resume"] = "safe"
        write(demo.config_path, JSON.json(cfg))
        Moose.MOOSE_from_config(demo.config_path; quiet = true)

        root = joinpath(demo.simulation_dir, "z", "Synchrotron")
        result_dir = joinpath(root, "WithFaraday")
        cube_relpaths = [
            "ne.fits",
            joinpath.(Ref("WithFaraday"),
                ["Qnu.fits", "Unu.fits", "Tnu.fits", "Pnu.fits", "polfrac.fits",
                 "FDF.fits", "realFDF.fits", "imagFDF.fits"])...,
        ]
        expected_cubes = Dict(rel => read(FITS(joinpath(root, rel))[1]) for rel in cube_relpaths)

        map_paths = Dict(
            name => joinpath(root, "$(name).fits") for name in
            ("intBtotal", "sigmaBtotal", "intne", "sigmane", "sigmaT", "intBLOS", "sigmaBLOS", "intBperp")
        )
        merge!(map_paths, Dict(
            name => joinpath(result_dir, "$(name).fits") for name in
            ("RMmap", "Pnumax", "polfracmax", "Pmax", "alpha", "alpha_err")
        ))
        maps2d = Dict(name => read(FITS(path)[1]) for (name, path) in map_paths)

        run_cfg, _ = build_config(cfg, demo.config_path)
        resume_cfg = Dict{String, Any}(String(k) => v for (k, v) in cfg)
        delete!(resume_cfg, "resume")
        resume_hash = Moose.moose_config_hash(resume_cfg)
        inputs = [Moose._fingerprint_file(path) for path in Moose._resume_input_paths(run_cfg, demo.simulation_dir)]
        canonical_inputs = [[input["path"], input["size"], input["mtime"]] for input in inputs]
        signature = bytes2hex(Moose.sha256(resume_hash * JSON.json(canonical_inputs)))

        rm(joinpath(root, Moose.COMPLETION_MANIFEST))
        for rel in cube_relpaths
            path = joinpath(root, rel)
            mv(path, path * ".part")
        end
        Moose._write_tile_checkpoint(joinpath(root, Moose.TILE_CHECKPOINT), signature, 1, 3, maps2d)

        @test_logs match_mode = :any (:info, "Resuming tiled processing from checkpoint")
            Moose.MOOSE_from_config(demo.config_path; quiet = true)

        @test !isfile(joinpath(root, Moose.TILE_CHECKPOINT))
        @test isempty(filter(path -> endswith(path, ".part"),
            [joinpath(dirpath, name) for (dirpath, _, names) in walkdir(root) for name in names]))
        for (rel, expected) in expected_cubes
            actual = read(FITS(joinpath(root, rel))[1])
            @test all(@. (isnan(actual) & isnan(expected)) | isapprox(actual, expected; rtol = 1e-10, atol = 1e-12))
        end

        # A checkpoint is only valid when every selected streamed output is
        # still present. If one .part file disappeared, discard the checkpoint
        # and recompute every band instead of creating a sparse/empty cube.
        rm(joinpath(root, Moose.COMPLETION_MANIFEST))
        for rel in cube_relpaths
            path = joinpath(root, rel)
            mv(path, path * ".part")
        end
        rm(joinpath(result_dir, "Qnu.fits.part"))
        Moose._write_tile_checkpoint(joinpath(root, Moose.TILE_CHECKPOINT), signature, 1, 3, maps2d)

        @test_logs match_mode = :any (:warn, r"missing partial output files")
            Moose.MOOSE_from_config(demo.config_path; quiet = true)

        @test !isfile(joinpath(root, Moose.TILE_CHECKPOINT))
        @test isempty(filter(path -> endswith(path, ".part"),
            [joinpath(dirpath, name) for (dirpath, _, names) in walkdir(root) for name in names]))
        for (rel, expected) in expected_cubes
            actual = read(FITS(joinpath(root, rel))[1])
            @test all(@. (isnan(actual) & isnan(expected)) | isapprox(actual, expected; rtol = 1e-10, atol = 1e-12))
        end
    end

    # Selective output groups are honored in tiled cartesian mode.
    mktempdir() do dir
        demo = make_demo_data(dir; npix = 5)
        cfg = JSON.parsefile(demo.config_path)
        cfg["tile_size"] = 2
        cfg["outputs"] = ["rm"]
        write(demo.config_path, JSON.json(cfg))

        Moose.MOOSE_from_config(demo.config_path; quiet = true)

        root = joinpath(demo.simulation_dir, "z", "Synchrotron")
        @test Moose._output_files(root) == [joinpath("WithFaraday", "RMmap.fits")]
        rmmap = read(FITS(joinpath(root, "WithFaraday", "RMmap.fits"))[1])
        @test all(isapprox.(rmmap, demo.expected.rm; rtol = 1e-10))
    end

    # float32 + tiling combine.
    mktempdir() do dir
        demo = make_demo_data(dir; npix = 4)
        cfg = JSON.parsefile(demo.config_path)
        cfg["tile_size"] = 3
        cfg["precision"] = "float32"
        write(demo.config_path, JSON.json(cfg))
        Moose.MOOSE_from_config(demo.config_path; quiet = true)

        result_dir = joinpath(demo.simulation_dir, "z", "Synchrotron", "WithFaraday")
        tnu = read(FITS(joinpath(result_dir, "Tnu.fits"))[1])
        @test eltype(tnu) == Float32
        for k in axes(tnu, 3)
            @test all(isapprox.(tnu[:, :, k], demo.expected.Tnu[k]; rtol = 1e-5))
        end
    end

    # Config-time incompatibilities.
    mktempdir() do dir
        demo = make_demo_data(dir; npix = 4)
        base = JSON.parsefile(demo.config_path)

        noisy = copy(base)
        noisy["tile_size"] = 2
        noisy["add_noise"] = "Y"
        noisy["SNR_nu"] = 5.0
        write(demo.config_path, JSON.json(noisy))
        @test_throws Exception Moose.MOOSE_from_config(demo.config_path; quiet = true)

        filtered = copy(base)
        filtered["tile_size"] = 2
        filtered["responseSynchrotron"] = "Y"
        filtered["kernel_size_synchrotron"] = 2.0
        write(demo.config_path, JSON.json(filtered))
        @test_throws Exception Moose.MOOSE_from_config(demo.config_path; quiet = true)

        cleaned = copy(base)
        cleaned["tile_size"] = 2
        cleaned["rm_clean"] = Dict("enabled" => true, "gain" => 0.1, "niter" => 5, "threshold" => 0.0)
        write(demo.config_path, JSON.json(cleaned))
        @test_throws Exception Moose.MOOSE_from_config(demo.config_path; quiet = true)

        invalid = copy(base)
        invalid["tile_size"] = 0
        write(demo.config_path, JSON.json(invalid))
        @test_throws Exception Moose.MOOSE_from_config(demo.config_path; quiet = true)
    end
end

@testset "CLI reproducibility options" begin
    _, _, _, overrides = parse_cli_args(["--rng-seed", "42"])
    @test overrides["rng_seed"] == 42
end

@testset "HEALPix support" begin
    q_maps = [
        Moose.healpix_map(fill(1.0, 12); order=:ring),
        Moose.healpix_map(fill(0.5, 12); order=:ring),
    ]
    u_maps = [
        Moose.healpix_map(fill(0.0, 12); order=:ring),
        Moose.healpix_map(fill(0.25, 12); order=:ring),
    ]

    q_stack = Moose.HealpixStack(q_maps)
    @test size(q_stack) == (12, 2)
    @test q_stack.nside == 1
    @test q_stack.order == :ring

    result = Moose.RMSynthesisHealpix(q_stack, u_maps, [1.0e9, 1.1e9], [-10.0, 0.0, 10.0])
    @test size(result.fdf) == (12, 3)
    @test result.nside == 1
    @test result.order == :ring
    @test result.phi == [-10.0, 0.0, 10.0]

    mktempdir() do dir
        map_path = joinpath(dir, "q.fits")
        Moose.write_healpix_map(map_path, q_stack[:, 1]; nside=q_stack.nside, order=q_stack.order)
        @test Moose.detect_fits_grid(map_path) == :healpix
        @test Moose.is_healpix_fits(map_path)
        reread = Moose.read_healpix_map(map_path)
        @test collect(reread) == q_stack[:, 1]

        paths = Moose.write_healpix_rm_result(joinpath(dir, "rm"), result; prefix="test")
        @test length(paths.fdf) == 3
        @test all(isfile, paths.fdf)
        @test length(Moose.read_healpix_map(first(paths.fdf))) == 12

        q_paths = [joinpath(dir, "q$(i).fits") for i in eachindex(q_maps)]
        u_paths = [joinpath(dir, "u$(i).fits") for i in eachindex(u_maps)]
        for i in eachindex(q_maps)
            Moose.write_healpix_map(q_paths[i], collect(q_maps[i]); nside=1, order=:ring)
            Moose.write_healpix_map(u_paths[i], collect(u_maps[i]); nside=1, order=:ring)
        end
        auto_hp = Moose.RMSynthesisAuto(q_paths, u_paths, [1.0e9, 1.1e9], [-10.0, 0.0, 10.0])
        @test auto_hp isa Moose.HealpixRMResult
        @test auto_hp.fdf ≈ result.fdf

        q_cube_path = joinpath(dir, "q_cube.fits")
        u_cube_path = joinpath(dir, "u_cube.fits")
        q_cube = reshape([1.0, 0.5, 1.0, 0.5], 2, 1, 2)
        u_cube = reshape([0.0, 0.25, 0.0, 0.25], 2, 1, 2)
        write_test_fits(q_cube_path, q_cube)
        write_test_fits(u_cube_path, u_cube)
        @test Moose.detect_fits_grid(q_cube_path) == :image
        @test Moose.is_image_fits(q_cube_path)
        auto_cube = Moose.RMSynthesisAuto(q_cube_path, u_cube_path, [1.0e9, 1.1e9], [-10.0, 0.0, 10.0])
        direct_cube = Moose.RMSynthesis(q_cube, u_cube, [1.0e9, 1.1e9], [-10.0, 0.0, 10.0])
        @test auto_cube[1] ≈ direct_cube[1]
        @test auto_cube[2] ≈ direct_cube[2]
        @test auto_cube[3] ≈ direct_cube[3]
    end
end

@testset "HEALPix ordering dispatch" begin
    @test Moose.healpix_order(Moose.healpix_map(fill(0.0, 12); order=:ring)) == :ring
    @test Moose.healpix_order(Moose.healpix_map(fill(0.0, 12); order=:nested)) == :nested
    @test Moose.normalize_healpix_order(" RING ") == :ring

    # Mixed Q/U ordering must be rejected, not silently combined.
    q_ring = [Moose.healpix_map(fill(1.0, 12); order=:ring)]
    u_nested = [Moose.healpix_map(fill(0.0, 12); order=:nested)]
    @test_throws ErrorException Moose.RMSynthesisHealpix(q_ring, u_nested, [1.0e9], [-10.0, 0.0, 10.0])
    @test_throws ErrorException Moose.RMCleanHealpix(q_ring, u_nested, [1.0e9], [-10.0, 0.0, 10.0])
end

@testset "HEALPix UNSEEN handling" begin
    mktempdir() do dir
        # NaN pixels are stored as the UNSEEN sentinel...
        map_path = joinpath(dir, "masked.fits")
        values = fill(1.0, 12)
        values[3] = NaN
        Moose.write_healpix_map(map_path, values; nside=1, order=:ring)
        raw = Moose.read_healpix_map(map_path)
        @test raw[3] <= -1.6e30

        # ...and converted back to NaN when read as a stack.
        stack = Moose.read_healpix_stack(map_path)
        @test isnan(stack.pixels[3, 1])
        @test stack.pixels[1, 1] == 1.0

        # Opt-out flags on both sides.
        raw_stack = Moose.read_healpix_stack(map_path; unseen_to_nan=false)
        @test raw_stack.pixels[3, 1] <= -1.6e30
        plain_path = joinpath(dir, "plain.fits")
        Moose.write_healpix_map(plain_path, values; nside=1, order=:ring, nan_to_unseen=false)
        @test isnan(Moose.read_healpix_map(plain_path)[3])
    end
end

@testset "HEALPix cube IO & COORDSYS" begin
    pixels = hcat(fill(1.0, 12), fill(2.0, 12))
    pixels[1, 1] = Moose.HEALPIX_UNSEEN

    mktempdir() do dir
        cube_path = joinpath(dir, "cube.fits")
        Moose.write_healpix_cube(cube_path, pixels, [10.0, 20.0];
            order=:ring, coordsys="G", unit="Jy", coordname="PHI", nan_to_unseen=false)
        @test Moose.detect_fits_grid(cube_path) == :healpix

        stack, coords = Moose.read_healpix_cube(cube_path)
        @test stack isa Moose.HealpixStack
        @test size(stack.pixels) == (12, 2)
        @test stack.nside == 1
        @test stack.order == :ring
        @test stack.coordsys == "G"
        @test coords == [10.0, 20.0]
        @test isnan(stack.pixels[1, 1])   # UNSEEN masked on read
        @test stack.pixels[2, 1] == 1.0
        @test stack.pixels[1, 2] == 2.0

        # COORDSYS survives the single-map writer too.
        map_path = joinpath(dir, "map.fits")
        Moose.write_healpix_map(map_path, fill(1.0, 12); nside=1, order=:ring, coordsys="g")
        @test Moose.read_healpix_stack(map_path).coordsys == "G"

        # format=:cube writes a single file per stack.
        clean = copy(pixels)
        clean[1, 1] = 1.0
        paths = Moose.write_healpix_stack(dir, clean, "stackcube", [1.0, 2.0];
            order=:ring, format=:cube, coordsys="C")
        @test length(paths) == 1
        st, c = Moose.read_healpix_cube(only(paths))
        @test st.pixels == clean
        @test st.coordsys == "C"
        @test c == [1.0, 2.0]

        # write_healpix_rm_result honours format=:cube and propagates COORDSYS.
        q = Moose.HealpixStack(hcat(fill(1.0, 12), fill(0.5, 12)); order=:ring, coordsys="G")
        u = Moose.HealpixStack(hcat(fill(0.0, 12), fill(0.25, 12)); order=:ring)
        result = Moose.RMSynthesisHealpix(q, u, [1.0e9, 1.1e9], [-10.0, 0.0, 10.0])
        @test result.coordsys == "G"
        rm_paths = Moose.write_healpix_rm_result(joinpath(dir, "rm_cube"), result; prefix="t", format=:cube)
        @test length(rm_paths.fdf) == 1
        fdf_stack, fdf_phi = Moose.read_healpix_cube(only(rm_paths.fdf))
        @test size(fdf_stack.pixels) == (12, 3)
        @test fdf_stack.coordsys == "G"
        @test fdf_phi == [-10.0, 0.0, 10.0]
    end
end

@testset "HEALPix reorder & udgrade" begin
    stack = Moose.HealpixStack(reshape(collect(1.0:12.0), 12, 1); order=:ring, coordsys="G")

    nested = Moose.healpix_reorder(stack, :nested)
    @test nested.order == :nested
    @test sort(vec(nested.pixels)) == collect(1.0:12.0)
    back = Moose.healpix_reorder(nested, :ring)
    @test back.pixels == stack.pixels
    @test back.coordsys == "G"
    @test Moose.healpix_reorder(stack, :ring) === stack

    # Upgrade replicates parents (sum scales by 4); degrading back is exact.
    up = Moose.healpix_udgrade(stack, 2)
    @test up.nside == 2
    @test size(up.pixels) == (48, 1)
    @test sum(up.pixels) ≈ 4 * sum(stack.pixels)
    @test up.coordsys == "G"
    down = Moose.healpix_udgrade(up, 1)
    @test down.pixels ≈ stack.pixels

    # Degrading averages children and ignores NaN; all-NaN children stay NaN.
    v = fill(1.0, 48)
    v[1] = NaN
    v[2], v[3], v[4] = 2.0, 4.0, 6.0     # children of nested parent 0
    v[5:8] .= NaN                          # children of nested parent 1
    st = Moose.HealpixStack(reshape(v, 48, 1); order=:nested)
    dg = Moose.healpix_udgrade(st, 1)
    @test dg.pixels[1, 1] ≈ 4.0
    @test isnan(dg.pixels[2, 1])
    @test all(dg.pixels[3:end, 1] .≈ 1.0)

    # Vector variant and invalid NSIDE.
    @test Moose.healpix_udgrade(fill(1.0, 12), 2; order=:ring) ≈ fill(1.0, 48)
    @test_throws ErrorException Moose.healpix_udgrade(stack, 3)
end

@testset "HEALPix smoothing" begin
    npix8 = 768  # NSIDE = 8

    # A constant map is a pure monopole: invariant under beam smoothing.
    stack = Moose.HealpixStack(fill(2.5, npix8, 1); order=:ring, coordsys="G")
    sm = Moose.healpix_smooth(stack; fwhm_deg=10.0)
    @test sm.nside == 8
    @test sm.order == :ring
    @test sm.coordsys == "G"
    @test all(isapprox.(sm.pixels, 2.5; rtol=1e-4))

    # Masked smoothing: NaN pixels stay NaN, the rest stays exactly constant
    # (numerator and denominator see the same window).
    v = fill(1.0, npix8)
    v[10] = NaN
    sm2 = Moose.healpix_smooth(v; fwhm_deg=10.0, order=:ring)
    @test isnan(sm2[10])
    valid = .!isnan.(sm2)
    @test count(valid) == npix8 - 1
    @test all(isapprox.(sm2[valid], 1.0; rtol=1e-4))

    # Nested input comes back nested.
    nstack = Moose.healpix_reorder(stack, :nested)
    smn = Moose.healpix_smooth(nstack; fwhm_deg=10.0)
    @test smn.order == :nested
    @test all(isapprox.(smn.pixels, 2.5; rtol=1e-4))

    # Exactly one FWHM keyword is required.
    @test_throws ErrorException Moose.healpix_smooth(stack)
    @test_throws ErrorException Moose.healpix_smooth(stack; fwhm_deg=1.0, fwhm_arcmin=60.0)
end

@testset "RM synthesis & RM-CLEAN masked pixels" begin
    nu = [1.0e9, 1.1e9]
    phi = [-10.0, 0.0, 10.0]

    q = fill(1.0, 4, 2)
    u = fill(0.5, 4, 2)
    q[2, :] .= NaN

    fdf, re, im_ = Moose.RMSynthesis(q, u, nu, phi)
    @test all(isnan, fdf[2, :])
    @test all(isnan, re[2, :])
    @test all(isnan, im_[2, :])
    @test all(isfinite, fdf[1, :])

    # Valid rows are identical to a fully unmasked computation.
    fdf0, _, _ = Moose.RMSynthesis(fill(1.0, 4, 2), fill(0.5, 4, 2), nu, phi)
    @test fdf[1, :] ≈ fdf0[1, :]
    @test fdf[3, :] ≈ fdf0[3, :]

    @test_throws ArgumentError Moose.RMSynthesis(q, u[:, 1:1], nu, phi)
    @test_throws ArgumentError Moose.RMSynthesis(q, u, [1.0e9, NaN], phi)
    @test_throws ArgumentError Moose.RMSynthesis(q, u, [1.0e9, 0.0], phi)
    @test_throws ArgumentError Moose.RMSynthesis(q, u, nu, [0.0, NaN])
    @test_throws ArgumentError Moose.RMSynthesis(q, u, [1.0e9], phi)

    clean = Moose.RMClean(q, u, nu, phi; gain=0.1, niter=10)
    @test all(isnan, clean.cleanFDF[2, :])
    @test all(isnan, clean.residual[2, :])
    @test all(isfinite, clean.cleanFDF[1, :])
end

@testset "HEALPix simulation NSIDE unification" begin
    mktempdir() do base_dir
        sim_dir = joinpath(base_dir, "mixed_hp")
        mkdir(sim_dir)

        npix = 12
        write_test_healpix_cube(joinpath(sim_dir, "Bx.fits"), hcat(fill(1.0, npix), fill(1.2, npix)))
        # By is provided at NSIDE=2 and must be degraded to Bx's NSIDE=1.
        write_test_healpix_cube(joinpath(sim_dir, "By.fits"), hcat(fill(0.5, 48), fill(0.4, 48)); nside=2)
        write_test_healpix_cube(joinpath(sim_dir, "Bz.fits"), hcat(fill(0.25, npix), fill(0.3, npix)))
        write_test_healpix_cube(joinpath(sim_dir, "density.fits"), hcat(fill(1.0, npix), fill(1.1, npix)))
        write_test_healpix_cube(joinpath(sim_dir, "temperature.fits"), hcat(fill(100.0, npix), fill(110.0, npix)))

        B1, B2, BLOS, T, n, nH2, nHp = Moose.ReadSimulation(sim_dir, "z", 1.0, 1.0, 1.0)
        @test size(B1) == (npix, 1, 2)
        @test size(B2) == (npix, 1, 2)   # degraded from 48 to 12 pixels
        @test all(B2[:, 1, 1] .≈ 0.5)
        @test all(B2[:, 1, 2] .≈ 0.4)

        # HEALPix lines of sight are radial: only LOS="z" is meaningful, and
        # Bx/By/Bz must be per-pixel (e_theta, e_phi, e_r) tangent components.
        @test_throws Exception Moose.ReadSimulation(sim_dir, "x", 1.0, 1.0, 1.0)
        @test_throws Exception Moose.ReadSimulation(sim_dir, "y", 1.0, 1.0, 1.0)
    end
end

@testset "Tiled HEALPix end-to-end run" begin
    mktempdir() do base_dir
        function make_hp_sim(dir)
            mkdir(dir)
            npix = 12
            write_test_healpix_cube(joinpath(dir, "Bx.fits"), hcat(fill(1.0, npix), fill(1.2, npix)))
            write_test_healpix_cube(joinpath(dir, "By.fits"), hcat(fill(0.5, npix), fill(0.4, npix)))
            write_test_healpix_cube(joinpath(dir, "Bz.fits"), hcat(fill(0.25, npix), fill(0.3, npix)))
            write_test_healpix_cube(joinpath(dir, "density.fits"), hcat(fill(1.0, npix), fill(1.1, npix)))
            write_test_healpix_cube(joinpath(dir, "temperature.fits"), hcat(fill(100.0, npix), fill(110.0, npix)))
            return dir
        end

        tiled_dir = make_hp_sim(joinpath(base_dir, "hp_tiled"))
        plain_dir = make_hp_sim(joinpath(base_dir, "hp_plain"))

        interpolation_path = joinpath(base_dir, "emissivity.csv")
        write_test_emissivity(interpolation_path)

        base_cfg = Dict(
            "base_dir" => base_dir,
            "chosen_LOS" => ["z"],
            "interpolation_file_path" => interpolation_path,
            "faraday" => Dict("enabled" => "Y", "phimin" => -2.0, "phimax" => 2.0, "dphi" => 2.0),
            "responseSynchrotron" => "N",
            "add_noise" => "N",
            "ne_option" => "2",
            "IonizationFraction" => 0.1,
            "freq" => Dict("start" => 100.0, "end" => 101.0, "step" => 1.0),
            "BoxLength_pc" => 2.0,
            "BoxLength_pix" => 2,
            "log_progress" => false,
        )

        tiled_cfg = copy(base_cfg)
        tiled_cfg["simulations"] = [tiled_dir]
        tiled_cfg["tile_size"] = 5   # 12 pixels -> bands of 5, 5, 2
        tiled_config_path = joinpath(base_dir, "moose_config_tiled.json")
        write(tiled_config_path, JSON.json(tiled_cfg))

        plain_cfg = copy(base_cfg)
        plain_cfg["simulations"] = [plain_dir]
        plain_config_path = joinpath(base_dir, "moose_config_plain.json")
        write(plain_config_path, JSON.json(plain_cfg))

        Moose.MOOSE_from_config(tiled_config_path; quiet = true)
        Moose.MOOSE_from_config(plain_config_path; quiet = true)

        tiled_root = joinpath(tiled_dir, "z", "Synchrotron")
        tiled_result = joinpath(tiled_root, "WithFaraday")
        plain_result = joinpath(plain_dir, "z", "Synchrotron", "WithFaraday")

        # Tiled 3D products are single-file HEALPix cubes.
        for name in ("ne",)
            @test isfile(joinpath(tiled_root, "$(name).fits"))
        end
        for name in ("Qnu", "Unu", "Tnu", "Pnu", "polfrac", "FDF", "realFDF", "imagFDF")
            @test isfile(joinpath(tiled_result, "$(name).fits"))
        end
        for name in ("intBtotal", "intne", "sigmaT")
            @test isfile(joinpath(tiled_root, "$(name).fits"))
        end
        for name in ("RMmap", "Pnumax", "polfracmax", "Pmax", "alpha", "alpha_err", "RMSF")
            @test isfile(joinpath(tiled_result, "$(name).fits"))
        end

        # Tiled values match the non-tiled reference run.
        q_cube, q_freqs = Moose.read_healpix_cube(joinpath(tiled_result, "Qnu.fits"))
        @test size(q_cube.pixels) == (12, 2)
        @test q_freqs ≈ [1.0e8, 1.01e8]
        q_plain_1 = Moose.read_healpix_map(joinpath(plain_result, "Qnu_0001_p1p0e8.fits"))
        q_plain_2 = Moose.read_healpix_map(joinpath(plain_result, "Qnu_0002_p1p01e8.fits"))
        @test q_cube.pixels[:, 1] ≈ collect(q_plain_1) rtol = 1e-10
        @test q_cube.pixels[:, 2] ≈ collect(q_plain_2) rtol = 1e-10

        fdf_cube, fdf_phi = Moose.read_healpix_cube(joinpath(tiled_result, "FDF.fits"))
        @test size(fdf_cube.pixels) == (12, 3)
        @test fdf_phi ≈ [-2.0, 0.0, 2.0]
        fdf_plain_2 = Moose.read_healpix_map(joinpath(plain_result, "FDF_0002_p0p0.fits"))
        @test fdf_cube.pixels[:, 2] ≈ collect(fdf_plain_2) rtol = 1e-10

        int_tiled = Moose.read_healpix_stack(joinpath(tiled_root, "intBtotal.fits"))
        int_plain = Moose.read_healpix_stack(joinpath(plain_dir, "z", "Synchrotron", "intBtotal.fits"))
        @test int_tiled.pixels ≈ int_plain.pixels rtol = 1e-10

        # Selective output groups also apply to the tiled HEALPix path.
        selective_dir = make_hp_sim(joinpath(base_dir, "hp_selective"))
        selective_cfg = copy(base_cfg)
        selective_cfg["simulations"] = [selective_dir]
        selective_cfg["tile_size"] = 5
        selective_cfg["outputs"] = ["spectral_index"]
        selective_config_path = joinpath(base_dir, "moose_config_selective.json")
        write(selective_config_path, JSON.json(selective_cfg))

        Moose.MOOSE_from_config(selective_config_path; quiet = true)

        selective_root = joinpath(selective_dir, "z", "Synchrotron")
        @test Moose._output_files(selective_root) == [
            joinpath("WithFaraday", "alpha.fits"),
            joinpath("WithFaraday", "alpha_err.fits"),
        ]
        alpha = Moose.read_healpix_stack(joinpath(selective_root, "WithFaraday", "alpha.fits"))
        @test size(alpha.pixels) == (12, 1)
    end
end

@testset "RMSF diagnostics" begin
    nu = collect(range(1.0e9, 1.5e9, length = 64))    # Hz
    phi = collect(range(-100.0, 100.0, length = 201)) # rad/m^2, dphi = 1.0
    diag = Moose.rmsf_diagnostics(nu, phi)

    # Symmetric lag grid matched to the Faraday-depth spacing, centred on zero.
    @test length(diag.phi) == 2 * length(phi) - 1
    @test diag.phi[argmin(abs.(diag.phi))] == 0.0
    @test isapprox(diag.phi[2] - diag.phi[1], 1.0; atol = 1e-9)

    # RMSF is normalised to unit peak at zero lag.
    @test isapprox(maximum(abs.(diag.rmsf)), 1.0; atol = 1e-6)
    @test isapprox(abs(diag.rmsf[argmin(abs.(diag.phi))]), 1.0; atol = 1e-6)

    # Resolution metrics are positive and finite.
    @test isfinite(diag.fwhm) && diag.fwhm > 0
    @test diag.phi_max > 0
    @test diag.max_scale > 0

    # Theoretical resolution matches 2*sqrt(3)/Δλ²; measured agrees to tens of %.
    lambda2 = (Moose.C_m ./ nu) .^ 2
    dlambda2 = maximum(lambda2) - minimum(lambda2)
    @test isapprox(diag.fwhm_theoretical, 2 * sqrt(3) / dlambda2; rtol = 1e-6)
    @test isapprox(diag.fwhm, diag.fwhm_theoretical; rtol = 0.5)
    @test isapprox(Moose.rmsf_diagnostics(reverse(nu), phi).fwhm_theoretical,
                   diag.fwhm_theoretical; rtol = 1e-12)
    @test_throws ErrorException Moose.rmsf_diagnostics(nu, [-100.0, -50.0, 10.0, 100.0])
    @test_throws ArgumentError Moose.rmsf_diagnostics([1.0e9, NaN], phi)
    @test_throws ArgumentError Moose.rmsf_diagnostics(fill(1.0e9, 2), phi)
    @test_throws ArgumentError Moose.getRMSF([1.0e9, 0.0], phi)

    mktempdir() do dir
        path = Moose.write_rmsf(dir, diag)
        @test isfile(path)
        data, hdr = FITS(path) do f
            (read(f[1]), read_header(f[1]))
        end
        @test size(data) == (length(diag.phi), 3)
        @test isapprox(hdr["RMSFFWHM"], diag.fwhm; rtol = 1e-6)
        @test isapprox(hdr["PHIMAX"], diag.phi_max; rtol = 1e-6)
    end
end

@testset "RM-CLEAN" begin
    nu = collect(range(1.0e9, 1.5e9, length = 64))     # Hz
    phi = collect(range(-100.0, 100.0, length = 201))  # rad/m^2, dphi = 1.0
    phi0 = 15.0
    lambda2 = (Moose.C_m ./ nu) .^ 2
    P = @. 1.0 * exp(2im * (0.3 + phi0 * lambda2))     # one Faraday-thin source
    Q = real.(P)
    U = imag.(P)

    absF, realF, imagF = Moose.RMSynthesis(Q, U, nu, phi)
    diag = Moose.rmsf_diagnostics(nu, phi)
    result = Moose.RMClean(Q, U, nu, phi; gain = 0.1, niter = 2000,
                           threshold = 1e-3, diagnostics = diag)

    # Output shapes match the (1D) input FDF, and the RMSF is carried through.
    @test size(result.cleanFDF) == size(absF)
    @test length(result.phi) == length(phi)
    @test result.rmsf === diag

    # The restored peak and the dominant clean component land on the source.
    @test isapprox(phi[argmax(result.cleanFDF)], phi0; atol = 1.5)
    @test isapprox(phi[argmax(abs.(result.model))], phi0; atol = 1.5)

    # Cleaning removes the RMSF sidelobe power.
    dirty_energy = sum(abs2, complex.(realF, imagF))
    residual_energy = sum(abs2, result.residual)
    @test residual_energy < 0.25 * dirty_energy
end

@testset "QU fitting" begin
    nu = collect(range(1.0e9, 2.0e9, length = 64))  # Hz
    lambda2 = (Moose.C_m ./ nu) .^ 2

    # --- External Faraday screen: exact recovery without noise.
    p0, chi0, rm = 0.7, 0.3, 15.0
    P = @. p0 * exp(2im * (chi0 + rm * lambda2))
    fit = Moose.QUFit(real.(P), imag.(P), nu; model = :screen)
    @test fit.converged
    @test isapprox(fit.params[1], p0; atol = 1e-6)
    @test isapprox(fit.params[2], chi0; atol = 1e-6)
    @test isapprox(fit.params[3], rm; atol = 1e-6)
    @test fit.chi2 < 1e-10
    @test fit.dof == 2 * length(nu) - 3

    # --- Noisy external dispersion: parameters recovered within tolerance
    # and chi2_red ≈ 1 for correctly declared uncertainties.
    rng = Moose.Random.MersenneTwister(42)
    sigma_rm, noise = 5.0, 0.005
    Pd = @. p0 * exp(-2 * sigma_rm^2 * lambda2^2) * exp(2im * (chi0 + rm * lambda2))
    qn = real.(Pd) .+ noise .* randn(rng, length(nu))
    un = imag.(Pd) .+ noise .* randn(rng, length(nu))
    fitd = Moose.QUFit(qn, un, nu; model = :external_dispersion,
                       sigma_q = noise, sigma_u = noise)
    @test isapprox(fitd.params[1], p0; atol = 0.05)
    @test isapprox(fitd.params[2], chi0; atol = 0.05)
    @test isapprox(fitd.params[3], rm; atol = 0.5)
    @test isapprox(fitd.params[4], sigma_rm; atol = 0.5)
    @test 0.5 < fitd.chi2_red < 2.0
    @test all(isfinite, fitd.stderr)

    # --- Burn slab: the sinc-like depolarization is recovered.
    phi_slab = 40.0
    x = phi_slab .* lambda2
    Ps = @. p0 * (sin(x) / x) * exp(2im * (chi0 + x / 2))
    fits = Moose.QUFit(real.(Ps), imag.(Ps), nu; model = :burn_slab)
    @test isapprox(fits.params[1], p0; atol = 1e-4)
    @test isapprox(fits.params[3], phi_slab; atol = 1e-3)

    # --- Model comparison: BIC prefers the true (simpler) screen model.
    qn2 = real.(P) .+ noise .* randn(rng, length(nu))
    un2 = imag.(P) .+ noise .* randn(rng, length(nu))
    best, results = Moose.QUFitCompare(qn2, un2, nu;
                                       sigma_q = noise, sigma_u = noise)
    @test best.model == :screen
    @test length(results) == length(Moose.QU_FIT_MODELS)
    @test all(r -> r.bic >= best.bic, values(results))

    # --- Cube fitting: per-pixel maps, NaN pixels stay masked.
    nx, ny = 2, 2
    Qc = Array{Float64}(undef, nx, ny, length(nu))
    Uc = similar(Qc)
    rms = [10.0 -20.0; 5.0 0.0]
    for j in 1:ny, i in 1:nx
        Pij = @. p0 * exp(2im * (chi0 + rms[i, j] * lambda2))
        Qc[i, j, :] .= real.(Pij)
        Uc[i, j, :] .= imag.(Pij)
    end
    Qc[2, 2, :] .= NaN  # masked pixel (e.g. HEALPix UNSEEN)
    Uc[2, 2, :] .= NaN
    params, perr, chi2map = Moose.QUFitCube(Qc, Uc, nu; model = :screen)
    @test size(params) == (nx, ny, 3)
    @test isapprox(params[1, 1, 3], 10.0; atol = 1e-4)
    @test isapprox(params[1, 2, 3], -20.0; atol = 1e-4)
    @test isapprox(params[2, 1, 3], 5.0; atol = 1e-4)
    @test all(isnan, params[2, 2, :])
    @test isnan(chi2map[2, 2])
    @test all(isfinite, chi2map[1:2, 1] )

    # --- Validation errors.
    @test_throws ArgumentError Moose.qu_model(:nope, [1.0], lambda2)
    @test_throws ArgumentError Moose.QUFit(real.(P), imag.(P)[1:10], nu)
    @test_throws ArgumentError Moose.QUFit(real.(P), imag.(P), nu; model = :nope)
end

@testset "Polarization gradient map" begin
    nx, ny = 16, 16
    Q = [3.0 * i for i in 1:nx, j in 1:ny]
    U = [4.0 * j for i in 1:nx, j in 1:ny]

    # Linear ramps: |∇P| = √(3² + 4²) = 5 exactly, including the edges
    # (one-sided differences are exact on a linear field).
    grad = Moose.polarization_gradient_map(Q, U)
    @test all(isapprox.(grad, 5.0; atol = 1e-12))

    # Physical pixel size scales the gradient.
    grad_h = Moose.polarization_gradient_map(Q, U; pixel_size = 0.5)
    @test all(isapprox.(grad_h, 10.0; atol = 1e-12))

    # |∇P| is invariant under addition of a uniform polarized screen
    # (the key property from Gaensler et al. 2011).
    grad_shift = Moose.polarization_gradient_map(Q .+ 7.0, U .- 11.0)
    @test isapprox(grad_shift, grad; atol = 1e-12)

    # Masked pixels stay masked; neighbours fall back to one-sided
    # differences and keep a finite value.
    Qm = copy(Q); Um = copy(U)
    Qm[8, 8] = NaN
    gm = Moose.polarization_gradient_map(Qm, Um)
    @test isnan(gm[8, 8])
    @test isfinite(gm[7, 8]) && isfinite(gm[9, 8]) && isfinite(gm[8, 7])
    @test count(isnan, gm) == 1

    # Normalized variant: |∇P| / |P|.
    gn = Moose.polarization_gradient_map(Q, U; normalized = true)
    @test isapprox(gn[5, 5], 5.0 / hypot(Q[5, 5], U[5, 5]); atol = 1e-12)

    # Cube input: channel-by-channel, same conventions as RMSynthesis.
    Qc = cat(Q, 2 .* Q; dims = 3)
    Uc = cat(U, 2 .* U; dims = 3)
    gc = Moose.polarization_gradient_map(Qc, Uc)
    @test size(gc) == (nx, ny, 2)
    @test all(isapprox.(gc[:, :, 2], 10.0; atol = 1e-12))

    @test_throws ArgumentError Moose.polarization_gradient_map(Q, U[:, 1:8])
    @test_throws ArgumentError Moose.polarization_gradient_map(Q, U; pixel_size = 0.0)
end

@testset "Structure function" begin
    rng = Moose.Random.MersenneTwister(7)
    nx, ny = 64, 64

    # White noise with variance σ²: SF(r) = 2σ² at every separation.
    sigma = 1.5
    noise = sigma .* randn(rng, nx, ny)
    sf = Moose.structure_function(noise; npairs = 400_000, rng = rng)
    good = sf.counts .> 2000
    @test any(good)
    @test all(abs.(sf.sf[good] .- 2 * sigma^2) .< 0.4 * sigma^2)
    @test sum(sf.counts) <= sf.npairs
    @test length(sf.separation) == length(sf.sf) == length(sf.counts) == 20

    # Smooth ramp: SF grows with separation (∝ r² for a linear field).
    ramp = [Float64(i) for i in 1:nx, j in 1:ny]
    sfr = Moose.structure_function(ramp; max_sep = 16.0, npairs = 400_000, rng = rng)
    populated = findall(sfr.counts .> 2000)
    @test length(populated) >= 2
    k1, k2 = populated[1], populated[end]
    @test sfr.separation[k2] > 2 * sfr.separation[k1]
    @test sfr.sf[k2] > 2 * sfr.sf[k1]

    # Polarization angles: ±(π/2 − δ) are nearly parallel, so the wrapped
    # SF must be tiny while the naive SF is of order π².
    delta = 0.05
    psi = [rand(rng, Bool) ? (pi / 2 - delta) : (-pi / 2 + delta) for i in 1:nx, j in 1:ny]
    sfw = Moose.structure_function(psi; angle = true, npairs = 200_000, rng = rng)
    sfn = Moose.structure_function(psi; angle = false, npairs = 200_000, rng = rng)
    goodw = sfw.counts .> 2000
    @test any(goodw)
    @test all(sfw.sf[goodw] .< (2 * delta)^2 + 1e-6)
    @test maximum(filter(isfinite, sfn.sf)) > 1.0

    # Masked pixels are ignored; a fully masked map raises.
    masked = copy(noise)
    masked[:, 1:32] .= NaN
    sfm = Moose.structure_function(masked; npairs = 100_000, rng = rng)
    @test sum(sfm.counts) > 0
    @test_throws ArgumentError Moose.structure_function(fill(NaN, 8, 8))
    @test_throws ArgumentError Moose.structure_function(noise; nbins = 0)
    @test_throws ArgumentError Moose.structure_function(noise; min_sep = 10.0, max_sep = 5.0)
end

@testset "RM-CLEAN cube and HEALPix" begin
    nu = collect(range(1.0e9, 1.5e9, length = 32))   # Hz
    phi = collect(range(-50.0, 50.0, length = 101))  # rad/m^2
    lambda2 = (Moose.C_m ./ nu) .^ 2

    # 3D cube: two spatial pixels with different Faraday depths.
    Q = Array{Float64}(undef, 2, 1, length(nu))
    U = Array{Float64}(undef, 2, 1, length(nu))
    for (j, phi0) in enumerate((-10.0, 12.0))
        P = @. exp(2im * (0.1 + phi0 * lambda2))
        Q[j, 1, :] .= real.(P)
        U[j, 1, :] .= imag.(P)
    end

    result = Moose.RMClean(Q, U, nu, phi; gain = 0.1, niter = 1500, threshold = 1e-3)
    @test size(result.cleanFDF) == (2, 1, length(phi))
    @test isapprox(phi[argmax(result.cleanFDF[1, 1, :])], -10.0; atol = 2.0)
    @test isapprox(phi[argmax(result.cleanFDF[2, 1, :])], 12.0; atol = 2.0)

    # HEALPix wrapper returns a writable HealpixRMResult of the restored FDF.
    q_maps = [Moose.healpix_map(fill(0.3, 12); order = :ring) for _ in 1:length(nu)]
    u_maps = [Moose.healpix_map(fill(0.1, 12); order = :ring) for _ in 1:length(nu)]
    hp = Moose.RMCleanHealpix(q_maps, u_maps, nu, phi; gain = 0.1, niter = 500, threshold = 1e-3)
    @test hp isa Moose.HealpixRMResult
    @test size(hp.fdf) == (12, length(phi))
    @test hp.nside == 1
    @test hp.order == :ring
    @test hp.phi == phi
    hp_auto = Moose.RMCleanAuto(q_maps, u_maps, nu, phi; gain = 0.1, niter = 500, threshold = 1e-3)
    @test hp_auto isa Moose.HealpixRMResult
    @test hp_auto.fdf ≈ hp.fdf

    mktempdir() do dir
        paths = Moose.write_healpix_rm_result(joinpath(dir, "clean"), hp; prefix = "clean")
        @test length(paths.fdf) == length(phi)
        @test all(isfile, paths.fdf)
    end
end
