"""
    run_moose()

Interactive tool to process synchrotron data for a set of simulations.

# Description
`MOOSE` is an interactive Julia-based tool designed to guide the user through processing synchrotron data for a set of simulations. The function allows the user to configure various parameters, including:
- The base directory of simulations.
- Specific simulations and lines of sight to process.
- Synchrotron data processing options, such as Faraday rotation and filtering.
- Conversion factors for physical units.

It processes the chosen simulations and lines of sight, performing calculations and saving results for later analysis.

# Pre-requisites
- The function assumes the presence of certain dependencies such as `CSV` and `DataFrames` packages.
- Ensure the interpolation file path points to a valid file, typically located in the Synchrotron data folder.

# Interactive Prompts
The function interacts with the user via the terminal to gather the following information:
### 1. Base Directory
- Prompt: `Enter the base directory for simulations`
- Default: Current working directory (`pwd()`)

### 2. Simulation Selection
- Option to process all simulations or select specific ones by their indices.
- Prompt: `Do you want to process all simulations or choose specific ones?`
- Input format for specific simulations: Comma-separated indices, e.g., `1,3,5`.

### 3. Unit Conversions
- For Magnetic Field (`B`): Default unit is μG. Provide a conversion factor if in a different unit.
- For Number Density (`n`): Default unit is cm^-3. Provide a conversion factor if in a different unit.
- For Temperature (`T`): Default unit is K. Provide a conversion factor if in a different unit.

### 4. Synchrotron Data Processing
- Option to include Faraday rotation.
  - Prompt: `Do you want to include Faraday rotation in the computation of Q and U?`
- Option to deconvolve the Faraday dispersion function with RM-CLEAN when
  Faraday rotation is enabled.
  - Prompt: `Do you want to run RM-CLEAN on the Faraday dispersion function?`
- Option to perform filtering.
  - Prompt: `Do you want to perform filtering for Synchrotron data?`
  - Specify the largest Fourier scale retained by the instrumental 0/1 mask.

### 5. Lines of Sight
- Option to process all lines of sight (`x`, `y`, `z`) or choose specific ones.
- Input format for specific lines: Comma-separated values, e.g., `x,y`.

### 6. Interpolation File Path
- Prompt: `Enter the path to the interpolation file`
- Default: `Synchrotron/emissivity.dat` in the user's home directory.

# Returns
The function outputs the following computed data:
- **RM Map**: Rotation Measure map.
- **I, Q, U**: Stokes parameters for synchrotron data.
- **FDF**: Faraday Dispersion Function.

# Example
```julia
# To run the interactive workflow, simply call:
run_moose()

# Sample interaction
Enter the base directory for simulations: /path/to/simulations
Do you want to process all simulations or choose specific ones? (Enter 'all' or 'choose'): choose
Enter the indices of the simulations you want to process, separated by commas (e.g., 1,3,5): 1,2
Do you want to include Faraday rotation in the computation of Q and U? (Y/N): Y
Enter the lines of sight you want to process, separated by commas (e.g., x,y): x,z
```
"""


function print_logo()
   cols = displaysize(stdout)[2]
   rainbow = [:red, :yellow, :green, :cyan, :blue, :magenta]

   rainbow_text_line(s::AbstractString) = begin
       colors = Iterators.cycle(rainbow)
       join([Crayon(foreground=color, bold=true)(string(c)) for (c, color) in zip(collect(s), colors)], "")
   end

   if cols < 60
       colors = Iterators.cycle(rainbow)
       for c in collect("MOOSE")
           color = iterate(colors)[1]
           print(Crayon(foreground=color, bold=true)(string(c)))
           sleep(0.1)
       end
       println("│")
       println(Crayon(foreground = :light_green, bold = true)("Synchrotron Data Tool -- dev. by Jack Berat"))
   else
       println("\n" ^ 2)
       logo_text = raw"""
        ____    ____   ___      ___     ______   ________ 
       |_   \  /   _|.'   `.  .'   `. .' ____ \ |_   __  |
         |   \/   | /  .-.  \/  .-.  \| (___ \_|  | |_ \_|
         | |\  /| | | |   | || |   | | _.____`.   |  _| _ 
        _| |_\/_| |_\  `-'  /\  `-'  /| \____) | _| |__/ |
       |_____||_____|`.___.'  `.___.'  \______.'|________|
       """
       logo_lines = filter(line -> !isempty(strip(line)), split(logo_text, "\n"))
       max_len = maximum(length.(logo_lines))
       pad_left = max(0, (cols - max_len) ÷ 2)
       border = repeat("─", max_len)
       border = repeat("─", max_len)
       println(" "^pad_left * "╭" * border * "╮")
       colors = Iterators.cycle(rainbow)
       color_state = iterate(colors)
       for line in logo_lines
           print(" "^pad_left * "│")
           for c in collect(line)
               color = color_state[1]
               print(Crayon(foreground=color, bold=true)(string(c)))
               sleep(0.0001)
               color_state = iterate(colors, color_state[2])
           end
           println("│")
       end

       println(" "^pad_left * "╰" * border * "╯")

       elk_lines = [
           (:light_yellow, raw"     \\_//"),
           (:white,    raw"   __/ \" ."),
           (:white,    raw"  /__ |"),
           (:white,    raw"  || ||")
       ]
       for (color, line) in elk_lines
           for c in collect(line)
               print(Crayon(foreground=color)(string(c)))
               sleep(0.01)
           end
           println()
       end
       println(Crayon(foreground = :light_green, bold = true)("Mock Observation Of Synchrotron Emission -- dev. by Jack Berat"))
       println(Crayon(foreground = :light_red, bold = true)("Version 1.0"))
   end
end

function load_previous_config(config_path="moose_config.json")
    if isfile(config_path)
        return JSON.parsefile(config_path)
    else
        println("[Info] No existing config found at $(config_path). Starting with defaults.")
        return Dict{String, Any}()
    end
end

function save_config(config::AbstractDict, config_path="moose_config.json")
   atomic_write_text(config_path, JSON.json(config))
end

function format_duration(elapsed)
    total_ms = Dates.value(elapsed)
    hours = total_ms ÷ 3_600_000
    minutes = (total_ms % 3_600_000) ÷ 60_000
    seconds = (total_ms % 60_000) ÷ 1_000
    milliseconds = total_ms % 1_000
    padded_hours = lpad(hours, 2, "0")
    padded_minutes = lpad(minutes, 2, "0")
    padded_seconds = lpad(seconds, 2, "0")
    padded_ms = lpad(milliseconds, 3, "0")
    return string(padded_hours, ":", padded_minutes, ":", padded_seconds, ".", padded_ms)
end

function write_summary_log(base_dir, chosen_simu, chosen_LOS, elapsed; config_path="moose_config.json", config_source_path=nothing, config_saved_path=config_path, faraday="N", rm_clean=false, responseSynchrotron="N", add_noise="N", interpolation_file_path=nothing, conversionB=nothing, conversionn=nothing, conversionT=nothing, ne_option=nothing, rng_seed=nothing, config_hash=nothing, density_kind=nothing, mean_molecular_weight=nothing, hydrogen_mass_g=nothing)
    log_path = joinpath(base_dir, "MOOSE_summary.log")
    open(log_path, "a") do io
        println(io)
        println(io, "MOOSE Summary Log")
        println(io, "=================")
        println(io, "Run completed at: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(io, "Simulations processed:")
        for simu in chosen_simu
            println(io, simu)
        end
        println(io, "Lines of sight: $(join(chosen_LOS, ", "))")
        println(io, "Output directory: $base_dir")
        println(io, "Total execution time: $(format_duration(elapsed))")
        println(io, "Config read: $(config_source_path === nothing ? "<none>" : config_source_path)")
        println(io, "Config effective: $(config_path)")
        println(io, "Config saved: $(config_saved_path === nothing ? "<not written>" : config_saved_path)")
        config_hash !== nothing && println(io, "Config hash: $(config_hash)")
        conversionB !== nothing && println(io, "Conversion B (to μG): $(conversionB)")
        if conversionn !== nothing
            density_unit = density_kind == "mass_density" ? "rho to g cm^-3" : "n to cm^-3"
            println(io, "Conversion density ($(density_unit)): $(conversionn)")
        end
        conversionT !== nothing && println(io, "Conversion T (to K): $(conversionT)")
        density_kind !== nothing && println(io, "Density kind: $(density_kind)")
        mean_molecular_weight !== nothing && println(io, "Mean molecular weight: $(mean_molecular_weight)")
        hydrogen_mass_g !== nothing && println(io, "Hydrogen mass (g): $(hydrogen_mass_g)")
        println(io, "Faraday rotation: $(faraday)")
        println(io, "RM-CLEAN: $(rm_clean)")
        println(io, "Synchrotron filtering: $(responseSynchrotron)")
        println(io, "Noise added: $(add_noise)")
        interpolation_file_path !== nothing && println(io, "Interpolation file: $(interpolation_file_path)")
        ne_option !== nothing && println(io, "Electron density option: $(ne_option)")
        rng_seed !== nothing && println(io, "Random seed: $(rng_seed)")
    end
end

const LOSFloatTuple = NamedTuple{(:x, :y, :z), Tuple{Float64, Float64, Float64}}
const LOSIntTuple = NamedTuple{(:x, :y, :z), Tuple{Int, Int, Int}}

function _get_axis_value(values::AbstractDict, names, default)
    for name in names
        string_name = String(name)
        try
            haskey(values, string_name) && return values[string_name]
        catch
        end

        try
            haskey(values, name) && return values[name]
        catch
        end
    end

    return default
end

function normalize_los_float_values(value; default::Real = 50.0, fallback_keys = (:size_pc,))
    if value isa NamedTuple
        fallback = hasproperty(value, fallback_keys[1]) ? getproperty(value, fallback_keys[1]) : default
        return (;
            x = Float64(hasproperty(value, :x) ? getproperty(value, :x) : fallback),
            y = Float64(hasproperty(value, :y) ? getproperty(value, :y) : fallback),
            z = Float64(hasproperty(value, :z) ? getproperty(value, :z) : fallback),
        )
    elseif value isa AbstractDict
        fallback = _get_axis_value(value, fallback_keys, default)
        return (;
            x = Float64(_get_axis_value(value, (:x, :X), fallback)),
            y = Float64(_get_axis_value(value, (:y, :Y), fallback)),
            z = Float64(_get_axis_value(value, (:z, :Z), fallback)),
        )
    elseif value isa AbstractVector
        length(value) == 3 || throw_config_error("LOS length arrays must have three elements (x, y, z).")
        return (; x = Float64(value[1]), y = Float64(value[2]), z = Float64(value[3]))
    else
        scalar = Float64(value)
        return (; x = scalar, y = scalar, z = scalar)
    end
end

function normalize_los_int_values(value; default::Integer = 256, fallback_keys = (:npix,))
    if value isa NamedTuple
        fallback = hasproperty(value, fallback_keys[1]) ? getproperty(value, fallback_keys[1]) : default
        return (;
            x = Int(hasproperty(value, :x) ? getproperty(value, :x) : fallback),
            y = Int(hasproperty(value, :y) ? getproperty(value, :y) : fallback),
            z = Int(hasproperty(value, :z) ? getproperty(value, :z) : fallback),
        )
    elseif value isa AbstractDict
        fallback = _get_axis_value(value, fallback_keys, default)
        return (;
            x = Int(_get_axis_value(value, (:x, :X), fallback)),
            y = Int(_get_axis_value(value, (:y, :Y), fallback)),
            z = Int(_get_axis_value(value, (:z, :Z), fallback)),
        )
    elseif value isa AbstractVector
        length(value) == 3 || throw_config_error("LOS pixel arrays must have three elements (x, y, z).")
        return (; x = Int(value[1]), y = Int(value[2]), z = Int(value[3]))
    else
        scalar = Int(value)
        return (; x = scalar, y = scalar, z = scalar)
    end
end

function los_axis_value(values::NamedTuple, los::AbstractString)
    axis = Symbol(lowercase(String(los)))
    axis in (:x, :y, :z) || throw_config_error("Invalid line of sight: $(los). Allowed values are x, y, or z."; code=:invalid_los)
    return getproperty(values, axis)
end

function los_cube_shape(values::NamedTuple, los::AbstractString)
    los_norm = lowercase(String(los))
    if los_norm == "z"
        return (values.x, values.y, values.z)
    elseif los_norm == "x"
        return (values.y, values.z, values.x)
    elseif los_norm == "y"
        return (values.z, values.x, values.y)
    end

    throw_config_error("Invalid line of sight: $(los). Allowed values are x, y, or z."; code=:invalid_los)
end

function axis_values_for_json(values::NamedTuple)
    values.x == values.y == values.z && return values.x
    return Dict("x" => values.x, "y" => values.y, "z" => values.z)
end

function validate_los_float_values(values::NamedTuple, field_name)
    for axis in (:x, :y, :z)
        value = Float64(getproperty(values, axis))
        isfinite(value) || error("`$(field_name).$(axis)` must be finite. Got: $(value)")
        value > 0 || error("`$(field_name).$(axis)` must be > 0. Got: $(value)")
    end

    return values
end

function validate_los_int_values(values::NamedTuple, field_name)
    for axis in (:x, :y, :z)
        value = Int(getproperty(values, axis))
        value > 0 || error("`$(field_name).$(axis)` must be > 0. Got: $(value)")
    end

    return values
end

normalize_rng_seed(value) =
    value === nothing ? nothing : Int(value)

function normalize_field_sources(value)
    value === nothing && return Dict{String, Any}()
    value isa AbstractDict || throw_config_error("`field_sources` must be an object mapping canonical field names to paths or HDF5 datasets."; code=:invalid_field_sources)

    allowed = Set((REQUIRED_SIMULATION_FIELDS..., OPTIONAL_SIMULATION_FIELDS..., "amr"))
    field_sources = Dict{String, Any}()
    for (key, spec) in value
        field = String(key)
        field in allowed || throw_config_error(
            "`field_sources` contains unsupported field `$(field)`. Allowed fields: $(join(sort!(collect(allowed)), ", ")).";
            code=:invalid_field_sources,
        )
        field_sources[field] = spec
    end

    return field_sources
end

function _mask_config_get(mask_cfg::AbstractDict, names)
    for name in names
        for key in (String(name), Symbol(name))
            try
                haskey(mask_cfg, key) && return mask_cfg[key]
            catch
            end
        end
    end
    return nothing
end

function _normalize_mask_threshold(mask_cfg::AbstractDict, out::Dict{String, Any}, canonical::String, names)
    value = _mask_config_get(mask_cfg, names)
    value === nothing && return nothing
    parsed = Float64(value)
    isfinite(parsed) || throw_config_error("`physical_mask.$(canonical)` must be finite. Got: $(value)"; code=:invalid_physical_mask)
    out[canonical] = parsed
    return nothing
end

function normalize_physical_mask(value)
    value === nothing && return Dict{String, Any}()
    value === false && return Dict{String, Any}()
    value isa AbstractDict || throw_config_error("`physical_mask` must be an object with optional T_min/T_max/n_min/n_max thresholds."; code=:invalid_physical_mask)

    enabled = _mask_config_get(value, ("enabled",))
    enabled === false && return Dict{String, Any}()

    mask = Dict{String, Any}()
    _normalize_mask_threshold(value, mask, "T_min", ("T_min", "temperature_min", "Tmin"))
    _normalize_mask_threshold(value, mask, "T_max", ("T_max", "temperature_max", "Tmax"))
    _normalize_mask_threshold(value, mask, "n_min", ("n_min", "density_min", "nH_min", "nmin"))
    _normalize_mask_threshold(value, mask, "n_max", ("n_max", "density_max", "nH_max", "nmax"))

    if get(mask, "T_min", -Inf) > get(mask, "T_max", Inf)
        throw_config_error("`physical_mask.T_min` must be <= `physical_mask.T_max`."; code=:invalid_physical_mask)
    end
    if get(mask, "n_min", -Inf) > get(mask, "n_max", Inf)
        throw_config_error("`physical_mask.n_min` must be <= `physical_mask.n_max`."; code=:invalid_physical_mask)
    end

    return mask
end

function normalize_density_kind(value)
    kind = lowercase(strip(String(value)))
    kind in ("number_density", "mass_density") ||
        throw_config_error("`density_kind` must be \"number_density\" or \"mass_density\". Got: $(value)"; code=:invalid_density_kind)
    return kind
end

function normalize_positive_config_float(value, field_name::AbstractString)
    parsed = Float64(value)
    isfinite(parsed) || throw_config_error("`$(field_name)` must be finite. Got: $(value)"; code=:invalid_density_kind)
    parsed > 0 || throw_config_error("`$(field_name)` must be > 0. Got: $(value)"; code=:invalid_density_kind)
    return parsed
end

function build_density_parameters(cfg)
    density_cfg = get(cfg, "density", nothing)
    kind = density_cfg isa AbstractDict ? get(density_cfg, "kind", get(cfg, "density_kind", "number_density")) : get(cfg, "density_kind", "number_density")
    mu = density_cfg isa AbstractDict ? get(density_cfg, "mean_molecular_weight", get(density_cfg, "mu", get(cfg, "mean_molecular_weight", get(cfg, "mu", 1.0)))) : get(cfg, "mean_molecular_weight", get(cfg, "mu", 1.0))
    mH = density_cfg isa AbstractDict ? get(density_cfg, "hydrogen_mass_g", get(density_cfg, "mH_g", get(cfg, "hydrogen_mass_g", M_p))) : get(cfg, "hydrogen_mass_g", M_p)

    return (
        normalize_density_kind(kind),
        normalize_positive_config_float(mu, "mean_molecular_weight"),
        normalize_positive_config_float(mH, "hydrogen_mass_g"),
    )
end

struct RunConfig
    base_dir::String
    simulations::Vector{String}
    chosen_LOS::Vector{String}
    conversionB::Float64
    conversionn::Float64
    conversionT::Float64
    faraday_rotation::String
    phimin::Float64
    phimax::Float64
    dphi::Float64
    rm_clean_enabled::Bool
    rm_clean_gain::Float64
    rm_clean_niter::Int
    rm_clean_threshold::Float64
    responseSynchrotron::String
    kernel_size_synchrotron::Union{Nothing, Float64}
    add_noise::String
    SNR_nu::Union{Nothing, Float64}
    interpolation_file_path::String
    ne_option::String
    IonizationFraction::Union{Nothing, Float64}
    wolfire_constants::Union{Nothing, NTuple{4, Float64}}
    nustart::Float64
    nuend::Float64
    dnu::Float64
    BoxLength_pc::LOSFloatTuple
    BoxLength_pix::LOSIntTuple
    config_path::String
    log_progress::Bool
    rng_seed::Union{Nothing, Int}
    precision::String
    tile_size::Union{Nothing, Int}
    field_sources::Dict{String, Any}
    physical_mask::Dict{String, Any}
    density_kind::String
    mean_molecular_weight::Float64
    hydrogen_mass_g::Float64
    resume::String
    outputs::Set{String}

    function RunConfig(
        base_dir,
        simulations,
        chosen_LOS,
        conversionB,
        conversionn,
        conversionT,
        faraday_rotation,
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
        log_progress=true,
        rng_seed=nothing,
        ;
        rm_clean_enabled=false,
        rm_clean_gain=0.1,
        rm_clean_niter=1000,
        rm_clean_threshold=0.0,
        precision="float64",
        tile_size=nothing,
        field_sources=nothing,
        physical_mask=nothing,
        density_kind="number_density",
        mean_molecular_weight=1.0,
        hydrogen_mass_g=M_p,
        resume="off",
        outputs=["all"],
    )
        precision_flag = lowercase(String(precision))
        precision_flag in ("float64", "float32") ||
            error("`precision` must be \"float64\" or \"float32\". Got: $(precision)")
        tile_size = tile_size === nothing ? nothing : Int(tile_size)
        tile_size === nothing || tile_size > 0 ||
            error("`tile_size` must be a positive integer. Got: $(tile_size)")
        resume_mode = lowercase(strip(String(resume)))
        resume_mode in ("off", "safe") ||
            error("`resume` must be \"off\" or \"safe\". Got: $(resume)")
        output_set = Set(lowercase.(String.(outputs)))
        "all" in output_set && (output_set = Set(["integrated", "stokes", "rm", "fdf", "spectral_index", "diagnostics"]))

        faraday_flag = uppercase(faraday_rotation)
        response_flag = uppercase(responseSynchrotron)
        noise_flag = uppercase(add_noise)
        sim_paths = map(String, simulations)
        los_list = map(String, chosen_LOS)

        kernel_size_synchrotron = kernel_size_synchrotron === nothing ? nothing : Float64(kernel_size_synchrotron)
        SNR_nu = SNR_nu === nothing ? nothing : Float64(SNR_nu)
        IonizationFraction = IonizationFraction === nothing ? nothing : Float64(IonizationFraction)
        wolfire_constants = wolfire_constants === nothing ? nothing : (
            Float64(wolfire_constants[1]),
            Float64(wolfire_constants[2]),
            Float64(wolfire_constants[3]),
            Float64(wolfire_constants[4]),
        )

        return new(
            String(base_dir),
            sim_paths,
            los_list,
            Float64(conversionB),
            Float64(conversionn),
            Float64(conversionT),
            faraday_flag,
            Float64(phimin),
            Float64(phimax),
            Float64(dphi),
            Bool(rm_clean_enabled),
            Float64(rm_clean_gain),
            Int(rm_clean_niter),
            Float64(rm_clean_threshold),
            response_flag,
            kernel_size_synchrotron,
            noise_flag,
            SNR_nu,
            String(interpolation_file_path),
            String(ne_option),
            IonizationFraction,
            wolfire_constants,
            Float64(nustart),
            Float64(nuend),
            Float64(dnu),
            normalize_los_float_values(BoxLength_pc),
            normalize_los_int_values(BoxLength_pix),
            String(config_path),
            Bool(log_progress),
            normalize_rng_seed(rng_seed),
            precision_flag,
            tile_size,
            normalize_field_sources(field_sources),
            normalize_physical_mask(physical_mask),
            normalize_density_kind(density_kind),
            normalize_positive_config_float(mean_molecular_weight, "mean_molecular_weight"),
            normalize_positive_config_float(hydrogen_mass_g, "hydrogen_mass_g"),
            resume_mode,
            output_set,
        )
    end
end

function config_dict_from_struct(cfg::RunConfig)
    ne_config = Dict{String, Any}("mode" => cfg.ne_option)
    cfg.IonizationFraction !== nothing && (ne_config["ion_fraction"] = cfg.IonizationFraction)
    if cfg.wolfire_constants !== nothing
        zeta, Geff, phiPAH, XC = cfg.wolfire_constants
        ne_config["zeta"] = zeta
        ne_config["Geff"] = Geff
        ne_config["phiPAH"] = phiPAH
        ne_config["XC"] = XC
    end

    return Dict(
        "base_dir" => cfg.base_dir,
        "chosen_simu" => cfg.simulations,
        "chosen_LOS" => cfg.chosen_LOS,
        "conversionB" => cfg.conversionB,
        "conversionn" => cfg.conversionn,
        "conversionT" => cfg.conversionT,
        "FaradayRotation" => cfg.faraday_rotation,
        "phimin" => cfg.phimin,
        "phimax" => cfg.phimax,
        "dphi" => cfg.dphi,
        "rm_clean" => Dict(
            "enabled" => cfg.rm_clean_enabled,
            "gain" => cfg.rm_clean_gain,
            "niter" => cfg.rm_clean_niter,
            "threshold" => cfg.rm_clean_threshold,
        ),
        "responseSynchrotron" => cfg.responseSynchrotron,
        "kernel_size_synchrotron" => cfg.kernel_size_synchrotron,
        "add_noise" => cfg.add_noise,
        "SNR_nu" => cfg.SNR_nu,
        "interpolation_file_path" => cfg.interpolation_file_path,
        "ne_option" => cfg.ne_option,
        "IonizationFraction" => cfg.IonizationFraction,
        "ne" => ne_config,
        "nustart" => cfg.nustart,
        "nuend" => cfg.nuend,
        "dnu" => cfg.dnu,
        "BoxLength_pc" => axis_values_for_json(cfg.BoxLength_pc),
        "BoxLength_pix" => axis_values_for_json(cfg.BoxLength_pix),
        "config_path" => cfg.config_path,
        "log_progress" => cfg.log_progress,
        "rng_seed" => cfg.rng_seed,
        "precision" => cfg.precision,
        "tile_size" => cfg.tile_size,
        "field_sources" => cfg.field_sources,
        "physical_mask" => cfg.physical_mask,
        "density_kind" => cfg.density_kind,
        "mean_molecular_weight" => cfg.mean_molecular_weight,
        "hydrogen_mass_g" => cfg.hydrogen_mass_g,
        "resume" => cfg.resume,
        "outputs" => sort!(collect(cfg.outputs)),
    )
end

const COMPLETION_MANIFEST = ".moose-complete.json"

function _fingerprint_file(path::AbstractString)
    info = stat(path)
    return Dict("path" => abspath(path), "size" => info.size, "mtime" => info.mtime)
end

function _resume_input_paths(cfg::RunConfig, simu::AbstractString)
    paths = String[cfg.interpolation_file_path]
    for field in ("Bx", "By", "Bz", "density", "temperature")
        source = simulation_field_source(simu, field, cfg.field_sources)
        if source isa AbstractVector
            append!(paths, source_path.(source))
        else
            push!(paths, source_path(source))
        end
    end
    if cfg.ne_option == "3"
        source = simulation_field_source(simu, "densityHp", cfg.field_sources)
        if source isa AbstractVector
            append!(paths, source_path.(source))
        else
            push!(paths, source_path(source))
        end
    end
    amr = amr_config(cfg.field_sources)
    if amr !== nothing
        reference = simulation_field_source(simu, "Bx", cfg.field_sources)
        push!(paths, _amr_geometry_file(simu, amr, source_path(reference)))
    end
    return sort!(unique!(abspath.(paths)))
end

_los_output_root(simu::AbstractString, los::AbstractString) = joinpath(simu, los, "Synchrotron")

function _output_files(root::AbstractString)
    isdir(root) || return String[]
    files = String[]
    for (dir, _, names) in walkdir(root), name in names
        endswith(lowercase(name), ".fits") && push!(files, relpath(joinpath(dir, name), root))
    end
    return sort!(files)
end

function _write_completion_manifest(cfg::RunConfig, simu, los, config_hash)
    root = _los_output_root(simu, los)
    outputs = _output_files(root)
    isempty(outputs) && error("Cannot mark LOS=$(los) complete: no FITS outputs were produced in $(root).")
    manifest = Dict(
        "status" => "complete", "config_hash" => config_hash,
        "simulation" => abspath(simu), "los" => los,
        "moose_version" => moose_version(),
        "inputs" => [_fingerprint_file(path) for path in _resume_input_paths(cfg, simu)],
        "outputs" => outputs, "completed_at" => string(now()),
    )
    atomic_write_text(joinpath(root, COMPLETION_MANIFEST), JSON.json(manifest))
    # The completion manifest is now the durable LOS-level checkpoint. Remove
    # the finer tiled checkpoint only after that manifest exists atomically.
    rm(joinpath(root, TILE_CHECKPOINT); force = true)
end

function _can_resume_los(cfg::RunConfig, simu, los, config_hash)
    cfg.resume == "safe" || return false
    root = _los_output_root(simu, los)
    path = joinpath(root, COMPLETION_MANIFEST)
    isfile(path) || return false
    manifest = try
        JSON.parsefile(path)
    catch err
        @warn "Ignoring unreadable completion manifest" path exception=err
        return false
    end
    get(manifest, "status", "") == "complete" || return false
    get(manifest, "config_hash", "") == config_hash || return false
    get(manifest, "simulation", "") == abspath(simu) || return false
    get(manifest, "los", "") == los || return false
    current_inputs = try
        [_fingerprint_file(p) for p in _resume_input_paths(cfg, simu)]
    catch
        return false
    end
    get(manifest, "inputs", Any[]) == current_inputs || return false
    outputs = get(manifest, "outputs", Any[])
    !isempty(outputs) || return false
    return all(rel -> begin
        candidate = joinpath(root, String(rel))
        isfile(candidate) && filesize(candidate) > 0
    end, outputs)
end

function _fits_image_shape(path::AbstractString)
    return FITS(path) do fits
        for hdu in fits
            hdu isa ImageHDU && ndims(hdu) > 0 && return Tuple(Int.(size(hdu)))
        end
        error("No image HDU found in $(path).")
    end
end

function _source_shape(source, grid_kind::Symbol)
    if grid_kind == :amr
        error("AMR source shapes require the AMR geometry configuration.")
    end
    if source isa HDF5DatasetSource
        return h5open(source.file, "r") do h5
            Tuple(Int.(size(h5[source.dataset])))
        end
    elseif source isa AbstractString && is_hdf5_path(source)
        dataset = _resolve_hdf5_dataset(source, nothing)
        return h5open(source, "r") do h5
            Tuple(Int.(size(h5[dataset])))
        end
    elseif grid_kind == :healpix
        files = source isa AbstractString ? [source] : source
        info = _healpix_table_hdu_info(first(files))
        npix = Healpix.nside2npix(info.nside)
        if length(files) > 1
            return (npix, 1, length(files))
        end
        header = FITS(first(files)) do fits
            read_header(fits[info.hdu_index])
        end
        form = String(_fits_header_value(header, "TFORM1", "1D"))
        matched = match(r"^\s*(\d+)", form)
        nshell = matched === nothing ? 1 : parse(Int, matched.captures[1])
        return (npix, 1, nshell)
    elseif source isa AbstractVector
        plane = _fits_image_shape(first(source))
        length(plane) == 2 || error("Stacked FITS source must contain 2D planes; got $(plane).")
        all(_fits_image_shape(path) == plane for path in source) || error("FITS planes in a field stack have inconsistent dimensions.")
        return (plane..., length(source))
    end
    return _fits_image_shape(source)
end

_human_bytes(bytes::Real) = bytes < 1024 ? "$(round(Int, bytes)) B" :
    bytes < 1024^2 ? "$(round(bytes / 1024; digits=1)) KiB" :
    bytes < 1024^3 ? "$(round(bytes / 1024^2; digits=1)) MiB" :
    "$(round(bytes / 1024^3; digits=2)) GiB"

"""Validate input metadata and estimate the resources required by a run."""
function preflight_plan(cfg::RunConfig; io::IO=stdout)
    nfreq = length(range(start=cfg.nustart, stop=cfg.nuend, step=cfg.dnu))
    nphi = cfg.faraday_rotation == "Y" ? length(range(start=cfg.phimin, stop=cfg.phimax, step=cfg.dphi)) : 0
    bytes_per = cfg.precision == "float32" ? 4 : 8
    entries = NamedTuple[]
    total_disk = 0.0
    total_voxel_channels = 0
    want(name) = name in cfg.outputs
    need_qu = want("stokes") || want("fdf") || want("diagnostics")
    need_t = want("stokes") || want("spectral_index") || want("diagnostics")

    println(io, "MOOSE preflight plan")
    println(io, "Frequency channels: $(nfreq)" * (nphi > 0 ? " | Faraday-depth channels: $(nphi)" : ""))
    for simu in cfg.simulations
        grid = simulation_grid_kind(simu, cfg.field_sources)
        grid == :amr && cfg.tile_size !== nothing && throw_config_error(
            "`tile_size` is not supported with AMR inputs because rasterization must validate the complete leaf-cell coverage.";
            code=:unsupported_grid_operation)
        amr_plan = if grid == :amr
            reference = simulation_field_source(simu, "Bx", cfg.field_sources)
            load_amr_raster_plan(simu, amr_config(cfg.field_sources), source_path(reference))
        else
            nothing
        end
        amr_geometry = amr_plan === nothing ? nothing : amr_plan.geometry
        shapes = Dict{String, Tuple}()
        for field in ("Bx", "By", "Bz", "density", "temperature")
            source = simulation_field_source(simu, field, cfg.field_sources)
            _validate_simulation_source(source)
            if grid == :amr
                field_shape = _source_shape(source, :image)
                prod(field_shape) == size(amr_geometry.centers, 1) || throw_config_error(
                    "AMR field $(field) has $(prod(field_shape)) values but the geometry has $(size(amr_geometry.centers, 1)) cells.";
                    code=:cube_shape_mismatch)
                shapes[field] = amr_geometry.shape
            else
                shapes[field] = _source_shape(source, grid)
            end
        end
        cfg.ne_option == "3" && begin
            source = simulation_field_source(simu, "densityHp", cfg.field_sources)
            _validate_simulation_source(source)
            if grid == :amr
                field_shape = _source_shape(source, :image)
                prod(field_shape) == size(amr_geometry.centers, 1) || throw_config_error(
                    "AMR field densityHp has $(prod(field_shape)) values but the geometry has $(size(amr_geometry.centers, 1)) cells.";
                    code=:cube_shape_mismatch)
                shapes["densityHp"] = amr_geometry.shape
            else
                shapes["densityHp"] = _source_shape(source, grid)
            end
        end
        reference = shapes["Bx"]
        length(reference) == 3 || error("Simulation field Bx must be 3D; got $(reference).")
        mismatches = ["$(field)=$(shape)" for (field, shape) in shapes if shape != reference]
        isempty(mismatches) || throw_config_error("Preflight shape mismatch in $(simu): expected $(reference); $(join(mismatches, ", "))."; code=:cube_shape_mismatch)

        for los in cfg.chosen_LOS
            grid == :healpix && _validate_healpix_los(los)
            processed = grid == :healpix ? reference : _permuted_shape(reference, los)
            expected = grid == :healpix ? nothing : Tuple(los_cube_shape(cfg.BoxLength_pix, los))
            expected === nothing || processed == expected || throw_config_error(
                "Preflight shape mismatch in $(simu) after LOS=$(los): expected $(expected), got $(processed)."; code=:cube_shape_mismatch)
            nsky = processed[1] * processed[2]
            nvox = prod(processed)
            cube_elements = (want("integrated") ? nvox : 0) + (want("stokes") ? 5 * nsky * nfreq : 0) + (want("fdf") ? 3 * nsky * nphi : 0)
            map_elements = nsky * ((want("integrated") ? 8 : 0) + (want("stokes") ? 2 : 0) + (want("rm") ? 1 : 0) + (want("fdf") ? 1 : 0) + (want("spectral_index") ? 2 : 0))
            disk = (cube_elements + map_elements) * bytes_per
            working_rows = cfg.tile_size === nothing ? processed[2] : min(cfg.tile_size, processed[2])
            working_sky = processed[1] * working_rows
            ram_elements = 7 * processed[1] * working_rows * processed[3] + (need_qu ? 2 * working_sky * nfreq : 0) + (need_t ? working_sky * nfreq : 0) + (want("fdf") ? 3 * working_sky * nphi : 0)
            amr_ram = if grid == :amr
                ncell = size(amr_geometry.centers, 1)
                6 * ncell * sizeof(Float64) + nvox * sizeof(Int) + ncell * bytes_per
            else
                0
            end
            ram = ram_elements * bytes_per + amr_ram
            raster_work = grid == :amr ? 5 * nvox : 0
            work = (need_qu || need_t ? nsky * processed[3] * nfreq : 0) + (want("fdf") ? nsky * nfreq * nphi : 0) + raster_work
            total_disk += disk
            total_voxel_channels += work
            push!(entries, (; simulation=simu, los, grid, shape=processed, ram_bytes=ram, disk_bytes=disk, workload=work))
            amr_note = grid == :amr ? ", including AMR geometry/lookup ≈ $(_human_bytes(amr_ram))" : ""
            println(io, "- $(basename(simu)) LOS=$(los) grid=$(grid) shape=$(processed): peak RAM ≈ $(_human_bytes(ram))$(amr_note), FITS data ≈ $(_human_bytes(disk)), workload=$(work) cell-channel ops")
        end
    end
    println(io, "Total estimated FITS data: $(_human_bytes(total_disk))")
    println(io, "Total workload: $(total_voxel_channels) cell-channel ops")
    println(io, "Estimates exclude FITS headers, plots, allocator overhead, filtering FFT work, and RM-CLEAN iterations.")
    return (; frequency_channels=nfreq, faraday_channels=nphi, entries, disk_bytes=total_disk, workload=total_voxel_channels)
end

const DEFAULT_WOLFIRE_CONSTANTS = (2.5e-16, 1.0, 0.5, 1.4e-4)

function run_moose_processing(cfg::RunConfig; quiet::Bool = false, persisted_config::Union{Nothing, AbstractDict} = nothing, source_config_path=nothing, saved_config_path=cfg.config_path, write_config_file::Bool=true)
    if quiet
        with_logger(NullLogger()) do
            _run_moose_processing(cfg; quiet = quiet, persisted_config = persisted_config, source_config_path = source_config_path, saved_config_path = saved_config_path, write_config_file = write_config_file)
        end
    else
        _run_moose_processing(cfg; quiet = quiet, persisted_config = persisted_config, source_config_path = source_config_path, saved_config_path = saved_config_path, write_config_file = write_config_file)
    end
end

function _run_moose_processing(cfg::RunConfig; quiet::Bool = false, persisted_config::Union{Nothing, AbstractDict} = nothing, source_config_path=nothing, saved_config_path=cfg.config_path, write_config_file::Bool=true)
    nuArray = range(start = cfg.nustart, stop = cfg.nuend, step = cfg.dnu)
    PhiArray = cfg.faraday_rotation == "Y" ? range(start = cfg.phimin, stop = cfg.phimax, step = cfg.dphi) : nothing
    PixelLength_pc, PixelLength_cm, DistanceArray = los_pixel_scale(cfg.BoxLength_pc.x, cfg.BoxLength_pix.x)
    rng = cfg.rng_seed === nothing ? Random.default_rng() : Random.MersenneTwister(cfg.rng_seed)
    float_type = cfg.precision == "float32" ? Float32 : Float64

    isfile(cfg.interpolation_file_path) ||
        throw_config_error("The interpolation file $(cfg.interpolation_file_path) was not found."; code=:missing_interpolation_file)
    df = CSV.File(cfg.interpolation_file_path) |> DataFrame

    if !quiet
        print_logo()
    end

    start_time = now()
    config_to_save = persisted_config === nothing ? config_dict_from_struct(cfg) : persisted_config
    config_hash = moose_config_hash(config_to_save)
    resume_hash_config = Dict{String, Any}(String(k) => v for (k, v) in config_to_save)
    delete!(resume_hash_config, "resume")
    resume_hash = moose_config_hash(resume_hash_config)
    base_metadata = Dict{String, Any}(
        "MOOSEV" => moose_version(),
        "GITHASH" => moose_git_hash(),
        "CFGHASH" => config_hash,
        "CFGSRC" => source_config_path === nothing ? "" : String(source_config_path),
        "CFGSAVE" => saved_config_path === nothing ? "" : String(saved_config_path),
        "FARADAY" => cfg.faraday_rotation,
        "FILTER" => cfg.responseSynchrotron,
        "NOISE" => cfg.add_noise,
        "SNRNU" => cfg.SNR_nu,
        "RNGSEED" => cfg.rng_seed,
        "NSTART" => cfg.nustart,
        "NUEND" => cfg.nuend,
        "DNU" => cfg.dnu,
        "NUNIT" => "MHz input; Hz FITS",
        "PHIMIN" => cfg.phimin,
        "PHIMAX" => cfg.phimax,
        "DPHI" => cfg.dphi,
        "PHIUNIT" => "rad/m^2",
        "RMCLEAN" => cfg.rm_clean_enabled,
        "RMCGAIN" => cfg.rm_clean_gain,
        "RMCNITER" => cfg.rm_clean_niter,
        "RMCTHRES" => cfg.rm_clean_threshold,
        "CONVB" => cfg.conversionB,
        "CONVN" => cfg.conversionn,
        "CONVT" => cfg.conversionT,
        "NEOPT" => cfg.ne_option,
        "INTFILE" => cfg.interpolation_file_path,
        "PRECIS" => cfg.precision,
        "TILESIZE" => cfg.tile_size,
        "PHYMASK" => _has_physical_mask(cfg.physical_mask),
        "MSKTMIN" => get(cfg.physical_mask, "T_min", nothing),
        "MSKTMAX" => get(cfg.physical_mask, "T_max", nothing),
        "MSKNMIN" => get(cfg.physical_mask, "n_min", nothing),
        "MSKNMAX" => get(cfg.physical_mask, "n_max", nothing),
        "DENSKIND" => cfg.density_kind,
        "DENSMU" => cfg.mean_molecular_weight,
        "MHG" => cfg.hydrogen_mass_g,
    )
    if cfg.ne_option == "3"
        missing_cubes = String[]
        for simu in cfg.simulations
            source = simulation_field_source(simu, "densityHp", cfg.field_sources)
            try
                _validate_simulation_source(source)
            catch
                push!(missing_cubes, simu)
            end
        end
        !isempty(missing_cubes) && throw_config_error(
            "Electron density cube `densityHp` is missing for: $(join(missing_cubes, ", ")). Configure `field_sources.densityHp` if it is named differently.";
            code=:missing_density_cube,
        )
    end

    wolfire_constants = nothing
    ion_fraction = nothing
    if cfg.ne_option == "1"
        wolfire_constants = cfg.wolfire_constants === nothing ? DEFAULT_WOLFIRE_CONSTANTS : cfg.wolfire_constants
    elseif cfg.ne_option == "2"
        ion_fraction = cfg.IonizationFraction === nothing ? 0.01 : cfg.IonizationFraction
    end

    for (i, simu) in enumerate(cfg.simulations)
        @info "Processing simulation" simulation = simu

        checkpoint_signature = if cfg.resume == "safe" && cfg.tile_size !== nothing
            inputs = [_fingerprint_file(path) for path in _resume_input_paths(cfg, simu)]
            canonical_inputs = [[input["path"], input["size"], input["mtime"]] for input in inputs]
            bytes2hex(sha256(resume_hash * JSON.json(canonical_inputs)))
        else
            nothing
        end

        for LOS in cfg.chosen_LOS
            if _can_resume_los(cfg, simu, LOS, resume_hash)
                @info "Skipping completed line of sight" simulation=simu los=LOS manifest=joinpath(_los_output_root(simu, LOS), COMPLETION_MANIFEST)
                continue
            end
            @info "Processing line of sight" los = LOS
            box_length_pc = los_axis_value(cfg.BoxLength_pc, LOS)
            box_length_pix = los_axis_value(cfg.BoxLength_pix, LOS)
            grid_kind = simulation_grid_kind(simu, cfg.field_sources)
            expected_shape = grid_kind == :healpix ? nothing : los_cube_shape(cfg.BoxLength_pix, LOS)
            PixelLength_pc, PixelLength_cm, DistanceArray = los_pixel_scale(box_length_pc, box_length_pix)

            if cfg.ne_option == "1"
                zeta, Geff, omegaPAH, XC = wolfire_constants
                ProcessSynchrotron(simu, LOS, cfg.faraday_rotation, cfg.responseSynchrotron, df, cfg.add_noise, cfg.SNR_nu,
                    cfg.kernel_size_synchrotron, zeta, Geff, omegaPAH, XC, nuArray, PhiArray, PixelLength_pc, PixelLength_cm,
                    box_length_pc, DistanceArray, cfg.conversionn, cfg.conversionT, cfg.conversionB;
                    log_progress = cfg.log_progress, rng = rng, expected_shape = expected_shape, metadata = base_metadata,
                    rm_clean_enabled = cfg.rm_clean_enabled, rm_clean_gain = cfg.rm_clean_gain,
                    rm_clean_niter = cfg.rm_clean_niter, rm_clean_threshold = cfg.rm_clean_threshold,
                    float_type = float_type, tile_rows = cfg.tile_size, field_sources = cfg.field_sources,
                    physical_mask = cfg.physical_mask, density_kind = cfg.density_kind,
                    mean_molecular_weight = cfg.mean_molecular_weight, hydrogen_mass_g = cfg.hydrogen_mass_g,
                    outputs = cfg.outputs, checkpoint_signature = checkpoint_signature)
            elseif cfg.ne_option == "2"
                ProcessSynchrotron(simu, LOS, cfg.faraday_rotation, cfg.responseSynchrotron, df, cfg.add_noise, cfg.SNR_nu,
                    cfg.kernel_size_synchrotron, ion_fraction, nuArray, PhiArray, PixelLength_pc, PixelLength_cm,
                    box_length_pc, DistanceArray, cfg.conversionn, cfg.conversionT, cfg.conversionB;
                    log_progress = cfg.log_progress, rng = rng, expected_shape = expected_shape, metadata = base_metadata,
                    rm_clean_enabled = cfg.rm_clean_enabled, rm_clean_gain = cfg.rm_clean_gain,
                    rm_clean_niter = cfg.rm_clean_niter, rm_clean_threshold = cfg.rm_clean_threshold,
                    float_type = float_type, tile_rows = cfg.tile_size, field_sources = cfg.field_sources,
                    physical_mask = cfg.physical_mask, density_kind = cfg.density_kind,
                    mean_molecular_weight = cfg.mean_molecular_weight, hydrogen_mass_g = cfg.hydrogen_mass_g,
                    outputs = cfg.outputs, checkpoint_signature = checkpoint_signature)
            else
                ProcessSynchrotron(simu, LOS, cfg.faraday_rotation, cfg.responseSynchrotron, df, cfg.add_noise, cfg.SNR_nu,
                    cfg.kernel_size_synchrotron, nuArray, PhiArray, PixelLength_pc, PixelLength_cm,
                    box_length_pc, DistanceArray, cfg.conversionn, cfg.conversionT, cfg.conversionB;
                    log_progress = cfg.log_progress, rng = rng, expected_shape = expected_shape, metadata = base_metadata,
                    rm_clean_enabled = cfg.rm_clean_enabled, rm_clean_gain = cfg.rm_clean_gain,
                    rm_clean_niter = cfg.rm_clean_niter, rm_clean_threshold = cfg.rm_clean_threshold,
                    float_type = float_type, tile_rows = cfg.tile_size, field_sources = cfg.field_sources,
                    physical_mask = cfg.physical_mask, density_kind = cfg.density_kind,
                    mean_molecular_weight = cfg.mean_molecular_weight, hydrogen_mass_g = cfg.hydrogen_mass_g,
                    outputs = cfg.outputs, checkpoint_signature = checkpoint_signature)
            end
            _write_completion_manifest(cfg, simu, LOS, resume_hash)
        end

        if length(cfg.simulations) > 1
            @info "Finished processing all chosen LOS" simulation = simu
            print_progress(i, length(cfg.simulations); label="Completed simulations")
        end
    end

    @info "Finished processing all simulations"

    elapsed = now() - start_time
    @info "Summary" simulations = join(cfg.simulations, ", ") los = join(cfg.chosen_LOS, ", ") output_directory = cfg.base_dir elapsed = format_duration(elapsed)

    write_config_file && save_config(config_to_save, cfg.config_path)
    write_summary_log(cfg.base_dir, map(basename, cfg.simulations), cfg.chosen_LOS, elapsed; config_path=cfg.config_path,
        config_source_path=source_config_path, config_saved_path=saved_config_path, config_hash=config_hash,
        faraday=cfg.faraday_rotation, rm_clean=cfg.rm_clean_enabled, responseSynchrotron=cfg.responseSynchrotron, add_noise=cfg.add_noise,
        interpolation_file_path=cfg.interpolation_file_path, conversionB=cfg.conversionB, conversionn=cfg.conversionn,
        conversionT=cfg.conversionT, ne_option=cfg.ne_option, rng_seed=cfg.rng_seed,
        density_kind=cfg.density_kind, mean_molecular_weight=cfg.mean_molecular_weight, hydrogen_mass_g=cfg.hydrogen_mass_g)
end


function run_moose(; quiet::Bool = false, reset_config::Bool = true, help::Bool = false)
   if help
       println("""
MOOSE v$(moose_version()) — Mock Observation Of Synchrotron Emission

Usage:
 run_moose(; quiet=false, reset_config=true, help=false)

Options:
 --quiet             Disable the rainbow logo at startup.
 --reset-config      Ignore previous config and prompt user again.
 --help              Show this help message and exit.

Description:
 Author: Jack Berat
 MOOSE is an interactive tool to process synchrotron mock observation data.

 Non-interactive CLI:
  For scripted or assistive use, run `julia src/MOOSE_cli.jl --help` to see the
  full set of command-line flags (e.g., --config, --base-dir, --simu, --los,
  --rng-seed, --quiet) for configuring MOOSE without prompts.

 It computes the synchrotron Stokes parameters Q and U for a set of simulation outputs,
 optionally applying Faraday rotation, noise and an interferometric Fourier mask.

 The tool interacts with the user to define simulation directories, select simulations and
 lines of sight, and configure physical unit conversions and data processing options.

 Previous configuration is loaded from `moose_config.json`, unless reset.
Outputs are saved to the base simulation directory, with logs in `MOOSE_summary.log`.

Flow:
 ┌─────────────────────────────┐
 │ Interactive config (or JSON)│
 └────────────┬────────────────┘
              │
              ▼
   Read simulation + unit setup
              │
              ▼
   Select LOS + process Q/U/I
              │
              ▼
   Apply noise + interferometric Fourier mask
              │
              ▼
   Process RM synthesis on Q and U
              │
              ▼
   Save FITS in the simulation file + logs

""")
       return
   end
    cfg, persisted_config = run_moose_interactive(; quiet = quiet, reset_config = reset_config)
    run_moose_processing(cfg; quiet = true, persisted_config = persisted_config)
end

function run_moose_interactive(; quiet::Bool = false, reset_config::Bool = true)
    if !quiet
        print_logo()
    end

    default_config_path = joinpath(pwd(), "moose_config.json")
    prompt_section("Configuration source")
    if reset_config
        println("[Info] Previous configuration ignored (reset_config=true)")
        config = Dict{String, Any}()
        config_path = default_config_path
    else
        config_path = ask_user("Enter the path to the configuration file to load", default_config_path)
        config = load_previous_config(config_path)
    end

    prompt_section("Simulation directory")
    base_dir = ""
    simu_list = String[]
    while true
        candidate_dir = ask_user("Enter the base directory for simulations", get(config, "base_dir", pwd()))
        validation_error = ensure_directory_access(candidate_dir)
        if validation_error !== nothing
            warn_user(validation_error)
            continue
        end

        candidate_list = get_simulation_list(candidate_dir)
        if isempty(candidate_list)
            warn_user("No simulations containing FITS files were found in $(candidate_dir). Please enter another path.")
            continue
        end

        base_dir = candidate_dir
        simu_list = candidate_list
        break
    end
    config["base_dir"] = base_dir

    prompt_section("Simulation selection")
    display_simulations(simu_list)

    simu_prompt = "Enter 'all' to process all simulations or provide comma-separated indices (e.g., 1,3,5)"
    simu_choice = ask_user(simu_prompt, get(config, "simu_choice", "all"))
    chosen_simu = begin
        parsed_indices = Int[]
        while true
            if uppercase(strip(simu_choice)) == "ALL"
                empty!(parsed_indices)
                append!(parsed_indices, 1:length(simu_list))
                break
            end

            raw_indices = split(simu_choice, ",")
            empty!(parsed_indices)
            for raw_idx in raw_indices
                candidate = strip(raw_idx)
                isempty(candidate) && continue
                parsed = tryparse(Int, candidate)
                if parsed === nothing
                    warn_user("Ignoring invalid simulation index: $(candidate). Enter comma-separated integers like 1,3,5.")
                elseif parsed < 1 || parsed > length(simu_list)
                    warn_user("Simulation index $(parsed) is out of range (1-$(length(simu_list))).")
                else
                    push!(parsed_indices, parsed)
                end
            end

            if !isempty(parsed_indices)
                break
            end

            error_user("No valid simulation indices provided. Use 'all' or comma-separated integers like 1,3,5 and try again.")
            simu_choice = ask_user(simu_prompt, simu_choice)
        end

        config["simu_choice"] = simu_choice
        config["simu_indices"] = uppercase(strip(simu_choice)) == "ALL" ? "" : join(parsed_indices, ",")
        map(i -> simu_list[i], parsed_indices)
    end
    config["chosen_simu"] = chosen_simu

    prompt_section("Unit conversions (to μG, cm⁻³, K)")
    conversionB = ask_user("Enter the conversion factor for magnetic field B to μG (microGauss)", Float64(get(config, "conversionB", 1.0)))
    conversionn = ask_user("Enter the conversion factor for number density n to cm^-3", Float64(get(config, "conversionn", 1.0)))
    conversionT = ask_user("Enter the conversion factor for temperature T to K", Float64(get(config, "conversionT", 1.0)))
    config["conversionB"] = conversionB
    config["conversionn"] = conversionn
    config["conversionT"] = conversionT

    prompt_section("Simulation geometry")
    BoxLength_pc = ask_user("Side of the Box size (pc), please give a Float", Float64(get(config, "BoxLength_pc", 50.0)))
    BoxLength_pix = ask_user("Number of pixels along the line of sight", Int(get(config, "BoxLength_pix", 256)))
    config["BoxLength_pc"] = BoxLength_pc
    config["BoxLength_pix"] = BoxLength_pix

    prompt_section("Frequency setup")
    nustart = ask_user("Frequency range start (MHz)", Float64(get(config, "nustart", 120)))
    nuend = ask_user("Frequency range end (MHz)", Float64(get(config, "nuend", 167)))
    dnu = ask_user("Frequency resolution (MHz)", Float64(get(config, "dnu", 0.098)))
    config["nustart"] = nustart
    config["nuend"] = nuend
    config["dnu"] = dnu

    prompt_section("Faraday rotation")
    FaradayRotation = ask_user("Do you want to include Faraday rotation in the computation of Q and U? (Y/N)", get(config,"FaradayRotation", "Y");
        validate = is_yes_no, error_message = "Please answer Y or N.")
    faraday_flag = uppercase(FaradayRotation)
    config["FaradayRotation"] = FaradayRotation
    phimin = get(config, "phimin", -10.0)
    phimax = get(config, "phimax", 10.0)
    dphi = get(config, "dphi", 0.25)
    if faraday_flag == "Y"
        phimin = ask_user("Faraday depth range start (rad/m^2)", Float64(phimin))
        phimax = ask_user("Faraday depth range end (rad/m^2)", Float64(phimax))
        dphi = ask_user("Faraday depth resolution (rad/m^2)", Float64(dphi))
    end
    config["phimin"] = phimin
    config["phimax"] = phimax
    config["dphi"] = dphi

    rm_clean_cfg = get(config, "rm_clean", Dict{String, Any}())
    rm_clean_default_enabled = if rm_clean_cfg isa AbstractDict
        get(rm_clean_cfg, "enabled", false)
    else
        get(config, "do_rm_clean", false)
    end
    rm_clean_default_flag = rm_clean_default_enabled isa Bool ?
        (rm_clean_default_enabled ? "Y" : "N") :
        (uppercase(strip(String(rm_clean_default_enabled))) in ("Y", "YES", "TRUE", "1") ? "Y" : "N")

    rm_clean_enabled = false
    rm_clean_gain = Float64(rm_clean_cfg isa AbstractDict ? get(rm_clean_cfg, "gain", get(config, "rm_clean_gain", 0.1)) : get(config, "rm_clean_gain", 0.1))
    rm_clean_niter = Int(rm_clean_cfg isa AbstractDict ? get(rm_clean_cfg, "niter", get(config, "rm_clean_niter", 1000)) : get(config, "rm_clean_niter", 1000))
    rm_clean_threshold = Float64(rm_clean_cfg isa AbstractDict ? get(rm_clean_cfg, "threshold", get(config, "rm_clean_threshold", 0.0)) : get(config, "rm_clean_threshold", 0.0))

    if faraday_flag == "Y"
        rm_clean_choice = ask_user("Do you want to run RM-CLEAN on the Faraday dispersion function? (Y/N)", rm_clean_default_flag;
            validate = is_yes_no, error_message = "Please answer Y or N.")
        rm_clean_enabled = uppercase(rm_clean_choice) == "Y"
        if rm_clean_enabled
            rm_clean_gain = ask_user("RM-CLEAN loop gain (0 < gain <= 1)", rm_clean_gain)
            rm_clean_niter = ask_user("RM-CLEAN maximum number of iterations", rm_clean_niter)
            rm_clean_threshold = ask_user("RM-CLEAN absolute stopping threshold", rm_clean_threshold)
        end
    end
    config["rm_clean"] = Dict(
        "enabled" => rm_clean_enabled,
        "gain" => rm_clean_gain,
        "niter" => rm_clean_niter,
        "threshold" => rm_clean_threshold,
    )
    config["do_rm_clean"] = rm_clean_enabled

    prompt_section("Instrumental effects")
    responseSynchrotron = ask_user("Do you want to perform interferometric Fourier filtering for Synchrotron data? (Y/N)", get(config, "responseSynchrotron", "N");
        validate = is_yes_no, error_message = "Please answer Y or N.")
    kernel_size_synchrotron = uppercase(responseSynchrotron) == "Y" ? ask_user("Largest Fourier scale to keep for Synchrotron filtering (in pixels, e.g. 154)", get(config, "kernel_size_synchrotron", 154.0)) : nothing
    config["responseSynchrotron"] = responseSynchrotron
    config["kernel_size_synchrotron"] = kernel_size_synchrotron

    add_noise = ask_user("Do you want to add noise to Q and U? (Y/N)", get(config, "add_noise", "N");
        validate = is_yes_no, error_message = "Please answer Y or N.")
    config["add_noise"] = add_noise

    SNR_nu = uppercase(add_noise) == "Y" ? ask_user("Enter the desired SNR in the frequency space", get(config, "SNR_nu", 0.9)) : nothing
    config["SNR_nu"] = SNR_nu

    rng_seed = get(config, "rng_seed", nothing)
    if uppercase(add_noise) == "Y"
        default_seed = rng_seed === nothing ? 0 : Int(rng_seed)
        raw_seed = ask_user("Random seed for reproducible noise (integer; 0 for a random seed)", default_seed)
        rng_seed = Int(raw_seed) == 0 ? nothing : Int(raw_seed)
    end
    config["rng_seed"] = rng_seed

    prompt_section("Lines of sight")
    list_LOS = ["x", "y", "z"]
    los_prompt = "Enter 'all' to process all lines of sight (x, y, z) or specify comma-separated ones (e.g., x,y):"

    function parse_los_choice(los_choice)
        valid_los = String[]
        if uppercase(strip(los_choice)) == "ALL"
            append!(valid_los, list_LOS)
            return valid_los
        end

        for los in split(los_choice, ",")
            candidate = lowercase(strip(los))
            isempty(candidate) && continue
            if candidate in list_LOS
                push!(valid_los, candidate)
            else
                warn_user("Ignoring invalid line of sight: $(candidate). Valid options are x, y, z.")
            end
        end

        return valid_los
    end

    los_choice = get(config, "chosen_LOS_input", "all")
    valid_los = String[]
    while isempty(valid_los)
        los_choice = ask_user(
            los_prompt,
            los_choice;
            validate = choice -> !isempty(parse_los_choice(choice)),
            error_message = "[Error] No valid lines of sight provided. Use 'all' or comma-separated values like x,y and try again.",
        )
        valid_los = parse_los_choice(los_choice)
    end

    config["chosen_LOS_input"] = uppercase(strip(los_choice)) == "ALL" ? "all" : join(valid_los, ",")
    chosen_LOS = valid_los
    config["chosen_LOS"] = chosen_LOS
    field_sources = normalize_field_sources(get(config, "field_sources", get(config, "input_fields", nothing)))
    physical_mask = normalize_physical_mask(get(config, "physical_mask", get(config, "mask", nothing)))
    density_kind, mean_molecular_weight, hydrogen_mass_g = build_density_parameters(config)

    prompt_section("Synchrotron emissivity")
    interpolation_default = get(config, "interpolation_file_path", joinpath(homedir(), "emissivity.dat"))
    interpolation_file_path = interpolation_default
    while true
        interpolation_file_path = ask_user("Enter the path to the interpolation file", interpolation_file_path)
        validation_error = ensure_readable_file(interpolation_file_path; expected_exts=[".dat"])
        if validation_error === nothing
            break
        else
            warn_user(validation_error)
        end
    end
    config["interpolation_file_path"] = interpolation_file_path

    prompt_section("Electron density")
    ne_option = ""
    while true
        ne_option = ask_user("Choose electron density prescription: (1) Wolfire et al. 2003, (2) Proportional to nH, (3) Provide ne cube", get(config, "ne_option", "1"))
        ne_option in ("1", "2", "3") && break
        warn_user("Please choose 1, 2, or 3 for the electron density prescription.")
    end
    config["ne_option"] = ne_option

    IonizationFraction = get(config, "IonizationFraction", nothing)
    wolfire_constants = get(config, "wolfire_constants", nothing)
    if wolfire_constants === nothing && haskey(config, "ne")
        ne_config = config["ne"]
        if ne_config isa AbstractDict &&
                all(haskey(ne_config, key) for key in ("zeta", "Geff", "phiPAH", "XC"))
            wolfire_constants = (
                Float64(ne_config["zeta"]),
                Float64(ne_config["Geff"]),
                Float64(ne_config["phiPAH"]),
                Float64(ne_config["XC"]),
            )
        end
    end
    if ne_option == "2"
        IonizationFraction = ask_user("Enter the ionization fraction for the alternative prescription", get(config, "IonizationFraction", 0.01))
        config["IonizationFraction"] = IonizationFraction
    elseif ne_option == "3"
        missing_cubes = String[]
        for simu in chosen_simu
            source = simulation_field_source(simu, "densityHp", field_sources)
            try
                _validate_simulation_source(source)
            catch
                push!(missing_cubes, simu)
            end
        end
        if !isempty(missing_cubes)
            throw_config_error("Electron density cube `densityHp` is missing for: $(join(missing_cubes, ", ")). Configure `field_sources.densityHp` if it is named differently.";
                code=:missing_density_cube)
        end
    else
        zeta = ask_user("zeta (ionization rate by Cosmic Rays)", get(config, "zeta", 2.5e-16))
        Geff = ask_user("Geff (effective radiation field)", get(config, "Geff", 1.0))
        phiPAH = ask_user("phiPAH (collision rate parameter for PAH)", get(config, "phiPAH", 0.5))
        XC = ask_user("XC (Conversion factor of H into C)", get(config, "XC", 1.4e-4))
        wolfire_constants = (Float64(zeta), Float64(Geff), Float64(phiPAH), Float64(XC))
        config["zeta"] = wolfire_constants[1]
        config["Geff"] = wolfire_constants[2]
        config["phiPAH"] = wolfire_constants[3]
        config["XC"] = wolfire_constants[4]
        config["ne"] = merge(get(config, "ne", Dict{String, Any}()), Dict(
            "mode" => "1",
            "zeta" => wolfire_constants[1],
            "Geff" => wolfire_constants[2],
            "phiPAH" => wolfire_constants[3],
            "XC" => wolfire_constants[4],
        ))
    end

    prompt_section("Save configuration")
    save_path_default = get(config, "config_path", joinpath(base_dir, "moose_config.json"))
    config_path = ask_user("Enter the path where the configuration should be saved", save_path_default)
    config["config_path"] = config_path

    cfg = RunConfig(
        base_dir,
        chosen_simu,
        chosen_LOS,
        conversionB,
        conversionn,
        conversionT,
        faraday_flag,
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
        get(config, "log_progress", true),
        rng_seed;
        rm_clean_enabled = rm_clean_enabled,
        rm_clean_gain = rm_clean_gain,
        rm_clean_niter = rm_clean_niter,
        rm_clean_threshold = rm_clean_threshold,
        field_sources = field_sources,
        physical_mask = physical_mask,
        density_kind = density_kind,
        mean_molecular_weight = mean_molecular_weight,
        hydrogen_mass_g = hydrogen_mass_g,
    )

    return cfg, config
end
