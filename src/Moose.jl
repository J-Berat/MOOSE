module Moose

using Crayons
using Dates
using JSON
using CSV
using DataFrames
using Distributions
using FITSIO
using HDF5
using Healpix
using StatsBase
using Interpolations
using Dierckx
using FFTW
using LinearAlgebra
using Logging
using Random
using SHA
using Serialization
using TOML

include(joinpath("Utils", "ArrayMath.jl"))
include(joinpath("Utils", "Prompts.jl"))
include(joinpath("Utils", "Progress.jl"))
include(joinpath("Utils", "Errors.jl"))
include(joinpath("Utils", "InputValidation.jl"))
include(joinpath("Utils", "AtomicWrite.jl"))

include(joinpath("FileIO", "FITSUtils.jl"))
include(joinpath("FileIO", "HDF5Utils.jl"))
include(joinpath("FileIO", "AMRGrid.jl"))
include(joinpath("FileIO", "Header.jl"))
include(joinpath("FileIO", "ReadSimulation.jl"))
include(joinpath("FileIO", "WriteDataOnDisk.jl"))
include(joinpath("FileIO", "HealpixIO.jl"))

include(joinpath("PhysicalParameters", "ConversionJyBeamtoK.jl"))
include(joinpath("PhysicalParameters", "BrightnessTemperature.jl"))
include(joinpath("PhysicalParameters", "Borientation.jl"))
include(joinpath("PhysicalParameters", "ElectronDensity.jl"))
include(joinpath("PhysicalParameters", "IntrinsicAngle.jl"))
include(joinpath("PhysicalParameters", "MagneticField.jl"))
include(joinpath("PhysicalParameters", "PolarizationAngle.jl"))
include(joinpath("PhysicalParameters", "PolarizationFraction.jl"))
include(joinpath("PhysicalParameters", "Pressure.jl"))
include(joinpath("PhysicalParameters", "RM.jl"))
include(joinpath("PhysicalParameters", "Constants.jl"))

include(joinpath("Frequencies", "FreqFile.jl"))

include(joinpath("Synchrotron", "EmissInterp.jl"))
include(joinpath("Synchrotron", "Pnu.jl"))
include(joinpath("Synchrotron", "ProcessSynchrotron.jl"))
include(joinpath("Synchrotron", "TiledProcessing.jl"))
include(joinpath("Synchrotron", "QUnu.jl"))
include(joinpath("Synchrotron", "Tnu.jl"))

include(joinpath("Faraday", "FaradayParameters.jl"))
include(joinpath("Faraday", "RMSynthesis.jl"))
include(joinpath("Faraday", "RMClean.jl"))
include(joinpath("Faraday", "QUFitting.jl"))

include(joinpath("Filtering", "Filter.jl"))

include(joinpath("Statistics", "EffectiveWidth.jl"))
include(joinpath("Statistics", "Moments.jl"))
include(joinpath("Statistics", "RMS.jl"))
include(joinpath("Statistics", "PowerSpectrum.jl"))
include(joinpath("Statistics", "PolarizationDiagnostics.jl"))
include(joinpath("Statistics", "PolarizationGradient.jl"))
include(joinpath("Statistics", "StructureFunction.jl"))
include(joinpath("Statistics", "SpectralIndex.jl"))
include(joinpath("Statistics", "Statistics.jl"))

include(joinpath("SyntheticObservations", "SimulationDiscovery.jl"))
include(joinpath("SyntheticObservations", "DictHeaderParameters.jl"))
include(joinpath("SyntheticObservations", "InstrumentalParameters.jl"))
include(joinpath("SyntheticObservations", "Moose.jl"))
include(joinpath("SyntheticObservations", "MOOSE_from_config.jl"))
include(joinpath("SyntheticObservations", "DemoData.jl"))

using .MooseFromConfig: MOOSE_from_config

# Stable public API: names exported below are the compatibility surface. Other
# `Moose.foo` bindings are implementation details, even when regression tests
# exercise them through qualified access.
export run_moose, MOOSE_from_config, preflight_plan, MooseError, cli_error, config_error,
       HealpixStack, HealpixRMResult, RMSynthesisHealpix, healpix_map,
       healpix_maps_from_stack, read_healpix_map, read_healpix_stack,
       detect_fits_grid, is_healpix_fits, is_image_fits,
       read_fits_grid, read_fits_grid_stack,
       write_healpix_map, write_healpix_stack, write_healpix_rm_result,
       write_healpix_cube, read_healpix_cube, HEALPIX_UNSEEN,
       healpix_udgrade, healpix_reorder, healpix_smooth,
       RMSynthesisAuto,
       rmsf_diagnostics, RMSFDiagnostics, write_rmsf,
       RMClean, RMCleanHealpix, RMCleanAuto, RMCleanResult,
       QUFit, QUFitCompare, QUFitCube, QUFitResult, qu_model, QU_FIT_MODELS,
       polarization_diagnostic_spectra, write_polarization_diagnostic_plots,
       polarization_gradient_map, structure_function, StructureFunctionResult,
       spectral_index_map, make_demo_data

const MOOSE_PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))

function moose_version()
    project_path = joinpath(MOOSE_PROJECT_ROOT, "Project.toml")
    try
        return String(get(TOML.parsefile(project_path), "version", "unknown"))
    catch
        return "unknown"
    end
end

function moose_git_hash()
    git = Sys.which("git")
    git === nothing && return "unknown"

    try
        revision = readchomp(`$git -C $MOOSE_PROJECT_ROOT rev-parse --short=12 HEAD`)
        dirty = success(`$git -C $MOOSE_PROJECT_ROOT diff --quiet --ignore-submodules HEAD`) ? "" : "+dirty"
        return revision * dirty
    catch
        return "unknown"
    end
end

function _canonical_json(value)
    if value isa AbstractDict
        parts = String[]
        for key in sort(collect(keys(value)); by = string)
            push!(parts, JSON.json(string(key)) * ":" * _canonical_json(value[key]))
        end
        return "{" * join(parts, ",") * "}"
    elseif value isa NamedTuple
        parts = String[]
        for key in sort(collect(keys(value)); by = string)
            push!(parts, JSON.json(string(key)) * ":" * _canonical_json(value[key]))
        end
        return "{" * join(parts, ",") * "}"
    elseif value isa Tuple
        return "[" * join((_canonical_json(item) for item in value), ",") * "]"
    elseif value isa AbstractVector
        return "[" * join((_canonical_json(item) for item in value), ",") * "]"
    else
        return JSON.json(value)
    end
end

function moose_config_hash(config::AbstractDict)
    compact_json = _canonical_json(config)
    return bytes2hex(sha256(Vector{UInt8}(codeunits(compact_json))))
end

end
