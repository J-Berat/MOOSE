# Moose

**Moose** (*Mock Observation Of Synchrotron Emission*) is a Julia toolkit for
turning magnetohydrodynamic (MHD) simulation data into synthetic radio
observations. It reads cartesian or HEALPix simulation fields and produces
reproducible FITS products for synchrotron and Faraday analysis.

## Main features

- Compute synchrotron Stokes **I**, **Q**, and **U** cubes.
- Include Faraday rotation, instrumental filtering, and noise.
- Perform RM synthesis, RMSF diagnostics, and RM-CLEAN deconvolution.
- Fit physical Faraday models to Q/U spectra and compare them with AIC/BIC.
- Produce rotation-measure, spectral-index, polarization-fraction, and
  polarization-gradient maps.
- Compute structure functions and polarization diagnostics.
- Process cartesian FITS/HDF5, leaf-cell AMR HDF5, and HEALPix FITS data.
- Process large datasets in tiles to reduce memory use.
- Run interactively, from a JSON configuration, or through the Python wrapper.

Moose makes it possible to build mock radio-observation pipelines from MHD
simulations, compare simulated polarization with observations, study Faraday
depth structure and turbulence, and export analysis-ready products with
configuration and provenance metadata.

## Installation

Moose requires [Julia 1.10 or later](https://julialang.org/downloads/).
Clone the repository, enter its directory, and install the dependencies:

```bash
julia --startup-file=no --project -e 'using Pkg; Pkg.instantiate()'
```

Check the installation without input data:

```bash
julia --startup-file=no --project -e 'using Moose; run_moose(help=true)'
```

## Usage

Start the interactive workflow:

```julia
using Moose
run_moose()
```

Run an existing JSON configuration non-interactively:

```bash
julia --startup-file=no --project src/MOOSE_cli.jl /path/to/config.json --quiet
```

A template is available at `config/default_config.json`. The Python wrapper
provides the same command-line workflow:

```bash
python3 python/moose_frontend.py --config /path/to/config.json --quiet
```

To validate the complete pipeline with analytically known data:

```julia
using Moose
demo = make_demo_data("moose_demo")
MOOSE_from_config(demo.config_path; quiet=true)
```

Results are written as FITS files alongside the selected simulation. Each run
also records its configuration, provenance, and timing in `MOOSE_summary.log`.

### Optional all-sky example data

The Pluto tutorial can analyze the public Galactic Faraday rotation sky 2020
map by Hutschenreuter et al. The 24 MB normalized FITS file is intentionally not
stored in this repository. Download the current `faradaysky2020v2` release from
the [official MPA data page](https://wwwmpa.mpa-garching.mpg.de/~ensslin/research/data/faraday2020.html)
and convert it to MOOSE's HEALPix convention with:

```bash
julia --startup-file=no --project=. scripts/download_faraday2020.jl
```

The script verifies the upstream SHA-256 checksum, extracts
`faraday_sky_mean`, and writes `data/faraday2020v2.fits` as an NSIDE 512,
RING-ordered HEALPix map in Galactic coordinates. The data are distributed by
their authors under the ODC-By 1.0 license.

### Reprise après interruption

Pour les grands calculs tuilés, activez ensemble `tile_size` et
`"resume": "safe"`. MOOSE enregistre alors un checkpoint atomique après
chaque bande terminée. Si le processus est interrompu, relancer exactement la
même commande reprend à la bande suivante. Le checkpoint est ignoré et le
calcul redémarre proprement si la configuration, les entrées ou le découpage
en tuiles ont changé. Avec `"resume": "off"`, les fichiers partiels sont
supprimés comme auparavant.

```json
{
  "tile_size": 256,
  "resume": "safe"
}
```

### AMR inputs

MOOSE accepts AMR leaf cells stored in HDF5 and conservatively rasterizes each
intensive field onto the regular output grid. Configure the physical fields as
HDF5 datasets and add a shared `amr` geometry entry:

```json
"field_sources": {
  "Bx":          {"path": "amr.h5", "dataset": "cells/Bx"},
  "By":          {"path": "amr.h5", "dataset": "cells/By"},
  "Bz":          {"path": "amr.h5", "dataset": "cells/Bz"},
  "density":     {"path": "amr.h5", "dataset": "cells/density"},
  "temperature": {"path": "amr.h5", "dataset": "cells/temperature"},
  "amr": {
    "file": "amr.h5",
    "x": "cells/x", "y": "cells/y", "z": "cells/z",
    "size": "cells/dx",
    "bounds": [[0, 1], [0, 1], [0, 1]],
    "shape": [256, 256, 256]
  }
}
```

`size` may contain one width per cell or three axis widths. Alternatively use
`"level": "cells/level"`; a level `l` has width `domain_size / 2^l` (use
`level_offset` when levels are numbered relative to another root). Inputs must
contain leaf cells only. By default MOOSE rejects gaps and overlaps; set
`"strict": false` only for intentionally partial domains. AMR currently uses
the in-memory path and is therefore incompatible with `tile_size`.

For physical correctness, every leaf-cell boundary must align with the target
grid and `shape` must resolve the finest AMR level. MOOSE rejects coarser or
misaligned targets instead of averaging magnetic field, density, and
temperature before the nonlinear emissivity calculation. The cell-to-voxel
assignment is cached and reused for all fields and lines of sight.

## Citation

If you use Moose in scientific work, please cite the associated paper:
[Berat et al. (2026), A&A 708, A245](https://ui.adsabs.harvard.edu/abs/2026A%26A...708A.245B/abstract).

## License

Moose is distributed under the [MIT License](LICENSE).

## Author

**Jack Berat** — main developer
