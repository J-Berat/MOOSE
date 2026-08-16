# MOOSE tutorial notebook

`MOOSE_tutorial.jl` is a self-contained [Pluto.jl](https://plutojl.org/) notebook: a hands-on,
reactive tutorial covering the full MOOSE pipeline (test cube → line-of-sight geometry → Faraday
depth → synthetic Stokes Q/U → instrumental effects → RM synthesis / RM-CLEAN → validation against
`make_demo_data` → a full end-to-end run through `MOOSE_from_config`).

## Running it

```bash
julia --project=notebooks -e 'using Pkg; Pkg.instantiate(); using Pluto; Pluto.run()'
```

Then open `notebooks/MOOSE_tutorial.jl` from the Pluto launcher. The reproducible notebook
environment is stored in `notebooks/Project.toml` and `notebooks/Manifest.toml`; the local Moose
package is developed through the relative path `..`.

## Reproducible HTML export

```bash
julia --project=notebooks notebooks/export.jl
```

The export aborts if any Pluto cell fails. The generated
`notebooks/MOOSE_tutorial.html` file is intentionally ignored by Git.

## Automated notebook validation

```bash
julia --project=notebooks notebooks/check_notebook.jl
```

This executes both the one-screen and two-screen mock scenarios.

Figures saved from the notebook (§16, "Save the summary figure") are written to
`outputs/` at the repository root. Generated PDF and PNG files in that
directory are intentionally ignored by Git; `outputs/.gitkeep` only preserves
the directory.

The optional real-data section uses the Galactic Faraday rotation sky 2020
map. Follow the download and checksum instructions in the repository
`README.md`; the notebook expects the file at `data/faraday2020v2.fits`.
