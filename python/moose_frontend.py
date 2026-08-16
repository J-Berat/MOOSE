"""
Python front-end to invoke the Julia-based MOOSE CLI.

The script mirrors the flags available in `src/MOOSE_cli.jl` and forwards them
to a Julia process. It allows users who prefer Python tooling to prepare a
configuration, tweak overrides, and launch the MOOSE pipeline without having
to memorize the Julia command invocation.
"""

from __future__ import annotations

import argparse
import json
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Iterable, List

REPO_ROOT = Path(__file__).resolve().parent.parent
JULIA_ENTRYPOINT = REPO_ROOT / "src" / "MOOSE_cli.jl"


class JuliaInvocationError(RuntimeError):
    """Raised when the Julia subprocess returns a non-zero exit code."""


@dataclass(frozen=True)
class FrontendOptions:
    """Validated, normalized state needed to compose and run the Julia CLI."""

    config_path: Path | None
    julia_binary: str
    quiet: bool
    print_command: bool
    dry_run: bool
    write_back: bool
    plan: bool
    base_dir: str | None
    simu: tuple[str, ...]
    los: tuple[str, ...]
    interpolation: str | None
    conversionB: float | None
    conversionn: float | None
    conversionT: float | None
    density_kind: str | None
    mean_molecular_weight: float | None
    hydrogen_mass_g: float | None
    faraday: str | None
    phimin: float | None
    phimax: float | None
    dphi: float | None
    filtering: str | None
    kernel_size: float | None
    noise: str | None
    snr: float | None
    rng_seed: int | None
    precision: str | None
    tile_size: int | None
    resume: str | None
    outputs: tuple[str, ...]
    ne_option: str | None
    zeta: float | None
    Geff: float | None
    phiPAH: float | None
    XC: float | None
    log_file: Path | None


def resolve_config_path(parsed: argparse.Namespace) -> Path | None:
    """Return the selected config path, favoring the positional argument.

    The path is made absolute because the Julia subprocess is launched with
    ``cwd=REPO_ROOT``: a relative path is validated here against the caller's
    working directory but would be resolved by Julia against the repository
    root, so it has to be pinned down before it crosses that boundary.
    """

    config_value = parsed.config or parsed.config_path
    return _absolute_path(config_value) if config_value else None


def _absolute_path(value: str) -> Path:
    """Expand ``~`` and make ``value`` absolute against the caller's cwd."""

    return Path(value).expanduser().resolve()


def build_julia_args(options: FrontendOptions) -> List[str]:
    args: List[str] = []

    if options.config_path:
        args.append(str(options.config_path))

    if options.base_dir:
        args.extend(["--base-dir", str(options.base_dir)])

    for simu in options.simu:
        args.extend(["--simu", str(simu)])

    for los in options.los:
        args.extend(["--los", los])

    if options.interpolation:
        args.extend(["--interpolation", str(options.interpolation)])

    if options.conversionB is not None:
        args.extend(["--conversionB", str(options.conversionB)])

    if options.conversionn is not None:
        args.extend(["--conversionn", str(options.conversionn)])

    if options.conversionT is not None:
        args.extend(["--conversionT", str(options.conversionT)])

    if options.density_kind:
        args.extend(["--density-kind", options.density_kind])

    if options.mean_molecular_weight is not None:
        args.extend(["--mean-molecular-weight", str(options.mean_molecular_weight)])

    if options.hydrogen_mass_g is not None:
        args.extend(["--hydrogen-mass-g", str(options.hydrogen_mass_g)])

    if options.faraday:
        args.extend(["--faraday", options.faraday])

    if options.phimin is not None:
        args.extend(["--phimin", str(options.phimin)])

    if options.phimax is not None:
        args.extend(["--phimax", str(options.phimax)])

    if options.dphi is not None:
        args.extend(["--dphi", str(options.dphi)])

    if options.filtering:
        args.extend(["--filtering", options.filtering])

    if options.kernel_size is not None:
        args.extend(["--kernel-size", str(options.kernel_size)])

    if options.noise:
        args.extend(["--noise", options.noise])

    if options.snr is not None:
        args.extend(["--snr", str(options.snr)])

    if options.rng_seed is not None:
        args.extend(["--rng-seed", str(options.rng_seed)])

    if options.precision:
        args.extend(["--precision", options.precision])

    if options.tile_size is not None:
        args.extend(["--tile-size", str(options.tile_size)])

    if options.resume:
        args.extend(["--resume", options.resume])

    if options.outputs:
        args.extend(["--outputs", ",".join(options.outputs)])

    if options.ne_option:
        args.extend(["--ne-option", str(options.ne_option)])

    if options.zeta is not None:
        args.extend(["--zeta", str(options.zeta)])

    if options.Geff is not None:
        args.extend(["--Geff", str(options.Geff)])

    if options.phiPAH is not None:
        args.extend(["--phiPAH", str(options.phiPAH)])

    if options.XC is not None:
        args.extend(["--XC", str(options.XC)])

    if options.quiet:
        args.append("--quiet")

    if options.write_back:
        args.append("--write-back")

    if options.plan:
        args.append("--plan")

    return args


def format_command(command: Iterable[str]) -> str:
    """Return a shell-safe string representation of the command."""

    return shlex.join(list(command))


def _normalize_los_values(values: List[str] | None, parser: argparse.ArgumentParser) -> list[str]:
    """Return a normalized list of LOS choices, mirroring the Julia CLI semantics."""

    if not values:
        return []

    normalized: list[str] = []
    seen: set[str] = set()
    for raw_value in values:
        parts = [value.strip().lower() for value in raw_value.split(",") if value.strip()]
        for value in parts:
            if value == "all":
                for expanded in ("x", "y", "z"):
                    if expanded not in seen:
                        normalized.append(expanded)
                        seen.add(expanded)
                continue

            if value not in {"x", "y", "z"}:
                parser.error("--los must be x, y, z, or 'all'.")

            if value not in seen:
                normalized.append(value)
                seen.add(value)

    return normalized


def _normalize_flag(value: str | None) -> str | None:
    """Return normalized Y/N flags while preserving unspecified values."""

    return value.upper() if value else None


def validate_args(parser: argparse.ArgumentParser, parsed: argparse.Namespace) -> FrontendOptions:
    """Validate CLI arguments and return normalized front-end state."""

    if parsed.config and parsed.config_path:
        parser.error("Provide either a positional config path or --config, not both.")

    los = tuple(_normalize_los_values(parsed.los, parser))
    config_path = resolve_config_path(parsed)
    if config_path and not config_path.exists():
        parser.error(f"Config file not found: {config_path}")
    if parsed.write_back and config_path is None:
        parser.error("--write-back requires a config file path (positional or --config).")

    if parsed.filtering and parsed.filtering.upper() == "Y" and parsed.kernel_size is None:
        parser.error("--kernel-size is required when enabling filtering.")

    if parsed.noise and parsed.noise.upper() == "Y" and parsed.snr is None:
        parser.error("--snr is required when enabling noise injection.")

    if parsed.faraday and parsed.faraday.upper() == "Y":
        missing_faraday = [
            flag
            for flag, value in {
                "--phimin": parsed.phimin,
                "--phimax": parsed.phimax,
                "--dphi": parsed.dphi,
            }.items()
            if value is None
        ]
        if missing_faraday:
            parser.error(
                "Faraday rotation enabled; please provide values for "
                + ", ".join(missing_faraday)
            )

    # `base_dir` goes through `abspath` in normalize_base_dir, which resolves
    # against the Julia process cwd (REPO_ROOT), so it needs the same treatment
    # as the config path. `simu` and `interpolation` are deliberately left
    # alone: collect_simulations and the interpolation lookup join them onto
    # `base_dir` when they are relative, so absolutizing them here against the
    # caller's cwd would break the documented "base_dir + folder name" usage.
    base_dir = str(_absolute_path(parsed.base_dir)) if parsed.base_dir else None

    return FrontendOptions(
        config_path=config_path,
        julia_binary=parsed.julia_binary,
        quiet=parsed.quiet,
        print_command=parsed.print_command,
        dry_run=parsed.dry_run,
        write_back=parsed.write_back,
        plan=parsed.plan,
        base_dir=base_dir,
        simu=tuple(parsed.simu or ()),
        los=los,
        interpolation=parsed.interpolation,
        conversionB=parsed.conversionB,
        conversionn=parsed.conversionn,
        conversionT=parsed.conversionT,
        density_kind=parsed.density_kind,
        mean_molecular_weight=parsed.mean_molecular_weight,
        hydrogen_mass_g=parsed.hydrogen_mass_g,
        faraday=_normalize_flag(parsed.faraday),
        phimin=parsed.phimin,
        phimax=parsed.phimax,
        dphi=parsed.dphi,
        filtering=_normalize_flag(parsed.filtering),
        kernel_size=parsed.kernel_size,
        noise=_normalize_flag(parsed.noise),
        snr=parsed.snr,
        rng_seed=parsed.rng_seed,
        precision=parsed.precision,
        tile_size=parsed.tile_size,
        resume=parsed.resume,
        outputs=tuple(parsed.outputs.split(",")) if parsed.outputs else (),
        ne_option=parsed.ne_option,
        zeta=parsed.zeta,
        Geff=parsed.Geff,
        phiPAH=parsed.phiPAH,
        XC=parsed.XC,
        log_file=Path(parsed.log_file) if parsed.log_file else None,
    )


def _log_invocation(log_file: Path, command: Iterable[str], status: int, message: str | None) -> None:
    """Append a JSON log entry describing the invocation outcome."""

    log_file.parent.mkdir(parents=True, exist_ok=True)
    entry = {
        "timestamp": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
        "command": list(command),
        "command_string": format_command(command),
        "status": status,
        "message": message,
    }
    with log_file.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(entry))
        handle.write("\n")


def run_julia_frontend(options: FrontendOptions) -> None:
    cmd = [options.julia_binary, f"--project={REPO_ROOT}", str(JULIA_ENTRYPOINT)]
    cmd.extend(build_julia_args(options))
    command_string = format_command(cmd)

    if options.print_command or options.dry_run:
        print(f"Julia command: {command_string}")

    if options.dry_run:
        message = "Dry run: Julia command not executed."
        if options.log_file:
            _log_invocation(options.log_file, cmd, status=0, message=message)
        return

    julia_path = shutil.which(options.julia_binary)
    if julia_path is None:
        raise FileNotFoundError(options.julia_binary)

    process = subprocess.run(cmd, cwd=REPO_ROOT)

    failure_message = None
    if process.returncode != 0:
        failure_message = (
            "Julia front-end exited with status "
            f"{process.returncode}. Command: {command_string}"
        )

    if options.log_file:
        _log_invocation(options.log_file, cmd, process.returncode, failure_message)

    if failure_message:
        raise JuliaInvocationError(failure_message)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Python wrapper for the Julia MOOSE CLI. If no config is provided, "
            "the Julia tool will prompt for missing values and save them."
        )
    )

    parser.add_argument(
        "config",
        nargs="?",
        help="Optional path to a JSON config; passed as the first positional argument.",
    )
    parser.add_argument(
        "--config",
        dest="config_path",
        help="Explicit path to a JSON config file (equivalent to the positional argument).",
    )
    parser.add_argument(
        "--julia-binary",
        default="julia",
        help="Julia executable to invoke (default: julia in PATH).",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress status output from the Julia pipeline.",
    )
    parser.add_argument(
        "--print-command",
        action="store_true",
        help="Display the composed Julia command before running it.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the Julia command and exit without invoking Julia.",
    )
    parser.add_argument(
        "--write-back",
        action="store_true",
        help="Persist merged CLI overrides into the provided config file.",
    )
    parser.add_argument(
        "--plan",
        action="store_true",
        help="Validate input metadata and print channel, workload, RAM, and disk estimates without processing.",
    )

    parser.add_argument("--base-dir", help="Base directory containing the simulations.")
    parser.add_argument(
        "--simu",
        action="append",
        help="Simulation path (may be provided multiple times).",
    )
    parser.add_argument(
        "--los",
        action="append",
        help=(
            "Line of sight options; may be repeated. Accepts x, y, z, 'all', "
            "or comma-separated combinations such as x,y,z."
        ),
    )
    parser.add_argument(
        "--interpolation",
        help="Path to an emissivity interpolation file.",
    )
    parser.add_argument("--conversionB", type=float, help="Magnetic field conversion factor.")
    parser.add_argument("--conversionn", type=float, help="Density conversion factor.")
    parser.add_argument("--conversionT", type=float, help="Temperature conversion factor.")
    parser.add_argument(
        "--density-kind",
        choices=["number_density", "mass_density"],
        help="Interpret density as number density nH or mass density rho.",
    )
    parser.add_argument(
        "--mean-molecular-weight",
        type=float,
        help="Mean molecular weight used when --density-kind mass_density.",
    )
    parser.add_argument(
        "--hydrogen-mass-g",
        type=float,
        help="Hydrogen mass in grams used when --density-kind mass_density.",
    )
    parser.add_argument(
        "--faraday",
        choices=["Y", "N", "y", "n"],
        help="Enable or disable Faraday rotation.",
    )
    parser.add_argument("--phimin", type=float, help="Minimum Faraday depth.")
    parser.add_argument("--phimax", type=float, help="Maximum Faraday depth.")
    parser.add_argument("--dphi", type=float, help="Faraday depth resolution.")
    parser.add_argument(
        "--filtering",
        choices=["Y", "N", "y", "n"],
        help="Enable or disable synchrotron filtering.",
    )
    parser.add_argument("--kernel-size", dest="kernel_size", type=float, help="Filtering kernel size.")
    parser.add_argument(
        "--noise",
        choices=["Y", "N", "y", "n"],
        help="Enable or disable noise injection.",
    )
    parser.add_argument("--snr", type=float, help="Signal-to-noise ratio for noise injection.")
    parser.add_argument("--rng-seed", type=int, help="Random seed for reproducible noise injection.")
    parser.add_argument(
        "--precision",
        choices=["float64", "float32"],
        help="Working precision of the processed cubes (float32 halves the memory footprint).",
    )
    parser.add_argument(
        "--tile-size",
        dest="tile_size",
        type=int,
        help="Process the sky plane in bands of this many rows (for cubes larger than RAM). "
        "Incompatible with filtering, noise injection, and RM-CLEAN.",
    )
    parser.add_argument(
        "--resume",
        choices=["off", "safe"],
        help="Safely skip simulation/LOS work whose completion manifest still matches its inputs and config.",
    )
    parser.add_argument(
        "--outputs",
        help="Comma-separated output groups: integrated,stokes,rm,fdf,spectral_index,diagnostics (or all).",
    )
    parser.add_argument(
        "--ne-option",
        choices=["1", "2", "3"],
        help="Electron density prescription option (1, 2, or 3).",
    )
    parser.add_argument("--zeta", type=float, help="Wolfire model zeta constant.")
    parser.add_argument("--Geff", type=float, help="Wolfire model Geff constant.")
    parser.add_argument("--phiPAH", type=float, help="Wolfire model phiPAH constant.")
    parser.add_argument("--XC", type=float, help="Wolfire model XC constant.")
    parser.add_argument(
        "--log-file",
        help="Optional path to write a JSONL log entry for each invocation.",
    )

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    parsed = parser.parse_args(argv)
    options = validate_args(parser, parsed)
    run_julia_frontend(options)


if __name__ == "__main__":
    try:
        main()
    except JuliaInvocationError as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError as exc:
        missing = exc.filename if exc.filename else "julia"
        print(f"Unable to find required executable: {missing}", file=sys.stderr)
        sys.exit(1)
