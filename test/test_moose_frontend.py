from __future__ import annotations

import io
import json
import os
import subprocess
import tempfile
import unittest
from contextlib import contextmanager, redirect_stdout
from datetime import datetime, timedelta
from pathlib import Path
from unittest.mock import patch

from python import moose_frontend


@contextmanager
def working_directory(path: str):
    """Run the block with ``path`` as the process working directory."""

    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


class FrontendTestCase(unittest.TestCase):
    """Shared helpers; holds no tests of its own."""

    def parse_options(self, argv: list[str]) -> moose_frontend.FrontendOptions:
        parser = moose_frontend.build_parser()
        parsed = parser.parse_args(argv)
        return moose_frontend.validate_args(parser, parsed)

    def run_with_stub_julia(self, options, returncode: int = 0):
        """Invoke run_julia_frontend with the Julia subprocess stubbed out."""

        with patch.object(moose_frontend.shutil, "which", return_value="/usr/bin/julia"):
            with patch.object(moose_frontend.subprocess, "run") as run_mock:
                run_mock.return_value = subprocess.CompletedProcess([], returncode)
                moose_frontend.run_julia_frontend(options)
        return run_mock


class MooseFrontendStateTests(FrontendTestCase):
    def test_validation_normalizes_without_mutating_raw_namespace(self) -> None:
        parser = moose_frontend.build_parser()
        parsed = parser.parse_args(["--los", "Z,all", "--dry-run"])

        options = moose_frontend.validate_args(parser, parsed)

        self.assertEqual(parsed.los, ["Z,all"])
        self.assertEqual(options.los, ("z", "x", "y"))

    def test_build_julia_args_uses_validated_state(self) -> None:
        options = self.parse_options(
            [
                "--simu",
                "/data/sim-a",
                "--los",
                "x,y",
                "--faraday",
                "y",
                "--phimin",
                "-10",
                "--phimax",
                "10",
                "--dphi",
                "0.5",
                "--precision",
                "float32",
                "--tile-size",
                "64",
                "--resume",
                "safe",
                "--plan",
                "--density-kind",
                "mass_density",
                "--mean-molecular-weight",
                "1.4",
                "--hydrogen-mass-g",
                "1.6726231e-24",
                "--quiet",
            ]
        )

        self.assertEqual(
            moose_frontend.build_julia_args(options),
            [
                "--simu",
                "/data/sim-a",
                "--los",
                "x",
                "--los",
                "y",
                "--density-kind",
                "mass_density",
                "--mean-molecular-weight",
                "1.4",
                "--hydrogen-mass-g",
                "1.6726231e-24",
                "--faraday",
                "Y",
                "--phimin",
                "-10.0",
                "--phimax",
                "10.0",
                "--dphi",
                "0.5",
                "--precision",
                "float32",
                "--tile-size",
                "64",
                "--resume",
                "safe",
                "--quiet",
                "--plan",
            ],
        )

    def test_dry_run_prints_command_without_requiring_julia(self) -> None:
        options = self.parse_options(
            [
                "--julia-binary",
                "definitely-not-installed-julia",
                "--simu",
                "/data/sim-a",
                "--los",
                "z",
                "--dry-run",
            ]
        )

        stdout = io.StringIO()
        with patch.object(moose_frontend.shutil, "which", side_effect=AssertionError):
            with redirect_stdout(stdout):
                moose_frontend.run_julia_frontend(options)

        self.assertIn("Julia command:", stdout.getvalue())
        self.assertIn("definitely-not-installed-julia", stdout.getvalue())
        self.assertIn(f"--project={moose_frontend.REPO_ROOT}", stdout.getvalue())
        self.assertIn(str(moose_frontend.JULIA_ENTRYPOINT), stdout.getvalue())

    def test_dry_run_logs_composed_command(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "invocations.jsonl"
            options = self.parse_options(
                [
                    "--julia-binary",
                    "definitely-not-installed-julia",
                    "--los",
                    "all",
                    "--dry-run",
                    "--log-file",
                    str(log_path),
                ]
            )

            with redirect_stdout(io.StringIO()):
                moose_frontend.run_julia_frontend(options)

            log_text = log_path.read_text(encoding="utf-8")
            self.assertIn('"status": 0', log_text)
            self.assertIn('"message": "Dry run: Julia command not executed."', log_text)
            self.assertIn('"--los", "x", "--los", "y", "--los", "z"', log_text)

            # UTC instant, "Z" suffix rather than "+00:00". Guards the switch
            # from datetime.UTC (3.11+) to timezone.utc.
            timestamp = json.loads(log_text)["timestamp"]
            self.assertTrue(timestamp.endswith("Z"), timestamp)
            parsed = datetime.fromisoformat(timestamp.replace("Z", "+00:00"))
            self.assertEqual(parsed.utcoffset(), timedelta(0))


class MooseFrontendPathResolutionTests(FrontendTestCase):
    """The Julia subprocess runs with cwd=REPO_ROOT, not the caller's cwd.

    Anything the front-end validates against the caller's working directory has
    to be absolutized before it is handed over, or Julia resolves it somewhere
    else entirely.
    """

    def test_relative_config_path_is_absolutized_before_reaching_julia(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            config = Path(tmpdir) / "myconfig.json"
            config.write_text("{}", encoding="utf-8")

            with working_directory(tmpdir):
                options = self.parse_options(["myconfig.json"])

            self.assertEqual(options.config_path, config.resolve())

            run_mock = self.run_with_stub_julia(options)
            command = run_mock.call_args[0][0]

            self.assertEqual(run_mock.call_args[1]["cwd"], moose_frontend.REPO_ROOT)
            self.assertIn(str(config.resolve()), command)
            # The bare name would be looked up inside the repository instead.
            self.assertNotIn("myconfig.json", command)

    def test_relative_base_dir_is_absolutized(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            sims = Path(tmpdir) / "sims"
            sims.mkdir()

            with working_directory(tmpdir):
                options = self.parse_options(["--base-dir", "sims"])

            self.assertEqual(options.base_dir, str(sims.resolve()))

            args = moose_frontend.build_julia_args(options)
            self.assertIn(str(sims.resolve()), args)
            self.assertNotIn("sims", args)

    def test_user_home_shorthand_is_expanded(self) -> None:
        options = self.parse_options(["--base-dir", "~/simulations"])

        self.assertEqual(options.base_dir, str(Path.home() / "simulations"))

    def test_relative_simu_stays_relative_to_base_dir(self) -> None:
        # collect_simulations joins a relative `simu` onto `base_dir`, and the
        # interpolation lookup does the same, so neither may be absolutized
        # against the caller's cwd. This guards against over-correcting the
        # config/base_dir fix.
        options = self.parse_options(
            [
                "--base-dir",
                "/data",
                "--simu",
                "simulation_001",
                "--interpolation",
                "emissivity.dat",
            ]
        )

        args = moose_frontend.build_julia_args(options)
        self.assertEqual(options.simu, ("simulation_001",))
        self.assertIn("simulation_001", args)
        self.assertIn("emissivity.dat", args)

    def test_nonzero_exit_raises_and_is_logged(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "invocations.jsonl"
            options = self.parse_options(["--log-file", str(log_path)])

            with self.assertRaises(moose_frontend.JuliaInvocationError):
                self.run_with_stub_julia(options, returncode=2)

            log_text = log_path.read_text(encoding="utf-8")
            self.assertIn('"status": 2', log_text)
            self.assertIn("exited with status 2", log_text)


if __name__ == "__main__":
    unittest.main()
