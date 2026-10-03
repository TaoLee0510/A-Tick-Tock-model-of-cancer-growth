#!/usr/bin/env python3
"""Measure production prefixes and verify a restart in a separate process."""

import argparse
import json
from pathlib import Path
import platform
import subprocess


STATE_KEYS = ("model", "grid_shape", "configured_hours", "time_hours", "mass",
              "active_mass", "abm_fraction", "radius_99", "vascular_roots",
              "vascular_length", "state_checksum", "field_checksum")


def run(executable, config, model, stop, threads, extra):
    command = [str(executable.resolve()), "--model", model, "--config", str(config.resolve()),
               "--stop-hours", str(stop), "--threads", str(threads), *extra]
    process = subprocess.run(command, capture_output=True, text=True, check=False)
    if process.returncode:
        raise RuntimeError(f"{model} production run failed: {process.stderr.strip()}")
    return json.loads(process.stdout)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--stop-hours", type=float, default=1.0)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--resume-threads", type=int, default=8)
    parser.add_argument("--models", nargs="+", choices=("abm", "pde", "hybrid"),
                        default=("abm", "pde", "hybrid"))
    args = parser.parse_args()
    if args.stop_hours <= 0.0 or min(args.threads, args.resume_threads) < 1:
        parser.error("positive duration and thread counts are required")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    configs = {
        "abm": root / "ATCG3D_SharedRules/config/production_abm_2d_2000_r200_v17.yaml",
        "pde": root / "ATCG3D_SharedRules/config/production_pde_2d_2000_r200_v17.yaml",
        "hybrid": root / "ATCG3D_Hybrid/config/production_2d_2000_r200_v5.yaml"}
    report = {"protocol": "production_benchmark_contract.md", "passed": False,
              "stop_hours": args.stop_hours, "configured_hours": 2160.0,
              "full_duration": args.stop_hours == 2160.0, "models": {},
              "system": platform.system(), "machine": platform.machine(),
              "timing_scope": "construction, initialization or restore, operators, diagnostics, checkpoint and teardown"}
    for model in args.models:
        checkpoint = args.output.parent / f"production_{model}.checkpoint"
        for suffix in ("", ".pde.bin", ".resource.bin"):
            Path(str(checkpoint) + suffix).unlink(missing_ok=True)
        reference = run(args.exe, configs[model], model, args.stop_hours, args.threads, [])
        prefix = run(args.exe, configs[model], model, args.stop_hours / 2.0, args.threads,
                     ["--checkpoint", str(checkpoint.resolve())])
        resumed = run(args.exe, configs[model], model, args.stop_hours, args.resume_threads,
                      ["--resume", str(checkpoint.resolve())])
        exact = all(reference[key] == resumed[key] for key in STATE_KEYS)
        interior = max(reference["radius_99_half_width_ratio"],
                       prefix["radius_99_half_width_ratio"]) < 0.8
        report["models"][model] = {"config": configs[model].name,
                                    "reference": reference, "checkpoint_prefix": prefix,
                                    "resumed": resumed, "restart_exact": exact,
                                    "interior_at_observations": interior,
                                    "checkpoint_bytes": prefix["checkpoint_bytes"]}
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        if not exact:
            raise RuntimeError(f"{model} production restart changed biological state")
        if reference["grid_shape"] != [2000, 2000, 1] or reference["configured_hours"] != 2160.0:
            raise RuntimeError("unexpected production grid or configured duration")
        print(f"{model}: {reference['wall_seconds']:.17g}s, "
              f"{reference['peak_resident_bytes']} bytes, restart exact, "
              f"r99/half_width={reference['radius_99_half_width_ratio']:.17g}", flush=True)
    report["passed"] = all(record["restart_exact"] and record["interior_at_observations"]
                           for record in report["models"].values())
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
