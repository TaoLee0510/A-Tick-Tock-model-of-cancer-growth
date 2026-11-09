#!/usr/bin/env python3
"""Measure complete 3D shared models and require thread-independent states."""

import argparse
import json
from pathlib import Path
import platform
import statistics
import subprocess
import time


def biological_state(record):
    return {key: value for key, value in record.items()
            if key not in ("wall_seconds", "peak_resident_bytes", "threads")}


def wait_for_test_suites(manifest):
    if not manifest.exists():
        return
    dependencies = json.loads(manifest.read_text())
    paths = [Path(value) for value in dependencies["logs"]]
    deadline = time.monotonic() + dependencies.get("timeout_seconds", 14400)
    print("Waiting for declared test suites before isolated timing.", flush=True)
    while True:
        outputs = [path.read_text() if path.exists() else "" for path in paths]
        if all("Total Test time (real)" in output for output in outputs):
            if not all("100% tests passed" in output for output in outputs):
                raise RuntimeError("a prerequisite test suite failed")
            return
        if time.monotonic() >= deadline:
            raise RuntimeError("benchmark prerequisites did not finish")
        time.sleep(5)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=root / "build-codex/atcg3d_shared_3d_benchmark")
    parser.add_argument("--config", type=Path,
                        default=root / "ATCG3D_SharedRules/config/three_dimensional_r20_v16.yaml")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--test-suites", type=Path,
                        help="optional JSON manifest of test logs to await before timing")
    args = parser.parse_args()
    if args.repeat < 1:
        parser.error("use at least one repetition")
    manifest = args.test_suites or args.output.parent / "benchmark_dependencies.json"
    wait_for_test_suites(manifest)
    records = []
    groups = {}
    for model in ("abm", "pde"):
        reference = None
        realizations = {threads: [] for threads in (1, 4, 8)}
        for repetition in range(1, args.repeat + 1):
            for threads in (1, 4, 8):
                command = [str(args.exe.resolve()), "--config", str(args.config.resolve()),
                           "--model", model, "--threads", str(threads)]
                result = subprocess.run(command, capture_output=True, text=True, check=True)
                record = json.loads(result.stdout)
                if record["grid_shape"] != [256, 256, 256] or record["hours"] != 1.0:
                    raise RuntimeError("unexpected 3D benchmark dimensions or duration")
                if min(record[key] for key in ("active_mass", "vascular_roots", "vascular_length")) <= 0.0:
                    raise RuntimeError("benchmark did not exercise activation and angiogenesis")
                state = biological_state(record)
                if reference is None:
                    reference = state
                if state != reference:
                    raise RuntimeError("3D biological state changed with threads or repetition")
                record["repetition"] = repetition
                records.append(record)
                realizations[threads].append(record)
                print(f"{model} threads={threads} repetition={repetition} "
                      f"seconds={record['wall_seconds']:.17g}", flush=True)
        for threads, group in realizations.items():
            groups[f"{model}_{threads}"] = {
                "model": model, "threads": threads,
                "median_wall_seconds": statistics.median(record["wall_seconds"] for record in group),
                "maximum_peak_resident_bytes": max(record["peak_resident_bytes"] for record in group)}
        for threads in (1, 4, 8):
            groups[f"{model}_{threads}"]["speedup"] = (
                groups[f"{model}_1"]["median_wall_seconds"] /
                groups[f"{model}_{threads}"]["median_wall_seconds"])
    report = {"protocol": "three_dimensional_sector_contract.md", "config": args.config.name,
              "passed": True, "repetitions": args.repeat, "groups": groups, "records": records,
              "system": platform.system(), "release": platform.release(), "machine": platform.machine(),
              "timing_scope": "construction, coupled initialization, all operators, diagnostics and teardown"}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
