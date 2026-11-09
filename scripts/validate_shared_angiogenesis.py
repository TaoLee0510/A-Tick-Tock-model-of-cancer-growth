#!/usr/bin/env python3
"""Compare unique vascular centerlines and event counts under the shared law."""

import argparse
import concurrent.futures
import json
import math
from pathlib import Path
import subprocess

from validate_abm_vs_pde import boundary_summary, run_pair, statistical_comparison


MARGINS = {
    "vascular_length": 0.75,
    "perfused_volume": 0.75,
    "lesion_perfused_fraction": 0.15,
    "vascular_branches": 0.75,
    "vascular_anastomoses": 0.75,
    "vascular_roots": 0.75,
}


def analyze(samples):
    times = [snapshot["time_hours"] for snapshot in samples[0]["abm"]["time_series"]]
    if times[0] != 0.0 or times[-1] < 240.0:
        raise RuntimeError("vascular validation requires initialization and at least 240 hours")
    tables = []
    for index, time in enumerate(times):
        records = [{model: pair[model]["time_series"][index] for model in ("abm", "pde", "baseline")}
                   for pair in samples]
        for pair in samples:
            for model in ("abm", "pde", "baseline"):
                actual = [snapshot["time_hours"] for snapshot in pair[model]["time_series"]]
                if actual != times:
                    raise RuntimeError("vascular ensembles have different sampling times")
        for record in records:
            for snapshot in record.values():
                if snapshot.get("vascular_length_model") != "unique_lattice_centerline_growth_v2":
                    raise RuntimeError("vascular length must be a unique centerline, not a volume proxy")
                if min(snapshot["grid_shape"][:2]) < 256:
                    raise RuntimeError("vascular validation requires at least 256-square voxels")
                for metric in MARGINS:
                    if not math.isfinite(snapshot[metric]) or snapshot[metric] < 0.0:
                        raise RuntimeError("invalid vascular metric")
        metrics = statistical_comparison(records, MARGINS)
        if time == 0.0:
            for metric, result in metrics.items():
                if result.get("undefined_reason") and all(
                        snapshot[metric] == 0.0 for record in records for snapshot in record.values()):
                    result.update(passed=True, inference_performed=False,
                                  initialization_identity=True,
                                  interpretation="exact structural initialization check; no TOST inference")
        tables.append({"time_hours": time, "metrics": metrics,
                       "passed": all(result["passed"] for result in metrics.values())})
    return {"samples": tables, "passed": all(table["passed"] for table in tables)}


def markdown(report):
    lines = ["# Shared VEGF angiogenesis validation", "",
             "Result: " + ("PASS" if report["passed"] else "FAIL") + ".",
             "Initialization zeros are exact setup checks, not equivalence inference.",
             "Post-initialization criteria use paired TOST at unchanged margins.",
             "The configuration is a controlled frozen-lesion mechanism experiment.", "",
             "| Hours | Metric | ABM mean | PDE mean | ABM baseline difference | Model/baseline RMS | Error | Margin | Paired t | Paired p | TOST lower p | TOST upper p | Result |",
             "| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |"]
    for table in report["trajectory"]["samples"]:
        for metric, result in table["metrics"].items():
            zero, intrinsic, tost = result["paired_zero_test"], result["baseline"], result["tost"]
            values = [result["abm_mean"], result["pde_mean"], intrinsic["difference"],
                      intrinsic["model_to_baseline_rms_ratio"], result["error"], result["margin"],
                      zero["t"], zero["p"], tost["lower"]["p"], tost["upper"]["p"]]
            rendered = ["null" if value is None else format(value, ".17g") for value in values]
            status = "INITIALIZATION" if result.get("initialization_identity") else (
                "PASS" if result["passed"] else "FAIL")
            lines.append(f"| {table['time_hours']:.17g} | {metric} | " + " | ".join(rendered) +
                         " | " + status + " |")
    return "\n".join(lines) + "\n"


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=root / "build-codex/atcg3d_shared_abm")
    parser.add_argument("--config", type=Path,
                        default=root / "ATCG3D_SharedRules/config/angiogenesis_controlled_v15.yaml")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", type=int, default=16)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--analysis-only", action="store_true")
    args = parser.parse_args()
    if args.seeds < 16 or args.jobs < 1 or args.threads < 1:
        parser.error("at least 16 seeds and positive jobs/threads are required")
    args.output.mkdir(parents=True, exist_ok=True)
    description = subprocess.run([str(args.exe.resolve()), "--config", str(args.config.resolve()), "--dry-run"],
                                 text=True, capture_output=True, check=True)
    effective = json.loads(description.stdout)
    continuum = effective["continuum_config"]
    if continuum["angiogenesis"]["model"] != "shared_vegf_lattice_v2":
        parser.error("use the shared vascular model, not a published legacy model")
    if min(continuum["grid"]["shape"][:2]) < 256 or continuum["numerics"]["end_time_hours"] < 240.0:
        parser.error("vascular validation requires at least 256-square voxels and 240 hours")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_pair, args.exe.resolve(), args.config.resolve(), seed,
                               args.output, args.threads, 24.0, args.seeds, args.analysis_only)
                   for seed in range(1, args.seeds + 1)]
        samples = [future.result() for future in futures]
    trajectory = analyze(samples)
    boundary = boundary_summary(samples)
    report = {"schema_version": 1, "config": args.config.name, "seeds": args.seeds,
              "seed_groups": {"abm_a": list(range(1, args.seeds + 1)),
                              "abm_b": list(range(args.seeds + 1, 2 * args.seeds + 1))},
              "sample_every_hours": 24.0, "margins": MARGINS,
              "protocol": "shared_angiogenesis_contract.md", "trajectory": trajectory,
              "boundary": boundary, "passed": trajectory["passed"] and all(
                  not value.get("boundary_influenced", False) for value in boundary.values()),
              "effective_vascular_law": continuum["angiogenesis"]}
    (args.output / "validation_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    (args.output / "validation_report.md").write_text(markdown(report))
    print("Shared vascular equivalence: " + ("PASS" if report["passed"] else "FAIL"))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
