#!/usr/bin/env python3
"""Compare mixed populations with native limits and refine splitting/exchange intervals."""

import argparse
import concurrent.futures
import json
import math
from pathlib import Path
import statistics
import subprocess

from validate_abm_vs_pde import TOLERANCES, compare


SCENARIOS = {
    "abm": ("all_abm", 0.25, 1.0),
    "pde": ("all_pde", 0.25, 1.0),
    "coarse_step": ("adaptive", 0.25, 1.0),
    "middle_step": ("adaptive", 0.125, 1.0),
    "fine_step": ("adaptive", 0.0625, 1.0),
    "coarse_exchange": ("adaptive", 0.125, 2.0),
    "fine_exchange": ("adaptive", 0.125, 0.5),
}
REFINEMENT_TOLERANCES = {key: 0.10 for key in TOLERANCES}
REFINEMENT_TOLERANCES.update(active_fraction=0.02, radial_profile_L2=0.15)


def run_seed(executable, config, seed, output, threads):
    result = {}
    for name, (mode, step, exchange) in SCENARIOS.items():
        path = output / f"seed_{seed:04d}_{name}.json"
        command = [str(executable), "--config", str(config), "--mode", mode,
                   "--seed", str(seed), "--threads", str(threads), "--no-output",
                   "--validation-report", "--step-hours", str(step),
                   "--exchange-hours", str(exchange), "--report", str(path)]
        process = subprocess.run(command, text=True, capture_output=True, check=False)
        if process.returncode:
            raise RuntimeError(f"{name} seed {seed} failed: {process.stderr.strip()}")
        result[name] = json.loads(path.read_text())
        if result[name]["K_mass"] <= 0:
            raise RuntimeError("hybrid comparison requires positive K mass")
        for key in TOLERANCES:
            if key != "radial_profile_L2" and not math.isfinite(result[name][key]):
                raise RuntimeError(f"nonfinite {key} in {name}")
    times = [r["time_hours"] for r in result.values()]
    if min(times) <= 0 or max(times) - min(times) > 1e-9:
        raise RuntimeError("hybrid scenario endpoints differ")
    return result


def comparison(samples, reference, candidate, tolerances):
    return compare([{"abm": s[reference], "pde": s[candidate]} for s in samples], tolerances)


def refinement(samples, coarse, middle, fine):
    first = comparison(samples, fine, coarse, REFINEMENT_TOLERANCES)
    second = comparison(samples, fine, middle, REFINEMENT_TOLERANCES)
    order = {}
    for metric in REFINEMENT_TOLERANCES:
        if metric == "radial_profile_L2":
            jackknife = []
            for i in range(len(samples)):
                remaining = samples[:i] + samples[i + 1:]
                a = comparison(remaining, fine, coarse, {metric: 1})[metric]["error"]
                b = comparison(remaining, fine, middle, {metric: 1})[metric]["error"]
                jackknife.append(b - a)
            mean = statistics.mean(jackknife)
            standard_error = math.sqrt((len(samples) - 1) / len(samples) *
                sum((value - mean) ** 2 for value in jackknife))
        else:
            values = [s[middle][metric] - s[coarse][metric] for s in samples]
            scale = 1.0 if metric == "active_fraction" else max(1e-12, first[metric]["abm_mean"])
            standard_error = statistics.stdev(values) / (math.sqrt(len(samples)) * scale)
        allowance = 2.131 * standard_error
        order[metric] = {"coarse_error": first[metric]["error"],
                         "middle_error": second[metric]["error"],
                         "sampling_allowance": allowance,
                         "passed": second[metric]["error"] <= first[metric]["error"] + allowance + 1e-12}
    passed = all(v["passed"] for v in second.values()) and all(v["passed"] for v in order.values())
    return {"coarse_vs_fine": first, "middle_vs_fine": second,
            "refinement_order": order, "passed": passed}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=root / "build-codex/atcg3d_hybrid")
    parser.add_argument("--config", type=Path, default=root / "ATCG3D_Hybrid/config/hybrid_regular_cycle_v2.yaml")
    parser.add_argument("--seeds", type=int, default=16)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output", type=Path, default=root / "hybrid-validation")
    args = parser.parse_args()
    if args.seeds < 16 or args.jobs < 1 or args.threads < 1:
        parser.error("use at least 16 seeds and positive jobs/threads")
    args.output.mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_seed, args.exe.resolve(), args.config.resolve(), seed,
                               args.output, args.threads) for seed in range(1, args.seeds + 1)]
        samples = [future.result() for future in futures]
    limits = {name: comparison(samples, name, "coarse_step", TOLERANCES) for name in ("abm", "pde")}
    refinements = {"step": refinement(samples, "coarse_step", "middle_step", "fine_step"),
                   "exchange": refinement(samples, "coarse_exchange", "middle_step", "fine_exchange")}
    coverage = {key: statistics.mean(s["coarse_step"][key] for s in samples)
                for key in ("abm_mass", "pde_mass", "to_abm", "to_pde")}
    passed = all(v > 0 for v in coverage.values()) and all(
        metric["passed"] for case in limits.values() for metric in case.values()) and all(
        result["passed"] for result in refinements.values())
    report = {"schema_version": 1, "seeds": args.seeds, "config": args.config.name,
              "time_hours": samples[0]["abm"]["time_hours"], "passed": passed,
              "scenarios": SCENARIOS, "mixed_coverage": coverage,
              "native_comparisons": limits, "refinements": refinements,
              "limits": "regular-cycle growing mixed case; spatial and distribution correlations remain approximate",
              "refinement_definition": "middle/fine drift passes fixed tolerances; drift decreases within paired t(15) uncertainty, with profile jackknife uncertainty"}
    (args.output / "validation_report.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = ["# Hybrid ensemble and refinement", "",
             f"{args.seeds} paired seeds, {report['time_hours']:g} hours. Result: {'PASS' if passed else 'FAIL'}.", "",
             "| Comparison | Metric | Error | Tolerance | Sampling allowance | Result |",
             "| --- | --- | ---: | ---: | ---: | --- |"]
    cases = {f"adaptive/{name}": metrics for name, metrics in limits.items()}
    cases.update({f"{name} middle/fine": result["middle_vs_fine"] for name, result in refinements.items()})
    for name, metrics in cases.items():
        for metric, value in metrics.items():
            lines.append(f"| {name} | {metric} | {value['error']:.6f} | {value['tolerance']:.4f} | "
                         f"{value['sampling_allowance']:.6f} | {'PASS' if value['passed'] else 'FAIL'} |")
    lines += ["", "Coverage: " + json.dumps(coverage, sort_keys=True), "",
              "Refinement order: " + json.dumps({key: value["passed"] for key, value in refinements.items()}),
              "", report["refinement_definition"], report["limits"]]
    (args.output / "validation_report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
