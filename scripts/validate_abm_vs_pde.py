#!/usr/bin/env python3
"""Run paired realizations of the shared ABM/PDE rule contract and report errors."""

import argparse
import concurrent.futures
import json
import math
from pathlib import Path
import statistics
import subprocess


TOLERANCES = {
    "total_mass": 0.35,
    "r_K_ratio": 0.35,
    "active_fraction": 0.05,
    "r50": 0.25,
    "r90": 0.25,
    "r99": 0.30,
    "radial_profile_L2": 0.35,
}
VASCULAR_TOLERANCES = {
    "vascular_length": 0.75,
    "perfused_volume": 0.75,
    "lesion_perfused_fraction": 0.15,
}


def run_pair(executable, config, seed, output, threads):
    result = {}
    for model in ("abm", "pde"):
        path = output / f"seed_{seed:04d}_{model}.json"
        command = [str(executable), "--config", str(config), "--model", model,
                   "--seed", str(seed), "--threads", str(threads), "--report", str(path)]
        process = subprocess.run(command, capture_output=True, text=True, check=False)
        if process.returncode:
            raise RuntimeError(f"{model} seed {seed} failed: {process.stderr.strip()}")
        result[model] = json.loads(path.read_text())
        if not math.isfinite(result[model]["time_hours"]) or result[model]["time_hours"] <= 0.0:
            raise RuntimeError("invalid validation endpoint")
        if result[model]["K_mass"] <= 0.0:
            raise RuntimeError("r/K ratio requires positive K mass")
        for metric in TOLERANCES:
            if metric == "radial_profile_L2":
                continue
            if not math.isfinite(result[model][metric]):
                raise RuntimeError(f"nonfinite {metric} in {model} seed {seed}")
    if not math.isclose(result["abm"]["time_hours"], result["pde"]["time_hours"], abs_tol=1e-9):
        raise RuntimeError("ABM/PDE realization endpoints differ")
    return result


def radial_profile(samples, model):
    """Normalized cell-number density in width-four annuli, area-weighted in L2."""
    size = max(len(sample[model]["radial_mass"]) for sample in samples)
    bins = (size + 3) // 4
    mass = [0.0] * bins
    for sample in samples:
        for i, value in enumerate(sample[model]["radial_mass"]):
            mass[i // 4] += value / len(samples)
    total = sum(mass)
    areas = [math.pi * ((4 * i + 4) ** 2 - (4 * i) ** 2) for i in range(bins)]
    return [value / (max(total, 1e-15) * area) for value, area in zip(mass, areas)], areas


def compare(samples, tolerances):
    metrics = {}
    for metric, tolerance in tolerances.items():
        if metric == "radial_profile_L2":
            abm, areas = radial_profile(samples, "abm")
            pde, _ = radial_profile(samples, "pde")
            error = math.sqrt(sum(area * (a - b) ** 2 for a, b, area in zip(abm, pde, areas)) /
                              max(1e-30, sum(area * a * a for a, area in zip(abm, areas))))
            metrics[metric] = {"error": error, "tolerance": tolerance,
                               "sampling_allowance": 0.0, "passed": error <= tolerance,
                               "abm_profile": abm, "pde_profile": pde}
            continue
        values = {model: [sample[model][metric] for sample in samples] for model in ("abm", "pde")}
        means = {model: statistics.mean(v) for model, v in values.items()}
        differences = [p - a for a, p in zip(values["abm"], values["pde"])]
        standard_error = statistics.stdev(differences) / math.sqrt(len(samples))
        scale = 1.0 if metric in ("active_fraction", "lesion_perfused_fraction") else max(1e-12, abs(means["abm"]))
        error = abs(means["pde"] - means["abm"]) / scale
        allowance = 2.131 * standard_error / scale  # paired t(15), two-sided 95% for N=16
        metrics[metric] = {"abm_mean": means["abm"], "pde_mean": means["pde"],
                           "paired_standard_error": standard_error, "error": error,
                           "tolerance": tolerance, "sampling_allowance": allowance,
                           "passed": error <= tolerance + allowance}
    return metrics


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=root / "build-codex/atcg3d_shared_abm")
    parser.add_argument("--config", type=Path, default=root / "ATCG3D_SharedRules/config/validation_v7.yaml")
    parser.add_argument("--seeds", type=int, default=16)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--output", type=Path, default=root / "validation-output")
    parser.add_argument("--tolerances", type=Path, help="JSON object overriding named tolerances")
    parser.add_argument("--vascular", action="store_true", help="compare vascular length, volume and lesion perfusion")
    args = parser.parse_args()
    if args.seeds < 16 or args.jobs < 1 or args.threads < 1:
        parser.error("use at least 16 seeds and positive jobs/threads")
    executable, config = args.exe.resolve(), args.config.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    tolerances = dict(VASCULAR_TOLERANCES if args.vascular else TOLERANCES)
    if args.tolerances:
        overrides = json.loads(args.tolerances.read_text())
        if set(overrides) - set(tolerances) or any(not math.isfinite(v) or v < 0 for v in overrides.values()):
            parser.error("unknown or invalid tolerance")
        tolerances.update(overrides)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_pair, executable, config, seed, args.output, args.threads)
                   for seed in range(1, args.seeds + 1)]
        samples = [future.result() for future in futures]
    metrics = compare(samples, tolerances)
    passed = all(metric["passed"] for metric in metrics.values())
    report = {"schema_version": 1, "seeds": args.seeds, "time_hours": samples[0]["abm"]["time_hours"],
              "config": config.name, "passed": passed, "metrics": metrics,
              "sampling_method": "paired seeds, conservative t(15) allowance; fixed closure tolerances",
              "limits": "stochastic roots and continuous vessel/tip densities differ; early vascular smoke tolerances are broad" if args.vascular else
                        "mean-field closure and initial age distribution differ; this smoke case has sparse activation"}
    (args.output / "validation_report.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = ["# Shared ABM/PDE validation", "", f"{args.seeds} paired seeds, {samples[0]['abm']['time_hours']:g} hours. Result: {'PASS' if passed else 'FAIL'}.", "",
             "| Metric | ABM mean | PDE mean | Error | Tolerance | Sampling allowance | Result |",
             "| --- | ---: | ---: | ---: | ---: | ---: | --- |"]
    for name, metric in metrics.items():
        lines.append(f"| {name} | {metric.get('abm_mean', 0):.6g} | {metric.get('pde_mean', 0):.6g} | "
                     f"{metric['error']:.4f} | {metric['tolerance']:.4f} | {metric['sampling_allowance']:.4f} | "
                     f"{'PASS' if metric['passed'] else 'FAIL'} |")
    lines += ["", "Relative errors use the ABM ensemble mean. Active and lesion perfusion fractions use absolute errors.",
              "Scalar errors include a declared paired sampling allowance."]
    if not args.vascular:
        lines.append("The radial error uses normalized cell-number density in four-voxel annuli and an area-weighted L2 norm, with its fixed tolerance.")
    lines.append(report["limits"])
    (args.output / "validation_report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
