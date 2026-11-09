#!/usr/bin/env python3
"""Run paired realizations of the shared ABM/PDE rule contract and report errors."""

import argparse
import concurrent.futures
import json
import math
from pathlib import Path
import statistics
import subprocess

from equivalence_stats import baseline, jackknife_distance_tost, paired_tost, paired_zero_test


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


def run_model(executable, config, model, seed, output, threads, interval, analysis_only):
    path = output / f"seed_{seed:04d}_{model}.json"
    if not analysis_only:
        command = [str(executable), "--config", str(config), "--model", model,
                   "--seed", str(seed), "--threads", str(threads), "--report", str(path)]
        if interval:
            command += ["--sample-every-hours", str(interval)]
        process = subprocess.run(command, capture_output=True, text=True, check=False)
        if process.returncode:
            raise RuntimeError(f"{model} seed {seed} failed: {process.stderr.strip()}")
    result = json.loads(path.read_text())
    if result["seed"] != seed or result["model"] != model:
        raise RuntimeError(f"wrong realization identity in {path}")
    snapshots = result.get("time_series", [result])
    if interval and "time_series" not in result:
        raise RuntimeError(f"missing time series in {path}")
    previous = -1.0
    for snapshot in snapshots:
        time = snapshot["time_hours"]
        if not math.isfinite(time) or time < 0.0 or time <= previous:
            raise RuntimeError(f"invalid or non-increasing sampling times in {path}")
        previous = time
        if snapshot["K_mass"] <= 0.0:
            raise RuntimeError("r/K ratio requires positive K mass")
        for metric in TOLERANCES:
            if metric != "radial_profile_L2" and not math.isfinite(snapshot[metric]):
                raise RuntimeError(f"nonfinite {metric} in {model} seed {seed}")
        if any(not math.isfinite(value) or value < 0.0 for value in snapshot["radial_mass"]):
            raise RuntimeError("invalid radial mass")
    if result["time_hours"] <= 0.0 or not math.isclose(previous, result["time_hours"], abs_tol=1e-9):
        raise RuntimeError("invalid validation endpoint")
    return result


def run_pair(executable, config, seed, output, threads, interval=0.0,
             baseline_offset=0, analysis_only=False):
    result = {model: run_model(executable, config, model, seed, output, threads,
                              interval, analysis_only) for model in ("abm", "pde")}
    if baseline_offset:
        result["baseline"] = run_model(executable, config, "abm", seed + baseline_offset,
                                       output, threads, interval, analysis_only)
    times = [record["time_hours"] for record in result.values()]
    if max(times) - min(times) > 1e-9:
        raise RuntimeError("ensemble realization endpoints differ")
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
    dimensions = samples[0][model].get("spatial_dimensions", 2)
    if any(sample[model].get("spatial_dimensions", 2) != dimensions for sample in samples):
        raise RuntimeError("mixed spatial dimensions")
    coefficient = math.pi if dimensions == 2 else 4.0 * math.pi / 3.0
    areas = [coefficient * ((4 * i + 4) ** dimensions - (4 * i) ** dimensions)
             for i in range(bins)]
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


def profile_distance(samples, candidate="pde"):
    reference, measures = radial_profile(samples, "abm")
    profile, candidate_measures = radial_profile(samples, candidate)
    if measures != candidate_measures:
        raise RuntimeError("profile domains or spatial dimensions differ")
    distance = math.sqrt(sum(weight * (a - b) ** 2 for a, b, weight in
                             zip(reference, profile, measures)) /
                         max(1e-30, sum(weight * a * a for a, weight in zip(reference, measures))))
    return distance, reference, profile


def statistical_comparison(samples, tolerances):
    metrics = {}
    for metric, margin in tolerances.items():
        if metric == "radial_profile_L2":
            distance, reference, candidate = profile_distance(samples)
            deleted = [profile_distance(samples[:i] + samples[i + 1:])[0]
                       for i in range(len(samples))]
            result = jackknife_distance_tost(distance, deleted, margin)
            baseline_distance, _, independent = profile_distance(samples, "baseline")
            individual_baseline = [profile_distance([sample], "baseline")[0] for sample in samples]
            individual_model = [profile_distance([sample])[0] for sample in samples]
            excess = [candidate - intrinsic for candidate, intrinsic in
                      zip(individual_model, individual_baseline)]
            rms = math.sqrt(statistics.mean(value * value for value in individual_baseline))
            result.update(error=distance, abm_profile=reference, pde_profile=candidate,
                          baseline_profile=independent,
                          paired_zero_test=paired_zero_test(excess),
                          paired_zero_test_contrast="individual L2(P,A) minus L2(B,A), exploratory",
                          individual_model_distances=individual_model,
                          baseline={"error": baseline_distance, "rms_distance": rms,
                                    "model_to_baseline_mean_ratio": distance / baseline_distance if baseline_distance else None,
                                    "model_to_baseline_rms_ratio": distance / rms if rms else None,
                                    "individual_distances": individual_baseline})
        else:
            values = {model: [sample[model][metric] for sample in samples]
                      for model in ("abm", "pde", "baseline")}
            relative = metric not in ("active_fraction", "lesion_perfused_fraction")
            result = paired_tost(values["abm"], values["pde"], margin, relative)
            scale = statistics.mean(values["abm"]) if relative else 1.0
            result.update(abm_mean=statistics.mean(values["abm"]),
                          pde_mean=statistics.mean(values["pde"]),
                          error=abs(result["difference"]) / scale if scale > 0.0 else None,
                          realizations=values,
                          baseline=baseline(values["abm"], values["baseline"], values["pde"], relative))
        metrics[metric] = result
    return metrics


def trajectory_comparison(samples, tolerances):
    times = [record["time_hours"] for record in samples[0]["abm"]["time_series"]]
    for sample in samples:
        for model in ("abm", "pde", "baseline"):
            actual = [record["time_hours"] for record in sample[model]["time_series"]]
            if len(actual) != len(times) or any(abs(a - b) > 1e-9 for a, b in zip(actual, times)):
                raise RuntimeError("ensembles have different sampling times")
    tables = []
    for i, time in enumerate(times):
        snapshots = [{model: sample[model]["time_series"][i]
                      for model in ("abm", "pde", "baseline")} for sample in samples]
        metrics = statistical_comparison(snapshots, tolerances)
        tables.append({"time_hours": time, "metrics": metrics,
                       "passed": all(result["passed"] for result in metrics.values())})
    maximum = {metric: max(({"time_hours": table["time_hours"],
                            "error": table["metrics"][metric]["error"]} for table in tables),
                           key=lambda result: result["error"] if result["error"] is not None else -math.inf)
               for metric in tolerances}
    return {"samples": tables, "maximum_errors": maximum,
            "passed": all(table["passed"] for table in tables)}


def boundary_summary(samples, geometry=None):
    result = {}
    for model in samples[0]:
        records = [record for sample in samples for record in
                   sample[model].get("time_series", [sample[model]])]
        if geometry:
            dimensions = 2 if geometry["shape"][2] == 1 else 3
            half_width = min(min(-geometry["origin"][axis], geometry["origin"][axis] +
                                 geometry["shape"][axis] * geometry["spacing_voxels"])
                             for axis in range(dimensions))
            records = [dict(record, r99_to_half_width=record["r99"] / half_width) for record in records]
        annotated = [record for record in records if "r99_to_half_width" in record]
        result[model] = {
            "annotated": len(annotated) == len(records),
            "maximum_r99_to_half_width": max((record["r99_to_half_width"] for record in annotated), default=None),
            "maximum_boundary_mass": max((record["boundary_mass"] for record in annotated
                                           if "boundary_mass" in record), default=None),
            "outer_face_mass_present": any(record.get("boundary_mass", 0.0) > 0.0 for record in annotated),
            "boundary_influenced": any(record["r99_to_half_width"] >= 0.8 for record in annotated),
        }
    return result


def statistical_markdown(report):
    lines = ["# Shared ABM/PDE statistical validation", "",
             f"{report['seeds']} paired realizations plus {report['seeds']} disjoint ABM baseline seeds. "
             f"Result: {'PASS' if report['passed'] else 'FAIL'}.", ""]
    if report.get("short_endpoint_case"):
        lines += ["Duration is shorter than the sampling interval: endpoint statistical "
                  "diagnostic only, without a resolved time-series claim.", ""]
    lines += [
             "| Hours | Metric | ABM mean | PDE mean | Error | Baseline error | Model/baseline | Paired t | Paired p | TOST lower p | TOST upper p | Equivalent |",
             "| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |"]
    def number(value):
        return "null" if value is None else format(value, ".17g")
    for table in report["trajectory"]["samples"]:
        for name, metric in table["metrics"].items():
            zero = metric.get("paired_zero_test", {})
            values = [table["time_hours"], name, metric.get("abm_mean"), metric.get("pde_mean"),
                      metric["error"], metric["baseline"]["error"],
                      metric["baseline"]["model_to_baseline_mean_ratio"], zero.get("t"), zero.get("p"),
                      metric["tost"]["lower"]["p"], metric["tost"]["upper"]["p"]]
            cells = [str(value) if isinstance(value, str) else number(value) for value in values]
            lines.append("| " + " | ".join(cells) + f" | {'PASS' if metric['passed'] else 'FAIL'} |")
    lines += ["", "Boundary diagnostics: " + json.dumps(report["boundary"], sort_keys=True), "",
              "A failed equivalence test is not proof of inequivalence. Review the fixed margin, "
              "confidence interval and coverage. The profile test is a nonlinear jackknife approximation.",
              "Raw arrays, baseline RMS comparisons and every test contrast are in validation_report.json."]
    return lines


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=root / "build-codex/atcg3d_shared_abm")
    parser.add_argument("--config", type=Path, default=root / "ATCG3D_SharedRules/config/validation_v7.yaml")
    parser.add_argument("--seeds", type=int, default=16)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--output", type=Path, default=root / "validation-output")
    parser.add_argument("--tolerances", type=Path, help="JSON object tightening named tolerances")
    parser.add_argument("--vascular", action="store_true")
    parser.add_argument("--minimum-active-fraction", type=float, default=0.0)
    parser.add_argument("--require-vascular-activity", action="store_true")
    parser.add_argument("--smoke-only", action="store_true", help="retain the published endpoint smoke criterion")
    parser.add_argument("--analysis-only", action="store_true", help="analyze existing realization reports")
    parser.add_argument("--sample-every-hours", type=float, default=4.0)
    parser.add_argument("--require-interior", action="store_true", help="reject boundary-influenced ensembles")
    args = parser.parse_args()
    if args.seeds < 16 or args.jobs < 1 or args.threads < 1:
        parser.error("use at least 16 seeds and positive jobs/threads")
    if not math.isfinite(args.sample_every_hours) or not 4.0 <= args.sample_every_hours <= 24.0:
        parser.error("sampling interval must be between four and 24 hours")
    if not 0.0 <= args.minimum_active_fraction <= 1.0:
        parser.error("minimum active fraction must be between zero and one")
    if args.smoke_only and args.require_interior:
        parser.error("interior checks require statistical time-series reports")
    executable, config = args.exe.resolve(), args.config.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    tolerances = dict(VASCULAR_TOLERANCES if args.vascular else TOLERANCES)
    if args.tolerances:
        overrides = json.loads(args.tolerances.read_text())
        if not isinstance(overrides, dict) or set(overrides) - set(tolerances) or any(
                not isinstance(value, (int, float)) or not math.isfinite(value) or
                value < 0.0 or value > tolerances[name] for name, value in overrides.items()):
            parser.error("tolerances may only tighten the declared finite nonnegative margins")
        tolerances.update(overrides)
    interval = 0.0 if args.smoke_only else args.sample_every_hours
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_pair, executable, config, seed, args.output, args.threads,
                               interval, 0 if args.smoke_only else args.seeds, args.analysis_only)
                   for seed in range(1, args.seeds + 1)]
        samples = [future.result() for future in futures]
    smoke = compare(samples, tolerances)
    coverage = {"active_fraction": {model: statistics.mean(sample[model]["active_fraction"] for sample in samples)
                                    for model in ("abm", "pde")}}
    coverage_passed = all(value >= args.minimum_active_fraction
                          for value in coverage["active_fraction"].values())
    if args.require_vascular_activity:
        coverage["vascular_activity"] = {
            "abm_roots": sum(sample["abm"]["vascular_roots"] for sample in samples),
            "abm_anastomoses": sum(sample["abm"]["vascular_anastomoses"] for sample in samples),
            "pde_branching_rate": statistics.mean(sample["pde"]["tip_branching_rate"] for sample in samples),
            "pde_anastomosis_rate": statistics.mean(sample["pde"]["tip_anastomosis_rate"] for sample in samples)}
        coverage_passed = coverage_passed and all(value > 0 for value in coverage["vascular_activity"].values())
    smoke_passed = all(metric["passed"] for metric in smoke.values()) and coverage_passed
    report = {"schema_version": 2, "seeds": args.seeds, "time_hours": samples[0]["abm"]["time_hours"],
              "config": config.name, "metrics": smoke, "coverage": coverage,
              "coverage_passed": coverage_passed, "smoke_passed": smoke_passed,
              "equivalence_margins": tolerances, "statistical_equivalence": None,
              "mode": "legacy_smoke" if args.smoke_only else "statistical_equivalence"}
    if args.smoke_only:
        passed = smoke_passed
        report["passed"] = passed
        lines = ["# Shared ABM/PDE legacy smoke check", "",
                 f"{args.seeds} paired seeds, endpoint {report['time_hours']:.17g} hours. "
                 f"Result: {'PASS' if passed else 'FAIL'}. No statistical equivalence claim.", "",
                 "| Metric | ABM mean | PDE mean | Error | Tolerance | Sampling allowance | Smoke |",
                 "| --- | ---: | ---: | ---: | ---: | ---: | --- |"]
        for name, metric in smoke.items():
            lines.append(f"| {name} | {metric.get('abm_mean', 0):.17g} | {metric.get('pde_mean', 0):.17g} | "
                         f"{metric['error']:.17g} | {metric['tolerance']:.17g} | "
                         f"{metric['sampling_allowance']:.17g} | {'PASS' if metric['passed'] else 'FAIL'} |")
        description = subprocess.run([str(executable), "--config", str(config), "--dry-run"],
                                     capture_output=True, text=True, check=False)
        if description.returncode:
            raise RuntimeError(description.stderr.strip())
        geometry = json.loads(description.stdout)["continuum_config"]["grid"]
        report["grid"] = geometry
        report["boundary"] = boundary_summary(samples, geometry)
        lines += ["", "Endpoint boundary diagnostics (face mass unavailable in legacy reports): " +
                  json.dumps(report["boundary"], sort_keys=True)]
    else:
        report["trajectory"] = trajectory_comparison(samples, tolerances)
        report["boundary"] = boundary_summary(samples)
        interior = all(value["annotated"] and not value["boundary_influenced"]
                       for value in report["boundary"].values())
        report["interior_passed"] = interior
        report["statistical_equivalence"] = report["trajectory"]["passed"]
        report["seed_groups"] = {"abm_a": list(range(1, args.seeds + 1)),
                                 "abm_b": list(range(args.seeds + 1, 2 * args.seeds + 1))}
        report["sample_every_hours"] = interval
        report["short_endpoint_case"] = report["time_hours"] < interval
        passed = report["statistical_equivalence"] and coverage_passed and (interior or not args.require_interior)
        report["passed"] = passed
        lines = statistical_markdown(report)
    (args.output / "validation_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    (args.output / "validation_report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
