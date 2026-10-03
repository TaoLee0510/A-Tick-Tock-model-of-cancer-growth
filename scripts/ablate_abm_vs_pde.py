#!/usr/bin/env python3
"""Attribute paired model differences to prespecified operator interventions."""

import argparse
import concurrent.futures
import json
import os
from pathlib import Path
import re
import statistics
import subprocess

from equivalence_stats import paired_zero_test, standard_error
from validate_abm_vs_pde import (TOLERANCES, boundary_summary, profile_distance,
                                run_pair, trajectory_comparison)


CASES = ("full", "no_migration", "no_activation", "no_division", "no_exchange", "enlarged_boundary")
CORRECTION_CASES = ("legacy_growth_mean", "axial_normal_transport", "unconditional_lattice")


def configurations(executable, source, output, cases=CASES):
    process = subprocess.run([str(executable), "--config", str(source), "--dry-run"],
                             capture_output=True, text=True, check=False)
    if process.returncode:
        raise RuntimeError(process.stderr.strip())
    config = json.loads(process.stdout)
    if config["version"] < 14:
        raise ValueError("operator interventions require schema v14")
    structured = source.read_text()
    if re.search(r"^operators:", structured, re.MULTILINE):
        raise ValueError("intervention source must use all published operators")
    match = re.search(r"^continuum_config:\s*(\S+)\s*$", structured, re.MULTILINE)
    if not match:
        raise ValueError("a scalar continuum_config reference is required")
    continuum_source = (source.parent / match.group(1)).resolve()
    continuum = continuum_source.read_text()
    base_match = re.search(r"^base_config:\s*(\S+)\s*$", continuum, re.MULTILINE)
    if not base_match:
        raise ValueError("a scalar base_config reference is required")
    base_source = (continuum_source.parent / base_match.group(1)).resolve()
    result = {}
    for case in cases:
        directory = output / case
        directory.mkdir(parents=True, exist_ok=True)
        field = continuum.replace(base_match.group(1), os.path.relpath(base_source, directory), 1)
        if case == "enlarged_boundary":
            geometry = config["continuum_config"]["grid"]
            shape = geometry["shape"]
            origin = geometry["origin"]
            dimensions = 2 if shape[2] == 1 else 3
            enlarged = [2 * shape[axis] if axis < dimensions else shape[axis] for axis in range(3)]
            shifted = [origin[axis] - shape[axis] * geometry["spacing_voxels"] / 2.0
                       if axis < dimensions else origin[axis] for axis in range(3)]
            field = re.sub(r"shape:\s*\[[^]]+\]", "shape: " + json.dumps(enlarged), field, count=1)
            field = re.sub(r"origin:\s*\[[^]]+\]", "origin: " + json.dumps(shifted), field, count=1)
        (directory / "continuum.yaml").write_text(field)
        rules = structured.replace(match.group(1), "continuum.yaml", 1)
        switches = {name: case != "no_" + name for name in ("migration", "activation", "division", "exchange")}
        if case == "legacy_growth_mean":
            rules = rules.replace("truncated_normal_expectation_v2", "clipped_location_mean_v1")
        elif case == "axial_normal_transport":
            rules = rules.replace("feasible_fixed_lattice_jump_v3", "axial_diffusion_v1")
        elif case == "unconditional_lattice":
            rules = rules.replace("feasible_fixed_lattice_jump_v3", "fixed_lattice_jump_v2")
        rules += "\noperators:\n  model: shared_operator_switches_v1\n"
        rules += "".join(f"  {name}: {'true' if enabled else 'false'}\n" for name, enabled in switches.items())
        target = directory / "structured.yaml"
        target.write_text(rules)
        result[case] = target
    return result


def interaction_tables(cases):
    full = cases["full"]
    times = [snapshot["time_hours"] for snapshot in full[0]["abm"]["time_series"]]
    report = {}
    for name, samples in cases.items():
        if len(samples) != len(full):
            raise RuntimeError("interventions have different ensemble sizes")
        for sample in samples:
            for model in ("abm", "pde"):
                actual = [snapshot["time_hours"] for snapshot in sample[model]["time_series"]]
                if len(actual) != len(times) or any(abs(a - b) > 1e-9 for a, b in zip(actual, times)):
                    raise RuntimeError("interventions have different sampling times")
        if name == "full":
            continue
        tables = []
        for index, time in enumerate(times):
            current = [{model: sample[model]["time_series"][index] for model in ("abm", "pde")}
                       for sample in samples]
            reference = [{model: sample[model]["time_series"][index] for model in ("abm", "pde")}
                         for sample in full]
            metrics = {}
            for metric in TOLERANCES:
                if metric == "radial_profile_L2":
                    before = profile_distance(reference)[0]
                    after = profile_distance(current)[0]
                    metrics[metric] = {"full_error": before, "intervention_error": after,
                                       "error_change": after - before}
                    continue
                effects = {model: [current[seed][model][metric] - reference[seed][model][metric]
                                   for seed in range(len(samples))] for model in ("abm", "pde")}
                interactions = [p - a for a, p in zip(effects["abm"], effects["pde"])]
                metrics[metric] = {"abm_effect": statistics.mean(effects["abm"]),
                                   "pde_effect": statistics.mean(effects["pde"]),
                                   "interaction": statistics.mean(interactions),
                                   "interaction_standard_error": standard_error(interactions),
                                   "paired_zero_test": paired_zero_test(interactions),
                                   "individual_interactions": interactions}
            tables.append({"time_hours": time, "metrics": metrics})
        report[name] = tables
    return report


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=root / "build-codex/atcg3d_shared_abm")
    parser.add_argument("--config", type=Path, default=root / "ATCG3D_SharedRules/config/activation_ablation_v14.yaml")
    parser.add_argument("--output", type=Path, default=root / "ablation-output")
    parser.add_argument("--seeds", type=int, default=16)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--analysis-only", action="store_true")
    parser.add_argument("--correction-study", action="store_true",
                        help="run closure reversions using a previously completed full intervention")
    parser.add_argument("--interventions", nargs="+", choices=CASES,
                        help="run a declared subset including the full reference")
    args = parser.parse_args()
    if args.seeds < 16 or args.jobs < 1 or args.threads < 1:
        parser.error("use at least 16 seeds and positive jobs/threads")
    if args.interventions and ("full" not in args.interventions or args.correction_study):
        parser.error("an intervention subset must include full and cannot be a correction study")
    args.output.mkdir(parents=True, exist_ok=True)
    selected = ("full", *CORRECTION_CASES) if args.correction_study else args.interventions or CASES
    configs = configurations(args.exe.resolve(), args.config.resolve(), args.output, selected)
    cases = {}
    for name, config in configs.items():
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(run_pair, args.exe.resolve(), config.resolve(), seed,
                                   args.output / name, args.threads, 4.0, args.seeds,
                                   args.analysis_only or (args.correction_study and name == "full"))
                       for seed in range(1, args.seeds + 1)]
            cases[name] = [future.result() for future in futures]
        print(f"Completed intervention {name}", flush=True)
    interactions = interaction_tables(cases)
    equivalence = {name: trajectory_comparison(samples, TOLERANCES) for name, samples in cases.items()}
    report = {"schema_version": 1, "seeds": args.seeds, "config": args.config.name,
              "interventions": interactions, "equivalence": equivalence,
              "boundary": {name: boundary_summary(samples) for name, samples in cases.items()},
              "seed_groups": {"abm_a": list(range(1, args.seeds + 1)),
                              "abm_b": list(range(args.seeds + 1, 2 * args.seeds + 1))},
              "sample_every_hours": 4.0, "equivalence_margins": TOLERANCES,
              "interpretation": "matched-seed one-at-a-time effects; interactions need not add; completion is not equivalence"}
    (args.output / ("correction_report.json" if args.correction_study else "ablation_report.json")).write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    lines = ["# Operator interventions", "", report["interpretation"], "",
             "| Intervention | Hours | Metric | ABM effect | PDE effect | Interaction/error change | SE | Paired p |",
             "| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |"]
    for name, tables in interactions.items():
        for table in tables:
            for metric, values in table["metrics"].items():
                cells = [values.get("abm_effect"), values.get("pde_effect"),
                         values.get("interaction", values.get("error_change")),
                         values.get("interaction_standard_error"), values.get("paired_zero_test", {}).get("p")]
                numbers = ["null" if value is None else format(value, ".17g") for value in cells]
                lines.append(f"| {name} | {table['time_hours']:.17g} | {metric} | " + " | ".join(numbers) + " |")
    (args.output / ("correction_report.md" if args.correction_study else "ablation_report.md")).write_text("\n".join(lines) + "\n")
    print("Equivalence by intervention: " + json.dumps({name: value["passed"] for name, value in equivalence.items()}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
