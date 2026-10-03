#!/usr/bin/env python3
"""Test invasive mixed representation with prespecified paired equivalence bounds."""

import argparse
import concurrent.futures
import copy
import json
import math
from pathlib import Path
import subprocess

from validate_abm_vs_pde import (TOLERANCES, boundary_summary, run_model,
                                statistical_markdown, trajectory_comparison)
from validate_hybrid import REFINEMENT_TOLERANCES, SCENARIOS


ADAPTIVE = {name: values for name, values in SCENARIOS.items() if values[0] == "adaptive"}


def validate_record(record, seed, expected_model="hybrid_invasion_front_v4"):
    if record.get("seed") != seed or record.get("model") != expected_model:
        raise RuntimeError("wrong hybrid realization identity")
    if record.get("mode") != "adaptive":
        raise RuntimeError("invasion validation requires adaptive mode")
    snapshots = record.get("time_series", [])
    if len(snapshots) < 3 or snapshots[0]["time_hours"] != 0.0:
        raise RuntimeError("invasion validation requires an initialization and resolved time series")
    previous = -1.0
    for snapshot in snapshots:
        time = snapshot["time_hours"]
        if not math.isfinite(time) or time <= previous:
            raise RuntimeError("invalid hybrid sampling times")
        previous = time
        if snapshot["K_mass"] <= 0.0:
            raise RuntimeError("hybrid r/K ratio requires positive K mass")
        for metric in TOLERANCES:
            if metric != "radial_profile_L2" and not math.isfinite(snapshot[metric]):
                raise RuntimeError("nonfinite hybrid metric")
        if any(not math.isfinite(value) or value < 0.0 for value in snapshot["radial_mass"]):
            raise RuntimeError("invalid hybrid radial profile")
        if not math.isclose(snapshot["abm_mass"] + snapshot["pde_mass"],
                            snapshot["total_mass"], abs_tol=1e-9):
            raise RuntimeError("hybrid representation masses do not sum")
    if snapshots[-1]["time_hours"] != record["time_hours"]:
        raise RuntimeError("hybrid trace does not reach its endpoint")


def coverage(record):
    extrema = {key: record[key] for key in ("minimum_abm_fraction", "minimum_active_fraction",
                                           "maximum_pde_fraction", "to_pde", "to_abm")}
    extrema["passed"] = (extrema["minimum_abm_fraction"] >= 0.20 and
                         extrema["minimum_active_fraction"] >= 0.10 and
                         extrema["maximum_pde_fraction"] >= 0.10 and
                         extrema["to_pde"] > 0 and extrema["to_abm"] > 0)
    extrema["time_series"] = [{key: snapshot[key] for key in (
        "time_hours", "abm_fraction", "pde_fraction", "active_fraction", "to_pde", "to_abm")}
        for snapshot in record["time_series"]]
    return extrema


def run_hybrid(executable, config, name, seed, output, threads, analysis_only):
    path = output / f"seed_{seed:04d}_{name}.json"
    if not analysis_only:
        _, step, exchange = ADAPTIVE[name]
        command = [str(executable), "--config", str(config), "--seed", str(seed),
                   "--threads", str(threads), "--no-output", "--validation-report",
                   "--step-hours", str(step), "--exchange-hours", str(exchange),
                   "--sample-every-hours", "4", "--report", str(path)]
        process = subprocess.run(command, capture_output=True, text=True, check=False)
        if process.returncode:
            raise RuntimeError(f"hybrid {name} seed {seed} failed: {process.stderr.strip()}")
    record = json.loads(path.read_text())
    validate_record(record, seed)
    return record


def run_seed(args, seed):
    sample = {}
    for name, model, offset in (("abm", "abm", 0), ("baseline", "abm", args.seeds),
                                ("pde", "pde", 0), ("pde_baseline", "pde", args.seeds)):
        sample[name] = run_model(args.native_exe, args.rules_config, model, seed + offset,
                                 args.output, args.threads, 4.0, args.analysis_only)
    for name in ADAPTIVE:
        sample[name] = run_hybrid(args.exe, args.config, name, seed, args.output,
                                  args.threads, args.analysis_only)
    for name in ("fine_step", "fine_exchange"):
        sample[name + "_baseline"] = run_hybrid(args.exe, args.config, name, seed + args.seeds,
                                                args.output, args.threads, args.analysis_only)
    return sample


def comparison(samples, reference, candidate, baseline, margins):
    records = [{"abm": sample[reference], "pde": sample[candidate],
                "baseline": sample[baseline]} for sample in samples]
    return trajectory_comparison(records, margins)


def analyze(samples):
    primary = comparison(samples, "abm", "coarse_step", "baseline", TOLERANCES)
    secondary = comparison(samples, "pde", "coarse_step", "pde_baseline", TOLERANCES)
    native = comparison(samples, "abm", "pde", "baseline", TOLERANCES)
    refinements = {}
    for name, reference, candidate, independent in (
            ("coarse_step", "fine_step", "coarse_step", "fine_step_baseline"),
            ("middle_step", "fine_step", "middle_step", "fine_step_baseline"),
            ("coarse_exchange", "fine_exchange", "coarse_exchange", "fine_exchange_baseline"),
            ("middle_exchange", "fine_exchange", "middle_step", "fine_exchange_baseline")):
        refinements[name] = comparison(samples, reference, candidate, independent, REFINEMENT_TOLERANCES)
    representation = []
    for sample in samples:
        for name in (*ADAPTIVE, "fine_step_baseline", "fine_exchange_baseline"):
            record = sample[name]
            representation.append(dict(coverage(record), seed=record["seed"], scenario=name))
    boundary = boundary_summary(samples)
    interior = all(value["annotated"] and not value["boundary_influenced"] for value in boundary.values())
    return {"primary_hybrid_vs_abm": primary, "secondary_hybrid_vs_pde": secondary,
            "native_pde_vs_abm": native, "refinements": refinements,
            "representation_coverage": representation, "boundary": boundary,
            "passed": primary["passed"] and secondary["passed"] and native["passed"] and
                      all(value["passed"] for value in refinements.values()) and interior and
                      all(value["passed"] for value in representation)}


def markdown(report):
    lines = ["# Invasion-front hybrid equivalence", "",
             "Result: " + ("PASS" if report["passed"] else "FAIL") + ".",
             "This 8-hour case measures early r20 invasion, not long-time production accuracy.", ""]
    tables = {name: report[name] for name in ("primary_hybrid_vs_abm", "secondary_hybrid_vs_pde",
                                              "native_pde_vs_abm")}
    tables.update({"refinement_" + name: table for name, table in report["refinements"].items()})
    for name, table in tables.items():
        lines += ["## " + name, "",
                  "Column ABM denotes the declared reference; PDE denotes the candidate in this table.", ""]
        wrapper = {"seeds": report["seeds"], "passed": table["passed"], "trajectory": table,
                   "boundary": report["boundary"]}
        lines += statistical_markdown(wrapper)[4:]
    lines += ["", "## Mixed coverage", "",
              "| Seed | Scenario | Min ABM fraction | Min active fraction | Max PDE fraction | To PDE | To ABM | Result |",
              "| ---: | --- | ---: | ---: | ---: | ---: | ---: | --- |"]
    for record in report["representation_coverage"]:
        values = [str(record["seed"]), record["scenario"]] + [format(record[key], ".17g") for key in (
            "minimum_abm_fraction", "minimum_active_fraction", "maximum_pde_fraction", "to_pde", "to_abm")]
        lines.append("| " + " | ".join(values) + " | " + ("PASS" if record["passed"] else "FAIL") + " |")
    return "\n".join(lines) + "\n"


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=root / "build-codex/atcg3d_hybrid")
    parser.add_argument("--native-exe", type=Path, default=root / "build-codex/atcg3d_shared_abm")
    parser.add_argument("--config", type=Path, default=root / "ATCG3D_Hybrid/config/hybrid_invasion_r20_v4.yaml")
    parser.add_argument("--rules-config", type=Path, default=root / "ATCG3D_SharedRules/config/active_r20_v14.yaml")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", type=int, default=16)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--analysis-only", action="store_true")
    args = parser.parse_args()
    if args.seeds < 16 or args.jobs < 1 or args.threads < 1:
        parser.error("use at least 16 seeds and positive jobs/threads")
    descriptions = []
    for executable, config in ((args.exe, args.config), (args.native_exe, args.rules_config)):
        process = subprocess.run([str(executable.resolve()), "--config", str(config.resolve()),
                                  "--threads", str(args.threads), "--dry-run"],
                                 text=True, capture_output=True, check=True)
        descriptions.append(json.loads(process.stdout))
    if descriptions[0]["model"] != "hybrid_invasion_front_v4":
        parser.error("use the new invasion-front policy")
    mixed, native = copy.deepcopy(descriptions[0]["rules_config"]), copy.deepcopy(descriptions[1])
    for description in (mixed, native):
        description["continuum_config"].pop("output")
    if mixed != native:
        parser.error("hybrid and native rule configurations differ")
    if min(mixed["continuum_config"]["grid"]["shape"][:2]) < 256:
        parser.error("invasion validation requires at least a 256-square grid")
    args.output.mkdir(parents=True, exist_ok=True)
    args.exe, args.native_exe = args.exe.resolve(), args.native_exe.resolve()
    args.config, args.rules_config = args.config.resolve(), args.rules_config.resolve()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        samples = list(pool.map(lambda seed: run_seed(args, seed), range(1, args.seeds + 1)))
    report = analyze(samples)
    report.update(schema_version=1, seeds=args.seeds, protocol="hybrid_invasion_contract.md",
                  config=args.config.name, rules_config=args.rules_config.name,
                  seed_groups={"a": list(range(1, args.seeds + 1)),
                               "b": list(range(args.seeds + 1, 2 * args.seeds + 1))},
                  margins=TOLERANCES, refinement_margins=REFINEMENT_TOLERANCES,
                  sample_every_hours=4.0)
    (args.output / "validation_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    (args.output / "validation_report.md").write_text(markdown(report))
    print("Hybrid invasion equivalence: " + ("PASS" if report["passed"] else "FAIL"))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
