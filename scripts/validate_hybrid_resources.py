#!/usr/bin/env python3
"""Validate the current resource-initialized hybrid against the invasive ABM."""

import argparse
import concurrent.futures
import copy
import json
from pathlib import Path
import subprocess

from validate_abm_vs_pde import (boundary_summary, run_model, statistical_markdown,
                                TOLERANCES, trajectory_comparison)
from validate_hybrid_invasion import coverage, validate_record


def realization(args, seed):
    reference = run_model(args.native_exe, args.reference_config, "abm", seed,
                          args.reference_directory, 1, 4.0,
                          (args.reference_directory / f"seed_{seed:04d}_abm.json").exists())
    other = seed + args.seeds
    baseline = run_model(args.native_exe, args.reference_config, "abm", other,
                         args.reference_directory, 1, 4.0,
                         (args.reference_directory / f"seed_{other:04d}_abm.json").exists())
    output = args.output / f"seed_{seed:04d}_hybrid_v5.json"
    subprocess.run([str(args.exe.resolve()), "--config", str(args.config.resolve()),
                    "--seed", str(seed), "--threads", "1", "--no-output", "--validation-report",
                    "--step-hours", "0.25", "--exchange-hours", "1", "--sample-every-hours", "4",
                    "--report", str(output.resolve())], check=True, capture_output=True, text=True)
    candidate = json.loads(output.read_text())
    validate_record(candidate, seed, "hybrid_resource_restart_v5")
    return {"abm": reference, "baseline": baseline, "pde": candidate}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--native-exe", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=root / "docs/recommended/hybrid/model.yaml")
    parser.add_argument("--reference-config", type=Path,
                        default=root / "ATCG3D_SharedRules/config/active_r20_v14.yaml")
    parser.add_argument("--reference-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", type=int, default=16)
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()
    if args.seeds < 16 or args.jobs < 1:
        parser.error("use two groups of at least 16 seeds and a positive worker count")
    descriptions = []
    for executable, config in ((args.exe, args.config), (args.native_exe, args.reference_config)):
        process = subprocess.run([str(executable.resolve()), "--config", str(config.resolve()),
                                  "--threads", "1", "--dry-run"], capture_output=True, text=True, check=True)
        descriptions.append(json.loads(process.stdout))
    if descriptions[0]["model"] != "hybrid_resource_restart_v5":
        parser.error("the resource-initialized v5 model is required")
    candidate = copy.deepcopy(descriptions[0]["rules_config"])
    reference = copy.deepcopy(descriptions[1])
    for description in (candidate, reference):
        description.pop("version")
        sector = description.pop("sector_mean_model", "published_sector_sums_v1")
        if sector != "published_sector_sums_v1":
            parser.error("the reference comparison requires published sector sums")
        description["continuum_config"].pop("output")
        description["continuum_config"].pop("base_config_path")
        description["continuum_config"]["base"].pop("output_directory")
    if candidate != reference:
        parser.error("hybrid and ABM must select the same live numerical operators and parameters")
    args.output.mkdir(parents=True, exist_ok=True)
    args.reference_directory.mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        samples = list(pool.map(lambda seed: realization(args, seed), range(1, args.seeds + 1)))
    comparison = trajectory_comparison(samples, TOLERANCES)
    guards = {str(seed): coverage(sample["pde"]) for seed, sample in enumerate(samples, 1)}
    boundaries = boundary_summary(samples)
    passed = (comparison["passed"] and all(value["passed"] for value in guards.values()) and
              all(value["annotated"] and not value["boundary_influenced"]
                  for value in boundaries.values()))
    report = {"protocol": "statistical_validation.md", "model": "hybrid_resource_restart_v5",
              "seeds": args.seeds, "seed_groups": {"a": list(range(1, args.seeds + 1)),
                                                     "b": list(range(args.seeds + 1, 2 * args.seeds + 1))},
              "primary_hybrid_vs_abm": comparison, "coverage": guards,
              "boundaries": boundaries, "passed": passed,
              "native_records": samples}
    (args.output / "validation_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    markdown = statistical_markdown({"seeds": args.seeds, "passed": passed,
                                     "trajectory": comparison, "boundary": boundaries})
    (args.output / "validation_report.md").write_text(
        "# Resource-initialized hybrid v5\n\n" + "\n".join(markdown) + "\n")
    print("Resource hybrid equivalence: " + ("PASS" if passed else "FAIL"), flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
