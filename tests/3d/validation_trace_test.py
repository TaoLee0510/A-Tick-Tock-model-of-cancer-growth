"""Sampling must expose exact times without changing either model trajectory."""

import argparse
import json
from pathlib import Path
import subprocess
import tempfile


def run(executable, config, output, model, threads, sampled, extra=()):
    command = [str(executable), "--config", str(config), "--model", model,
               "--threads", str(threads), "--report", str(output)]
    if sampled:
        command += ["--sample-every-hours", "4"]
    command += list(extra)
    process = subprocess.run(command, capture_output=True, text=True, check=False)
    assert process.returncode == 0, process.stderr
    return json.loads(output.read_text())


parser = argparse.ArgumentParser()
parser.add_argument("--exe", type=Path, required=True)
parser.add_argument("--root", type=Path, required=True)
parser.add_argument("--checkpoint", action="store_true")
args = parser.parse_args()
source = args.root / "ATCG3D_SharedRules/config"
with tempfile.TemporaryDirectory(prefix="atcg-trace-") as temporary:
    directory = Path(temporary)
    continuum = (source / "continuum_regular_cycle_v10.yaml").read_text()
    continuum = continuum.replace("end_time_hours: 48.0", "end_time_hours: 8.0")
    continuum = continuum.replace("../../configs/shared_resource_regular_cycle_v1.yaml",
                                  str(args.root / "configs/shared_resource_regular_cycle_v1.yaml"))
    (directory / "continuum.yaml").write_text(continuum)
    structured = (source / "regular_cycle_v13.yaml").read_text()
    structured = structured.replace("continuum_regular_cycle_v10.yaml", "continuum.yaml")
    config = directory / "structured.yaml"
    config.write_text(structured)
    for model in ("abm", "pde"):
        plain = run(args.exe, config, directory / f"{model}-plain.json", model, 1, False)
        sampled = run(args.exe, config, directory / f"{model}-sampled.json", model, 1, True)
        parallel = run(args.exe, config, directory / f"{model}-parallel.json", model, 4, True)
        assert all(sampled[key] == value for key, value in plain.items())
        assert sampled == parallel
        assert [record["time_hours"] for record in sampled["time_series"]] == [0.0, 4.0, 8.0]
        assert sampled["spatial_dimensions"] == 2
        assert sampled["domain_half_width"] == 128.0
        assert sampled["boundary_mass"] == 0.0
    rejected = subprocess.run([str(args.exe), "--config", str(config), "--sample-every-hours", "2"],
                              capture_output=True, text=True, check=False)
    assert rejected.returncode != 0
    if args.checkpoint:
        config.write_text((source / "regular_cycle_v14.yaml").read_text().replace(
            "continuum_regular_cycle_v10.yaml", "continuum.yaml"))
        checkpoint = directory / "cells.h5"
        stopped = run(args.exe, config, directory / "stopped.json", "abm", 1, True,
                      ["--checkpoint", str(checkpoint)])
        (directory / "continuum.yaml").write_text(continuum.replace("end_time_hours: 8.0", "end_time_hours: 12.0"))
        resumed = run(args.exe, config, directory / "resumed.json", "abm", 4, True,
                      ["--resume-checkpoint", str(checkpoint)])
        full = run(args.exe, config, directory / "full.json", "abm", 1, True)
        assert [record["time_hours"] for record in resumed["time_series"]] == [8.0, 12.0]
        assert resumed["time_series"][0] == stopped["time_series"][-1]
        assert resumed["time_series"] == full["time_series"][-2:]
        assert all(value == full[key] for key, value in resumed.items() if key != "time_series")
