"""Hybrid tracing retains native limits and reproduces adaptive restart samples."""

import argparse
import json
from pathlib import Path
import subprocess
import tempfile


def run(executable, config, target, arguments):
    process = subprocess.run([str(executable), "--config", str(config), "--threads", "1",
                              "--report", str(target), *arguments], text=True, capture_output=True)
    assert process.returncode == 0, process.stderr
    return json.loads(target.read_text())


parser = argparse.ArgumentParser()
parser.add_argument("--exe", type=Path, required=True)
parser.add_argument("--native-exe", type=Path, required=True)
parser.add_argument("--root", type=Path, required=True)
args = parser.parse_args()
with tempfile.TemporaryDirectory(prefix="atcg-hybrid-invasion-") as temporary:
    folder = Path(temporary)
    source = args.root / "ATCG3D_SharedRules/config"
    base = (args.root / "configs/shared_activation_r20_v1.yaml").read_text()
    base = base.replace("r_cells: 500", "r_cells: 32").replace("K_cells: 500", "K_cells: 32")
    base = base.replace("outer_radius: 24", "outer_radius: 9").replace("shell_inner_radius: 18", "shell_inner_radius: 6")
    base = base.replace("inner_radius: 18", "inner_radius: 6")
    base = base.replace("inner_small_radius: 22", "inner_small_radius: 5")
    (folder / "base.yaml").write_text(base)
    continuum = (source / "continuum_active_r20_v8.yaml").read_text().replace(
        "../../configs/shared_activation_r20_v1.yaml", "base.yaml")
    (folder / "continuum.yaml").write_text(continuum)
    structured = (source / "active_r20_v14.yaml").read_text().replace("continuum_active_r20_v8.yaml", "continuum.yaml")
    (folder / "structured.yaml").write_text(structured)
    hybrid = (args.root / "ATCG3D_Hybrid/config/hybrid_invasion_r20_v4.yaml").read_text().replace(
        "../../ATCG3D_SharedRules/config/active_r20_v14.yaml", "structured.yaml")
    config = folder / "hybrid.yaml"
    config.write_text(hybrid)
    for mode, native_model in (("all_abm", "abm"), ("all_pde", "pde")):
        native = run(args.native_exe, folder / "structured.yaml", folder / f"native-{mode}.json",
                     ["--model", native_model, "--sample-every-hours", "4"])
        mixed = run(args.exe, config, folder / f"limit-{mode}.json",
                    ["--mode", mode, "--no-output", "--sample-every-hours", "4"])
        assert [value["time_hours"] for value in mixed["time_series"]] == [0.0, 4.0, 8.0]
        assert [value["state_checksum"] for value in mixed["time_series"]] == [
            value["state_checksum"] for value in native["time_series"]]
        for candidate, reference in zip(mixed["time_series"], native["time_series"]):
            for key in ("total_mass", "r_mass", "K_mass", "active_fraction", "radial_mass"):
                assert candidate[key] == reference[key], (mode, key)
    sampled = run(args.exe, config, folder / "sampled.json", ["--no-output", "--sample-every-hours", "4"])
    plain = run(args.exe, config, folder / "plain.json", ["--no-output"])
    assert all(sampled[key] == value for key, value in plain.items())
    (folder / "continuum.yaml").write_text(continuum.replace("end_time_hours: 8.0", "end_time_hours: 4.0"))
    checkpoint = folder / "saved.bin"
    stopped = run(args.exe, config, folder / "stopped.json",
                  ["--no-output", "--sample-every-hours", "4", "--checkpoint", str(checkpoint)])
    (folder / "continuum.yaml").write_text(continuum)
    resumed = run(args.exe, config, folder / "resumed.json",
                  ["--no-output", "--sample-every-hours", "4", "--threads", "4",
                   "--resume-checkpoint", str(checkpoint)])
    assert resumed["time_series"][0] == stopped["time_series"][-1]
    assert resumed["time_series"] == sampled["time_series"][-2:]
