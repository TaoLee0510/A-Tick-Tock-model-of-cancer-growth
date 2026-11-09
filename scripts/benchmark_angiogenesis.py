#!/usr/bin/env python3
"""Measure the bounded vascular field at one/eight threads in separate processes."""

import argparse
import json
from pathlib import Path
import platform
import statistics
import subprocess


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=root / "build-codex/atcg3d_angiogenesis_benchmark")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--edge", type=int, default=2000)
    parser.add_argument("--source-edge", type=int, default=128)
    parser.add_argument("--hours", type=float, default=24.0)
    parser.add_argument("--repetitions", type=int, default=3)
    args = parser.parse_args()
    if args.repetitions < 1:
        parser.error("use at least one repetition")
    records = []
    for repetition in range(args.repetitions):
        for threads in (1, 8):
            command = [str(args.exe.resolve()), "--edge", str(args.edge), "--threads", str(threads),
                       "--source-edge", str(args.source_edge), "--hours", str(args.hours)]
            result = subprocess.run(command, text=True, capture_output=True, check=True)
            record = json.loads(result.stdout)
            record["repetition"] = repetition + 1
            records.append(record)
            print(f"Completed {threads}-thread repetition {repetition + 1}", flush=True)
    deterministic = len({record["state_checksum"] for record in records}) == 1
    medians = {str(threads): statistics.median(record["wall_seconds"] for record in records
                                             if record["threads"] == threads) for threads in (1, 8)}
    report = {"schema_version": 1, "system": platform.system(), "release": platform.release(),
              "machine": platform.machine(), "records": records, "deterministic": deterministic,
              "median_wall_seconds": medians, "speedup_1_to_8": medians["1"] / medians["8"],
              "timing_scope": "field advance only; peak resident memory includes all buffers and inputs"}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return 0 if deterministic else 1


if __name__ == "__main__":
    raise SystemExit(main())
