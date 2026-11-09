"""Check structural initialization separately from post-initialization inference."""

import copy
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from validate_shared_angiogenesis import MARGINS, analyze


def samples():
    result = []
    for seed in range(16):
        pair = {}
        for model in ("abm", "pde", "baseline"):
            trace = []
            for time in (0.0, 24.0, 240.0):
                snapshot = {"time_hours": time, "grid_shape": [256, 256, 1],
                            "vascular_length_model": "unique_lattice_centerline_growth_v2"}
                for metric in MARGINS:
                    snapshot[metric] = 0.0 if time == 0.0 else (
                        0.3 if metric == "lesion_perfused_fraction" else 4.0 + 0.01 * seed)
                trace.append(snapshot)
            pair[model] = {"time_series": trace}
        result.append(pair)
    return result


base = samples()
report = analyze(base)
assert report["passed"]
assert report["samples"][0]["metrics"]["vascular_branches"]["initialization_identity"]
assert not report["samples"][0]["metrics"]["vascular_branches"]["inference_performed"]

post_zero = copy.deepcopy(base)
for pair in post_zero:
    for model in pair.values():
        model["time_series"][-1]["vascular_branches"] = 0.0
report = analyze(post_zero)
assert not report["passed"]
assert report["samples"][-1]["metrics"]["vascular_branches"]["undefined_reason"]

unequal_initial = copy.deepcopy(base)
unequal_initial[0]["pde"]["time_series"][0]["vascular_length"] = 1.0
assert not analyze(unequal_initial)["passed"]

for fault in ("time", "length", "grid"):
    broken = copy.deepcopy(base)
    snapshot = broken[0]["pde"]["time_series"][-1]
    if fault == "time":
        snapshot["time_hours"] = 239.0
    elif fault == "length":
        snapshot["vascular_length_model"] = "volume_proxy"
    else:
        snapshot["grid_shape"] = [128, 128, 1]
    rejected = False
    try:
        analyze(broken)
    except RuntimeError:
        rejected = True
    assert rejected
