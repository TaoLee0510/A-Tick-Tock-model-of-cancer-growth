"""Invasion evidence must reject missing trajectories and endpoint-only coverage."""

import copy
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from validate_hybrid_invasion import ADAPTIVE, analyze, coverage, validate_record


def record(seed):
    snapshot = {"time_hours": 0.0, "total_mass": 100.0, "abm_mass": 80.0, "pde_mass": 20.0,
                "r_mass": 50.0, "K_mass": 50.0, "r_K_ratio": 1.0, "active_fraction": 0.3,
                "r50": 5.0, "r90": 10.0, "r99": 12.0, "radial_mass": [0.0, 50.0, 50.0],
                "abm_fraction": 0.8, "pde_fraction": 0.2, "minimum_abm_fraction": 0.8,
                "minimum_active_fraction": 0.3, "maximum_pde_fraction": 0.2,
                "to_pde": 10, "to_abm": 5, "spatial_dimensions": 2, "grid_shape": [256, 256, 1],
                "domain_half_width": 128.0, "r99_to_half_width": 12.0 / 128.0,
                "boundary_mass": 0.0, "model": "hybrid_invasion_front_v4", "mode": "adaptive",
                "seed": seed}
    trace = [dict(snapshot, time_hours=time) for time in (0.0, 4.0, 8.0)]
    return dict(trace[-1], time_series=trace)


samples = []
for seed in range(1, 17):
    value = record(seed)
    validate_record(value, seed)
    assert coverage(value)["passed"]
    sample = {name: copy.deepcopy(value) for name in ("abm", "baseline", "pde", "pde_baseline", *ADAPTIVE,
                                                      "fine_step_baseline", "fine_exchange_baseline")}
    samples.append(sample)
assert analyze(samples)["passed"]
for key, invalid in (("minimum_abm_fraction", 0.19), ("minimum_active_fraction", 0.09),
                     ("maximum_pde_fraction", 0.09), ("to_pde", 0), ("to_abm", 0)):
    altered = record(1)
    altered[key] = invalid
    assert not coverage(altered)["passed"]
    bad = copy.deepcopy(samples)
    bad[0]["coarse_step"] = altered
    assert not analyze(bad)["passed"]
for change in (lambda value: value.pop("time_series"),
               lambda value: value.update(seed=17),
               lambda value: value.update(mode="all_abm"),
               lambda value: value["time_series"][1].update(time_hours=0.0)):
    altered = record(1)
    change(altered)
    try:
        validate_record(altered, 1)
    except RuntimeError:
        pass
    else:
        raise AssertionError("invalid invasion record accepted")
bad = copy.deepcopy(samples)
bad[0]["coarse_step"]["time_series"][1]["r99_to_half_width"] = 0.8
assert not analyze(bad)["passed"]
