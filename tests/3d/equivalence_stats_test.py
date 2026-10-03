"""Analytic distribution fixtures and adversarial equivalence decisions."""

import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from equivalence_stats import (baseline, jackknife_distance_tost, paired_tost,
                               student_quantile, student_survival)
from validate_abm_vs_pde import TOLERANCES, compare, trajectory_comparison
from ablate_abm_vs_pde import interaction_tables


def distribution_fixtures():
    # Cauchy (df=1) and the analytic df=2 CDF independently fix both tails.
    for value in (0.0, 0.1, 1.0, 3.0, 100.0):
        expected = 0.5 - math.atan(value) / math.pi
        assert abs(student_survival(value, 1) - expected) < 2.0e-14
        expected = 0.5 - value / (2.0 * math.sqrt(2.0 + value * value))
        assert abs(student_survival(value, 2) - expected) < 2.0e-14
        for df in (1, 2, 15, 31, 1000):
            assert abs(student_survival(value, df) + student_survival(-value, df) - 1.0) < 1.0e-14
    # NIST's tabulated df=15 and df=30 quantiles are rounded to three decimals.
    for probability, df, expected in ((0.95, 15, 1.753), (0.975, 15, 2.131),
                                      (0.95, 30, 1.697), (0.975, 30, 2.042)):
        assert abs(student_quantile(probability, df) - expected) < 0.0005
    assert abs(student_quantile(0.95, 1) - math.tan(0.45 * math.pi)) < 1.0e-12
    assert student_survival(20.0, 15) < 1.0e-11


def decisions():
    reference = [100.0] * 16
    assert paired_tost(reference, [105.0] * 16, 0.10, True)["passed"]
    at_bound = paired_tost(reference, [110.0] * 16, 0.10, True)
    assert not at_bound["passed"]
    assert at_bound["tost"]["upper"]["p"] == 1.0
    assert paired_tost([0.0] * 16, [0.0] * 16, 0.05, False)["passed"]
    zero_relative = paired_tost([0.0] * 16, [0.0] * 16, 0.75, True)
    assert not zero_relative["passed"]
    assert zero_relative["tost"]["lower"]["p"] is None
    assert "undefined_reason" in zero_relative
    assert baseline([0.0] * 16, [0.0] * 16, [0.0] * 16, True)["error"] is None
    outside = paired_tost(reference, [140.0 + (-1) ** i * 0.1 for i in range(16)], 0.35, True)
    assert not outside["passed"]
    candidate = [100.0 + (-1) ** i * 100.0 for i in range(16)]
    noisy = paired_tost(reference, candidate, 0.10, True)
    assert noisy["paired_zero_test"]["p"] == 1.0
    assert not noisy["passed"]
    smoke = compare([{"abm": {"total_mass": a}, "pde": {"total_mass": p}}
                     for a, p in zip(reference, candidate)], {"total_mass": 0.10})
    assert smoke["total_mass"]["passed"]
    significant = paired_tost(reference, [105.0 + (-1) ** i * 0.5 for i in range(16)], 0.10, True)
    assert significant["passed"]
    assert significant["paired_zero_test"]["p"] < 0.00001
    zero = baseline(reference, reference, [105.0] * 16, True)
    assert zero["model_to_baseline_mean_ratio"] is None
    assert zero["model_to_baseline_rms_ratio"] is None
    cancel = baseline(reference, [99.0, 101.0] * 8, [105.0] * 16, True)
    assert cancel["model_to_baseline_mean_ratio"] is None
    assert cancel["model_to_baseline_rms_ratio"] == 5.0
    json.dumps([noisy, zero, cancel], allow_nan=False)
    assert jackknife_distance_tost(0.0, [0.0] * 16, 0.35)["passed"]
    assert not jackknife_distance_tost(0.35, [0.35] * 16, 0.35)["passed"]
    assert not jackknife_distance_tost(0.30, [0.15, 0.45] * 8, 0.35)["passed"]


def trajectories():
    def snapshot(time, value):
        result = {name: value for name in TOLERANCES if name != "radial_profile_L2"}
        result.update(time_hours=time, active_fraction=0.0, radial_mass=[1.0, 0.0, 0.0, 0.0])
        return result
    samples = []
    for seed in range(16):
        samples.append({model: {"time_series": [snapshot(0.0, 100.0),
                         snapshot(4.0, 160.0 if model == "pde" else 100.0),
                         snapshot(8.0, 100.0)]} for model in ("abm", "pde", "baseline")})
    report = trajectory_comparison(samples, TOLERANCES)
    assert report["samples"][-1]["passed"]
    assert not report["passed"]
    assert report["maximum_errors"]["total_mass"]["time_hours"] == 4.0
    samples[0]["baseline"]["time_series"][1]["time_hours"] = 5.0
    try:
        trajectory_comparison(samples, TOLERANCES)
    except RuntimeError:
        pass
    else:
        assert False, "different sample times must be rejected"


def interventions():
    def realization(abm, pde):
        return {model: {"time_series": [{"time_hours": 4.0, "radial_mass": [1.0],
                 **{name: value for name in TOLERANCES if name != "radial_profile_L2"}}]}
                for model, value in (("abm", abm), ("pde", pde))}
    full = [realization(100.0, 110.0) for _ in range(16)]
    changed = [realization(95.0, 102.0) for _ in range(16)]
    values = interaction_tables({"full": full, "no_migration": changed})
    metric = values["no_migration"][0]["metrics"]["total_mass"]
    assert metric["abm_effect"] == -5.0
    assert metric["pde_effect"] == -8.0
    assert metric["interaction"] == -3.0
    assert metric["interaction_standard_error"] == 0.0
    assert metric["paired_zero_test"]["p"] == 0.0
    changed[0]["pde"]["time_series"][0]["time_hours"] = 5.0
    try:
        interaction_tables({"full": full, "no_migration": changed})
    except RuntimeError:
        pass
    else:
        assert False, "misaligned intervention times must be rejected"


distribution_fixtures()
decisions()
trajectories()
interventions()
