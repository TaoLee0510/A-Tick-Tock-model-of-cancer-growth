"""Paired Student tests with fixed equivalence margins and no optional dependency."""

import math
import statistics


ALPHA = 0.05


def beta_fraction(a, b, x):
    """Evaluate the incomplete-beta continued fraction using modified Lentz."""
    tiny = 1.0e-300
    c = 1.0
    d = 1.0 - (a + b) * x / (a + 1.0)
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    result = d
    for m in range(1, 401):
        for coefficient in (
            m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m)),
            -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1)),
        ):
            d = 1.0 + coefficient * d
            if abs(d) < tiny:
                d = tiny
            c = 1.0 + coefficient / c
            if abs(c) < tiny:
                c = tiny
            d = 1.0 / d
            change = d * c
            result *= change
        if abs(change - 1.0) < 3.0e-14:
            return result
    raise ArithmeticError("incomplete-beta continued fraction did not converge")


def regularized_beta(x, a, b):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    factor = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) +
                      a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return factor * beta_fraction(a, b, x) / a
    return 1.0 - factor * beta_fraction(b, a, 1.0 - x) / b


def student_survival(t, degrees_of_freedom):
    if degrees_of_freedom < 1 or not math.isfinite(t):
        raise ValueError("finite t and positive degrees of freedom are required")
    if t == 0.0:
        return 0.5
    df = float(degrees_of_freedom)
    tail = 0.5 * regularized_beta(df / (df + t * t), df / 2.0, 0.5)
    return tail if t > 0.0 else 1.0 - tail


def student_quantile(probability, degrees_of_freedom):
    if not 0.5 <= probability < 1.0:
        raise ValueError("positive Student quantiles require probability in [0.5, 1)")
    lower, upper = 0.0, 1.0
    while student_survival(upper, degrees_of_freedom) > 1.0 - probability:
        upper *= 2.0
    for _ in range(80):
        middle = (lower + upper) / 2.0
        if student_survival(middle, degrees_of_freedom) > 1.0 - probability:
            lower = middle
        else:
            upper = middle
    return (lower + upper) / 2.0


def standard_error(values):
    if len(values) < 2 or any(not math.isfinite(value) for value in values):
        raise ValueError("at least two finite observations are required")
    return statistics.stdev(values) / math.sqrt(len(values))


def one_sided(mean, error, df, alternative):
    if error == 0.0:
        favorable = mean > 0.0 if alternative == "greater" else mean < 0.0
        return {"mean_contrast": mean, "standard_error": error, "t": None,
                "p": 0.0 if favorable else 1.0, "deterministic": True}
    t = mean / error
    p = student_survival(t if alternative == "greater" else -t, df)
    return {"mean_contrast": mean, "standard_error": error, "t": t,
            "p": p, "deterministic": False}


def paired_zero_test(differences):
    mean = statistics.mean(differences)
    error = standard_error(differences)
    if error == 0.0:
        return {"t": None, "p": 1.0 if mean == 0.0 else 0.0,
                "deterministic": True}
    t = mean / error
    return {"t": t, "p": min(1.0, 2.0 * student_survival(abs(t), len(differences) - 1)),
            "deterministic": False}


def paired_tost(reference, candidate, margin, relative):
    if len(reference) != len(candidate) or margin < 0.0 or not math.isfinite(margin):
        raise ValueError("paired observations and a finite nonnegative margin are required")
    differences = [p - a for a, p in zip(reference, candidate)]
    error = standard_error(differences)
    df = len(differences) - 1
    if relative:
        if statistics.mean(reference) <= 0.0:
            return {"degrees_of_freedom": df, "alpha": ALPHA, "margin": margin,
                    "relative": relative, "difference": statistics.mean(differences),
                    "standard_error": error, "paired_zero_test": paired_zero_test(differences),
                    "tost": {"lower": {"p": None}, "upper": {"p": None}},
                    "undefined_reason": "relative equivalence requires a positive reference mean",
                    "passed": False}
        lower = [difference + margin * a for a, difference in zip(reference, differences)]
        upper = [difference - margin * a for a, difference in zip(reference, differences)]
    else:
        lower = [value + margin for value in differences]
        upper = [value - margin for value in differences]
    tests = {
        "lower": one_sided(statistics.mean(lower), standard_error(lower), df, "greater"),
        "upper": one_sided(statistics.mean(upper), standard_error(upper), df, "less"),
    }
    critical = student_quantile(1.0 - ALPHA, df)
    mean = statistics.mean(differences)
    return {"degrees_of_freedom": df, "alpha": ALPHA, "margin": margin,
            "relative": relative, "difference": mean, "standard_error": error,
            "difference_ci90": [mean - critical * error, mean + critical * error],
            "paired_zero_test": paired_zero_test(differences), "tost": tests,
            "passed": all(test["p"] < ALPHA for test in tests.values())}


def jackknife_distance_tost(estimate, delete_one, margin):
    n = len(delete_one)
    if n < 3 or any(not math.isfinite(value) for value in [estimate, *delete_one]):
        raise ValueError("a finite distance and at least three jackknife estimates are required")
    center = statistics.mean(delete_one)
    error = math.sqrt((n - 1) / n * sum((value - center) ** 2 for value in delete_one))
    critical = student_quantile(1.0 - ALPHA, n - 1)
    tests = {"lower": one_sided(estimate + margin, error, n - 1, "greater"),
             "upper": one_sided(estimate - margin, error, n - 1, "less")}
    return {"degrees_of_freedom": n - 1, "alpha": ALPHA, "margin": margin,
            "distance": estimate, "standard_error": error,
            "distance_ci90": [estimate - critical * error, estimate + critical * error],
            "delete_one_distances": delete_one, "tost": tests,
            "method": "nonlinear delete-one-pair jackknife t approximation",
            "passed": all(test["p"] < ALPHA for test in tests.values())}


def baseline(reference, independent, candidate, relative):
    if len(reference) != len(independent) or len(reference) != len(candidate):
        raise ValueError("baseline and candidate ensembles must have the same size")
    differences = [b - a for a, b in zip(reference, independent)]
    candidate_difference = statistics.mean(p - a for a, p in zip(reference, candidate))
    mean = statistics.mean(differences)
    rms = math.sqrt(statistics.mean(value * value for value in differences))
    scale = statistics.mean(reference) if relative else 1.0
    return {"group_a_mean": statistics.mean(reference),
            "group_b_mean": statistics.mean(independent),
            "difference": mean, "error": abs(mean) / scale if scale > 0.0 else None,
            "standard_error": standard_error(differences),
            "rms_difference": rms, "paired_zero_test": paired_zero_test(differences),
            "model_to_baseline_mean_ratio": abs(candidate_difference / mean) if mean else None,
            "model_to_baseline_rms_ratio": abs(candidate_difference / rms) if rms else None,
            "zero_baseline_mean": mean == 0.0, "zero_baseline_rms": rms == 0.0}
