"""
pps/residualize.py

Handling correlation and avoiding double counting. Two residualization cases are named
explicitly:

    Z_TX,resid(t)            = Z_TX(t) - E[Z_TX(t) | Z_IMM(t)]
    Z_TX,regulatory_resid(t) = Z_TX,regulatory(t) - E[Z_TX,regulatory(t) | Z_G3]

E[.|.] is a conditional expectation to be estimated from data. This
module provides a simple OLS estimator for the common single-upstream-
predictor case, plus the residualization step itself once that
conditional-expectation model is available.

prespecified module exclusion, and penalized/hierarchical modelling of correlated
predictors  are modelling choices made at the pps.score call site
(by simply omitting a module) or in an external calibration pipeline,
not represented as separate functions here.
"""
from typing import List, Sequence, Tuple


def fit_conditional_mean(y: Sequence[float], given: Sequence[float]) -> Tuple[float, float]:
    """Estimate E[y | given] as a simple OLS line: y ~ intercept + slope * given.

    Returns (slope, intercept).
    """
    n = len(y)
    if n != len(given) or n < 2:
        raise ValueError("y and given must be equal-length sequences with n >= 2.")
    mean_x = sum(given) / n
    mean_y = sum(y) / n
    cov = sum((g - mean_x) * (v - mean_y) for g, v in zip(given, y))
    var = sum((g - mean_x) ** 2 for g in given)
    if var == 0:
        raise ValueError("given has zero variance; cannot estimate a slope.")
    slope = cov / var
    intercept = mean_y - slope * mean_x
    return slope, intercept


def residualize(y: Sequence[float], given: Sequence[float]) -> Tuple[float, float, List[float]]:
    """Fit E[y | given] by OLS and return (slope, intercept, residuals),
    where residuals[i] = y[i] - (intercept + slope * given[i]).
    """
    slope, intercept = fit_conditional_mean(y, given)
    residuals = [v - (intercept + slope * g) for g, v in zip(given, y)]
    return slope, intercept, residuals


def residualize_single(y: float, given: float, slope: float, intercept: float) -> float:
    """Apply an already-estimated conditional-mean model to one new observation.

    e.g. Z_TX_resid = residualize_single(Z_TX, Z_IMM, slope, intercept)
    """
    return y - (intercept + slope * given)