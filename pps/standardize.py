"""
pps/standardize.py

Standardization of module inputs

    Z_k,i(t) = f_k( (X_k,i(t) - mu_k) / sigma_k )

f_k is a prespecified transformation (identity by default), mu_k is
the reference mean or median, and sigma_k is the reference standard
deviation or robust scale parameter. Skewed biomarkers (CRP, IL-6,
TNF-alpha, ferritin, etc.) may be log-transformed *before*
standardization -- i.e. before the (X - mu) / sigma step, on a mu/sigma
already estimated on the log scale.
"""
from dataclasses import dataclass
from typing import Callable, Optional
import math


@dataclass
class Standardizer:
    """One reusable standardization rule for a single raw biomarker.

    mu / sigma should be estimated on whatever scale the biomarker is
    standardized on (raw or log), i.e. after log_transform is applied.
    """
    mu: float
    sigma: float
    log_transform: bool = False
    f: Optional[Callable[[float], float]] = None  # optional f_k, applied to the z-score

    def __call__(self, x: float) -> float:
        if self.sigma == 0:
            raise ValueError("sigma must be non-zero.")
        value = math.log(x) if self.log_transform else x
        z = (value - self.mu) / self.sigma
        return self.f(z) if self.f is not None else z


def code_categorical(value, mapping: dict) -> float:
    """Prespecified coding for categorical/ordinal inputs (Section 3).

    e.g. inhibitor/inducer exposure coded by strength and expected
    direction of effect before being scaled onto the module axis:

        code_categorical("strong_inhibitor", {
            "none": 0.0, "weak_inhibitor": -0.3,
            "moderate_inhibitor": -0.6, "strong_inhibitor": -1.0,
        })
    """
    try:
        return mapping[value]
    except KeyError as e:
        raise ValueError(f"No prespecified code for {value!r}; expected one of {list(mapping)}") from e