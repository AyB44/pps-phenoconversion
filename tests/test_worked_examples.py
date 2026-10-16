"""
Tests for the corrected PPS definition: PPS = 2 * sigmoid(L).

Sections:
  1. Reference, inhibition and induction behaviour of the corrected PPS
  2. Bounds and numerical safety
  3. Optional ceiling generalization
  4. PK step (Vmax scaling)
  5. Supplementary Methods 2 worked examples
  6. A documented limitation of the default ceiling
"""
import math
import pytest

from pps.score import (
    PPSCoefficients, latent_score, pps_from_latent, compute_pps,
    PPS_CEILING_DEFAULT,
    pps_from_latent_logistic_deprecated, compute_pps_logistic_deprecated,
)
from pps.pk import MichaelisMentenParams, mm_rate


# --- 1. Reference, inhibition, induction ----------------------------------

def test_reference_patient_has_pps_one():
    # The canary for the original bug: sigmoid(0) was 0.5 here.
    assert compute_pps({}) == pytest.approx(1.0)
    assert compute_pps({"G1": 0.0, "IMM": 0.0}) == pytest.approx(1.0)
    assert pps_from_latent(0.0) == pytest.approx(1.0)


def test_inhibition_gives_pps_below_one():
    assert compute_pps({"G1": -0.5, "XEN": -0.6}) < 1.0


def test_induction_gives_pps_above_one():
    theta = PPSCoefficients(weights={"XEN": 1.0})
    assert compute_pps({"XEN": 0.8}, theta) > 1.0


@pytest.mark.parametrize("L, expected", [
    (-2.0, 0.238),   # values from the correction report's table
    (0.0, 1.000),
    (1.0, 1.462),
    (2.0, 1.762),
])
def test_report_table_values(L, expected):
    assert pps_from_latent(L) == pytest.approx(expected, abs=1e-3)


def test_pps_equals_two_sigmoid():
    for L in (-3.0, -0.7, 0.0, 0.4, 2.5):
        assert pps_from_latent(L) == pytest.approx(2.0 / (1.0 + math.exp(-L)))


def test_unrecognized_module_name_rejected():
    with pytest.raises(ValueError):
        latent_score({"NOT_A_MODULE": 1.0})


# --- 2. Bounds and numerical safety ---------------------------------------

def test_pps_within_zero_two_and_monotonic():
    Ls = [-6.0, -2.0, -0.5, 0.0, 0.5, 2.0, 6.0]
    vals = [pps_from_latent(L) for L in Ls]
    assert all(0.0 < v < 2.0 for v in vals)
    assert vals == sorted(vals)


def test_extreme_latent_scores_do_not_overflow():
    assert pps_from_latent(-800.0) == pytest.approx(0.0, abs=1e-12)
    assert pps_from_latent(800.0) == pytest.approx(2.0)
    assert 0.0 <= pps_from_latent(-800.0) <= 2.0


# --- 3. Optional ceiling ---------------------------------------------------

def test_default_ceiling_is_two():
    assert PPS_CEILING_DEFAULT == 2.0


@pytest.mark.parametrize("c", [1.5, 2.0, 3.0, 4.0])
def test_reference_is_one_for_any_ceiling(c):
    assert pps_from_latent(0.0, ceiling=c) == pytest.approx(1.0)


def test_ceiling_two_matches_default():
    assert pps_from_latent(0.9, ceiling=2.0) == pytest.approx(pps_from_latent(0.9))


def test_higher_ceiling_allows_larger_values():
    assert pps_from_latent(6.0, ceiling=4.0) > 2.0


@pytest.mark.parametrize("bad", [1.0, 0.5, 0.0, -2.0])
def test_invalid_ceiling_rejected(bad):
    with pytest.raises(ValueError):
        pps_from_latent(0.0, ceiling=bad)


# --- 4. PK step ------------------------------------------------------------

def test_pk_reference_preserves_vmax_and_rate():
    mm = MichaelisMentenParams(Vmax_ref=10.0, Km_ref=2.0)
    scaled = mm.scaled(compute_pps({}))
    assert scaled.Vmax_ref == pytest.approx(10.0)
    assert mm_rate(3.0, scaled) == pytest.approx(mm_rate(3.0, mm))


def test_pk_uses_pps_as_vmax_multiplier_and_keeps_km():
    mm = MichaelisMentenParams(Vmax_ref=10.0, Km_ref=2.0)
    low = mm.scaled(0.55)
    high = mm.scaled(1.462)
    assert low.Vmax_ref == pytest.approx(5.5)
    assert high.Vmax_ref == pytest.approx(14.62)
    assert low.Km_ref == high.Km_ref == 2.0


def test_pk_accepts_agreed_range_rejects_invalid():
    mm = MichaelisMentenParams(Vmax_ref=10.0, Km_ref=2.0)
    mm.scaled(0.0)
    mm.scaled(2.0)
    with pytest.raises(ValueError):
        mm.scaled(-0.01)
    with pytest.raises(ValueError):
        mm.scaled(2.01)


def test_pk_custom_ceiling_range():
    mm = MichaelisMentenParams(Vmax_ref=10.0, Km_ref=2.0)
    assert mm.scaled(3.0, pps_max=4.0).Vmax_ref == pytest.approx(30.0)
    with pytest.raises(ValueError):
        mm.scaled(4.5, pps_max=4.0)


# --- 5. Supplementary Methods 2 worked examples:  ------

EXAMPLES = [
    ("example_1", {"G1": -0.50, "IMM": -0.80}, None, -1.30, 0.214, 0.428),
    ("example_2", {"G1": -0.50, "IMM": -0.80}, PPSCoefficients(weights={"IMM": 0.65}), -1.02, 0.265, 0.530),
    ("example_3", {"G1": -0.50, "XEN": 0.30, "IMM": -0.10}, None, -0.30, 0.426, 0.851),
]


@pytest.mark.parametrize("name, z, theta, L, old_pps, new_pps", EXAMPLES)
def test_worked_examples_old_and_new(name, z, theta, L, old_pps, new_pps):
    assert latent_score(z, theta) == pytest.approx(L, abs=1e-9)
    assert compute_pps_logistic_deprecated(z, theta) == pytest.approx(old_pps, abs=1e-3)
    assert compute_pps(z, theta) == pytest.approx(new_pps, abs=1e-3)


def test_deprecated_reference_point_gives_point_five():
    # Permanent record of the original inconsistency.
    assert compute_pps_logistic_deprecated({}) == pytest.approx(0.5)


# --- 6. Documented limitation of the default ceiling ------------------------

def test_strongest_observed_inducer_exceeds_default_ceiling():
    # CYP3A5 cohort Group E, patient P22: CLind 4.3 L/h vs Group A mean
    # 2.12 L/h is a ratio of about 2.03, which is above the 2.0 limit.
    # With the default ceiling this patient cannot be matched exactly,
    # and ratios near 2 sit on the flat part of the curve.
    ratio = 4.3 / 2.12
    assert ratio > PPS_CEILING_DEFAULT
    assert pps_from_latent(40.0) <= PPS_CEILING_DEFAULT