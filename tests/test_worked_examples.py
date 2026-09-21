"""
Regression tests reproducing, digit-for-digit, the three worked
examples in Section 8 ("Worked example") of Supplementary Methods 2.

These pin down the two-step Section 2 calculation (latent score ->
logistic transform) against numbers already published in the
manuscript, so any future refactor that silently changes the math
will fail loudly here.
"""
import pytest

from pps.score import PPSCoefficients, latent_score, pps_from_latent, compute_pps


def test_example_1_unit_coefficients():
    # "a standardized primary genetic module score of -0.50 and a
    # standardized inflammatory module score of -0.80. With unit
    # coefficients used only for illustration:
    #   L = (1.0 x -0.50) + (1.0 x -0.80) = -1.30 ; PPS = 0.214"
    z = {"G1": -0.50, "IMM": -0.80}
    L = latent_score(z)
    assert L == pytest.approx(-1.30, abs=1e-9)
    assert pps_from_latent(L) == pytest.approx(0.214, abs=1e-3)
    assert compute_pps(z) == pytest.approx(0.214, abs=1e-3)


def test_example_2_calibrated_inflammatory_coefficient():
    # "If subsequent clinical calibration estimated theta_IMM = 0.65:
    #   L = -0.50 + (0.65 x -0.80) = -1.02 ; PPS = 0.265"
    z = {"G1": -0.50, "IMM": -0.80}
    theta = PPSCoefficients(weights={"IMM": 0.65})  # G1 stays at the default unit weight
    L = latent_score(z, theta)
    assert L == pytest.approx(-1.02, abs=1e-9)
    assert compute_pps(z, theta) == pytest.approx(0.265, abs=1e-3)


def test_example_3_second_patient_xenobiotic_and_mild_inflammation():
    # "the same standardized primary genetic module score of -0.50, a
    # standardized xenobiotic module score of +0.30, and mild
    # inflammation of -0.10:
    #   L = -0.50 + 0.30 - 0.10 = -0.30 ; PPS = 0.426"
    z = {"G1": -0.50, "XEN": 0.30, "IMM": -0.10}
    L = latent_score(z)
    assert L == pytest.approx(-0.30, abs=1e-9)
    assert compute_pps(z) == pytest.approx(0.426, abs=1e-3)


def test_pps_bounded_and_monotonic_in_L():
    # Sanity check on Eq. 2 independent of the worked examples: PPS is
    # strictly in (0, 1) and increasing in L.
    Ls = [-5.0, -1.0, 0.0, 1.0, 5.0]
    values = [pps_from_latent(L) for L in Ls]
    assert all(0.0 < v < 1.0 for v in values)
    assert values == sorted(values)
    assert pps_from_latent(0.0) == pytest.approx(0.5)


def test_unrecognized_module_name_rejected():
    with pytest.raises(ValueError):
        latent_score({"NOT_A_MODULE": 1.0})