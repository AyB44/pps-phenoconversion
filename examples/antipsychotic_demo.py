"""
examples/antipsychotic_demo.py

Illustration of the PPS calculation on the CYP3A5*3/*3 antipsychotic
(quetiapine) cohort, using ONLY values published in the open-access
paper (Kohalmy, Borthakur, Porrogi, Pharmaceutics 2026, 18(7), 782).

WHAT THIS IS: a worked example showing that the calculation can
represent both inhibition (PPS < 1) and induction (PPS > 1).

WHAT THIS IS NOT: a test of predictive validity.
  - The codings below were chosen with the observed outcomes known, so
    agreement in ordering is built in and is not evidence.
  - Anchors mix group means, subset means and single patients (see the
    'Unit' column). No per-patient data are used; those are in the
    paper's Supplementary File.
  - Module weights are unit weights. Nothing is fitted.

Run from the repository root:  python examples/antipsychotic_demo.py
"""
import math
from pps import PPSCoefficients, compute_pps, PPS_CEILING_DEFAULT

CL_REF = 2.12      # Group A mean CLind, L/h (cohort paper)
CRP_REF = 3.0      # mg/L, reference CRP in the cohort (cohort paper)

# ASSUMED codings (illustrative)
Z_GENOTYPE = {"gNM": 0.0, "gPM": -1.0}
Z_DDI = {"none": 0.0, "crowding": -0.5, "strong_inhibitor": -1.0, "strong_inducer": 1.0}


def z_imm(crp):
    return -math.log10(crp / CRP_REF) if crp > CRP_REF else 0.0


# label, unit, genotype, DDI class, CRP (mg/L), observed CLind (L/h)
# CRP for B and E was not reported for the anchor; it is set to the
# reference value (assumed). D uses the mean of the two reported CRPs.
ANCHORS = [
    ("A  concordant reference",    "group mean, n=6",     "gNM", "none",             3.0,  2.12),
    ("B  crowding, *3/*3",         "subset mean, n=3 of 8", "gPM", "crowding",       3.0,  0.53),
    ("C  inflammation",            "single patient P7",   "gNM", "none",            55.0,  0.70),
    ("D  strong inhibitor, *3/*3", "mean of 2 (P9, P18)", "gPM", "strong_inhibitor", 5.5,  0.475),
    ("E  strong inducer",          "single patient P22",  "gNM", "strong_inducer",   3.0,  4.30),
]


def parameter_table():
    rows = [
        ("Reference clearance", "2.12 L/h", "cohort paper (Group A mean)"),
        ("Reference CRP", "3 mg/L", "cohort paper"),
        ("D: CRP 5.5 mg/L", "mean of 6 and 5", "cohort paper"),
        ("B, E: CRP", "set to reference (not reported)", "assumed"),
        ("Genotype coding", "*3/*3 = -1", "assumed"),
        ("DDI scale", "-1 to +1", "range from cohort paper; placement of strata assumed"),
        ("Crowding", "-0.5", "assumed"),
        ("Z_IMM", "-log10(CRP / 3 mg/L), above 3 only", "assumed"),
        ("Weights", "all 1.0", "assumed"),
    ]
    print("| Parameter | Value | Status |")
    print("|---|---|---|")
    for r in rows:
        print("| " + " | ".join(r) + " |")


def main():
    parameter_table()
    theta = PPSCoefficients()  # unit weights, intercept 0

    print("\n| Anchor | Unit | L | PPS | Predicted CLind (L/h) | Observed CLind (L/h) | "
          "Predicted - observed (L/h) | Obs/ref | PPS - obs/ref |")
    print("|---|---|---|---|---|---|---|---|---|")
    for label, unit, g, ddi, crp, cl in ANCHORS:
        z = {"G1": Z_GENOTYPE[g], "XEN": Z_DDI[ddi], "IMM": z_imm(crp)}
        L = sum(z.values())
        p = compute_pps(z, theta)
        ratio = cl / CL_REF
        pred = p * CL_REF
        print(f"| {label} | {unit} | {L:.2f} | {p:.3f} | {pred:.2f} | {cl:g} | "
              f"{pred - cl:+.2f} | {ratio:.2f} | {p - ratio:+.2f} |")

    ratio_e = 4.30 / CL_REF
    p_e = compute_pps({"G1": 0.0, "XEN": Z_DDI["strong_inducer"], "IMM": 0.0}, theta)
    print(f"\nAnchor E: predicted ratio {p_e:.2f} vs observed {ratio_e:.2f}. The shortfall "
          f"reflects the assumed coding (inducer = +1, unit weight). Separately, observed "
          f"ratios above the PPS ceiling of {PPS_CEILING_DEFAULT} are outside the model's range.")


if __name__ == "__main__":
    main()