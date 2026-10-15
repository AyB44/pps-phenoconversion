"""
examples/thiopurine_simulation.py

Thiopurine sensitivity analyses: TPMT and NUDT15 allele variants x
simulated clinical conditions. Everything printed is computed from the
inputs below; no value is typed in by hand.

THESE ARE SENSITIVITY ANALYSES, NOT FINDINGS. The coding of genotype and
CRP into Z-scores, the module weights, and the pathway fractions f and g
are ASSUMPTIONS (see the parameter table printed at the top). None is
estimated from patient data, and the outputs are not predictions or
dosing guidance.

Definition used throughout: PPS_TPMT is the relative capacity of TPMT-
mediated inactivation of thiopurines. Lower PPS_TPMT means MORE drug
goes down the thioguanine-nucleotide route, i.e. HIGHER erythrocyte
6-TGN. (A lower PPS therefore does not mean "less active metabolite".)
PPS_NUDT15 is the relative capacity to degrade active thiopurine
triphosphates; lower means higher DNA-level exposure, which erythrocyte
6-TGN does not capture.

Run from the repository root:  python examples/thiopurine_simulation.py
"""
import math
from pps import PPSCoefficients, compute_pps

# ---------------------------------------------------------------------
# 1. Allele function (CPIC): 1 = normal function, 0 = no function.
#    Uncertain-function alleles (e.g. NUDT15*4, TPMT*6) are excluded.
# ---------------------------------------------------------------------
TPMT_ALLELES = {"*1": 1, "*2": 0, "*3A": 0, "*3B": 0, "*3C": 0}
NUDT15_ALLELES = {"*1": 1, "*2": 0, "*3": 0}


def activity_score(diplotype, table):
    a, b = diplotype.split("/")
    return table[a] + table[b]


def cpic_phenotype(score):
    return {2: "NM", 1: "IM", 0: "PM"}[score]


# ASSUMED coding of activity score -> Z_G1 (reference NM = 0), chosen so
# that PPS is about 1.0 / 0.5 / 0.05 for NM / IM / PM.
Z_GENOTYPE = {2: 0.0, 1: -1.1, 0: -3.7}

# ASSUMED coding of CRP -> Z_IMM: log10 fold above the reference value.
CRP_REF = 5.0  # mg/L


def z_imm(crp):
    return -math.log10(crp / CRP_REF) if crp > CRP_REF else 0.0


# ---------------------------------------------------------------------
# 2. PPS per pathway.
#    TPMT pathway  : Z_G1 = TPMT genotype; Z_IMM = CRP with weight
#                    theta_imm (direction UNPROVEN; theta_imm = 0 means
#                    CRP has no effect on TPMT capacity).
#    NUDT15 pathway: genotype only; no dynamic modifier is entered.
#    Co-medications are NOT entered (no validated effect size).
# ---------------------------------------------------------------------

def pps_tpmt(tpmt_diplotype, crp, theta_imm):
    z = {"G1": Z_GENOTYPE[activity_score(tpmt_diplotype, TPMT_ALLELES)],
         "IMM": z_imm(crp)}
    return compute_pps(z, PPSCoefficients(weights={"G1": 1.0, "IMM": theta_imm}))


def pps_nudt15(nudt15_diplotype):
    z = {"G1": Z_GENOTYPE[activity_score(nudt15_diplotype, NUDT15_ALLELES)]}
    return compute_pps(z, PPSCoefficients(weights={"G1": 1.0}))


# ---------------------------------------------------------------------
# 3. Readouts: PPS as a Vmax multiplier, Km constant, linear range,
#    constant dose and adherence.
#    f = share of competing thiopurine clearance via TPMT at reference.
#    g = share of active-nucleotide degradation via NUDT15 at reference.
# ---------------------------------------------------------------------

def rel_6tgn(pps, f):
    return 1.0 / ((1.0 - f) + f * pps)


def rel_dna_exposure(pps, g):
    return 1.0 / ((1.0 - g) + g * pps)


F_VALUES = (0.5, 0.8, 0.9)
THETAS = (0.0, 0.5, 1.0)


def table(title, header, rows):
    print(f"\n{title}")
    print("| " + " | ".join(header) + " |")
    print("|" + "|".join("---" for _ in header) + "|")
    for r in rows:
        print("| " + " | ".join(str(c) for c in r) + " |")


def parameter_table():
    table("Parameters and their status",
          ["Parameter", "Value", "Status"],
          [["No-function alleles", "TPMT *2, *3A, *3B, *3C; NUDT15 *2, *3", "CPIC"],
           ["Excluded alleles", "uncertain function (e.g. NUDT15*4, TPMT*6)", "CPIC classification"],
           ["Activity score -> Z_G1", "2: 0, 1: -1.1, 0: -3.7", "assumed (gives PPS ~ 1 / 0.5 / 0.05)"],
           ["Z_IMM", "-log10(CRP / 5 mg/L), above 5 mg/L only", "assumed"],
           ["theta_IMM", "0 / 0.5 / 1.0", "assumed; direction unproven (manuscript)"],
           ["f (TPMT share of competing clearance)", "0.5 / 0.8 / 0.9", "assumed, unknown"],
           ["g (NUDT15 share of nucleotide degradation)", "0.5 / 0.8 / 0.9", "assumed, unknown"],
           ["CRP for Table 1 scenarios B-E", "12 / 38 / 72 / 72 mg/L", "assumed"],
           ["Table 2 genotypes, CRP, co-medication", "as in the manuscript", "manuscript"],
           ["Co-medications", "not entered into PPS", "no validated effect size"]])


def main():
    parameter_table()

    # --- A. TPMT genotype grid, no clinical modifier -------------------
    rows = []
    for dip in ["*1/*1", "*1/*2", "*1/*3A", "*1/*3C", "*3A/*3A", "*3A/*3C", "*3C/*3C"]:
        s = activity_score(dip, TPMT_ALLELES)
        p = pps_tpmt(dip, crp=0.0, theta_imm=0.0)
        rows.append([dip, s, cpic_phenotype(s), f"{p:.3f}"] +
                    [f"{rel_6tgn(p, f):.2f}" for f in F_VALUES])
    table("A. TPMT allele variants (no clinical modifier; sensitivity to f)",
          ["TPMT diplotype", "Activity score", "CPIC phenotype", "PPS_TPMT"] +
          [f"rel. 6-TGN (f={f})" for f in F_VALUES], rows)

    # --- B. NUDT15 genotype grid ---------------------------------------
    rows = []
    for dip in ["*1/*1", "*1/*2", "*1/*3", "*2/*3", "*3/*3"]:
        s = activity_score(dip, NUDT15_ALLELES)
        p = pps_nudt15(dip)
        rows.append([dip, s, cpic_phenotype(s), f"{p:.3f}"] +
                    [f"{rel_dna_exposure(p, g):.2f}" for g in F_VALUES])
    table("B. NUDT15 allele variants (no clinical modifier; sensitivity to g)",
          ["NUDT15 diplotype", "Activity score", "CPIC phenotype", "PPS_NUDT15"] +
          [f"rel. DNA-level exposure (g={g})" for g in F_VALUES], rows)

    # --- C. Inflammation sensitivity -----------------------------------
    rows = []
    for dip in ["*1/*1", "*1/*3A"]:
        for crp in (2, 12, 38, 72):
            row = [dip, crp, f"{z_imm(crp):.2f}"]
            for th in THETAS:
                p = pps_tpmt(dip, crp, th)
                row.append(f"{p:.3f} ({rel_6tgn(p, 0.8):.2f}x)")
            rows.append(row)
    table("C. Inflammation sensitivity, TPMT pathway: PPS_TPMT (rel. 6-TGN at f=0.8)",
          ["TPMT diplotype", "CRP mg/L", "Z_IMM"] +
          [f"theta_IMM={th}" for th in THETAS], rows)

    # --- D. Manuscript Table 1 scenarios, all three theta_IMM ----------
    t1 = [
        ("A", "*1/*1", "*1/*1", 2, "none; stable disease"),
        ("B", "*1/*1", "*1/*1", 12, "mild flare"),
        ("C", "*1/*1", "*1/*1", 38, "moderate flare"),
        ("D", "*1/*3A", "*1/*1", 72, "severe flare, admission"),
        ("E", "*1/*3A", "*1/*3", 72, "severe flare, interacting co-medication"),
    ]
    rows = []
    for sc, t, n, crp, ctx in t1:
        row = [sc, t, n, crp, ctx]
        for th in THETAS:
            p = pps_tpmt(t, crp, th)
            row.append(f"{p:.3f} ({rel_6tgn(p, 0.8):.2f}x)")
        pn = pps_nudt15(n)
        row.append(f"{pn:.3f} ({rel_dna_exposure(pn, 0.8):.2f}x)")
        rows.append(row)
    table("D. Manuscript Table 1 scenarios: PPS_TPMT (rel. 6-TGN) at f=0.8, PPS_NUDT15 (rel. DNA exposure) at g=0.8",
          ["Scenario", "TPMT", "NUDT15", "CRP", "Context"] +
          [f"PPS_TPMT, theta_IMM={th}" for th in THETAS] + ["PPS_NUDT15"], rows)

    # --- E. Manuscript Table 2 profiles: primary theta_IMM = 0 ---------
    t2 = [
        (1, "*1/*1", "*1/*1", 2, "none"),
        (2, "*1/*1", "*1/*1", 12, "none"),
        (3, "*1/*1", "*1/*1", 38, "mesalazine"),
        (4, "*1/*3A", "*1/*1", 16, "mesalazine"),
        (5, "*1/*1", "*1/*3", 72, "trimethoprim-sulfamethoxazole"),
        (6, "*3A/*3A", "*1/*1", 3, "allopurinol"),
    ]
    rows = []
    for pid, t, n, crp, comed in t2:
        p0, p5 = pps_tpmt(t, crp, 0.0), pps_tpmt(t, crp, 0.5)
        pn = pps_nudt15(n)
        rows.append([pid, t, n, crp, comed,
                     f"{p0:.3f} ({rel_6tgn(p0, 0.8):.2f}x)",
                     f"{pn:.3f} ({rel_dna_exposure(pn, 0.8):.2f}x)",
                     f"{p5:.3f} ({rel_6tgn(p5, 0.8):.2f}x)"])
    table("E. Manuscript Table 2 profiles. Primary: theta_IMM=0 (CRP is not a measure of TPMT/NUDT15 "
          "function). Co-medication is an alert only and is not modelled.",
          ["Patient", "TPMT", "NUDT15", "CRP", "Co-medication alert",
           "PPS_TPMT (rel. 6-TGN, f=0.8), theta_IMM=0",
           "PPS_NUDT15 (rel. DNA exposure, g=0.8)",
           "Sensitivity: PPS_TPMT (rel. 6-TGN), theta_IMM=0.5"], rows)


if __name__ == "__main__":
    main()