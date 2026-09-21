"""
pps/modules.py

Modular architecture
static germline modules plus seven time-indexed dynamic modules.
"""
from dataclasses import dataclass
from typing import Dict, Optional


@dataclass
class GermlineProfile:
    """Section 4.1: the static, time-invariant baseline layer.

    Z_G1 -- primary pharmacogene module
            (e.g. CYP2D6, CYP2C19, CYP2C9, CYP3A5, TPMT, NUDT15, DPYD, UGT1A1)
    Z_G2 -- secondary pharmacogenomic module
            (transporters/targets, e.g. SLCO1B1, ABCB1, ABCC2, ABCC4, SLC22A-family)
    Z_G3 -- regulatory germline module
            (e.g. HNF4A, NR1I2/PXR, NR1I3/CAR, AHR, PPARA -- germline variation
            only; expression of these genes belongs to Z_TX, not here)
    """
    Z_G1: Optional[float] = None
    Z_G2: Optional[float] = None
    Z_G3: Optional[float] = None

    def as_dict(self) -> Dict[str, float]:
        out = {}
        if self.Z_G1 is not None:
            out["G1"] = self.Z_G1
        if self.Z_G2 is not None:
            out["G2"] = self.Z_G2
        if self.Z_G3 is not None:
            out["G3"] = self.Z_G3
        return out


@dataclass
class DynamicState:
    """Sections 4.2-4.3: the time-indexed layers at a single timepoint t.

    Z_XEN -- xenobiotic exposure (concomitant drugs, inhibitors/inducers, OTC,
             herbal products, alcohol, tobacco, dietary xenobiotics)
    Z_MIC -- microbiome (microbial activation/inactivation, beta-glucuronidase
             activity, bile acids, SCFAs, indoles)
    Z_IMM -- immune-inflammatory (CRP, IL-6, TNF-alpha, ESR, ferritin, albumin, NLR)
    Z_TX  -- transcriptomic (drug-metabolizing enzyme/transporter/target expression;
             residualize against Z_IMM and/or Z_G3 per Section 6 to avoid double counting)
    Z_PT  -- post-transcriptional/post-translational (miRNA, lncRNA, mRNA stability,
             protein turnover, PTMs -- may be exploratory/omitted)
    Z_HEP -- hepatic function (bilirubin, albumin, INR, ALT, AST, ALP, GGT, cholestasis)
    Z_REN -- renal function (eGFR, creatinine clearance, serum creatinine, cystatin C)

    A module left as None is treated as *not entered* for this
    timepoint (Section 6: prespecified module exclusion), not as zero
    -- silently defaulting an unmeasured module to zero would
    misrepresent missing data as "at reference."
    """
    Z_XEN: Optional[float] = None
    Z_MIC: Optional[float] = None
    Z_IMM: Optional[float] = None
    Z_TX: Optional[float] = None
    Z_PT: Optional[float] = None
    Z_HEP: Optional[float] = None
    Z_REN: Optional[float] = None

    def as_dict(self) -> Dict[str, float]:
        mapping = {
            "XEN": self.Z_XEN, "MIC": self.Z_MIC, "IMM": self.Z_IMM,
            "TX": self.Z_TX, "PT": self.Z_PT, "HEP": self.Z_HEP, "REN": self.Z_REN,
        }
        return {k: v for k, v in mapping.items() if v is not None}


def combine(germline: GermlineProfile, dynamic: DynamicState) -> Dict[str, float]:
    """Merge germline + dynamic modules into the z dict pps.score expects."""
    z = germline.as_dict()
    z.update(dynamic.as_dict())
    return z