# pps-phenoconversion

A Python reference implementation of the **Phenotype-Predictive Score (PPS)** framework, following the exact module structure and equations defined in *Supplementary Methods 2* of:


PPS is presented in the source manuscript as a **hypothesis-generating, not-yet-prospectively-validated modeling framework** for estimating a patient's time-indexed functional pharmacological capacity from germline pharmacogenomic and dynamic (inflammatory, xenobiotic, microbiome, transcriptomic, hepatic, renal) inputs. This repository implements the mathematics of that framework — it is **not** a validated clinical prediction tool, and should not be used to guide dosing decisions.

## Installation

```bash
git clone https://github.com/AyB44/pps-phenoconversion.git
cd pps-phenoconversion
pip install -e .
```

## Quick start

```python
import pps

z = {"G1": -0.50, "IMM": -0.80}   # primary pharmacogene, immune-inflammatory
pps.compute_pps(z)
# -> 0.214...
```

## What's implemented

| Module | File | Section (Supp. Methods 2) |
|---|---|---|
| Latent score + logistic transform (Eqs. 1–2) | `pps/score.py` | Section 2 |
| Standardization of raw biomarkers | `pps/standardize.py` | Section 3 |
| Germline (G1–G3) and dynamic (XEN, MIC, IMM, TX, PT, HEP, REN) modules | `pps/modules.py` | Section 4 |
| Optional temporal smoothing | `pps/temporal.py` | Section 5 |
| Double-counting / residualization | `pps/residualize.py` | Section 6 |
| PK application (Vmax scaling) | `pps/pk.py` | Section 9 |

## Testing

```bash
pip install pytest
pytest tests/ -v
```

The test suite reproduces the three worked numerical examples from Supplementary Methods Section 8 exactly, so the implementation's correctness can be checked independently of this repository.

## Status

Reference implementation of the model's *mathematics only*. Coefficients (θ) used in examples are unit values for illustration, as in the source manuscript — not calibrated or clinically validated.

## Citation

See [`CITATION.cff`](CITATION.cff).

## License

MIT — see [`LICENSE`](LICENSE).
