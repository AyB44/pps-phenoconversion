"""
pps/temporal.py

Time indexing and dynamic updating:

    PPS_i^updated(t) = lambda * PPS_i(t) + (1 - lambda) * PPS_i^updated(t-1)

Explicitly optional in the Supplementary Methods ("not required for
the core PPS definition") -- provided as an opt-in utility, never
applied automatically by pps.score or pps.pk.
"""
from typing import Optional


def update_pps(pps_t: float, previous_updated: Optional[float], lam: float) -> float:
    """One step of the optional temporal-smoothing rule.

    Parameters
    ----------
    pps_t : freshly computed PPS(t), e.g. from pps.score.compute_pps.
    previous_updated : previous PPS^updated(t-1), or None to initialize
        the recursion (returns pps_t unchanged, as there is nothing yet
        to smooth against).
    lam : smoothing weight in [0, 1]; higher = more weight on new data.
    """
    if not 0.0 <= lam <= 1.0:
        raise ValueError("lam (lambda) must be in [0, 1].")
    if previous_updated is None:
        return pps_t
    return lam * pps_t + (1.0 - lam) * previous_updated


class PPSTracker:
    """Stateful convenience wrapper around update_pps for one patient/pathway."""

    def __init__(self, lam: float):
        if not 0.0 <= lam <= 1.0:
            raise ValueError("lam (lambda) must be in [0, 1].")
        self.lam = lam
        self._updated: Optional[float] = None

    def push(self, pps_t: float) -> float:
        self._updated = update_pps(pps_t, self._updated, self.lam)
        return self._updated

    @property
    def current(self) -> Optional[float]:
        return self._updated