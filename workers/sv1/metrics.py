"""``sanntis.train.metrics`` plus the test-set extras SanntiS dropped before the alpha (27d2680).

Re-exports everything from the pinned ``sanntis.train.metrics`` (``Scored``, ``aggregate``, ``METRICS``,
``CLASSES``, ``threshold_grid``, ...) and restores, with the 842675a semantics:

- ``score_contig(..., orf_native=...)``: on shuffled backgrounds a prediction touching native-BGC ORFs and no
  true-BGC ORF is ignored (neither right nor wrong), and native ORFs do not count against Precision_T.
- ``bootstrap(..., gt_mask=...)``: resampling units stay those of the full ``Scored`` while only the masked
  ground-truth rows are scored.
- ``auprc``: step-wise area under a PR curve with interpolated precision.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
import pandas as pd
from sanntis.train import metrics as _M
from sanntis.train.metrics import *  # noqa: F401,F403  (re-export the public API)
from sanntis.train.metrics import METRICS, Region, Scored, _arrays, _metrics


def score_contig(pred: Sequence[Region], orf_bgc: np.ndarray, *, unit: object = 0,
                 orf_native: np.ndarray | None = None, bgc_classes: dict | None = None,
                 pred_unit: Callable[[int, int], object] | None = None) -> Scored:
    """``sanntis.train.metrics.score_contig`` with the native-ORF rule of 842675a."""
    if orf_native is None:
        return _M.score_contig(pred, orf_bgc, unit=unit, bgc_classes=bgc_classes, pred_unit=pred_unit)
    is_bgc = np.array([not (b is None or b == -1 or b == "" or (isinstance(b, float) and np.isnan(b)))
                       for b in orf_bgc], dtype=bool)
    pos_cum = np.concatenate([[0], np.cumsum(is_bgc)])
    nat_cum = np.concatenate([[0], np.cumsum(np.asarray(orf_native, bool))])
    kept, n_nat = [], []
    for a, b in pred:
        nat, true = int(nat_cum[b] - nat_cum[a]), int(pos_cum[b] - pos_cum[a])
        if nat and true == 0:
            continue  # a call about an unlabelled native BGC: neither right nor wrong
        kept.append((a, b))
        n_nat.append(nat)
    s = _M.score_contig(kept, orf_bgc, unit=unit, bgc_classes=bgc_classes, pred_unit=pred_unit)
    if len(s.pred):
        s.pred["denom"] = np.maximum(np.array([b - a for a, b in kept]) - np.array(n_nat), 1)
    return s


def bootstrap(s: Scored, *, n_boot: int = 1000, seed: int = 0,
              gt_mask: np.ndarray | None = None) -> dict[str, tuple[float, float]]:
    """95% percentile intervals over resampled units (``sanntis.train.metrics.bootstrap`` with ``gt_mask``)."""
    gt = s.gt if gt_mask is None else s.gt[gt_mask]
    units = pd.unique(pd.concat([s.gt["unit"], s.pred["unit"]], ignore_index=True))
    code = {u: i for i, u in enumerate(units)}
    gc = gt["unit"].map(code).to_numpy(int)
    pc = s.pred["unit"].map(code).to_numpy(int)
    arrays = _arrays(gt, s.pred)
    rng = np.random.default_rng(seed)
    vals: dict[str, list[float]] = {m: [] for m in METRICS}
    for _ in range(n_boot):
        counts = np.bincount(rng.integers(0, len(units), len(units)), minlength=len(units)).astype(float)
        a = _metrics(*arrays, counts[gc], counts[pc])
        for m in METRICS:
            vals[m].append(a[m])
    return {m: (float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))) for m, v in vals.items()}


def auprc(recall: np.ndarray, precision: np.ndarray) -> float:
    """Step-wise area under a PR curve with interpolated precision (running max from high recall down)."""
    o = np.argsort(recall)
    r = np.concatenate([[0.0], np.asarray(recall, float)[o]])
    p = np.asarray(precision, float)[o]
    p_interp = np.maximum.accumulate(p[::-1])[::-1]
    return float(np.sum(np.diff(r) * p_interp))
