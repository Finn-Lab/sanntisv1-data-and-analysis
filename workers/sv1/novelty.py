# Copied from SanntiS 842675a `src/sanntis/train/novelty.py`
# Test-set code leaves SanntiS in EMG-10233; only the imports were changed.
"""Novelty strata (spec §2.5): max Dice between a BGC's Pfam set and the Pfam set of any BGC the scored model was
trained on. Bins: novel < 0.5, related 0.5-0.9, near-duplicate >= 0.9."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import scipy.sparse as sp

BINS = ("novel", "related", "near-duplicate")


def _matrix(sets: Sequence[set], vocab: dict) -> sp.csr_matrix:
    r, c = [], []
    for i, s in enumerate(sets):
        for t in s:
            j = vocab.get(t)
            if j is not None:
                r.append(i)
                c.append(j)
    return sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(len(sets), len(vocab)))


def max_dice(query: Sequence[set], reference: Sequence[set]) -> np.ndarray:
    vocab = {t: i for i, t in enumerate(sorted(set().union(*query, *reference)))} if (query or reference) else {}
    Q, R = _matrix(query, vocab), _matrix(reference, vocab)
    inter = (Q @ R.T).toarray()
    den = np.asarray(Q.sum(1)).ravel()[:, None] + np.asarray(R.sum(1)).ravel()[None, :]
    d = np.divide(2 * inter, den, out=np.zeros_like(inter, dtype=float), where=den > 0)
    return d.max(1) if d.shape[1] else np.zeros(len(query))


def novelty_bin(x: float) -> str:
    return "novel" if x < 0.5 else ("related" if x < 0.9 else "near-duplicate")
