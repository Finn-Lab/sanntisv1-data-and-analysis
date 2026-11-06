"""Flat result tables (the ``expected/`` format) and their comparison.

    python -m sv1.expected compare <new dir> <expected dir> [--rtol 1e-5] [--files a.tsv b.tsv ...]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

FRAMES = ("strata", "radar", "curve", "npc", "npc_per_pathway")


def metric_rows(res: dict, **key) -> list[dict]:
    """Scalar metrics of one result (``metrics``, ``metrics_proka`` / ``_euk`` / ``_fungi``) with their CIs
    (``ci``, ``ci_proka`` / ``_euk`` / ``_fungi`` when present)."""
    rows = []
    for scope in ("metrics", "metrics_proka", "metrics_euk", "metrics_fungi"):
        for m, v in (res.get(scope) or {}).items():
            ci = (res.get(scope.replace("metrics", "ci")) or {}).get(m)
            rows.append({**key, "scope": scope.replace("metrics", "all").replace("all_", ""), "metric": m, "value": v,
                         "ci_lo": ci[0] if ci else None, "ci_hi": ci[1] if ci else None})
    return rows


def frames(res: dict, **key) -> dict[str, pd.DataFrame]:
    out = {name: res[name].assign(**key) for name in FRAMES if isinstance(res.get(name), pd.DataFrame)}
    if res.get("auprc"):
        out["auprc"] = pd.DataFrame([{**key, "curve": k, "auprc": v if not isinstance(v, tuple | list) else v[0]}
                                     for k, v in res["auprc"].items()])
    return out


def write(results: dict[str, dict[str, dict]], out: Path, prefix: str = "") -> list[Path]:
    """``results``: set → tool → result. Writes ``<prefix>metrics.tsv`` and one TSV per frame kind."""
    out.mkdir(parents=True, exist_ok=True)
    rows, extra = [], {}
    for s, tools in results.items():
        for tool, res in tools.items():
            rows += metric_rows(res, set=s, tool=tool)
            for k, df in frames(res, set=s, tool=tool).items():
                extra.setdefault(k, []).append(df)
    tables = {"metrics": pd.DataFrame(rows), **{k: pd.concat(v) for k, v in extra.items()}}
    paths = []
    for k, df in tables.items():
        paths.append(out / f"{prefix}{k}.tsv")
        df.to_csv(paths[-1], sep="\t", index=False, float_format="%.6g")
    return paths


def compare(new: Path, exp: Path, rtol: float = 1e-5, files: list[str] | None = None) -> bool:
    """Every expected table must match the new one: same rows (keyed on the non-numeric columns), numbers within
    ``rtol`` (the tables are written with 6 significant digits). Extra rows in ``new`` (e.g. new tools) are fine."""
    ok = True
    for f in files or sorted(p.name for p in exp.glob("*.tsv")):
        e, n = pd.read_csv(exp / f, sep="\t"), pd.read_csv(new / f, sep="\t") if (new / f).exists() else None
        if n is None:
            print(f"MISSING {f}")
            ok = False
            continue
        num = [c for c in e.columns if pd.api.types.is_numeric_dtype(e[c]) and c not in ("n", "n_bgc", "support")]
        keys = [c for c in e.columns if c not in num]
        if "threshold" in e.columns:
            keys = [*keys, "threshold"] if "threshold" not in keys else keys
            num = [c for c in num if c != "threshold"]
        m = e.merge(n, on=keys, how="left", suffixes=("", "_new"), indicator=True)
        missing = int((m["_merge"] != "both").sum())
        bad = 0
        for c in num:
            a, b = m[c].to_numpy(float), m[f"{c}_new"].to_numpy(float)
            diff = ~np.isclose(a, b, rtol=rtol, atol=1e-9, equal_nan=True)
            bad += int(diff.sum())
            if diff.any():
                i = int(np.flatnonzero(diff)[0])
                print(f"  {f}: {c} first mismatch {m.loc[i, keys].to_dict()}: expected {a[i]}, got {b[i]}")
        ok &= missing == 0 and bad == 0
        print(f"{'OK  ' if missing == 0 and bad == 0 else 'DIFF'} {f}: {len(e)} rows, {missing} missing, {bad} "
              f"values differ")
    return ok


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("cmd", choices=["compare"])
    p.add_argument("new", type=Path)
    p.add_argument("expected", type=Path)
    p.add_argument("--rtol", type=float, default=1e-5)
    p.add_argument("--files", nargs="*")
    a = p.parse_args()
    sys.exit(0 if compare(a.new, a.expected, a.rtol, a.files) else 1)
