"""MIBiG 3.1 validation report: markdown + figures from the v3-active training run and the Test A / B scoring.

    python report_v3active.py --run DIR --eval DIR [--expected DIR] --out DIR

``--run``: the ``sanntis train`` run (``collect/done.json``, ``tune/tune.json``, ``bundle/``); ``--eval``: the
``05_score.sh`` outputs (``score/{A,B}-{isolate,shuffled}.pkl``, ``test_bgcs.parquet``); ``--expected``: the frozen
tables of the pre-split release run (reference columns: its v3-active had 2,416 BGCs and its test backgrounds were
pool v1 genomes, so they are a reference, not a target). Figures: the threshold sweep per test set (reference
sweep dashed), Recall_T / BGC recall radars over class, lineage and novelty strata (one series per test set).
"""

import argparse
import json
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sanntis.train import metrics as M

from sv1.figures import ANTISMASH, GECCO, PROPHET, SANNTIS, fig_radar_sectors, fig_sweep, fmt, md_table

SETS = ("A-isolate", "B-isolate", "A-shuffled", "B-shuffled")
STYLES = {"A-isolate": ("o", SANNTIS), "B-isolate": ("s", GECCO), "A-shuffled": ("D", ANTISMASH),
          "B-shuffled": ("^", PROPHET)}
HEADLINE = ("precision_t", "recall_macro_t", "f1_macro_t", "bgc_precision", "bgc_macro_recall")
NPC_COLS = ("n", "pathway_macro_ap", "pathway_macro_f1", "pathway_top1", "superclass_top1", "superclass_top3",
            "coverage@0.9", "precision@0.9")


def tsv(d: Path | None, name: str) -> pd.DataFrame | None:
    f = d / name if d else None
    return pd.read_csv(f, sep="\t") if f and f.exists() else None


def training(run: Path, exp: Path | None) -> str:
    done = json.loads((run / "collect" / "done.json").read_text())["outputs"]
    tune = json.loads((run / "tune" / "tune.json").read_text())
    man = json.loads((run / "bundle" / "manifest.json").read_text())
    ref = tsv(exp, "tune.tsv")
    rows = [("Training BGCs (v3-active, collected loci)", f"{done['loci']:,} (excluded > 200 kb: {done['excluded']})",
             "2,416 BGCs"),
            ("by kingdom",
             f"proka {done['by_kingdom'].get('proka', 0):,} / euk {done['by_kingdom'].get('euk', 0):,}", ""),
            ("by primary class", ", ".join(f"{k} {v:,}" for k, v in sorted(done["by_primary_class"].items())), ""),
            ("NPC labels", ", ".join(f"{k} {v:,}" for k, v in sorted(done["by_npc_status"].items())), ""),
            ("Features", f"InterProScan {man['ips']['version']}, InterPro {man['ips']['interpro']}, "
                         f"{len(man['ips']['apps'])} apps, vocabulary {man['vocab_size']:,}", ""),
            ("Epochs (median of the CV best epochs)", f"{tune['epochs']} ({', '.join(map(str, tune['best_epochs']))})",
             f"{int(ref['epochs'][0])} {ref['best_epochs'][0]}" if ref is not None else ""),
            ("Operating point t*", f"{tune['t_star']:.2f} (low {tune['threshold_low']:.2f}; floor "
                                   f"{tune['floor_metric']} >= {tune['floor']}, met: {tune['floor_met']})",
             f"{ref['t_star'][0]:.2f}" if ref is not None else "")]
    return md_table(pd.DataFrame(rows, columns=["", "this run", "reference (release run)"]))


def validation(run: Path, exp: Path | None) -> str:
    v = pd.read_csv(run / "bundle" / "report" / "validation.tsv", sep="\t").set_index("stratum")
    ref = tsv(exp, "validation_oof.tsv")
    ref = ref.set_index(["stratum", "metric"])["value"] if ref is not None else None
    rows = []
    for m in M.METRICS:
        r = {"Metric": M.METRIC_LABELS.get(m, m)}
        for st, key in (("All", "all"), ("Proka", "proka"), ("Euk", "euk")):
            ci = v.at[st, f"{m}_ci"] if f"{m}_ci" in v and isinstance(v.at[st, f"{m}_ci"], str) else None
            nums = re.findall(r"\d*\.?\d+", ci) if ci else []  # "0.6220-0.6503"
            ci = (float(nums[0]), float(nums[1])) if len(nums) == 2 else None
            r[st] = fmt(float(v.at[st, m]), ci)
            if ref is not None and (key, m) in ref:
                r[f"{st} ref."] = fmt(float(ref[(key, m)]))
        rows.append(r)
    return md_table(pd.DataFrame(rows))


def test_tables(res: dict[str, dict], exp_metrics: pd.DataFrame | None) -> str:
    head = "| Metric | " + " | ".join(SETS) + " |\n|---|" + "---|" * len(SETS)
    lines = []
    for m in M.METRICS:
        cells = [fmt(res[s]["metrics"][m], res[s]["ci"].get(m)) for s in SETS]
        lines.append(f"| {M.METRIC_LABELS.get(m, m)} | " + " | ".join(cells) + " |")
    if exp_metrics is not None:
        e = exp_metrics[exp_metrics["scope"] == "all"].set_index(["set", "metric"])["value"]
        for m in HEADLINE:
            cells = [fmt(float(e[(s, m)])) if (s, m) in e else "–" for s in SETS]
            lines.append(f"| {M.METRIC_LABELS.get(m, m)} (reference) | " + " | ".join(cells) + " |")
    n = " · ".join(f"{s}: {int(res[s]['metrics']['n_gt'])} BGCs on {res[s]['n_contigs']} contigs" for s in SETS)
    king = []
    for s in SETS:
        for k, lab in (("proka", "prokaryotic"), ("euk", "eukaryotic")):
            if f"metrics_{k}" in res[s]:
                a = res[s][f"metrics_{k}"]
                king.append({"set": s, "contigs": lab, "n BGCs": int(a["n_gt"]), "Precision_T": fmt(a["precision_t"]),
                             "Recall_macro_T": fmt(a["recall_macro_t"]), "F1_macro_T": fmt(a["f1_macro_t"]),
                             "BGC precision": fmt(a["bgc_precision"]), "BGC macro recall": fmt(a["bgc_macro_recall"])})
    shuf = [f"{s}: {res[s]['ignored_native_predictions']} predictions over native BGCs ignored; natives as background: "
            f"Precision_T {res[s]['natives_as_background']['precision_t']:.3f}, BGC precision "
            f"{res[s]['natives_as_background']['bgc_precision']:.3f}"
            for s in SETS if "natives_as_background" in res[s]]
    return (head + "\n" + "\n".join(lines) + f"\n\n{n}.\n\n### By kingdom\n\n" + md_table(pd.DataFrame(king))
            + ("\n\nShuffled backgrounds: " + "; ".join(shuf) + ".\n" if shuf else "\n"))


def strata(res: dict[str, dict], exp: pd.DataFrame | None) -> str:
    frames = []
    for s in SETS:
        df = res[s]["strata"].copy()
        df.insert(0, "set", s)
        if exp is not None:
            e = exp[exp["set"] == s].set_index(["axis", "stratum"])
            df["BGC macro recall ref."] = [fmt(float(e.at[(a, b), "bgc_macro_recall"])) if (a, b) in e.index else "–"
                                           for a, b in zip(df["axis"], df["stratum"], strict=True)]
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    for c in ("bgc_macro_recall", "recall_macro_t"):
        df[c] = df[c].map(fmt)
    return md_table(df)


def npc(res: dict[str, dict], exp: pd.DataFrame | None) -> str:
    rows = []
    for s in SETS:
        t = res[s].get("npc")
        if t is None or not len(t):
            continue
        for _, r in t.iterrows():
            row = {"set": s, "input": r["input"], "stratum": r["stratum"]}
            row.update({c: (int(r[c]) if c == "n" else fmt(float(r[c]))) for c in NPC_COLS if c in r})
            if exp is not None:
                e = exp[(exp["set"] == s) & (exp["input"] == r["input"]) & (exp["stratum"] == r["stratum"])]
                row["pathway_top1 ref."] = fmt(float(e["pathway_top1"].iloc[0])) if len(e) else "–"
                row["superclass_top1 ref."] = fmt(float(e["superclass_top1"].iloc[0])) if len(e) else "–"
            rows.append(row)
    return md_table(pd.DataFrame(rows)) if rows else "(no NPC classifier in the bundle)"


def novelty(tb: pd.DataFrame) -> str:
    rows = []
    for test, col in (("Test A", "test_a"), ("Test B", "test_b")):
        sub = tb[tb[col]]
        c = sub["novelty"].value_counts()
        rows.append({"test": test, "BGCs": len(sub),
                     **{k: int(c.get(k, 0)) for k in ("near-duplicate", "related", "novel")},
                     "median max Dice": f"{sub['max_dice'].median():.2f}"})
    return md_table(pd.DataFrame(rows))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--eval", type=Path, required=True)
    p.add_argument("--expected", type=Path)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    figs = a.out / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    res = {s: pickle.loads((a.eval / "score" / f"{s}.pkl").read_bytes()) for s in SETS}
    summary = json.loads((a.eval / "summary.json").read_text())
    t_star = float(summary["t_star"])
    tb = pd.read_parquet(a.eval / "test_bgcs.parquet")
    exp = a.expected if a.expected and a.expected.is_dir() else None
    exp_metrics, exp_curve, exp_strata, exp_npc = (tsv(exp, f"test_{x}.tsv")
                                                   for x in ("metrics", "curve", "strata", "npc"))
    ref_t = float(tsv(exp, "tune.tsv")["t_star"][0]) if tsv(exp, "tune.tsv") is not None else t_star

    sweeps = []
    for s in SETS:
        ref = (("release run (pool v1 backgrounds)", exp_curve[exp_curve["set"] == s], ref_t)
               if exp_curve is not None else None)
        fig_sweep(res[s]["curve"], t_star, {}, s, figs / f"{s}_sweep.png", label="v3-active", reference=ref)
        sweeps.append(f"![{s} sweep](figures/{s}_sweep.png)")
    radars = []
    for metric, label in (("recall_t", "Recall_T"), ("bgc_recall", "BGC recall")):
        fig_radar_sectors({s: res[s]["radar"] for s in SETS}, metric, label, f"v3-active model at t* {t_star:.2f}",
                          figs / f"radar_{metric}.png", styles=STYLES, fill="B-isolate")
        radars.append(f"![{label} radar](figures/radar_{metric}.png)")

    md = ["# SanntiS 1.2.0: MIBiG 3.1 validation\n",
          "A SanntiS model trained on v3-active (MIBiG entries with accession < BGC0002750, active or pending in "
          "5.0rc1) and scored on the entries added since (>= BGC0002750): Test A = active, Test B = active + pending, "
          "each inserted into isolate and gene-shuffled nucleotide contigs whose genomes are not in the training "
          "negative pool. Metrics: the SanntiS 1.2.0 metric code at the model's operating threshold t*; 95% CIs by "
          "bootstrap over source contigs. Reference columns and dashed curves: the pre-split release run (v3-active "
          "of 2,416 BGCs, test backgrounds from pool v1 genomes) — a reference, not a target.\n",
          "## Training\n", training(a.run, exp),
          "\n## Validation (pooled out-of-fold CV on the v3-active pool)\n", validation(a.run, exp),
          "\n## Test sets\n", test_tables(res, exp_metrics), "\n### Threshold sweeps\n", "\n\n".join(sweeps),
          "\n### Strata\n", strata(res, exp_strata), "", " ".join(radars),
          "\n### Novelty of the test BGCs vs the training positives (max Dice over Pfam sets)\n", novelty(tb),
          "\n## Product classification (NPClassifier pathway / superclass)\n",
          "Test BGCs with an NPClassifier label: the NPC classifier on the true loci (ceiling) and on the predicted "
          "regions at t* (a region takes the labels of the BGC it overlaps most).\n", npc(res, exp_npc), ""]
    (a.out / "report.md").write_text("\n".join(md))
    pd.DataFrame([{"set": s, **{m: round(res[s]["metrics"][m], 4) for m in HEADLINE}} for s in SETS]) \
        .to_csv(a.out / "summary.tsv", sep="\t", index=False)
    print("wrote", a.out / "report.md")


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
