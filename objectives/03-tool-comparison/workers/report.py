"""Tool-comparison report: markdown + figures from ``score.py`` (and ``npc_compare.py``) outputs.

    python report.py --score DIR [--classification DIR] --out DIR

Figures and tables follow the SanntiS 1.0 validation report (``sv1.figures``): per test set the metric table with
95% CIs, the SanntiS threshold sweep with each tool as a point, and Recall_T / BGC recall radars over class, lineage
and novelty strata.
"""

import argparse
import pickle
from pathlib import Path

import pandas as pd
from sanntis.train import metrics as M

from sv1.figures import fig_radar_sectors, fig_sweep, fmt, md_table

RADAR_METHODS = ("SanntiS", "antiSMASH 7.1.0 relaxed", "GECCO 0.11", "SanntiS 0.9", "BGC-Prophet")
_fmt = fmt


def detection(results: dict, t_star: float, fig_dir: Path) -> str:
    head = "| Method | " + " | ".join(M.METRIC_LABELS[x] for x in M.METRICS) + " |\n|---|" + "---|" * len(M.METRICS)
    md = ""
    for name, tools in results.items():
        key = name.replace("-", "_").lower()
        rows, krows = [], []
        for tool, res in tools.items():
            label = f"SanntiS 1.0 (t* = {t_star:.2f})" if tool == "SanntiS" else tool
            rows.append(f"| {label} | " + " | ".join(_fmt(res["metrics"][x], res["ci"].get(x)) for x in M.METRICS)
                        + " |")
            kr = {"Method": label}
            for k, lab in (("proka", "proka"), ("euk", "euk"), ("fungi", "fungal BGCs")):
                if f"metrics_{k}" in res:
                    kr[f"F1_macro_T {lab}"] = f"{res[f'metrics_{k}']['f1_macro_t']:.3f}"
                    kr[f"Recall_macro_T {lab}"] = f"{res[f'metrics_{k}']['recall_macro_t']:.3f}"
            krows.append(kr)
        fig_sweep(tools["SanntiS"]["curve"], t_star, {t: r["metrics"] for t, r in tools.items() if t != "SanntiS"},
                  name, fig_dir / f"{key}_sweep.png")
        picked = [next((t for t in tools if t.startswith(p)), None) for p in RADAR_METHODS]
        rad = {(f"SanntiS 1.0 (t* = {t_star:.2f})" if t == "SanntiS" else t): tools[t]["radar"]
               for t in dict.fromkeys(picked) if t and "radar" in tools[t]}
        figs = " ".join(_radar(rad, metric, label, name, fig_dir / f"{key}_radar_{metric}.png")
                        for metric, label in (("recall_t", "Recall_T"), ("bgc_recall", "BGC recall")))
        md += (f"\n### {name}\n\n{head}\n" + "\n".join(rows)
               + f"\n\nn = {int(tools['SanntiS']['metrics']['n_gt'])} BGCs.\n\n" + md_table(pd.DataFrame(krows))
               + f"\n\n![{name} sweep](figures/{key}_sweep.png)\n\n{figs}\n")
    return md


EUK_METRICS = ("precision_t", "recall_t", "f1_macro_t", "bgc_precision", "bgc_recall")
# (set kind, scope, description): isolate contigs carry several BGCs of mixed lineages, so only the kingdom is
# attributable there; each shuffled contig carries one BGC, so the fungal scope is clean
EUK_SCOPES = (("isolate", "euk", "eukaryotic contigs (all their predictions and BGCs)"),
              ("shuffled", "fungi", "contigs carrying a fungal BGC (one BGC per contig)"))


def _lean(tools: dict, t_star: float) -> dict[str, dict]:
    """The report's reduced method list (RADAR_METHODS: one antiSMASH level, one GECCO version), labelled."""
    picked = [next((t for t in tools if t.startswith(p)), None) for p in RADAR_METHODS]
    return {(f"SanntiS 1.0 (t* = {t_star:.2f})" if t == "SanntiS" else t): tools[t] for t in dict.fromkeys(picked) if t}


def _scope_table(tools: dict[str, dict], scope: str) -> tuple[str, dict]:
    rows = [{"Method": lab, **{M.METRIC_LABELS[x]: _fmt(r[f"metrics_{scope}"][x], r[f"ci_{scope}"].get(x))
                               for x in EUK_METRICS}}
            for lab, r in tools.items() if f"metrics_{scope}" in r]
    first = next(r[f"metrics_{scope}"] for r in tools.values() if f"metrics_{scope}" in r)
    return md_table(pd.DataFrame(rows)), first


def _paired(path: Path) -> str:
    if not path.exists():
        return ""
    p = pd.read_csv(path, sep="\t")
    # the scopes of the tables above: eukaryotic contigs on isolate, fungal on shuffled
    p = p[p["set"].str.startswith("B-") & ((p["set"].str.endswith("isolate")) == (p["scope"] == "euk"))]
    rows = [{"Set, scope": f"{r.set}, {'fungal' if r.scope == 'fungi' else 'eukaryotic contigs'} "
                           f"({r.n_contigs} contigs)",
             "Metric": M.METRIC_LABELS[r.metric], "SanntiS 1.0": f"{r.a:.3f}", "fungiSMASH": f"{r.b:.3f}",
             "Difference (95% CI)": f"{r.diff:+.3f} ({r.ci_lo:+.3f} to {r.ci_hi:+.3f})",
             "p": "< 0.001" if r.p < 1 / r.n_boot else f"{r.p:.3f}"} for r in p.itertuples()]
    return ("\n#### SanntiS 1.0 vs fungiSMASH, paired (Test B)\n\nSanntiS 1.0 − antiSMASH 7.1.0 relaxed "
            "(lineage-matched). Paired bootstrap: each replicate resamples the scope's contigs once and scores both "
            "methods on it; 95% percentile interval, two-sided bootstrap p (resolution 1 / replicates; replicates "
            "with a difference of exactly 0 count in both tails). Test A and the shuffled eukaryotic scope in "
            "`score/paired.tsv`.\n\n" + md_table(pd.DataFrame(rows)) + "\n")


def eukaryotes(results: dict, t_star: float, fig_dir: Path, paired_tsv: Path) -> str:
    md = ("Two scopes, chosen so precision is attributable: "
          "on the **isolate** sets the eukaryotic contigs (fungal / protist backbones carrying several test BGCs of "
          "mixed lineages, so precision is per kingdom, not per lineage); on the **shuffled** sets the contigs "
          "whose one test BGC is fungal (precision and recall both about fungal BGCs). Precision does not depend "
          "on the test pool on the isolate sets (same contigs for A and B); on the shuffled fungal scope it does "
          "(the pool picks the contigs). 95% CIs: bootstrap over the contigs of the scope. Methods: SanntiS 1.0, "
          "antiSMASH 7.1.0 relaxed (lineage-matched: fungiSMASH, plantiSMASH over plant BGCs), GECCO 0.11.0 "
          "(regions identical to 0.9.5), SanntiS 0.9.5.0, BGC-Prophet 0.1.2.\n")
    for test, note in (("B", ""), ("A", " (fewer BGCs, wider CIs: supplementary to Test B)")):
        md += f"\n#### Test {test}{note}\n"
        for kind, scope, desc in EUK_SCOPES:
            tools = _lean(results[f"{test}-{kind}"], t_star)
            tab, m = _scope_table(tools, scope)
            md += (f"\n{test}-{kind}, {desc}: {int(m['n_contigs'])} contigs, {int(m['n_gt'])} test BGCs.\n\n"
                   f"{tab}\n")
        if test == "B":
            md += _paired(paired_tsv)
    md += ("\n#### Prokaryotic vs eukaryotic contigs (Test B)\n\nPrecision_T and Recall_T of each method on the "
           "prokaryotic (open marker) and eukaryotic (filled, 95% CI) contigs of each Test B set.\n\n"
           + _delta_fig(results, t_star, fig_dir / "b_proka_vs_euk.png") + "\n")
    tools = _lean(results["B-shuffled"], t_star)
    m = next(iter(tools.values()))["metrics_fungi"]
    classes = [c for c in M.CLASSES if m.get(f"n_gt[{c}]", 0) > 0]
    rows = [{"Class": f"{c} (n = {m[f'n_gt[{c}]']:.0f}){' †' if m[f'n_gt[{c}]'] < 10 else ''}",
             **{lab: f"{r['metrics_fungi'][f'recall_t[{c}]']:.3f}" for lab, r in tools.items()}} for c in classes]
    md += ("\n#### Fungal BGCs per class (B-shuffled)\n\nRecall_T per BGC class on the B-shuffled fungal scope "
           f"({int(m['n_gt'])} BGCs; n: BGCs of the class, a hybrid counts in each of its classes; no CIs). "
           "† n < 10: shown, not to be interpreted.\n\n" + md_table(pd.DataFrame(rows)) + "\n")
    return md


def _delta_fig(results: dict, t_star: float, path: Path) -> str:
    from sv1.figures import GRID, INK2, SURFACE, _style, plt, style

    p = plt()
    fig, axes = p.subplots(2, 2, figsize=(9, 5.2), sharex=True, sharey=True)
    for i, kind in enumerate(("isolate", "shuffled")):
        tools = _lean(results[f"B-{kind}"], t_star)
        labels = list(tools)[::-1]
        for j, x in enumerate(("precision_t", "recall_t")):
            ax = axes[i, j]
            _style(ax)
            ax.grid(axis="x", color=GRID, linewidth=0.6)
            for y, lab in enumerate(labels):
                r = tools[lab]
                mk, col = style("SanntiS" if lab.startswith("SanntiS 1.0") else lab)
                a, b = r["metrics_proka"][x], r["metrics_euk"][x]
                lo, hi = r["ci_euk"][x]
                ax.plot([a, b], [y, y], color=INK2, lw=1, alpha=0.5, zorder=1)
                ax.plot([lo, hi], [y, y], color=col, lw=2, zorder=2)
                ax.plot(a, y, mk, ms=8, mfc=SURFACE, mec=col, mew=1.5, zorder=3)
                ax.plot(b, y, mk, ms=8, color=col, mec=SURFACE, mew=1, zorder=4)
            ax.set_yticks(range(len(labels)), labels, fontsize=8)
            ax.set_xlim(0, 1.02)
            ax.set_title(f"B-{kind}: {M.METRIC_LABELS[x]}", fontsize=9, loc="left")
    h = [p.Line2D([], [], marker="o", ls="", ms=7, mfc=SURFACE, mec=INK2, mew=1.5, label="prokaryotic contigs"),
         p.Line2D([], [], marker="o", ls="", ms=7, color=INK2, label="eukaryotic contigs (95% CI)")]
    fig.legend(handles=h, loc="lower center", ncol=2, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(path, dpi=150)
    p.close(fig)
    return f"![Test B, prokaryotic vs eukaryotic contigs](figures/{path.name})"


def _radar(tables: dict, metric: str, label: str, title: str, path: Path) -> str:
    fig_radar_sectors(tables, metric, label, title, path)
    return f"![{title} {label} radar](figures/{path.name})"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--score", type=Path, required=True)
    p.add_argument("--classification", type=Path)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    fig_dir = a.out / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    r = pickle.loads((a.score / "results.pkl").read_bytes())
    md = ["# SanntiS 1.0: tool comparison\n",
          "Every method on the same contigs and CDS (objective 02's nucleotide test sets on backgrounds disjoint from "
          "the training negatives: isolate = MIBiG test BGCs inserted into BGC-free genome stretches; shuffled = into "
          "gene-shuffled genome contigs), scored with the SanntiS 1.0 metric code. Test A: MIBiG 5.0rc1 active "
          "entries new since v3 (>= BGC0002750); Test B: active + pending. "
          "SanntiS 1.0 is the model trained on v3-active (objective 02) at its operating threshold t*. antiSMASH is "
          "lineage-matched: bacteria on prokaryotic contigs, fungiSMASH on eukaryotic ones, plantiSMASH 1.0 for "
          "regions over plant BGCs. 95% CIs: bootstrap over source contigs.\n",
          "## Detection\n", detection(r["results"], r["t_star"], fig_dir),
          "\n### Eukaryotic and fungal BGCs\n", eukaryotes(r["results"], r["t_star"], fig_dir, a.score / "paired.tsv")]
    if a.classification and (a.classification / "classification.tsv").exists():
        c = pd.read_csv(a.classification / "classification.tsv", sep="\t").round(3)
        md += ["\n## Product classification (NPClassifier)\n",
               "Test B BGCs with an NPClassifier label; SanntiS 1.0 NPC on the true loci, BGCat on the biosynthetic "
               "genes antiSMASH 7.1.0 finds in them (its documented input), ranking rolled up from NPClassifier "
               "classes. Superclasses in the SanntiS label space (rare ones folded to `other <pathway>`).\n",
               md_table(c), ""]
    (a.out / "report.md").write_text("\n".join(md))
    print("wrote", a.out / "report.md")


if __name__ == "__main__":
    main()
