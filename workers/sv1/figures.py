"""Report figures and markdown helpers shared by the objectives' reports (functions adapted from SanntiS 842675a
``src/sanntis/train/report.py``): the threshold sweep with methods as points, the three-sector strata radar,
markdown tables. Colours: EMBL-EBI Visual Framework palette."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

# EMBL-EBI Visual Framework (vf-core tokens)
SANNTIS, ANTISMASH, GECCO = "#18974C", "#A6093D", "#193f90"
SANNTIS09, PROPHET, REFERENCE = "#6CC24A", "#8B4799", "#54585A"
INK, INK2, GRID, BAND, SURFACE = "#1a1c1a", "#54585A", "#E4E4E4", "#F3F3F3", "#FFFFFF"
# label prefix → (marker, colour); the first match wins, so the more specific prefixes come first
STYLE = {"SanntiS 0.9": ("P", SANNTIS09), "SanntiS": ("o", SANNTIS), "antiSMASH 7.1.0 relaxed": ("s", ANTISMASH),
         "antiSMASH 7.1.0 strict": ("^", ANTISMASH), "antiSMASH 7.1.0 loose": ("v", ANTISMASH),
         "GECCO 0.9": ("D", GECCO), "GECCO": ("d", GECCO), "BGC-Prophet": ("X", PROPHET)}
SECTOR_TITLES = {"class": "BGC class", "lineage": "Lineage", "novelty": "Novelty"}
Style = tuple[str, str]  # (marker, colour)


def style(name: str, styles: dict[str, Style] | None = None) -> Style:
    """``styles`` (exact names) first, then the prefix table."""
    if styles and name in styles:
        return styles[name]
    return next((v for k, v in STYLE.items() if name.startswith(k)), ("^", INK2))


def _style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "text.color": INK, "axes.labelcolor": INK2, "figure.facecolor": SURFACE})
    return plt


def fig_radar_sectors(tables: dict[str, pd.DataFrame], metric: str, label: str, title: str, path: Path, *,
                      styles: dict[str, Style] | None = None, fill: str | None = None) -> None:
    """One radar, three sectors (class, lineage, novelty): ``metric`` per stratum for each series (name → radar
    table with ``sector``, ``axis``, ``n`` and the metric). Axes come from the first table (n = its BGCs); a gap and
    a shaded band separate the sectors. ``fill``: the series drawn with a shaded area (default: the ``SanntiS (t*``
    one)."""
    p = plt()
    ref = next(iter(tables.values()))
    sectors = [x for x in SECTOR_TITLES if (ref["sector"] == x).any()]
    slots, ang_of = [], {}
    for sec in sectors:
        slots += [(sec, a) for a in ref.loc[ref["sector"] == sec, "axis"]] + [None]  # None = gap
    step = 2 * np.pi / len(slots)
    for i, sl in enumerate(slots):
        if sl:
            ang_of[sl] = (i + 0.5) * step
    fig = p.figure(figsize=(6.4, 7.0))
    ax = fig.add_subplot(projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)  # clockwise from the top
    i0 = 0
    for k, sec in enumerate(sectors):
        n = int((ref["sector"] == sec).sum())
        lo, hi = i0 * step, (i0 + n) * step
        if k % 2 == 0:
            ax.fill_between(np.linspace(lo, hi, 50), 0, 1, color=BAND, zorder=0)
        ax.text((lo + hi) / 2, 1.36, SECTOR_TITLES[sec], ha="center", va="center", fontsize=9, fontweight="bold",
                color=INK)
        i0 += n + 1
    keys = [sl for sl in slots if sl]
    angs = [ang_of[k] for k in keys]
    for name, t in tables.items():
        mk, col = style(name, styles)
        val = t.set_index(["sector", "axis"])[metric]
        v = [float(val.get(k, np.nan)) for k in keys]
        ax.plot(angs + angs[:1], v + v[:1], color=col, lw=2, marker=mk, ms=5, mec=SURFACE, mew=1, label=name,
                zorder=3)
        if name == fill if fill else name.startswith("SanntiS (t*"):
            ax.fill(angs + angs[:1], v + v[:1], color=col, alpha=0.10, zorder=2)
    n_of = ref.set_index(["sector", "axis"])["n"]
    ax.set_xticks(angs, [f"{a}\n(n={int(n_of[(s_, a)])})" for s_, a in keys], fontsize=7.5, color=INK)
    ax.tick_params(axis="x", pad=10)
    ax.set_rlabel_position(np.degrees((len(slots) - 0.5) * step))  # in the last gap
    ax.set_xlim(0, 2 * np.pi)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0], ["0.25", "0.5", "0.75", "1"], fontsize=7, color=INK2)
    ax.grid(color=GRID)
    ax.spines["polar"].set_color(GRID)
    ax.set_title(f"{title}: {label}", fontsize=9, loc="left", color=INK, pad=34)
    h, lab = ax.get_legend_handles_labels()
    fig.legend(h, lab, frameon=False, fontsize=8, loc="lower center", ncol=min(len(tables), 3))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(path, dpi=150)
    p.close(fig)


PANELS = (("recall_macro_t", "precision_t", "Recall_macro_T", "Precision_T"),
          ("bgc_macro_recall", "bgc_precision", "BGC macro recall", "BGC precision"))


def fig_sweep(curve: pd.DataFrame, t_star: float, points: dict[str, dict], title: str, path: Path, *,
              label: str = "SanntiS 1.2.0", reference: tuple[str, pd.DataFrame, float] | None = None,
              styles: dict[str, Style] | None = None) -> None:
    """SanntiS threshold sweep (Precision_T vs Recall_macro_T, BGC precision vs BGC macro recall) with its t*; each
    entry of ``points`` (name → metrics) a marker; ``reference``: (label, curve, t*) drawn dashed in grey."""
    p = plt()
    fig, axs = p.subplots(1, 2, figsize=(8.4, 3.8))
    at = curve[np.isclose(curve["threshold"], t_star)]
    for ax, (x, y, xl, yl) in zip(axs, PANELS, strict=True):
        if reference is not None:
            rl, rc, rt = reference
            ax.plot(rc[x], rc[y], color=REFERENCE, lw=1.5, ls="--", label=f"{rl} (threshold sweep)")
            rat = rc[np.isclose(rc["threshold"], rt)]
            if len(rat):
                ax.plot(rat[x], rat[y], "o", color=REFERENCE, ms=6, mec=SURFACE, mew=1.5, label=f"{rl} t* {rt:.2f}")
        ax.plot(curve[x], curve[y], color=SANNTIS, lw=2, label=f"{label} (threshold sweep)")
        if len(at):
            ax.plot(at[x], at[y], "o", color=SANNTIS, ms=8, mec=SURFACE, mew=2, label=f"{label} t* {t_star:.2f}")
        for name, m in points.items():
            mk, col = style(name, styles)
            ax.plot(m[x], m[y], mk, color=col, ms=7, mec=SURFACE, mew=1.5, label=name)
        ax.set_xlim(0, 1.02)
        ax.set_ylim(0, 1.02)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        _style(ax)
    axs[0].set_title(title, fontsize=9, loc="left")
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=3, fontsize=7, frameon=False)
    fig.tight_layout(rect=(0, 0.2, 1, 1))
    fig.savefig(path, dpi=150)
    p.close(fig)


def md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(str(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join(lines)


def fmt(v: float | None, ci: tuple[float, float] | None = None) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "–"
    return f"{v:.3f}" + (f" ({ci[0]:.3f}–{ci[1]:.3f})" if ci else "")
