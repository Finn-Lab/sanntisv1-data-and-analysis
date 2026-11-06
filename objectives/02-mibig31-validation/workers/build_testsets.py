"""The nucleotide test sets: the Test B BGCs (MIBiG 5.0rc1 active + pending, accession >= BGC0002750, locus <= 200 kb,
objective 01's archives) inserted into the test-background pool of 03_background.sh, written per contig as the tool
input GenBank (CDS + translations), the CDS truth table and the ground-truth BED that objective 03 consumes.
Ported from the sanntis_dev nucleotide builders (synthetic_isolate_builder / _shuffled) and the feature-level
builders of the release run (SanntiS 842675a ``train/testsets.py``); the pool tables replace their inputs.

- isolate: the pool's isolate test backbones (55 proka + 15 euk, longest BGC-free run >= 50 CDS of a genome),
  sliced first CDS start -> last CDS end; the kingdom's test BGCs dealt round-robin (sorted by accession) over the
  backbones and each inserted (+ strand) at a random intergenic position, never splitting a backbone CDS.
- shuffled: every background genome's contigs with >= --min-genes gene loci (overlapping CDS share a locus),
  cropped to a random window of <= --max-genes loci; the BGCs (a random order) dealt round-robin over the
  kingdom's genomes, then over each genome's eligible contigs (random order, cycled); the loci are permuted, each
  with its upstream spacer, and one BGC inserted at a locus boundary in the middle 40-60 %. Backbone CDS inside an
  antiSMASH 8.0.1 region of the source genome are flagged ``native_bgc`` (scored as background).
Draws: ``default_rng(seed + crc32(track) % 1000)`` per set and track; every draw happens in the planning pass, so
the output does not depend on the order genomes are read.

    python build_testsets.py --archives <01 inputs/mibig_5.0rc1-train> --background <02/background> --out DIR
        [--seed 0] [--min-genes 60] [--max-genes 400] [--limit N]

writes <out>/testsets/<set>/<track>/<track>_contig_<k>.{gbk,cds.parquet,ground_truth.bed} (+ README.txt),
<out>/test_bgcs.parquet (Test A / B membership, kingdom, lineage, classes; tokens and novelty are added at scoring),
<out>/contigs.tsv, <out>/manifest.json
"""

from __future__ import annotations

import argparse
import json
import logging
import warnings
import zlib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import BiopythonWarning, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from sanntis.io.genbank import Cds, _span, iter_cds, read_records
from sanntis.train import mibig
from sanntis.train.ncbi import LineageCache, lineage_group

from sv1.archives import sha256

log = logging.getLogger("build_testsets")

V3_MAX = "BGC0002750"  # first accession after MIBiG 3.1
MAX_LOCUS_BP = 200_000  # sanntis train's locus limit (recipe dedup.max_locus_bp)
TRACKS = ("proka", "euk")
RECORD_PREFIX = {"isolate": "synth_iso", "shuffled": "synth_shuf"}
COPIED_QUALIFIERS = ("transl_table", "codon_start", "product")
CDS_COLUMNS = ["contig_id", "start", "end", "strand", "parts", "cds_id", "source_cds_id", "source_kind", "entity_id",
               "source_genome_id", "class", "product", "translation", "native_bgc"]
BED_COLUMNS = ["contig_id", "start", "end", "entity_id", "class"]


def rng_for(seed: int, track: str) -> np.random.Generator:
    return np.random.default_rng(seed + zlib.crc32(track.encode()) % 1000)


# ── sources ──────────────────────────────────────────────────────────────────


@dataclass
class Positive:
    entity_id: str
    seq: str
    cds: list[tuple[Cds, SeqFeature]]
    bgc_class: str
    source_contig: str


@dataclass
class Plan:
    """One output contig: where its pieces come from (source coordinates) and which BGCs go where."""

    set_name: str
    track: str
    k: int
    genome: str
    contig: str
    unit_id: str  # run id (isolate) or ``<contig>_shuf<n>`` (shuffled)
    units: list[tuple[int, int, list[int]]]  # (src_start, src_end, cds_idx list) in output order
    inserts: dict[int, list[tuple[int, str]]] = field(default_factory=dict)  # unit index → [(offset, accession)]

    @property
    def contig_id(self) -> str:
        return f"{RECORD_PREFIX[self.set_name]}_{self.track}_{self.k:03d}"

    @property
    def prefix(self) -> str:
        return f"{self.track}_contig_{self.k:03d}"


def cds_features(rec: SeqRecord, **kw) -> list[tuple[Cds, SeqFeature]]:
    """(CDS, feature) pairs in SanntiS's CDS order (``iter_cds``: /pseudo dropped, sorted by span), so that the
    pair index equals the pool's ``cds_idx``."""
    feats = [f for f in rec.features if f.type == "CDS" and "pseudo" not in f.qualifiers]
    feats.sort(key=lambda f: _span(f)[:2])
    cds = list(iter_cds(rec, **kw))
    if len(cds) != len(feats) or any((c.start, c.end) != _span(f)[:2] for c, f in zip(cds, feats, strict=True)):
        raise RuntimeError(f"{rec.id}: CDS features and iter_cds disagree")
    return list(zip(cds, feats, strict=True))


def load_positives(archives: Path, limit: int | None) -> tuple[pd.DataFrame, dict[str, Positive], pd.DataFrame]:
    """Test B candidates of objective 01's manifest → (test-BGC table, loci, excluded rows)."""
    man = pd.read_csv(archives / "manifest.tsv", sep="\t", dtype=str, keep_default_na=False)
    man["length"] = man["length"].astype(int)
    cand = man[(man["accession"] >= V3_MAX) & man["status"].isin(["active", "pending"])].sort_values("accession")
    excluded = cand[cand["length"] > MAX_LOCUS_BP].assign(reason=f"locus > {MAX_LOCUS_BP} bp")
    sel = cand[cand["length"] <= MAX_LOCUS_BP].reset_index(drop=True)
    if limit:
        sel = sel.head(limit)
    gbk = mibig.read_gbk_archive(archives / "mibig_5.0rc1-train_gbk.tar.gz", set(sel["accession"]))
    entries = mibig.read_entries(archives / "mibig_5.0rc1-train_json.tar.gz")
    taxid = {a: mibig.taxid_of(entries[a]) for a in sel["accession"]}
    lineage = LineageCache().resolve(t for t in taxid.values() if t)
    rows, loci = [], {}
    for acc, status in zip(sel["accession"], sel["status"], strict=True):
        rec = mibig.archive_record(gbk[acc])
        pairs = cds_features(rec, missing_translation="error")  # what sanntis train reads
        classes = mibig.classes_of(entries[acc])
        lin = lineage.get(taxid[acc], {}) if taxid[acc] else {}
        # the release run's rules (sanntis collect / steps.lineage_of): Eukaryota → euk, anything else (incl. the
        # giant-virus entry BGC0002846 and taxon-less loci) → proka / prokaryotes; the NCBI lineage splits eukaryotes
        kingdom = "euk" if lin.get("superkingdom") == "Eukaryota" else "proka"
        loci[acc] = Positive(acc, str(rec.seq).upper(), pairs, ";".join(classes), rec.id)
        rows.append({"entity_id": acc, "status": status, "kingdom": kingdom,
                     "lineage": "prokaryotes" if kingdom == "proka" else lineage_group(lin),
                     "ncbi_tax_id": taxid[acc] or -1, "class_raw": ";".join(classes),
                     "classes": classes, "length": len(rec.seq), "n_cds": len(pairs), "test_a": status == "active",
                     "test_b": True})
    tb = pd.DataFrame(rows)
    missing = tb[tb["lineage"].isna()]
    if not missing.empty:
        raise SystemExit(f"{len(missing)} test BGCs without an NCBI lineage: {list(missing['entity_id'][:5])}")
    return tb, loci, excluded


def pool_tables(background: Path) -> dict[str, pd.DataFrame]:
    b = background / "background"
    return {n: pd.read_parquet(b / f"{n}.parquet") for n in ("genomes", "cds", "clean_runs", "native_mask",
                                                              "test_backbones")}


# ── planning (every random draw) ─────────────────────────────────────────────


def gene_loci(starts: np.ndarray, ends: np.ndarray) -> list[list[int]]:
    """Positions (sorted by start) grouped into loci: overlapping CDS share a locus (release-run rule)."""
    order = np.lexsort((ends, starts))
    loci: list[list[int]] = []
    reach = -1
    for i in order:
        if loci and starts[i] < reach:
            loci[-1].append(int(i))
        else:
            loci.append([int(i)])
        reach = max(reach, int(ends[i]))
    return loci


def plan_isolate(pool: dict[str, pd.DataFrame], tests: pd.DataFrame, seed: int) -> list[Plan]:
    by_contig = {k: g.sort_values("cds_idx") for k, g in pool["cds"].groupby(["genome", "contig"], sort=False)}
    runs = pool["clean_runs"].set_index(["genome", "run_id"])
    plans = []
    for track in TRACKS:
        bb = pool["test_backbones"][pool["test_backbones"]["kingdom"] == track].reset_index(drop=True)
        pos = sorted(tests.loc[tests["kingdom"] == track, "entity_id"])
        if bb.empty or not pos:
            continue
        rng = rng_for(seed, track)
        assign: list[list[str]] = [[] for _ in range(len(bb))]
        for i, acc in enumerate(pos):
            assign[i % len(bb)].append(acc)
        for b, row in bb.iterrows():
            idx = [int(i) for i in runs.loc[(row["genome"], row["run_id"]), "cds_idx"]]
            cds = by_contig[(row["genome"], row["contig"])].set_index("cds_idx").loc[idx]
            lo, hi = int(cds["start"].min()), int(cds["end"].max())
            # intergenic gaps of the slice (sanntis_dev rule): the two ends, then every gap between CDS reaches
            gaps, reach = [(0, 0), (hi - lo, hi - lo)], 0
            for s, e in zip(cds["start"] - lo, cds["end"] - lo, strict=True):
                if s >= reach:
                    gaps.append((reach, int(s)))
                reach = max(reach, int(e))
            picks = rng.integers(0, len(gaps), size=len(assign[b]))
            cuts = sorted(int(rng.integers(gaps[i][0], gaps[i][1] + 1)) for i in picks)
            plans.append(Plan("isolate", track, int(b), row["genome"], row["contig"], row["run_id"],
                              [(lo, hi, idx)], {0: list(zip(cuts, assign[b], strict=True))}))
    return plans


def plan_shuffled(pool: dict[str, pd.DataFrame], tests: pd.DataFrame, seed: int, min_genes: int,
                  max_genes: int, insert_range: tuple[float, float] = (0.4, 0.6)) -> list[Plan]:
    kingdom = dict(zip(pool["genomes"]["genome"], pool["genomes"]["kingdom"], strict=True))
    by_contig = {k: g.sort_values("cds_idx") for k, g in pool["cds"].groupby(["genome", "contig"], sort=False)}
    plans = []
    for track in TRACKS:
        rng = rng_for(seed, track)
        gs = sorted(g for g, k in kingdom.items() if k == track)
        pools: dict[str, list[str]] = {}
        loci_of: dict[tuple[str, str], list[list[int]]] = {}
        for (g, c), df in by_contig.items():
            if g not in gs:
                continue
            loci = gene_loci(df["start"].to_numpy(), df["end"].to_numpy())
            if len(loci) >= min_genes:
                pools.setdefault(g, []).append(c)
                loci_of[(g, c)] = loci
        used = [g for g in gs if g in pools]
        pos = sorted(tests.loc[tests["kingdom"] == track, "entity_id"])
        if not used or not pos:
            continue
        log.info("shuffled %s: %d / %d genomes with a >= %d-locus contig (%d contigs), %d BGCs", track, len(used),
                 len(gs), min_genes, sum(len(v) for v in pools.values()), len(pos))
        order = rng.permutation(len(pos))
        cycle = {g: [sorted(pools[g])[i] for i in rng.permutation(len(pools[g]))] for g in used}
        seen = dict.fromkeys(used, 0)
        for k, pi in enumerate(order):
            g = used[k % len(used)]
            c = cycle[g][seen[g] % len(cycle[g])]
            seen[g] += 1
            df = by_contig[(g, c)]
            loci = loci_of[(g, c)]
            if len(loci) > max_genes:
                s0 = int(rng.integers(0, len(loci) - max_genes + 1))
                loci = loci[s0:s0 + max_genes]
            starts, ends, idx = df["start"].to_numpy(), df["end"].to_numpy(), df["cds_idx"].to_numpy()
            units, prev = [], int(min(starts[i] for i in loci[0]))
            for locus in loci:  # [previous locus end, locus end): the locus with its upstream spacer
                end = int(max(ends[i] for i in locus))
                units.append((prev, end, [int(idx[i]) for i in locus]))
                prev = end
            shuffled = [units[i] for i in rng.permutation(len(units))]
            lo, hi = int(np.floor(insert_range[0] * len(units))), int(np.ceil(insert_range[1] * len(units)))
            cut = int(rng.integers(lo, hi + 1))
            plans.append(Plan("shuffled", track, k, g, c, f"{c}_shuf{seen[g] - 1}", shuffled,
                              {cut: [(0, pos[int(pi)])]}))
    return plans


# ── assembly ─────────────────────────────────────────────────────────────────


class Assembler:
    def __init__(self, plan: Plan, native: set[int]):
        self.plan, self.native = plan, native
        self.pieces: list[str] = []
        self.features: list[SeqFeature] = []
        self.rows: list[dict] = []
        self.bed: list[dict] = []
        self.cursor = 0

    def add_cds(self, c: Cds, f: SeqFeature, shift: int, kind: str, entity: str, genome: str, bgc_class: str,
                native: bool) -> None:
        loc = f.location + shift
        cds_id = f"{entity}|{c.cds_id}" if kind == "positive" else c.cds_id
        translation = (f.qualifiers.get("translation") or [c.sequence])[0]
        quals = {"locus_tag": [cds_id], "protein_id": [cds_id],
                 **{q: f.qualifiers[q] for q in COPIED_QUALIFIERS if q in f.qualifiers}, "translation": [translation]}
        self.features.append(SeqFeature(loc, type="CDS", qualifiers=quals))
        self.rows.append({"contig_id": self.plan.contig_id, "start": c.start + shift, "end": c.end + shift,
                          "strand": "-" if c.strand == -1 else "+",
                          "parts": [[int(p.start), int(p.end)] for p in loc.parts],
                          "cds_id": cds_id, "source_cds_id": c.cds_id, "source_kind": kind, "entity_id": entity,
                          "source_genome_id": genome, "class": bgc_class,
                          "product": (f.qualifiers.get("product") or [""])[0], "translation": translation,
                          "native_bgc": native})

    def add_positive(self, p: Positive) -> None:
        self.pieces.append(p.seq)
        for c, f in p.cds:
            self.add_cds(c, f, self.cursor, "positive", p.entity_id, p.entity_id, p.bgc_class, False)
        self.bed.append({"contig_id": self.plan.contig_id, "start": self.cursor, "end": self.cursor + len(p.seq),
                         "entity_id": p.entity_id, "class": p.bgc_class})
        self.cursor += len(p.seq)

    def add_backbone(self, seq: str, a: int, b: int, pairs: list[tuple[Cds, SeqFeature]], idx: list[int]) -> None:
        """Source interval [a, b) with the CDS ``idx`` that lie inside it."""
        if b <= a:
            return
        self.pieces.append(seq[a:b])
        shift = self.cursor - a
        for i in idx:
            c, f = pairs[i]
            if c.start >= a and c.end <= b:
                self.add_cds(c, f, shift, "negative", self.plan.unit_id, self.plan.genome, "", i in self.native)
        self.cursor += b - a

    def finish(self) -> tuple[SeqRecord, pd.DataFrame, pd.DataFrame]:
        cid = self.plan.contig_id
        rec = SeqRecord(Seq("".join(self.pieces)), id=cid, name=cid, description=f"{self.plan.prefix} {cid}",
                        annotations={"molecule_type": "DNA", "topology": "linear"})
        rec.features = [SeqFeature(FeatureLocation(0, len(rec.seq)), type="source",
                                   qualifiers={"note": ["SanntiS 1.0 nucleotide test contig (objective 02)"]}),
                        *sorted(self.features, key=lambda f: (int(f.location.start), int(f.location.end)))]
        cds = pd.DataFrame(self.rows, columns=CDS_COLUMNS).sort_values(["start", "end"]).reset_index(drop=True)
        return rec, cds, pd.DataFrame(self.bed, columns=BED_COLUMNS)


def build(plan: Plan, rec: SeqRecord, pairs: list[tuple[Cds, SeqFeature]], loci: dict[str, Positive],
          native: set[int]) -> tuple[SeqRecord, pd.DataFrame, pd.DataFrame]:
    seq = str(rec.seq).upper()
    asm = Assembler(plan, native)
    for u, (a, b, idx) in enumerate(plan.units):
        prev = a
        for offset, acc in plan.inserts.get(u, []):  # offsets are relative to the unit start (isolate: the slice)
            asm.add_backbone(seq, prev, a + offset, pairs, idx)
            asm.add_positive(loci[acc])
            prev = a + offset
        asm.add_backbone(seq, prev, b, pairs, idx)
    return asm.finish()


def write_contig(out: Path, plan: Plan, rec: SeqRecord, cds: pd.DataFrame, bed: pd.DataFrame) -> dict:
    d = out / plan.set_name / plan.track
    d.mkdir(parents=True, exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)  # long LOCUS name (the record id)
        SeqIO.write([rec], d / f"{plan.prefix}.gbk", "genbank")
    cds.to_parquet(d / f"{plan.prefix}.cds.parquet", index=False)
    bed.to_csv(d / f"{plan.prefix}.ground_truth.bed", sep="\t", index=False)
    neg = cds["source_kind"] == "negative"
    return {"set": plan.set_name, "track": plan.track, "contig": plan.prefix, "record_id": plan.contig_id,
            "genome": plan.genome, "source_contig": plan.contig, "unit_id": plan.unit_id, "bp": len(rec.seq),
            "n_units": len(plan.units), "n_bgc": len(bed), "bgcs": ";".join(bed["entity_id"]),
            "bgc_cds": int((~neg).sum()), "backbone_cds": int(neg.sum()),
            "native_bgc_cds": int(cds.loc[neg, "native_bgc"].sum())}


README = """SanntiS 1.0 nucleotide test sets (objective 02; the input of objective 03's tool comparison).
<set>/<track>/<prefix>.gbk: contig + CDS (locus_tag = protein_id = CDS id, /translation), the input every tool
runs on; .cds.parquet: one row per CDS with source_kind (positive = test BGC, negative = backbone), entity_id
(MIBiG accession for positives), native_bgc (backbone CDS inside an antiSMASH 8.0.1 region of its genome, shuffled
set); .ground_truth.bed: BGC spans (0-based half-open). isolate: {n_iso} contigs ({iso_proka} proka + {iso_euk} euk
BGC-free backbones, several BGCs each); shuffled: one BGC per gene-shuffled genome contig ({n_shuf}). BGCs: MIBiG
5.0rc1 active + pending, accession >= BGC0002750, locus <= 200 kb ({n_bgc}); backgrounds: {n_genomes} genomes not in
SanntiS's negative pool v1 (objectives/00-data-collection/accessions/testset_background_genomes.tsv).
"""


def main(a: argparse.Namespace) -> None:
    out = Path(a.out)
    ts = out / "testsets"
    out.mkdir(parents=True, exist_ok=True)
    tests, loci, excluded = load_positives(Path(a.archives), a.limit)
    log.info("test BGCs: %d (Test A %d), excluded %d", len(tests), int(tests["test_a"].sum()), len(excluded))
    pool = pool_tables(Path(a.background))
    plans = plan_isolate(pool, tests, a.seed) + plan_shuffled(pool, tests, a.seed, a.min_genes, a.max_genes)
    placed = sorted(acc for p in plans for ins in p.inserts.values() for _, acc in ins)
    if placed != sorted(list(tests["entity_id"]) * 2):
        raise RuntimeError("not every test BGC is placed exactly once per set")
    native = {(g, c): set(d["cds_idx"].astype(int)) for (g, c), d in pool["native_mask"].groupby(["genome", "contig"])}

    for s in ("isolate", "shuffled"):
        for t in TRACKS:
            d = ts / s / t
            if d.exists():
                for f in d.iterdir():
                    f.unlink()
    stats = []
    by_genome: dict[str, list[Plan]] = {}
    for p in plans:
        by_genome.setdefault(p.genome, []).append(p)
    for n, (genome, ps) in enumerate(sorted(by_genome.items()), 1):
        recs = {r.id: r for r in read_records(Path(a.background) / "prepared" / genome / "input.gbk")}
        pairs_of: dict[str, list] = {}
        for p in ps:
            if p.contig not in pairs_of:
                pairs_of[p.contig] = cds_features(recs[p.contig])
            rec, cds, bed = build(p, recs[p.contig], pairs_of[p.contig], loci, native.get((genome, p.contig), set()))
            stats.append(write_contig(ts, p, rec, cds, bed))
        log.info("[%d/%d] %s: %d contigs", n, len(by_genome), genome, len(ps))

    st = pd.DataFrame(stats).sort_values(["set", "track", "contig"]).reset_index(drop=True)
    st.to_csv(out / "contigs.tsv", sep="\t", index=False)
    tests.to_parquet(out / "test_bgcs.parquet", index=False)
    counts = {s: {t: {"contigs": int(len(g)), "bgcs": int(g["n_bgc"].sum()), "bp": int(g["bp"].sum()),
                      "cds": int((g["bgc_cds"] + g["backbone_cds"]).sum()),
                      "native_bgc_cds": int(g["native_bgc_cds"].sum()),
                      "genomes": int(g["genome"].nunique())}
                  for t, g in st[st["set"] == s].groupby("track")} for s in ("isolate", "shuffled")}
    n_iso = {t: counts["isolate"].get(t, {}).get("contigs", 0) for t in TRACKS}
    (ts / "README.txt").write_text(README.format(
        n_iso=sum(n_iso.values()), iso_proka=n_iso["proka"], iso_euk=n_iso["euk"],
        n_shuf=sum(v["contigs"] for v in counts["shuffled"].values()), n_bgc=len(tests),
        n_genomes=len(pool["genomes"])))
    arch = Path(a.archives)
    manifest = {
        "params": {"seed": a.seed, "min_genes": a.min_genes, "max_genes": a.max_genes, "insert_range": [0.4, 0.6],
                   "max_locus_bp": MAX_LOCUS_BP, "v3_max": V3_MAX, "limit": a.limit},
        "sources": {"archives": {p.name: sha256(p) for p in sorted(arch.glob("mibig_5.0rc1-train_*.tar.gz"))},
                    "background_manifest": sha256(Path(a.background) / "background" / "manifest.json")},
        "test_bgcs": {"test_a": int(tests["test_a"].sum()), "test_b": int(tests["test_b"].sum()),
                      "by_kingdom": tests["kingdom"].value_counts().to_dict(),
                      "by_lineage": tests["lineage"].value_counts().to_dict(),
                      "excluded": excluded[["accession", "status", "length", "reason"]].to_dict("records")},
        "counts": counts}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    log.info("wrote %s: %s", out, json.dumps(counts))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    p.add_argument("--archives", required=True, help="objective 01's inputs/mibig_5.0rc1-train directory")
    p.add_argument("--background", required=True, help="objective 02's background build directory")
    p.add_argument("--out", required=True)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--min-genes", type=int, default=60)
    p.add_argument("--max-genes", type=int, default=400)
    p.add_argument("--limit", type=int, default=None, help="first N test BGCs (a trial)")
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    main(p.parse_args())
