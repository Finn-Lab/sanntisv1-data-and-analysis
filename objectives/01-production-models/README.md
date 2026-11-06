# 01 · Production models

**Objective:** train the SanntiS 1.0 production bundle with `sanntis train` on the MIBiG 5.0rc1 active + pending
entries (GBKs from the MIBiG 4.0 archive, else NCBI nuccore) and the negative pool v1, with the default protocol (CV +
tuning on the training pool itself).

```bash
bash objectives/01-production-models/01_inputs.sh   # MIBiG-format archives, ~2 min (MIBIG_TRAIN=deposit: from Zenodo)
bash objectives/01-production-models/02_train.sh    # sanntis train (NEGATIVES=<pool dir | Hub revision>, default v1)
```

Status (2026-11-02): inputs and training done. `02_train.sh` ran on Codon on 2026-10-30 (SLURM head job 14471873,
06:03–07:20 UTC, SanntiS 1.0.0b1 = `env/sanntis.ref`): 2,777 loci (2,229 proka / 548 euk), 8 epochs, t* 0.30
(BGC-precision floor 0.8 met). Out-of-fold validation (2,647 BGCs): Precision_T 0.698, Recall_T 0.895,
F1_macro_T 0.776, BGC precision 0.882, BGC recall 0.944 (95% CIs in `expected/validation.tsv`); NPC pathway top-1
0.868, superclass top-1 0.732 (5-fold CV, 1,883 regions). Frozen in `expected/` (below). `03_compare.sh` is still planned.

## Inputs
MIBiG 5.0rc1 JSON + 4.0 JSON / GBK archives and the nuccore loci (objective 00), negative pool v1 (`--negatives`),
InterProScan 6.0.3.0 / InterPro 107.0 data (`sanntis setup`). NPC labels are computed by `sanntis train` from the
compounds of the JSON entries (NPClassifier API, cached in `SANNTIS_CACHE`).

## Steps

| Script | Does | Output (`$DATA_ROOT/01/…`) |
|---|---|---|
| `01_inputs.sh` | worker `mibig_train_archives.py`: 5.0rc1 active + pending as a MIBiG-format archive pair, checked with SanntiS's reader for external data and against `expected/`; packs the Zenodo file | `inputs/mibig_5.0rc1-train/`, `deposit/sanntisv1_mibig5.0rc1_train.tar.gz` |
| `02_train.sh` | `sanntis train --mibig none --external-gbk … --external-json …`; checks the collected pool (2,777 loci) | `run/` (`bundle/`, `report/report.md`, `tune/tune.json`) |
| `03_compare.sh` | planned: a re-run's validation vs `expected/` | |

## Training data (`01_inputs.sh`)

MIBiG 5.0rc1 ships JSON only, so `sanntis train --mibig 5.0rc1` has no GBK archive to read, and SanntiS no longer
falls back to NCBI nuccore. The 5.0rc1 entries therefore go in as external data in the MIBiG format
(`mibig_5.0rc1-train_{gbk,json}.tar.gz`, one `<accession>.gbk` / `<accession>.json` per BGC):

| | Entries | |
|---|---|---|
| 5.0rc1 active + pending | 2,797 (2,598 + 199) | the pool of the release bundle (pending = active in 4.0, under re-review in 5.0) |
| GBK from the 4.0 archive, verbatim | 2,621 | first locus unchanged since 4.0 |
| GBK from nuccore, first 5.0rc1 locus | 174 | 159 not in the 4.0 archive + 15 whose locus changed since 4.0; cleaned as the release run did (`/pseudo` CDS, CDS cut by the slice and repeated CDS locations dropped) |
| excluded (`excluded.tsv`) | 2 | BGC0003182 (WGS contig without CDS), BGC0003249 (8.4 Mb genome); both drop out of training anyway |
| **packaged** | **2,795** | `sanntis train` then excludes the 18 loci > 200 kb → **2,777 loci**, as the release run |

JSON entries are the 5.0rc1 files, verbatim. Against the release run's `collect`, 2,762 of the 2,777 loci carry the
same proteins; the 15 changed loci are the difference (10 with other proteins; 5 only changed notation, e.g.
`0..0` vs `1..<end>`), since the release run kept their 4.0 GBKs.

The archives are written deterministically. `expected/mibig_train_manifest.tsv` lists every BGC (status, GBK source,
locus, length, CDS, md5 of the sequence and of the proteins as SanntiS reads them) and is what a rebuild must
reproduce; `env/zenodo.tsv` (`mibig5_train`) pins the deposit file, reproduced by a build with the same SanntiS and
Biopython (`provenance.json` records both).

## Expected
`expected/mibig_train_{manifest,excluded}.tsv`: the training archives (above). From the 2026-10-30 run (frozen
2026-11-02, `PROVENANCE.json`: source paths + sha256): `tune.tsv`, `validation_oof.tsv` (OOF metrics per stratum),
`validation.tsv` (the same with 95% CIs, bundle report), `validation_strata.tsv`, `validation_class_curves.tsv`,
`sweep.tsv`, `npc_validation.tsv`, `collect_counts.tsv` (the collected pool) and `bundle_{manifest,config,
npc_thresholds}.json`. A re-run with the same SanntiS and inputs should reproduce them up to training
nondeterminism; `03_compare.sh` (planned) is to check that.

`expected/prerepo_release/` holds the bundle published before this repository existed (`current_release_*.json`). Its
epochs (9) and t* (0.35) were tuned on the v3-active CV, not on this pool, so it is a reference, not a target. Against
the pre-split release run's v3-active OOF validation (`../02-mibig31-validation/expected/release_v1sets/
validation_oof.tsv`, 2,297 BGCs, t* 0.35), the reference recall and BGC-level values fall within this run's 95%
CIs in every stratum. Precision_T is lower (all 0.698 vs 0.715, proka 0.733 vs 0.746: reference above the CI; euk
0.565 vs 0.593, within), and so is F1_macro_T (all 0.776 vs 0.788). That is expected at the lower t* (0.30 vs 0.35) on a larger pool (2,647 vs 2,297
OOF BGCs), which are the t* / epoch changes (9 → 8 epochs) to report.
Accepted on 2026-11-02: the Precision_T miss (all and proka strata, reference just above the 95% CI) is
attributed to the lower t* and does not block the release.
