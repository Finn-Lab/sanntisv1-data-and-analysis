"""Product classification of the test BGCs: SanntiS 1.0 NPC (v3-active model, true loci) vs BGCat.

    python npc_compare.py --test-bgcs PARQUET --labels DIR --model BUNDLE --bgcat predictions.json \
        --npc-index index_v1.json --out DIR

Ground truth: the NPClassifier labels of the MIBiG test BGCs (Test B, ``classified`` only). Both methods are scored
in the SanntiS NPC label space: pathways as they are, superclasses folded as in SanntiS (a superclass the model does
not predict becomes ``other <pathway>``) for the truth and for both rankings. BGCat ranks NPClassifier classes; its
ranking is rolled up to superclasses and pathways with the NPClassifier hierarchy (first occurrence). Scores: top-1
pathway; top-1 / top-3 superclass (BGCs with a superclass label); coverage (share of the labelled BGCs a method
classifies). Rows: each method on the BGCs it covers, and both on the BGCs BGCat covers; with ``--mibig4-json``
also on the covered BGCs that are not in MIBiG 4.0 (BGCat's training data may include MIBiG 4.0 entries).
"""

import argparse
import json
import re
import tarfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sanntis.model.bundle import SanntisModel
from sanntis.train.npc_train import LabelSpace, labelled


def rolled_up(classes: list[str], index: dict) -> tuple[list[str], list[str]]:
    """NPClassifier class ranking → (pathway ranking, superclass ranking), first occurrence kept."""
    cls_idx = index["Class"]
    sup_names = {v: k for k, v in index["Superclass"].items()}
    path_names = {v: k for k, v in index["Pathway"].items()}
    pw, sc = [], []
    for c in classes:
        h = index["Class_hierarchy"].get(str(cls_idx.get(c, -1)), {})
        pw += [path_names[i] for i in h.get("Pathway", []) if path_names[i] not in pw]
        sc += [sup_names[i] for i in h.get("Superclass", []) if sup_names[i] not in sc]
    return pw, sc


def unique(xs):
    return list(dict.fromkeys(xs))


def scores(ids, pw_rank, sc_rank, truth_pw, truth_sc) -> dict:
    top1 = [bool(pw_rank[i]) and pw_rank[i][0] in truth_pw[i] for i in ids]
    with_sc = [i for i in ids if truth_sc[i]]
    s1 = [bool(sc_rank[i]) and sc_rank[i][0] in truth_sc[i] for i in with_sc]
    s3 = [bool(set(sc_rank[i][:3]) & truth_sc[i]) for i in with_sc]
    return {"n": len(ids), "pathway_top1": float(np.mean(top1)) if ids else np.nan, "n_superclass": len(with_sc),
            "superclass_top1": float(np.mean(s1)) if with_sc else np.nan,
            "superclass_top3": float(np.mean(s3)) if with_sc else np.nan}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--test-bgcs", type=Path, required=True)
    p.add_argument("--labels", type=Path, required=True)
    p.add_argument("--model", type=Path, required=True)
    p.add_argument("--bgcat", type=Path, required=True)
    p.add_argument("--npc-index", type=Path, required=True)
    p.add_argument("--mibig4-json", type=Path, help="MIBiG 4.0 JSON archive: adds the subset of BGCs not in 4.0")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    tb = pd.read_parquet(a.test_bgcs).set_index("entity_id")
    labels = pd.read_parquet(a.labels / "npc_labels.parquet").set_index("accession")
    parents = json.loads((a.labels / "npc_superclass_parents.json").read_text())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        npc = SanntisModel.from_pretrained(a.model).npc
    if npc is None:  # sanntis train switches the NPC head off below 50 NPC-labelled BGCs (a LIMIT trial)
        print(f"{a.model}: no NPC model in the bundle; no classification comparison")
        return
    ls = LabelSpace(superclasses=list(npc.superclasses), parents=parents, pathways=list(npc.pathways))
    ids = labelled(labels, sorted(tb.index[tb["test_b"]]))
    truth_pw = {i: set(labels.at[i, "pathway"]) for i in ids}
    truth_sc = {i: set(ls.fold(labels.at[i, "superclass"], list(labels.at[i, "pathway"]))) for i in ids}

    pp, ps = npc.proba([list(tb.at[i, "cds_tokens"]) for i in ids], tb.loc[ids, "kingdom"].to_numpy())
    pws, scs = list(npc.pathways), list(npc.superclasses)
    s_pw = {i: [pws[j] for j in np.argsort(-pp[r])] for r, i in enumerate(ids)}
    s_sc = {i: [scs[j] for j in np.argsort(-ps[r])] for r, i in enumerate(ids)}

    pred = json.loads(a.bgcat.read_text())
    index = json.loads(a.npc_index.read_text())
    b_pw, b_sc = {}, {}
    for i in ids:
        if i in pred["bgcs"]:
            pw, sc = rolled_up([pred["labels"][k] for k in pred["bgcs"][i]], index)
            b_pw[i] = [x for x in pw if x in set(pws)]
            b_sc[i] = unique(s for x in sc for s in ls.fold([x], pw[:1]))
    common = [i for i in ids if i in b_pw]

    subsets = [("all labelled", ids, False), ("BGCat-covered", common, True)]
    if a.mibig4_json:
        with tarfile.open(a.mibig4_json) as t:
            v4 = {m.group(0) for n in t.getnames() if (m := re.search(r"BGC\d{7}", n))}
        subsets.append(("BGCat-covered, not in MIBiG 4.0", [i for i in common if i not in v4], True))
    rows = []
    for subset, sub, both in subsets:
        rows.append({"method": "SanntiS 1.0 NPC", "subset": subset, "coverage": len(sub) / len(ids),
                     **scores(sub, s_pw, s_sc, truth_pw, truth_sc)})
        if both:
            rows.append({"method": "BGCat", "subset": subset, "coverage": len(sub) / len(ids),
                         **scores(sub, b_pw, b_sc, truth_pw, truth_sc)})
    out = pd.DataFrame(rows)
    out.to_csv(a.out / "classification.tsv", sep="\t", index=False, float_format="%.6g")
    pd.DataFrame([{"entity_id": i, "truth_pathway": ";".join(sorted(truth_pw[i])),
                   "truth_superclass": ";".join(sorted(truth_sc[i])), "sanntis_pathway": s_pw[i][0],
                   "sanntis_superclass": s_sc[i][0], "bgcat_pathway": (b_pw.get(i) or [""])[0],
                   "bgcat_superclass": (b_sc.get(i) or [""])[0]} for i in ids]).to_csv(
        a.out / "classification_per_bgc.tsv", sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
