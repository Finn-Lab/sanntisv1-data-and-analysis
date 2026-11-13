"""Stage the Zenodo deposits of env/zenodo.tsv (maintainers only; reviewers download them).

    python docs/zenodo/stage.py <out dir> [deposit ...]     # default: every deposit of env/zenodo.tsv

Deposits (env/zenodo.tsv names), each a deterministic tarball checked against the sha256 recorded in env/zenodo.tsv:
- mibig5_train: objective 01's 01_inputs.sh, the training archives            ($DATA_ROOT/01/deposit/)
- nucleotide_testsets_v3: objective 02's 04_testsets.sh, the test sets         ($DATA_ROOT/02/deposit/)
- v3active_eval_v3: objective 02's 06_report.sh, the MIBiG 3.1 evaluation      ($DATA_ROOT/02/deposit/)
- sanntis_model: the SanntiS 1.2.0 model bundle, built here from its Hugging Face Hub release tag (MODEL_TAG below),
  every file checked against the Hub's sha256; needs huggingface_hub (in the SanntiS environment).
Writes <out>/<deposit file stem>/<file>.tar.gz + SHA256SUMS + NOTICE.md (notices.py; licence from env/zenodo.tsv;
upload all three to the record) and <out>/<deposit file stem>.zenodo.json, the record's metadata (metadata.py; not a
file of the record). ZENODO_MIRROR=<out> then stands in for the records.
"""

import hashlib
import os
import shutil
import sys
from pathlib import Path

from metadata import write as write_metadata  # docs/zenodo/metadata.py: <deposit dir>.zenodo.json
from notices import notice  # docs/zenodo/notices.py: NOTICE.md next to each tarball

REPO = Path(__file__).resolve().parents[2]
DATA_ROOT = Path(os.environ.get("DATA_ROOT",
                                "/nfs/production/rdf/metagenomics/projects/sanntisv1-data-and-analysis/data"))
MODEL_REPO, MODEL_TAG = "microbiome-informatics/sanntis", "v1.2.0"
BUILT = {  # env/zenodo.tsv deposit -> the tarball a runbook built
    "nucleotide_testsets_v3": DATA_ROOT / "02" / "deposit" / "sanntisv1_nucleotide_testsets_v3.tar.gz",
    "v3active_eval_v3": DATA_ROOT / "02" / "deposit" / "sanntisv1_v3active_eval_v3.tar.gz",
    "mibig5_train": DATA_ROOT / "01" / "deposit" / "sanntisv1_mibig5.0rc1_train.tar.gz",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def model_tarball(dest: Path) -> None:
    """The Hub snapshot of MODEL_REPO at MODEL_TAG as a deterministic tarball (sanntis_model/<path>), each file
    checked against the sha256 the Hub records for it (LFS files) before packing."""
    sys.path.insert(0, str(REPO / "workers"))
    from huggingface_hub import HfApi, snapshot_download

    from sv1.archives import write_tar

    snap = Path(snapshot_download(MODEL_REPO, revision=MODEL_TAG))
    for f in HfApi().list_repo_tree(MODEL_REPO, revision=MODEL_TAG, recursive=True):
        lfs = getattr(f, "lfs", None)
        if lfs is not None:
            assert sha256(snap / f.path) == lfs.sha256, f"{f.path}: sha256 differs from the Hub"
    members = {f"sanntis_model/{p.relative_to(snap)}": p.read_bytes() for p in sorted(snap.rglob("*"))
               if p.is_file() and ".cache" not in p.relative_to(snap).parts}
    members["sanntis_model/SOURCE.txt"] = (
        f"https://huggingface.co/{MODEL_REPO}/tree/{MODEL_TAG}\nSanntiS model bundle, release tag {MODEL_TAG}; load "
        "it with `sanntis predict --model <this directory>` or from the Hub at that revision.\n").encode()
    dest.parent.mkdir(parents=True, exist_ok=True)
    write_tar(dest, members)


def stage(out: Path, deps: list[str]) -> None:
    rows = [line.rstrip("\n").split("\t") for line in (REPO / "env" / "zenodo.tsv").read_text().splitlines()]
    table = {r[0]: r for r in rows if r[0] != "deposit"}
    for dep in deps or table:
        _, _, name, want = table[dep][:4]
        dest = out / name.removesuffix(".tar.gz") / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not (dest.exists() and sha256(dest) == want):
            if dep == "sanntis_model":
                model_tarball(dest)
            else:
                assert BUILT[dep].name == name, f"{dep}: {BUILT[dep].name} != {name} (env/zenodo.tsv)"
                shutil.copyfile(BUILT[dep], dest)
        got = sha256(dest)
        if got != want:
            sys.exit(f"{dep}: sha256 {got} != {want} (env/zenodo.tsv) for {dest}; if the build is meant to "
                     "change, record the new sha256 in env/zenodo.tsv")
        (dest.parent / "SHA256SUMS").write_text(f"{got}  {name}\n")
        (dest.parent / "NOTICE.md").write_text(notice(dep))
        write_metadata(dep, out / f"{dest.parent.name}.zenodo.json")
        print(dest, dest.stat().st_size >> 20, "MB", got)


if __name__ == "__main__":
    stage(Path(sys.argv[1]).resolve(), sys.argv[2:])
