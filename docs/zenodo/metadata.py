"""Zenodo metadata of each deposit, as the deposition API takes it (`PUT /api/deposit/depositions/<id>` with
`{"metadata": ...}`): title, creators, description, licence (the `licence` column of env/zenodo.tsv), keywords,
version, related identifiers, the EUREMAP community and the EUREMAP grant. stage.py writes it as
`<out>/<deposit dir>.zenodo.json`, beside the staged directory, not in it: it is metadata, not a file of the record.

    python docs/zenodo/metadata.py <deposit>     # print one deposit's metadata (env/zenodo.tsv name)

Links between deposits (evaluation -> test sets, model -> training archives) are taken from the `record` column of
env/zenodo.tsv (a DOI, 10.5281/zenodo.<n>, or the bare record number); while a record is `TBD` its links are
left out and listed on stderr, so re-run after minting it and update the other record's metadata.
"""

import html
import json
import sys
from pathlib import Path

from notices import DEPOSITS, LICENCE_URL  # docs/zenodo/notices.py: what each deposit is, its licence URL

REPO = Path(__file__).resolve().parents[2]
ANALYSIS = "https://github.com/Finn-Lab/sanntisv1-data-and-analysis"
SANNTIS = "https://github.com/Finn-Lab/SanntiS"
HF_DATA = "https://huggingface.co/datasets/microbiome-informatics/sanntis-training-data"
HF_MODEL, MODEL_TAG = "https://huggingface.co/microbiome-informatics/sanntis", "v1.0.0-beta.1"

CREATORS = [  # Zenodo "Family, Given"; add "orcid": "0000-..." per person
    {"name": "Sanchez, Santiago", "affiliation": "EMBL - European Bioinformatics Institute"},
]

FUNDING = (
    "This work is part of EUREMAP (European Research Infrastructure for Marine Bioprospecting, https://euremap.eu), "
    "funded by the European Union's Horizon Europe research and innovation programme under grant agreement "
    "No 101131663 (doi:10.3030/101131663). Views and opinions expressed are however those of the author(s) only and "
    "do not necessarily reflect those of the European Union. Neither the European Union nor the granting authority "
    "can be held responsible for them.")
GRANT = "10.13039/501100000780::101131663"  # funder DOI (European Commission) :: grant agreement
COMMUNITY = "euremap"  # https://zenodo.org/communities/euremap
ZENODO_LICENCE = {"CC-BY-4.0": "cc-by-4.0", "Apache-2.0": "apache-2.0"}

KEYWORDS = ["biosynthetic gene clusters", "BGC detection", "natural products", "genome mining", "SanntiS", "MIBiG",
            "machine learning", "EUREMAP"]

META = {  # env/zenodo.tsv name -> (title, version, extra keywords)
    "mibig5_train": (
        "SanntiS 1.0 training archives: MIBiG 5.0rc1 active and pending entries in MIBiG format",
        "5.0rc1", ["training data"]),
    "nucleotide_testsets_v2": (
        "SanntiS 1.0 benchmark test sets (v2): MIBiG BGCs in isolate and shuffled genome backgrounds",
        "v2", ["benchmark", "test set"]),
    "v3active_eval_v2": (
        "SanntiS 1.0 MIBiG 3.1 validation: v3-active model and its evaluation on the v2 test sets",
        "v2", ["benchmark", "model evaluation"]),
    "sanntis_model": (
        "SanntiS 1.0 model bundle (release v1.0.0-beta.1)",
        "v1.0.0-beta.1", ["model", "deep learning"]),
}

LINKS = {  # deposit -> [(relation, other deposit)], resolved through the `record` column
    "v3active_eval_v2": [("references", "nucleotide_testsets_v2")],
    "sanntis_model": [("isDerivedFrom", "mibig5_train")],
}


def table() -> dict[str, dict[str, str]]:
    rows = [r.split("\t") for r in (REPO / "env" / "zenodo.tsv").read_text().splitlines()]
    return {r[0]: dict(zip(rows[0], r, strict=True)) for r in rows[1:]}


def doi(record: str) -> str | None:
    record = record.strip()
    if record in ("", "TBD"):
        return None
    return record if record.startswith("10.") else f"10.5281/zenodo.{record}"


def related(deposit: str, rows: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    out = [
        {"identifier": ANALYSIS, "relation": "isSupplementTo", "resource_type": "software", "scheme": "url"},
        {"identifier": SANNTIS, "relation": "isSupplementTo", "resource_type": "software", "scheme": "url"},
        {"identifier": "10.1101/2023.05.23.540769", "relation": "references",
         "resource_type": "publication-preprint", "scheme": "doi"},
    ]
    if "mibig" in DEPOSITS[deposit][1]:
        out.append({"identifier": "https://mibig.secondarymetabolites.org", "relation": "isDerivedFrom",
                    "resource_type": "dataset", "scheme": "url"})
    if deposit == "sanntis_model":  # the same files as the Hub release tag
        out.append({"identifier": f"{HF_MODEL}/tree/{MODEL_TAG}", "relation": "isIdenticalTo",
                    "resource_type": "dataset", "scheme": "url"})
    if deposit == "mibig5_train":  # same entries, packaged as one GBK and one JSON tarball on the Hub
        out.append({"identifier": f"{HF_DATA}/tree/mibig-5.0rc1-train", "relation": "isVariantFormOf",
                    "resource_type": "dataset", "scheme": "url"})
    for relation, other in LINKS.get(deposit, []):
        if (d := doi(rows[other]["record"])) is None:
            print(f"{deposit}: {relation} {other} left out (record TBD in env/zenodo.tsv)", file=sys.stderr)
        else:
            out.append({"identifier": d, "relation": relation, "resource_type": "dataset", "scheme": "doi"})
    return out


def description(deposit: str, row: dict[str, str]) -> str:
    what = DEPOSITS[deposit][0]
    e = html.escape
    return "".join([
        f"<p>{e(what[0].upper() + what[1:])}, from <a href=\"{ANALYSIS}\">SanntiS 1.0: data and analyses</a>, "
        "the runbooks that reproduce the data, models and analyses behind "
        f"<a href=\"{SANNTIS}\">SanntiS</a> 1.0.</p>",
        f"<p>Content: {e(row['content'])}.</p>",
        f"<p>Files: <code>{e(row['file'])}</code>; <code>SHA256SUMS</code>, its checksum (also recorded in "
        "<code>env/zenodo.tsv</code> of the repository, which the runbooks check on download); "
        "<code>NOTICE.md</code>, the licence and the third-party sources it carries, with what to cite.</p>",
        f"<p>Licence: <a href=\"{LICENCE_URL[row['licence']]}\">{e(row['licence'])}</a> for this project's "
        "contribution; third-party content keeps the terms of its source (<code>NOTICE.md</code>). No third-party "
        "software, container image or model weights are included.</p>",
        f"<p>Acknowledgement: {e(FUNDING)}</p>",
    ])


def metadata(deposit: str) -> dict:
    rows = table()
    row = rows[deposit]
    title, version, keywords = META[deposit]
    return {
        "upload_type": "dataset",
        "title": title,
        "creators": CREATORS,
        "description": description(deposit, row),
        "access_right": "open",
        "license": ZENODO_LICENCE[row["licence"]],
        "keywords": KEYWORDS + keywords,
        "version": version,
        "language": "eng",
        "related_identifiers": related(deposit, rows),
        "communities": [{"identifier": COMMUNITY}],
        "grants": [{"id": GRANT}],
        "notes": FUNDING,
    }


def write(deposit: str, path: Path) -> None:
    path.write_text(json.dumps({"metadata": metadata(deposit)}, indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    sys.stdout.write(json.dumps({"metadata": metadata(sys.argv[1])}, indent=2, ensure_ascii=False) + "\n")
