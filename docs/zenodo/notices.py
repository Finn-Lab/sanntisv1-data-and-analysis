"""NOTICE.md of each Zenodo deposit: its licence (the `licence` column of env/zenodo.tsv), the upstream sources whose
content it carries, and what to cite. stage.py writes it next to the tarball, as a file of its own in the record,
so the tarballs (and their sha256 in env/zenodo.tsv) stay as the runbooks build them.

    python docs/zenodo/notices.py <deposit>     # print one NOTICE.md (env/zenodo.tsv name)

Decisions: docs/FINDINGS.md F-00-1. No deposit carries GTDB taxonomy (checked 2026-11-04), so none is CC BY-SA.
"""

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ANALYSIS = "https://github.com/Finn-Lab/sanntisv1-data-and-analysis"

SOURCES = {
    "mibig": (
        "MIBiG (https://mibig.secondarymetabolites.org), CC BY 4.0. BGC loci, accessions and entry metadata of MIBiG "
        "5.0rc1 and 4.0. Cite: Zdouc et al., Nucleic Acids Research 2025, doi:10.1093/nar/gkae1115."),
    "ncbi": (
        "NCBI / INSDC (GenBank, RefSeq, Assembly, nuccore, Taxonomy). Genome and locus sequences, CDS coordinates and "
        "translations, taxonomy ids. NCBI places no restriction on the use or distribution of these data, but some "
        "submitters may claim patent, copyright or other intellectual-property rights in what they submitted "
        "(https://www.ncbi.nlm.nih.gov/home/about/policies/). This deposit does not relicense any INSDC record; the "
        "source accessions are kept in the files."),
    "interpro": (
        "InterPro (https://www.ebi.ac.uk/interpro), CC0 1.0. Domain annotations of the proteins. Cite: Blum et al., "
        "Nucleic Acids Research 2025, doi:10.1093/nar/gkae1082."),
    "antismash8": (
        "antiSMASH 8.0.1 regions of the background genomes (the `native_bgc` flag), output of antiSMASH (AGPL-3.0; the "
        "software is not included). Cite: Blin et al., Nucleic Acids Research 2025, doi:10.1093/nar/gkaf334."),
    "npclassifier": (
        "NPClassifier (https://npclassifier.gnps2.org) pathway / superclass labels of the MIBiG compounds, from its "
        "public API (no data licence stated; factual classifications). Cite: Kim et al., Journal of Natural Products "
        "2021, doi:10.1021/acs.jnatprod.1c00399."),
    "sanntis_model": (
        "The SanntiS model bundle, Apache-2.0 (its README.md, the model card), trained on MIBiG (CC BY 4.0) positives. "
        "Cite SanntiS: Sanchez et al., bioRxiv 2023, doi:10.1101/2023.05.23.540769."),
}

DEPOSITS = {  # env/zenodo.tsv name -> (what it is, sources it carries)
    "mibig5_train": ("MIBiG 5.0rc1 active + pending entries as a MIBiG-format GenBank + JSON archive pair (the SanntiS "
                     "1.0 training positives)", ["mibig", "ncbi"]),
    "nucleotide_testsets_v2": ("the SanntiS 1.0 nucleotide test sets: MIBiG test BGCs inserted into NCBI genome "
                               "backgrounds", ["mibig", "ncbi", "antismash8"]),
    "v3active_eval_v2": ("the v3-active SanntiS model and its evaluation on the v2 test sets",
                         ["sanntis_model", "mibig", "ncbi", "interpro", "npclassifier"]),
    "sanntis_model": ("the SanntiS 1.0 model bundle: detection ensemble, filter and class model, natural product "
                      "classifier, model card and validation report (Hugging Face Hub microbiome-informatics/sanntis, "
                      "release tag v1.0.0-beta.1)", ["sanntis_model", "mibig", "interpro", "npclassifier"]),
}

LICENCE_URL = {"CC-BY-4.0": "https://creativecommons.org/licenses/by/4.0/",
               "Apache-2.0": "https://www.apache.org/licenses/LICENSE-2.0"}


def licences() -> dict[str, str]:
    rows = [r.split("\t") for r in (REPO / "env" / "zenodo.tsv").read_text().splitlines()]
    col = rows[0].index("licence")
    return {r[0]: r[col] for r in rows[1:]}


def notice(deposit: str) -> str:
    what, sources = DEPOSITS[deposit]
    lic = licences()[deposit]
    lines = [f"# NOTICE: sanntisv1 `{deposit}`", "",
             f"Content: {what}. Built by {ANALYSIS} (env/zenodo.tsv, docs/zenodo/).", "",
             f"Licence: {lic} ({LICENCE_URL[lic]}), Copyright 2026 EMBL - European Bioinformatics Institute, for this "
             "project's contribution; third-party content keeps the terms of its source:", ""]
    lines += [f"- {SOURCES[s]}" for s in sources]
    lines += ["", "No third-party software, container image or model weights are included. Please cite SanntiS and "
              "the sources above.", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    sys.stdout.write(notice(sys.argv[1]))
