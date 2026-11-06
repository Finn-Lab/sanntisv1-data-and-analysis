"""The nucleotide test sets of this objective (objective 02's deposit ``sanntisv1_nucleotide_testsets_v2``: isolate
and gene-shuffled contigs on backgrounds disjoint from the training negatives): set names and kinds; the loaders live
in ``sv1.testsets``."""

from sv1.testsets import TRACKS, contigs, load_cds  # noqa: F401  (re-exported for the workers here)

SETS = ("isolate", "shuffled")  # isolate, gene-shuffled backgrounds
KIND = {"isolate": "isolate", "shuffled": "shuffled"}
