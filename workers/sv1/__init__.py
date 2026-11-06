"""Shared workers of the SanntiS 1.0 analysis runbooks.

The metric definitions come from SanntiS (`sanntis.train.metrics`, re-exported by `sv1.metrics` together with the
test-set extras SanntiS 1.0.0b1 dropped: native-ORF rule, masked bootstrap, AUPRC); this package holds the
test-set code that is not part of operating SanntiS: test contigs, their scoring and strata, novelty,
baseline-tool parsers."""
