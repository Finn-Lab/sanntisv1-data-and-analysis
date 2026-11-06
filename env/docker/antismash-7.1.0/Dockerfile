# antiSMASH 7.1.0 with its databases: the biocontainer ships none and its read-only filesystem blocks the one-time
# data preparation, so this layer downloads the 7.1.0 database set (versioned URLs + sha256 checks hard-coded in
# antiSMASH 7.1.0: Pfam 35.0, ClusterBlast 4.0, KnownClusterBlast / ClusterCompare MIBiG 3.1, CompaRiPPson
# ASDB 4.0, Resfams, TIGRFam, Stachelhaus 1.1, NRPS-SVM 2.0, transATor 2023.02.23) and prepares them.
#   docker build -t sanntisv1/antismash:7.1.0 env/docker/antismash-7.1.0
FROM quay.io/biocontainers/antismash:7.1.0--pyh7e72e81_1@sha256:1f1db64317d9f08b1049c10d118925cddbabac3fbc3b6ddcfc6a3a87d5e6729f
RUN download-antismash-databases --database-dir /opt/antismash-databases \
 && antismash --databases /opt/antismash-databases --check-prereqs
ENV ANTISMASH_DATABASES=/opt/antismash-databases
