# plantiSMASH 1.0 (tag 1.0 = 869a16b; antiSMASH 3 based, Python 2.7): no published container. The environment copies
# the one the archived runs used (conda: python 2.7.15, hmmer 3.1b2, hmmer2 2.3.2, glimmer 3.02, glimmerhmm 3.0.4,
# biopython 1.76, diamond 0.8.36, prodigal 2.6.3, muscle 3.8.1551, ... + the tag's pip requirements).
#   docker build -t sanntisv1/plantismash:1.0 env/docker/plantismash-1.0
FROM mambaorg/micromamba:1.5.10
USER root
RUN micromamba install -y -n base -c conda-forge -c bioconda \
      python=2.7.15 hmmer=3.1b2 hmmer2=2.3.2 glimmer=3.02 glimmerhmm=3.0.4 biopython=1.76 diamond=0.8.36 \
      prodigal=2.6.3 muscle=3.8.1551 blast=2.16.0 numpy=1.16.5 scipy=1.2.1 networkx=1.11 lxml=4.5.0 fasttree \
      cd-hit git pip \
 && micromamba run -n base pip install "straight.plugin==1.4.0-post-1" "pyquery==1.2.9" cssselect==1.1.0 \
      helperlibs==0.2.2 pysvg==0.2.2 pyExcelerator==0.6.4.1 backports.lzma "python-louvain==0.5" "bcbio-gff==0.6.4" \
 && git clone https://github.com/plantismash/plantismash /opt/plantismash \
 && git -C /opt/plantismash checkout 869a16b5f839fefc212a04f77a68ed1205b9e1c8 \
 && micromamba clean -a -y
# press the bundled HMM profiles at build time (the runtime filesystem may be read-only)
RUN cd /opt/plantismash && micromamba run -n base python run_antismash.py --check-prereqs
ENV PATH=/opt/conda/bin:$PATH
WORKDIR /opt/plantismash
