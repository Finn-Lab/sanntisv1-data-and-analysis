# BGC-Prophet 0.1.2 from PyPI on Python 3.12. The biocontainer (0.1.2--pyhdfd78af_0) does not run the pipeline:
# - it pins fair-esm 1.0.2, which has no ESM2 models (Alphabet.get_batch_converter() takes 1 positional argument);
# - it runs Python 3.14, whose default multiprocessing start method (forkserver) cannot pickle the LMDB handle of
#   BGC-Prophet's dataset (TypeError: cannot pickle 'Environment' object). Python 3.12 still forks.
#   docker build -t sanntisv1/bgc-prophet:0.1.2 env/docker/bgc-prophet-0.1.2
FROM python:3.12-slim
RUN pip install --no-cache-dir --extra-index-url https://download.pytorch.org/whl/cpu \
      bgc-prophet==0.1.2 fair-esm==2.0.0 torch==2.5.1
