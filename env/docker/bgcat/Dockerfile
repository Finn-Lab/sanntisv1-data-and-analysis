# BGCat (HassounLab, HEAD 4f77517; no releases, no container). Python 3.12 + esm 3.2.3 (the last release before that
# commit) + rdkit; the ESM C 600M weights are mounted at /opt/bgcat/data/weights (Hugging Face
# biohub/esmc-600m-2024-12 @ e4d83bc7, fetched by 02_tools.sh). One patch: embed.py hard-codes `.to('cuda:1')`
# (a second GPU); it is changed to the first GPU, or the CPU when there is none.
#   docker build -t sanntisv1/bgcat:4f77517 env/docker/bgcat
FROM python:3.12-slim
RUN apt-get update && apt-get install -y --no-install-recommends git && rm -rf /var/lib/apt/lists/* \
 && pip install --no-cache-dir esm==3.2.3 httpx rdkit joblib scikit-learn \
 && git clone https://github.com/HassounLab/BGCat /opt/bgcat \
 && git -C /opt/bgcat checkout 4f775172703f5c70883d7e4da3e00d5e257660a6 \
 && sed -i "s/.to('cuda:1')/.to('cuda' if torch.cuda.is_available() else 'cpu')/" /opt/bgcat/embed.py
WORKDIR /opt/bgcat
