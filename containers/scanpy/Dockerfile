FROM python:3.12-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    XDG_CACHE_HOME=/tmp/.cache \
    MPLCONFIGDIR=/tmp/matplotlib \
    NUMBA_CACHE_DIR=/tmp/numba-cache

RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p \
        /tmp/.cache \
        /tmp/matplotlib \
        /tmp/numba-cache \
    && chmod 1777 \
        /tmp/.cache \
        /tmp/matplotlib \
        /tmp/numba-cache

RUN python -m pip install --no-cache-dir --upgrade pip \
    && python -m pip install --no-cache-dir \
        "scanpy==1.12.4"

# Fail the image build immediately if the environment is broken.
RUN python - <<'PY'
import scanpy as sc
import anndata
import numpy
import scipy
import numba

print("scanpy:", sc.__version__)
print("anndata:", anndata.__version__)
print("numpy:", numpy.__version__)
print("scipy:", scipy.__version__)
print("numba:", numba.__version__)
PY

CMD ["python3"]
