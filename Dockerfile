FROM python:3.11-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/opt/carvefungi/bin \
    XDG_CACHE_HOME=/tmp/carvefungi-cache \
    TF_NUM_INTRAOP_THREADS=1 \
    TF_NUM_INTEROP_THREADS=1 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends libexpat1 libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/carvefungi
COPY requirements-reconstruction.txt requirements-prediction.txt ./
RUN pip install --no-cache-dir --only-binary=:all: \
    -r requirements-reconstruction.txt -r requirements-prediction.txt
ADD --chmod=644 --checksum=sha256:2765ce4f7cbf3a0dac9bd67e00c013df2a792bf820f2a06fb5a534d2d46f0ca0 \
    https://zenodo.org/api/records/18953301/files/model_binary_output10_18.keras/content \
    ./data/localization/model_binary_output10_18.keras
COPY data/localization/README.md ./data/localization/README.md
COPY bin/ ./bin/
COPY data/reactionDatabase/ ./data/reactionDatabase/
COPY tests/ ./tests/
COPY LICENSE README.md ./
# TensorFlow's x86 wheels require AVX, often unavailable under ARM emulation.
ARG SKIP_MODEL_TEST=0
RUN chmod -R a+rX data \
    && CARVEFUNGI_SKIP_MODEL_TEST="$SKIP_MODEL_TEST" python -m unittest discover -s tests -v \
    && rm -rf "$XDG_CACHE_HOME"

WORKDIR /work
ENTRYPOINT ["python", "/opt/carvefungi/bin/carvefungi.py"]
CMD ["--help"]
