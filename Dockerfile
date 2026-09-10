FROM python:3.11-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/opt/carvefungi/bin \
    XDG_CACHE_HOME=/tmp/carvefungi-cache \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends libexpat1 libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/carvefungi
COPY requirements-reconstruction.txt ./
RUN pip install --no-cache-dir --only-binary=:all: -r requirements-reconstruction.txt
COPY bin/ ./bin/
COPY data/reactionDatabase/ ./data/reactionDatabase/
COPY tests/ ./tests/
COPY LICENSE README.md ./
RUN python -m unittest discover -s tests -v && rm -rf "$XDG_CACHE_HOME"

WORKDIR /work
ENTRYPOINT ["python", "/opt/carvefungi/bin/carvefungi.py"]
CMD ["--help"]
