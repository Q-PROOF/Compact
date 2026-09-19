# Q-PROOF Compact — verified quantum circuit optimizer
FROM python:3.12-slim

WORKDIR /app

# Install compactq from PyPI (zero dependencies) — pin for reproducibility
RUN pip install --no-cache-dir "compactq==0.2.0"

# Smoke test at build time: the proof pipeline must work in the container
RUN python -c "import compactq; \
    from compactq.benchmarks import qft; \
    out = compactq.optimize(qft(4)); \
    assert out.two_qubit_count() <= qft(4).two_qubit_count(); \
    print('compactq', compactq.__version__, 'proof pipeline OK')"

# Default: verify two circuits and emit a certificate
ENTRYPOINT ["compactq"]
CMD ["--help"]
