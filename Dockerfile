FROM python:3.12-slim
RUN apt-get update && apt-get install -y --no-install-recommends gfortran libfftw3-dev libblas-dev liblapack-dev && rm -rf /var/lib/apt/lists/*
WORKDIR /work
COPY pyproject.toml README.md LICENSE CHANGELOG.md /work/
COPY src /work/src
RUN pip install --no-cache-dir . pandas matplotlib pytest
COPY experiments /work/experiments
COPY figures /work/figures
COPY tests /work/tests
ENV RESULTS_DIR=/work/results PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1
CMD ["bash", "experiments/run_all.sh"]
