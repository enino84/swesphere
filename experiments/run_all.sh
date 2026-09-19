#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python exp01_stability.py
python exp02_regimes.py
python exp03_footprint.py
python ../figures/make_figures.py
