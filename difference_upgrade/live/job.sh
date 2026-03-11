#!/usr/bin/env bash
set -euo pipefail
cd /home/primary/SUPERNOVA
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
exec .venv/bin/supernova-mission run-difference-upgrade \
  --output-dir difference_upgrade/live \
  --max-pairs 60 \
  --per-galaxy 4 \
  --compatibilities exact \
  --same-collection-only \
  --sign-mode fade \
  --resume
