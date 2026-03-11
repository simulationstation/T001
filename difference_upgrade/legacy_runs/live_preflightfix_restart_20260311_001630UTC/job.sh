#!/usr/bin/env bash
set -euo pipefail
cd /home/primary/SUPERNOVA
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
exec .venv/bin/supernova-mission run-difference-upgrade \
  --pair-subset-path archive/pixel_pair_queue.parquet \
  --output-dir difference_upgrade/live \
  --sign-mode fade \
  --resume
