#!/usr/bin/env bash

set -eo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate triplicate-r

mkdir -p scripts/logs cache

Rscript scripts/run_triplicate_threshold_analysis.R \
  --prevalence 0.25 \
  --min-count 10 \
  --threshold-min 0 \
  --threshold-max 1 \
  --threshold-by 0.01 \
  --method pearson \
  2>&1 | tee scripts/logs/triplicate_threshold_analysis.log
