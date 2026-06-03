#!/usr/bin/env bash
set -euo pipefail

# Rebuild the FNIH pancreas API-minimal outputs from scratch after the pancreas
# map or LIGER inputs have changed.
#
# This removes only the pancreas pipeline output directories:
#   results/fnih/pancreas/tmp
#   results/fnih/pancreas/out
#
# It then delegates to run_fnih_pancreas_api_data_pipeline.sh. Override the same
# environment variables supported by that script, for example:
#   CELL_SAMPLE_FRACTION=1.0 bash cell_state_de/scripts/rebuild_fnih_pancreas_api_data_pipeline.sh
#   PIGEAN_ENABLE=0 bash cell_state_de/scripts/rebuild_fnih_pancreas_api_data_pipeline.sh

TISSUE_ROOT="${TISSUE_ROOT:-results/fnih/pancreas}"
CELL_STATE_DE_DIR="${CELL_STATE_DE_DIR:-cell_state_de}"

rm -rf "${TISSUE_ROOT}/tmp" "${TISSUE_ROOT}/out"
mkdir -p "${TISSUE_ROOT}"

exec "${CELL_STATE_DE_DIR}/scripts/run_fnih_pancreas_api_data_pipeline.sh"
