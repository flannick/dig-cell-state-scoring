#!/usr/bin/env bash
set -euo pipefail

# Cleanly rebuild liver API-minimal outputs after an updated map or updated
# LIGER results. This removes only liver pipeline intermediates/final outputs,
# not raw inputs under data/external.
#
# Override TISSUE_ROOT if needed:
#   TISSUE_ROOT=/path/to/liver bash cell_state_de/scripts/rebuild_fnih_liver_api_data_pipeline.sh

TISSUE_ROOT="${TISSUE_ROOT:-results/fnih/liver}"
rm -rf "${TISSUE_ROOT}/tmp" "${TISSUE_ROOT}/out"
bash cell_state_de/scripts/run_fnih_liver_api_data_pipeline.sh
