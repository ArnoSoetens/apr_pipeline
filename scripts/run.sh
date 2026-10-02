#!/bin/bash
set -euo pipefail

PMEMD="${PMEMD:-pmemd}"
MAX_JOBS="${SLURM_CPUS_PER_TASK:-4}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$(cd "$SCRIPT_DIR/../data" && pwd)"
TOP="BCD-UNL-dum"

# Run minimize -> equilibrate -> production for ONE window
run_window() {
    local window="$1"

    echo "[$(date '+%H:%M:%S')] START ${window}"
    cd "$window"

    "$PMEMD" -O \
        -i "$SCRIPT_DIR/minimize_nopbc.in" \
        -o minimize.out \
        -p ${TOP}.prmtop \
        -c ${TOP}.rst7 \
        -r minimize.rst7 \
        -ref ${TOP}.rst7 \
        -x minimize.nc \
        -inf minimize.mdinfo || return 1

    "$PMEMD" -O \
        -i "$SCRIPT_DIR/production_nopbc1.in" \
        -o production.out \
        -p ${TOP}.prmtop \
        -c minimize.rst7 \
        -r production.rst7 \
        -ref ${TOP}.rst7 \
        -x production.nc \
        -v mdvel \
        -e production.mden \
        -inf production.mdinfo \
        -frc mdfrc || return 1

    "$PMEMD" -O \
        -i "$SCRIPT_DIR/production_nopbc2.in" \
        -o production2.out \
        -p ${TOP}.prmtop \
        -c production.rst7 \
        -r production.rst7 \
        -ref ${TOP}.rst7 \
        -x production.nc \
        -v mdvel \
        -e production.mden \
        -inf production.mdinfo \
        -frc mdfrc || return 1

    echo "[$(date '+%H:%M:%S')] DONE  ${window}"
}

# One SMILES at a time, its windows in parallel
for datapoint in "$DATA_DIR"/*/; do

    [ -d "${datapoint}windows" ] || continue

    echo "[$(date '+%H:%M:%S')] START DATAPOINT ${datapoint}"

    running=0
    for window in "${datapoint}"windows/*/; do

        { run_window "$window" || echo "FAILED ${window}"; } &

        running=$((running + 1))
        if [ "$running" -ge "$MAX_JOBS" ]; then
            wait -n
            running=$((running - 1))
        fi
    done

    wait
    echo "[$(date '+%H:%M:%S')] FINISHED DATAPOINT ${datapoint}"
done