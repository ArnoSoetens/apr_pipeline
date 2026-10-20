#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Settings come from config.toml (see config_env.py): DATA_DIR, INPUT_DIR, COMPLEX_PREFIX,
# PMEMD, MAX_JOBS. Anything already exported in the environment takes priority.
CONFIG_ENV="$(python "$SCRIPT_DIR/config_env.py")"
eval "$CONFIG_ENV"

# Fill the Amber .in templates with the values from config.toml
python "$SCRIPT_DIR/make_inputs.py"

read -r -a PMEMD_CMD <<< "$PMEMD"     # "mpirun -np 2 pmemd.MPI" -> one word per array item
TOP="${COMPLEX_PREFIX}-dum"

# Run minimize -> equilibrate -> production for ONE window
run_window() {
    local window="$1"

    echo "[$(date '+%H:%M:%S')] START ${window}"
    cd "$window" || return 1     # never run pmemd in the wrong folder

    # a finished window is skipped, so a resubmitted array only redoes what is left
    if [ -f "$DONE_MARKER" ]; then
        echo "[$(date '+%H:%M:%S')] SKIP  ${window} (already done)"
        return 0
    fi

    "${PMEMD_CMD[@]}" -O \
        -i "$INPUT_DIR/minimize_nopbc.in" \
        -o minimize.out \
        -p ${TOP}.prmtop \
        -c ${TOP}.rst7 \
        -r minimize.rst7 \
        -ref ${TOP}.rst7 \
        -x minimize.nc \
        -inf minimize.mdinfo || return 1

    "${PMEMD_CMD[@]}" -O \
        -i "$INPUT_DIR/production_nopbc1.in" \
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

    "${PMEMD_CMD[@]}" -O \
        -i "$INPUT_DIR/production_nopbc2.in" \
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

    touch "$DONE_MARKER"      # all three steps finished
    echo "[$(date '+%H:%M:%S')] DONE  ${window}"
}

# Which windows to run. Each argument is a data point name (all its windows) or
# <data point>/<window> (just that window), e.g. iso_00001  or  iso_00001/a000.
# Without arguments: every window of every data point in DATA_DIR.
WINDOWS=()
if [ "$#" -eq 0 ]; then
    set -- $(cd "$DATA_DIR" && ls -d */ | tr -d /)
fi
for arg in "$@"; do
    if [[ "$arg" == */* ]]; then
        WINDOWS+=("$DATA_DIR/${arg%%/*}/windows/${arg#*/}/")     # iso_00001/a000 -> .../iso_00001/windows/a000/
    else
        for window in "$DATA_DIR/$arg"/windows/*/; do
            if [ -d "$window" ]; then
                WINDOWS+=("$window")
            fi
        done
    fi
done

echo "[$(date '+%H:%M:%S')] ${#WINDOWS[@]} window(s) to run, $MAX_JOBS at a time"

# All windows in one list, MAX_JOBS at a time: a free CPU takes the next window,
# also when that window belongs to the next data point.
running=0
for window in "${WINDOWS[@]}"; do

    { run_window "$window" || echo "FAILED ${window}"; } &

    running=$((running + 1))
    if [ "$running" -ge "$MAX_JOBS" ]; then
        wait -n
        running=$((running - 1))
    fi
done

wait
echo "[$(date '+%H:%M:%S')] FINISHED"
