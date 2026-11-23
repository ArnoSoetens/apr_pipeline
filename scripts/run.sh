#!/bin/bash
# Run APR windows: minimize -> production part 1 -> production part 2, a few windows at a time.
# Plain bash, no Python: the settings come from run.env, which scripts/run.py writes into
# <data>/inputs/ (from config.toml) before starting this script.
#
#   run.sh --env <data>/inputs/run.env [--cpus N] WINDOW_FOLDER ...
#
# Normally started by scripts/run.py (which picks the windows). --cpus N: CPUs to use; windows at
# the same time = N / RANKS (MPI ranks per window, 1 for serial pmemd and cuda). Without --cpus:
# SLURM_CPUS_PER_TASK, or all cores. A free CPU takes the next window, also from the next data point.
#
# Per window:  no topology -> SKIP;  done marker + production.nc -> SKIP (never rerun);
# otherwise the three pmemd steps, then the done marker. A failing step writes run_failed.txt in the
# window (date, host, step, exit code, last line of that step's output) and the other windows go on.
# Exit code 1 if a window failed.
set -uo pipefail

usage() { sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

ENV_FILE=""
CPUS=""
WINDOWS=()
while [ "$#" -gt 0 ]; do
    case "$1" in
        --env)  ENV_FILE="${2:?--env needs a file}"; shift 2 ;;
        --cpus) CPUS="${2:?--cpus needs a number}"; shift 2 ;;
        -h|--help) usage 0 ;;
        -*) echo "run.sh: unknown option $1" >&2; usage 2 ;;
        *)  WINDOWS+=("${1%/}"); shift ;;
    esac
done
if [ ! -f "$ENV_FILE" ]; then
    echo "run.sh: settings file not found: '${ENV_FILE}' (run scripts/run.py, or pass --env)" >&2
    exit 2
fi
# shellcheck source=/dev/null
source "$ENV_FILE"                     # PMEMD, INPUT_DIR, TOPOLOGY, DONE_MARKER, RANKS

read -r -a PMEMD_CMD <<< "$PMEMD"     # "mpirun -np 2 pmemd.MPI" -> one word per array item
CPUS="${CPUS:-${SLURM_CPUS_PER_TASK:-$(nproc)}}"
SLOTS=$(( CPUS / ${RANKS:-1} ))
[ "$SLOTS" -lt 1 ] && SLOTS=1

stamp() { date '+%H:%M:%S'; }

# Write run_failed.txt in the current window folder and return 1.
failed() {   # failed <step> <exit code> <output file>
    local last
    last="$(grep -v '^[[:space:]]*$' "$3" 2>/dev/null | tail -n 1)"
    {
        echo "date: $(date '+%F %T')"
        echo "host: $(hostname)"
        echo "step: $1"
        echo "exit code: $2"
        echo "last line of $3: ${last:-(no output)}"
    } > run_failed.txt
    echo "[$(stamp)] FAILED ${LABEL} in $1 (exit code $2; see ${WINDOW}/run_failed.txt)"
    return 1
}

# One window. Exit status: 0 done, 10 skipped (no topology), 11 skipped (already done), 1 failed.
run_window() {
    WINDOW="$1"
    LABEL="$(basename "$(dirname "$(dirname "$WINDOW")")")/$(basename "$WINDOW")"   # iso_00001/a000
    if [ ! -f "$WINDOW/${TOPOLOGY}.prmtop" ] || [ ! -f "$WINDOW/${TOPOLOGY}.rst7" ]; then
        echo "[$(stamp)] SKIP  ${LABEL} (no topology: setup not finished)"
        return 10
    fi
    if [ -f "$WINDOW/$DONE_MARKER" ] && [ -f "$WINDOW/production.nc" ]; then
        echo "[$(stamp)] SKIP  ${LABEL} (already done)"
        return 11
    fi
    cd "$WINDOW" || return 1          # never run pmemd in the wrong folder
    rm -f run_failed.txt
    echo "[$(stamp)] START ${LABEL}"

    "${PMEMD_CMD[@]}" -O \
        -i "$INPUT_DIR/minimize_nopbc.in" \
        -o minimize.out \
        -p ${TOPOLOGY}.prmtop \
        -c ${TOPOLOGY}.rst7 \
        -r minimize.rst7 \
        -ref ${TOPOLOGY}.rst7 \
        -x minimize.nc \
        -inf minimize.mdinfo || { failed minimize $? minimize.out; return 1; }

    "${PMEMD_CMD[@]}" -O \
        -i "$INPUT_DIR/production_nopbc1.in" \
        -o production.out \
        -p ${TOPOLOGY}.prmtop \
        -c minimize.rst7 \
        -r production.rst7 \
        -ref ${TOPOLOGY}.rst7 \
        -x production.nc \
        -v mdvel \
        -e production.mden \
        -inf production.mdinfo \
        -frc mdfrc || { failed production_part1 $? production.out; return 1; }

    "${PMEMD_CMD[@]}" -O \
        -i "$INPUT_DIR/production_nopbc2.in" \
        -o production2.out \
        -p ${TOPOLOGY}.prmtop \
        -c production.rst7 \
        -r production.rst7 \
        -ref ${TOPOLOGY}.rst7 \
        -x production.nc \
        -v mdvel \
        -e production.mden \
        -inf production.mdinfo \
        -frc mdfrc || { failed production_part2 $? production2.out; return 1; }

    touch "$DONE_MARKER"              # all three steps finished
    echo "[$(stamp)] DONE  ${LABEL}"
    return 0
}

echo "[$(stamp)] ${#WINDOWS[@]} window(s) to run, ${SLOTS} at a time (${CPUS} CPU(s), ${RANKS:-1} per window; ${PMEMD})"

n_done=0
n_failed=0
n_skipped=0
count() {   # count <exit status of run_window>
    case "$1" in
        0) n_done=$((n_done + 1)) ;;
        10|11) n_skipped=$((n_skipped + 1)) ;;
        *) n_failed=$((n_failed + 1)) ;;
    esac
}

running=0
for window in "${WINDOWS[@]}"; do
    run_window "$window" &
    running=$((running + 1))
    if [ "$running" -ge "$SLOTS" ]; then
        wait -n; count $?
        running=$((running - 1))
    fi
done
while [ "$running" -gt 0 ]; do
    wait -n; count $?
    running=$((running - 1))
done

echo "[$(stamp)] run finished: ${n_done} done, ${n_failed} failed, ${n_skipped} skipped"
[ "$n_failed" -eq 0 ]
