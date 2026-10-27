#!/bin/bash
# What ONE slurm array task does.
#   array_task.sh setup <worklist> <points_per_task> [offset]     worklist lines: iso_00001
#   array_task.sh sim <worklist> <windows_per_task> [offset]      worklist lines: iso_00001/a000
#   array_task.sh analysis <worklist> <points_per_task> [offset] [--clean-after-analysis]
#                                                               worklist lines: iso_00001
#   array_task.sh collect
# Task number k (array index + offset; indices start at 0) handles worklist lines k*N+1 ... k*N+N.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE="${1:?stage: sim, analysis or collect}"

CONFIG_ENV="$(python "$SCRIPT_DIR/config_env.py")"
eval "$CONFIG_ENV"                                  # DATA_DIR, DONE_MARKER, ...

if [ "$STAGE" = "collect" ]; then
    python "$SCRIPT_DIR/collect_results.py"
    exit 0
fi
if [ "$STAGE" = "cleanup_logs" ]; then
    shift
    python "$SCRIPT_DIR/cleanup_batch_logs.py" "$@"
    exit 0
fi

WORKLIST="${2:?worklist file}"
if [[ ! -f "$WORKLIST" || ! -r "$WORKLIST" ]]; then
    echo "worklist is missing or unreadable: $WORKLIST" >&2
    exit 1
fi
PER_TASK="${3:?items per task}"
OFFSET="${4:-0}"          # when the list needs several arrays, each one has its own offset
CLEAN_AFTER_ANALYSIS="${5:-}"
TASK_ID=$(( ${SLURM_ARRAY_TASK_ID:?this must run inside a slurm array} + OFFSET ))

START=$(( TASK_ID * PER_TASK + 1 ))
END=$(( START + PER_TASK - 1 ))
mapfile -t POINTS < <(sed -n "${START},${END}p" "$WORKLIST")     # my lines of the worklist

echo "task $TASK_ID ($STAGE): ${#POINTS[@]} item(s), worklist lines $START-$END"
if [ "${#POINTS[@]}" -eq 0 ]; then
    exit 0
fi

case "$STAGE" in
    setup)
        python "$SCRIPT_DIR/main.py" "${POINTS[@]}"
        ;;
    sim)
        READY_WINDOWS=()
        for point in "${POINTS[@]}"; do
            data_point="${point%%/*}"
            window="${point#*/}"
            window_dir="$DATA_DIR/$data_point/windows/$window"
            if [[ "$point" != */* || ! -d "$window_dir" \
                  || ! -f "$window_dir/${COMPLEX_PREFIX}-dum.prmtop" \
                  || ! -f "$window_dir/${COMPLEX_PREFIX}-dum.rst7" ]]; then
                echo "SKIP $point: window setup is incomplete"
                continue
            fi
            READY_WINDOWS+=("$point")
        done
        if [ "${#READY_WINDOWS[@]}" -gt 0 ]; then
            bash "$SCRIPT_DIR/run.sh" "${READY_WINDOWS[@]}"
        fi
        ;;
    analysis)
        # only analyze data points whose windows are ALL finished
        READY=()
        for point in "${POINTS[@]}"; do
            complete=1
            for window in "$DATA_DIR/$point"/windows/*/; do
                [ -f "${window}${DONE_MARKER}" ] || { complete=0; break; }
            done
            if [ "$complete" -eq 1 ]; then
                READY+=("$point")
            else
                echo "SKIP $point: simulation not complete"
            fi
        done
        if [ "${#READY[@]}" -gt 0 ]; then
            if [[ "$CLEAN_AFTER_ANALYSIS" == "--clean-after-analysis" ]]; then
                python "$SCRIPT_DIR/analysis.py" --clean-after-analysis "${READY[@]}"
            else
                python "$SCRIPT_DIR/analysis.py" "${READY[@]}"
            fi
        fi
        ;;
    *)
        echo "unknown stage: $STAGE" >&2
        exit 1
        ;;
esac
