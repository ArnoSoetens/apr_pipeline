#!/bin/bash
# What ONE slurm array task does (started by job_array.slurm, submitted by submit_arrays.py).
#   array_task.sh setup    <worklist> <points per task>  <offset>                    lines: iso_00001
#   array_task.sh sim      <worklist> <windows per task> <offset> <run.env>          lines: iso_00001/a000
#   array_task.sh analysis <worklist> <points per task>  <offset> [analysis.py flags]  lines: iso_00001
#   array_task.sh implicit <worklist> <points per task>  <offset> [run_implicit.py flags]
#   array_task.sh collect
#   array_task.sh cleanup_logs <cleanup_batch_logs.py arguments>
# Task number k (array index + offset; indices start at 0) takes worklist lines k*N+1 ... k*N+N and
# hands them to the script; the scripts themselves skip what is not ready and say why.
# The sim stage needs no Python: run.env (written by submit_arrays.py at submit time) holds the
# data folder and run.sh's settings.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE="${1:?stage: setup, sim, analysis, implicit, collect or cleanup_logs}"

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
shift 4 || shift "$#"
EXTRA=("$@")              # sim: the run.env file; analysis/implicit: flags for the script
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
        RUN_ENV="${EXTRA[0]:?sim needs the run.env file}"
        # shellcheck source=/dev/null
        source "$RUN_ENV"                              # DATA_DIR (and run.sh's settings)
        FOLDERS=()
        for point in "${POINTS[@]}"; do
            if [[ "$point" != */* ]]; then
                echo "SKIP $point: not a data_point/window line"
                continue
            fi
            FOLDERS+=("$DATA_DIR/${point%%/*}/windows/${point#*/}")
        done
        if [ "${#FOLDERS[@]}" -gt 0 ]; then
            bash "$SCRIPT_DIR/run.sh" --env "$RUN_ENV" --cpus "${SLURM_CPUS_PER_TASK:-1}" "${FOLDERS[@]}"
        fi
        ;;
    analysis)
        # --no-collect: the collect job builds results.csv once all analysis tasks have ended
        python "$SCRIPT_DIR/analysis.py" --no-collect "${EXTRA[@]}" "${POINTS[@]}"
        ;;
    implicit)
        python "$SCRIPT_DIR/run_implicit.py" --cpus "${SLURM_CPUS_PER_TASK:-1}" "${EXTRA[@]}" "${POINTS[@]}"
        ;;
    *)
        echo "unknown stage: $STAGE" >&2
        exit 1
        ;;
esac
