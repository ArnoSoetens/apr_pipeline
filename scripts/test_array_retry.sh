#!/bin/bash
# Tests array_task.sh's rerun-after-filesystem-error logic with a fake python and a fake sleep.
# Runs nothing real (no Python, no pmemd), so it is fine on a login node.
#   bash test_array_retry.sh [/path/to/apr_pipeline/scripts/array_task.sh]
# Without an argument it uses array_task.sh next to this script, or in the current folder.
set -u
TASK="${1:-}"
if [ -z "$TASK" ]; then
    for candidate in "$(dirname "$0")/array_task.sh" ./array_task.sh; do
        if [ -f "$candidate" ]; then TASK="$candidate"; break; fi
    done
fi
if [ -z "$TASK" ] || [ ! -f "$TASK" ]; then
    echo "array_task.sh not found${TASK:+ at $TASK}; usage: bash $0 /path/to/scripts/array_task.sh" >&2
    exit 1
fi
TASK="$(realpath "$TASK")"
echo "testing $TASK"
T="$(mktemp -d)"
trap 'rm -rf "$T"' EXIT
mkdir -p "$T/bin" "$T/tmp"

cat > "$T/bin/sleep" <<'EOF'
#!/bin/bash
echo "(would sleep $1 s)"
EOF
cat > "$T/bin/python" <<'EOF'
#!/bin/bash
n=$(( $(cat "$T/count" 2>/dev/null || echo 0) + 1 )); echo "$n" > "$T/count"
shift
echo "fake python, call $n, items/args: $*"
case "$MODE" in
  crash)  [ "$n" -eq 1 ] && { echo "FileNotFoundError: [Errno 2] No such file or directory: '/fake/env/lib/python3.11/site-packages/x.py'" >&2; exit 1; } ;;
  item)   [ "$n" -eq 1 ] && { echo "iso_00002 FAILED: [Errno 2] No such file or directory: '/fake/env/lib/x.py'"; echo "iso_00001 done"; } ;;
  other)  echo "ValueError: a real bug" >&2; exit 1 ;;
  always) echo "FileNotFoundError: [Errno 2] No such file or directory: '/fake/amber/lib/y.py'" >&2; exit 1 ;;
esac
exit 0
EOF
chmod +x "$T/bin/"*
printf 'iso_00001\niso_00002\n' > "$T/worklist.txt"
echo 'DATA_DIR=/nonexistent' > "$T/run.env"

run() {   # run MODE STAGE EXPECTED_CALLS EXPECTED_EXIT ARGS...
    local mode=$1 stage=$2 want_calls=$3 want_exit=$4; shift 4
    rm -f "$T/count"
    echo "===== $mode ($stage)"
    T="$T" MODE="$mode" PATH="$T/bin:$PATH" TMPDIR="$T/tmp" ENV_DIR=/fake/env AMBERTOOLS_HOME=/fake/amber \
        SLURM_ARRAY_TASK_ID=0 bash "$TASK" "$stage" "$@"
    local rc=$? calls; calls=$(cat "$T/count" 2>/dev/null || echo 0)
    local left; left=$(ls -A "$T/tmp" | wc -l)
    if [ "$calls" = "$want_calls" ] && [ "$rc" = "$want_exit" ] && [ "$left" = 0 ]; then
        echo "PASS: $calls python call(s), exit $rc, no temp file left"
    else
        echo "FAIL: $calls python call(s) (want $want_calls), exit $rc (want $want_exit), $left temp file(s) left"
        FAILED=1
    fi
}
FAILED=0
run crash  analysis 2 0 "$T/worklist.txt" 2 0 --enhanced
run item   analysis 2 0 "$T/worklist.txt" 2 0
run other  analysis 1 1 "$T/worklist.txt" 2 0
run always analysis 3 1 "$T/worklist.txt" 2 0
run crash  sim      0 0 "$T/worklist.txt" 2 0 "$T/run.env"   # lines without '/': skipped, no pmemd
echo
[ "$FAILED" = 0 ] && echo "ALL PASSED" || echo "SOME TESTS FAILED"
