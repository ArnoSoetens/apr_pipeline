#!/usr/bin/env bash
#
# open_windows.sh  (lives in scripts/)
#
# For the given windows of one data point: run cpptraj to make a processed
# trajectory, then open ALL of them in one VMD session (one molecule each).
#
# Usage (from anywhere):
#   ./open_windows.sh <datapoint> <windows...>
#
#   ./open_windows.sh iso_00001 a000-a005
#   ./open_windows.sh iso_00001 a000 a003 p010-p012
#   ./open_windows.sh --traj production.nc iso_00001 a000-a005
#   ./open_windows.sh --no-cpptraj iso_00001 p000-p003     # reuse earlier cpptraj output
#   ./open_windows.sh --verbose iso_00001 a000-a005        # show cpptraj / VMD output
#
# Windows are read from  <scripts>/../data/<datapoint>/windows/<window>
#
# Requires: cpptraj (AmberTools) and vmd on PATH.

set -euo pipefail

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$(cd "$SCRIPT_DIR/../data" && pwd)"

TOPOLOGY_GLOB="*.prmtop"

# Trajectory to look at. Each window has minimize.nc, production.nc and
# production2.nc, so we name the one we want (change it here or with --traj).
TRAJECTORY="production.nc"

# Where cpptraj writes its result, inside each window folder.
OUTPUT_SUBDIR="analysis"
OUTPUT_TRAJ_NAME="traj_processed.nc"

# cpptraj options. The system is implicit solvent (no box), so no autoimage.
STRIP_MASK=""                 # e.g. ":DM1,DM2,DM3" to remove the dummy atoms
DO_AUTOIMAGE=false
DO_RMS_FIT=false
RMS_FIT_MASK="@CA,C,N,O"      # only used if DO_RMS_FIT=true

# Optional VMD settings (Tcl). Ignored if the file does not exist.
VMD_SETTINGS_SCRIPT="$SCRIPT_DIR/vmd_settings_multi.vmd"

# ---------------------------------------------------------------------------
# ARGUMENTS
# ---------------------------------------------------------------------------

RUN_CPPTRAJ=true
RUN_VMD=true
VERBOSE=false
POSITIONAL=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-cpptraj) RUN_CPPTRAJ=false ;;
        --no-vmd)     RUN_VMD=false ;;
        --verbose|-v) VERBOSE=true ;;
        --traj)       TRAJECTORY="$2"; shift ;;
        *)            POSITIONAL+=("$1") ;;
    esac
    shift
done

if [[ ${#POSITIONAL[@]} -lt 2 ]]; then
    echo "Usage: $0 [--traj FILE] [--no-cpptraj] [--no-vmd] [-v] <datapoint> <windows...>" >&2
    echo "  e.g. $0 iso_00001 a000-a005" >&2
    exit 1
fi

DATAPOINT="$(basename "${POSITIONAL[0]}")"
RAW_ARGS=("${POSITIONAL[@]:1}")
WINDOWS_DIR="$DATA_DIR/$DATAPOINT/windows"

if [[ ! -d "$WINDOWS_DIR" ]]; then
    echo "Error: no windows folder at $WINDOWS_DIR" >&2
    exit 1
fi

log() {
    if $VERBOSE; then
        echo "$@"
    fi
}

# ---------------------------------------------------------------------------
# EXPAND RANGES  (a000-a005 -> a000 a001 ... a005)
# ---------------------------------------------------------------------------

expand_token() {
    local token="$1"

    if [[ "$token" =~ ^([a-zA-Z]+)([0-9]+)-([a-zA-Z]+)([0-9]+)$ ]]; then
        local prefix1="${BASH_REMATCH[1]}"
        local num1="${BASH_REMATCH[2]}"
        local prefix2="${BASH_REMATCH[3]}"
        local num2="${BASH_REMATCH[4]}"
        local width=${#num1}

        if [[ "$prefix1" != "$prefix2" ]]; then
            echo "Error: range prefixes differ ('$prefix1' vs '$prefix2') in '$token'" >&2
            exit 1
        fi
        if (( 10#$num1 > 10#$num2 )); then
            echo "Error: range start > end in '$token'" >&2
            exit 1
        fi

        local i
        for (( i=10#$num1; i<=10#$num2; i++ )); do
            printf "%s%0${width}d\n" "$prefix1" "$i"
        done
    else
        echo "$token"
    fi
}

WINDOW_NAMES=()
for token in "${RAW_ARGS[@]}"; do
    while IFS= read -r name; do
        WINDOW_NAMES+=("$name")
    done < <(expand_token "$token")
done

log "== $DATAPOINT: requested windows (${#WINDOW_NAMES[@]}) =="
if $VERBOSE; then
    printf '  %s\n' "${WINDOW_NAMES[@]}"
fi

# ---------------------------------------------------------------------------
# PROCESS EACH WINDOW WITH CPPTRAJ
# ---------------------------------------------------------------------------

VALID_NAMES=()
VALID_TOPS=()
VALID_TRAJS=()

for name in "${WINDOW_NAMES[@]}"; do
    WINDOW_DIR="$WINDOWS_DIR/$name"

    if [[ ! -d "$WINDOW_DIR" ]]; then
        echo "Warning: skipping '$name' -- $WINDOW_DIR not found" >&2
        continue
    fi

    TOPOLOGY="$(find "$WINDOW_DIR" -maxdepth 1 -type f -name "$TOPOLOGY_GLOB" | sort | head -n 1)"
    if [[ -z "$TOPOLOGY" ]]; then
        echo "Warning: skipping '$name' -- no file matching '$TOPOLOGY_GLOB'" >&2
        continue
    fi

    RAW_TRAJ="$WINDOW_DIR/$TRAJECTORY"
    OUTPUT_DIR="$WINDOW_DIR/$OUTPUT_SUBDIR"
    OUTPUT_TRAJ="$OUTPUT_DIR/$OUTPUT_TRAJ_NAME"
    mkdir -p "$OUTPUT_DIR"

    if $RUN_CPPTRAJ; then
        if [[ ! -s "$RAW_TRAJ" ]]; then
            echo "Warning: skipping '$name' -- $TRAJECTORY not found (simulation not run?)" >&2
            continue
        fi

        CPPTRAJ_IN="$OUTPUT_DIR/cpptraj.in"
        {
            echo "parm $TOPOLOGY"
            echo "trajin $RAW_TRAJ"
            if [[ -n "$STRIP_MASK" ]]; then
                echo "strip $STRIP_MASK"
            fi
            if $DO_AUTOIMAGE; then
                echo "autoimage"
            fi
            if $DO_RMS_FIT; then
                echo "rms first $RMS_FIT_MASK"
            fi
            echo "trajout $OUTPUT_TRAJ netcdf"
            echo "run"
        } > "$CPPTRAJ_IN"

        log "== [$name] Running cpptraj =="
        if $VERBOSE; then
            ( cd "$OUTPUT_DIR" && cpptraj -i "$CPPTRAJ_IN" )
        else
            ( cd "$OUTPUT_DIR" && cpptraj -i "$CPPTRAJ_IN" ) > /dev/null
        fi

        if [[ ! -s "$OUTPUT_TRAJ" ]]; then
            echo "Error: cpptraj did not produce $OUTPUT_TRAJ for window '$name'" >&2
            exit 1
        fi
    else
        if [[ ! -s "$OUTPUT_TRAJ" ]]; then
            echo "Warning: skipping '$name' -- --no-cpptraj given but $OUTPUT_TRAJ does not exist" >&2
            continue
        fi
        log "[$name] Skipping cpptraj, using existing $OUTPUT_TRAJ"
    fi

    VALID_NAMES+=("$name")
    VALID_TOPS+=("$TOPOLOGY")
    VALID_TRAJS+=("$OUTPUT_TRAJ")
done

if [[ ${#VALID_NAMES[@]} -eq 0 ]]; then
    echo "Error: no valid windows to load." >&2
    exit 1
fi

log "== Windows ready for VMD (${#VALID_NAMES[@]}) =="

# ---------------------------------------------------------------------------
# LAUNCH VMD WITH ALL WINDOWS AS SEPARATE MOLECULES
# ---------------------------------------------------------------------------

if $RUN_VMD; then
    LOAD_SCRIPT="$(mktemp /tmp/vmd_load_XXXXXX.vmd)"

    {
        # Source the settings first so apply_window_reps exists before we call it.
        if [[ -f "$VMD_SETTINGS_SCRIPT" ]]; then
            echo "source {$VMD_SETTINGS_SCRIPT}"
        fi

        for i in "${!VALID_NAMES[@]}"; do
            echo "set molid [mol new {${VALID_TOPS[$i]}} type parm7 waitfor all]"
            echo "mol addfile {${VALID_TRAJS[$i]}} type netcdf waitfor all molid \$molid"
            echo "mol rename \$molid {${VALID_NAMES[$i]}}"
            if [[ -f "$VMD_SETTINGS_SCRIPT" ]]; then
                echo "apply_window_reps \$molid"
            fi
        done

        # Focus the first molecule, hide the rest, and fit the view.
        echo "catch {mol top 0}"
        echo "catch {hide_all_but_first}"
        echo "display resetview"
    } > "$LOAD_SCRIPT"

    log "== Launching VMD with ${#VALID_NAMES[@]} molecule(s) =="
    if $VERBOSE; then
        vmd -e "$LOAD_SCRIPT"
    else
        vmd -e "$LOAD_SCRIPT" > /dev/null 2>&1
    fi
else
    log "Skipping VMD (--no-vmd given)."
fi