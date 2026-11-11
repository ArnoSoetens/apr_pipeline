#!/usr/bin/env bash
#
# open_windows.sh  (lives in scripts/)
#
# For the given windows of one data point: run cpptraj to make a processed
# trajectory, then open ALL of them in one VMD session (one molecule each).
#
# Usage (from anywhere):
#   ./open_windows.sh <datapoint> <windows...>
#   ./open_windows.sh --complex-only iso_00001 - iso_00010
#   ./open_windows.sh --complex-only
#   ./open_windows.sh --failed-clashes
#   ./open_windows.sh --failed-clashes iso_00001-iso_00010
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

# DATA_DIR and TRAJECTORY come from config.toml (see config_env.py)
CONFIG_ENV="$(python "$SCRIPT_DIR/config_env.py")"
eval "$CONFIG_ENV"

TOPOLOGY_NAME="${COMPLEX_PREFIX}-dum.prmtop"

# Trajectory to look at: analysis.trajectory in config.toml, or --traj FILE.

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
COMPLEX_ONLY=false
FAILED_CLASHES=false
POSITIONAL=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-cpptraj) RUN_CPPTRAJ=false ;;
        --no-vmd)     RUN_VMD=false ;;
        --verbose|-v) VERBOSE=true ;;
        --complex-only) COMPLEX_ONLY=true ;;
        --failed-clashes) FAILED_CLASHES=true ;;
        --traj)       TRAJECTORY="$2"; shift ;;
        *)            POSITIONAL+=("$1") ;;
    esac
    shift
done

if $COMPLEX_ONLY && $FAILED_CLASHES; then
    echo "Error: --complex-only and --failed-clashes cannot be used together." >&2
    exit 1
fi

if ! $COMPLEX_ONLY && ! $FAILED_CLASHES && [[ ${#POSITIONAL[@]} -lt 2 ]]; then
    echo "Usage: $0 [--traj FILE] [--no-cpptraj] [--no-vmd] [-v] <datapoint> <windows...>" >&2
    echo "       $0 --complex-only [datapoint/range...]" >&2
    echo "       $0 --failed-clashes [datapoint/range...]" >&2
    echo "  e.g. $0 iso_00001 a000-a005" >&2
    exit 1
fi

log() {
    if $VERBOSE; then
        echo "$@"
    fi
}

expand_datapoint_args() {
    EXPANDED_DATAPOINTS=()
    local -a ranges=()
    local i token start end prefix1 prefix2 start_num end_num width number

    i=0
    while (( i < ${#POSITIONAL[@]} )); do
        token="${POSITIONAL[$i]}"
        if [[ "$token" == "-" ]]; then
            if (( i == 0 || i + 1 >= ${#POSITIONAL[@]} )); then
                echo "Error: a spaced range must be written as START - END." >&2
                return 1
            fi
            token="${POSITIONAL[$((i - 1))]}-${POSITIONAL[$((i + 1))]}"
            ranges[${#ranges[@]}-1]="$token"
            ((i += 2))
            continue
        fi
        ranges+=("$token")
        ((i += 1))
    done

    for token in "${ranges[@]}"; do
        token="$(basename "$token")"
        if [[ "$token" =~ ^([[:alpha:]_]+)([0-9]+)-([[:alpha:]_]+)([0-9]+)$ ]]; then
            prefix1="${BASH_REMATCH[1]}"
            start_num="${BASH_REMATCH[2]}"
            prefix2="${BASH_REMATCH[3]}"
            end_num="${BASH_REMATCH[4]}"
            if [[ "$prefix1" != "$prefix2" ]]; then
                echo "Error: range prefixes differ in '$token'." >&2
                return 1
            fi
            start=$((10#$start_num))
            end=$((10#$end_num))
            if (( start > end )); then
                echo "Error: range start is greater than end in '$token'." >&2
                return 1
            fi
            width=${#start_num}
            for ((number = start; number <= end; number++)); do
                EXPANDED_DATAPOINTS+=("$(printf "%s%0*d" "$prefix1" "$width" "$number")")
            done
        elif [[ "$token" =~ ^[[:alpha:]_]+[0-9]+$ ]]; then
            EXPANDED_DATAPOINTS+=("$token")
        else
            EXPANDED_DATAPOINTS+=("$token")
        fi
    done
}

run_pdb_mode() {
    local mode="$1"
    local datapoint pdb_path report_path guest_path guest_atom_index
    local -a valid_names=() valid_pdbs=() valid_indices=()

    if [[ ${#POSITIONAL[@]} -gt 0 ]]; then
        expand_datapoint_args || return 1
        local -a selected_datapoints=("${EXPANDED_DATAPOINTS[@]}")
    elif [[ "$mode" == "complex" ]]; then
        local -a selected_datapoints=()
        while IFS= read -r -d '' pdb_path; do
            selected_datapoints+=("$(basename "$(dirname "$pdb_path")")")
        done < <(find "$DATA_DIR" -mindepth 2 -maxdepth 2 -type f -name "$COMPLEX_PREFIX.pdb" -print0 | sort -z)
    else
        local -a selected_datapoints=()
        while IFS= read -r -d '' report_path; do
            selected_datapoints+=("$(basename "$(dirname "$report_path")")")
        done < <(find "$DATA_DIR" -mindepth 2 -maxdepth 2 -type f -name setup_failed.txt -print0 | sort -z)
    fi

    if [[ "$mode" == "complex" ]]; then
        for datapoint in "${selected_datapoints[@]}"; do
            pdb_path="$DATA_DIR/$datapoint/$COMPLEX_PREFIX.pdb"
            if [[ ! -s "$pdb_path" ]]; then
                echo "Warning: skipping '$datapoint' -- complex PDB not found at $pdb_path" >&2
                continue
            fi
            valid_names+=("$datapoint")
            valid_pdbs+=("$pdb_path")
        done
    else
        if [[ ! -s "$HOST_PDB" ]]; then
            echo "Error: host PDB not found at $HOST_PDB" >&2
            return 1
        fi
        for datapoint in "${selected_datapoints[@]}"; do
            report_path="$DATA_DIR/$datapoint/setup_failed.txt"
            guest_path="$DATA_DIR/$datapoint/$GUEST_PREFIX.pdb"
            if [[ ! -s "$report_path" || ! -s "$guest_path" ]]; then
                echo "Warning: skipping '$datapoint' -- setup failure report or guest PDB is missing" >&2
                continue
            fi
            guest_atom_index="$(sed -nE 's/.*guest atom ([0-9]+).*/\1/p' "$report_path" | head -n 1)"
            if [[ ! "$guest_atom_index" =~ ^[0-9]+$ ]]; then
                echo "Warning: skipping '$datapoint' -- no guest atom index in $report_path" >&2
                continue
            fi
            valid_names+=("$datapoint")
            valid_pdbs+=("$guest_path")
            valid_indices+=("$guest_atom_index")
        done
    fi

    if [[ ${#valid_names[@]} -eq 0 ]]; then
        echo "Error: no valid PDBs to load for --$mode mode." >&2
        return 1
    fi

    log "== PDBs ready for VMD (${#valid_names[@]}) =="
    if ! $RUN_VMD; then
        log "Skipping VMD (--no-vmd given)."
        return 0
    fi

    local load_script
    load_script="$(mktemp /tmp/vmd_load_XXXXXX.vmd)"
    {
        if [[ -f "$VMD_SETTINGS_SCRIPT" ]]; then
            echo "source {$VMD_SETTINGS_SCRIPT}"
        fi
        if [[ "$mode" == "complex" ]]; then
            for i in "${!valid_names[@]}"; do
                echo "set molid [mol new {${valid_pdbs[$i]}} type pdb waitfor all]"
                echo "mol rename \$molid {${valid_names[$i]}}"
                if [[ -f "$VMD_SETTINGS_SCRIPT" ]]; then
                    echo "apply_complex_reps \$molid"
                fi
            done
        else
            echo "set_cycle_offset 1"
            echo "set host_molid [mol new {$HOST_PDB} type pdb waitfor all]"
            echo "mol rename \$host_molid {$HOST_PREFIX}"
            for i in "${!valid_names[@]}"; do
                echo "set molid [mol new {${valid_pdbs[$i]}} type pdb waitfor all]"
                echo "mol rename \$molid {${valid_names[$i]}}"
                if [[ -f "$VMD_SETTINGS_SCRIPT" ]]; then
                    echo "apply_clash_reps \$molid ${valid_indices[$i]}"
                fi
            done
        fi
        echo "catch {mol top 0}"
        echo "catch {hide_all_but_first}"
        echo "display resetview"
    } > "$load_script"

    log "== Launching VMD with ${#valid_names[@]} PDB molecule(s) =="
    if $VERBOSE; then
        vmd -e "$load_script"
    else
        vmd -e "$load_script" > /dev/null 2>&1
    fi
}

if $COMPLEX_ONLY; then
    run_pdb_mode complex
    exit $?
fi

if $FAILED_CLASHES; then
    run_pdb_mode failed-clashes
    exit $?
fi

DATAPOINT="$(basename "${POSITIONAL[0]}")"
RAW_ARGS=("${POSITIONAL[@]:1}")
WINDOWS_DIR="$DATA_DIR/$DATAPOINT/windows"

if [[ ! -d "$WINDOWS_DIR" ]]; then
    echo "Error: no windows folder at $WINDOWS_DIR" >&2
    exit 1
fi

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

    TOPOLOGY="$WINDOW_DIR/$TOPOLOGY_NAME"
    if [[ ! -f "$TOPOLOGY" ]]; then
        TOPOLOGY="$WINDOWS_DIR/a000/$TOPOLOGY_NAME"
    fi
    if [[ ! -f "$TOPOLOGY" ]]; then
        echo "Warning: skipping '$name' -- topology missing in this window and windows/a000" >&2
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