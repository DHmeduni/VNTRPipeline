#!/bin/bash

parse_vntr_config() {
    local target_vntr="$1"

    # Get absolute path of the script
    local script="${BASH_SOURCE[0]}"
    local script_path
    script_path=$(dirname "$(readlink -f "$script" 2>/dev/null || echo "$script")")

    # Fallback if readlink is unavailable
    if [[ -z $script_path ]]; then
        script_path="$(dirname "$0")"
    fi

    # Path to config file
    local lib_path
    lib_path="${LIB_PATH:-$(dirname "$script_path")/lib}"
    local CONFIG_FILE="${CONFIG_FILE:-"$lib_path/vntr_variables.cfg"}"

    # Check for VNTR argument
    if [[ -z "$target_vntr" ]]; then
        echo "❌ Please specify the VNTR name as an argument."
        echo "Usage: parse_vntr_config <VNTR_NAME>"
        return 1
    fi

    target_vntr="${target_vntr^^}"  # Convert to uppercase
    local found_vntr=0

    # Parse config file
    while read -r line || [[ -n "$line" ]]; do
        [[ -z "$line" ]] && continue

        # Section header: e.g., "#MUC1"
        if [[ "$line" =~ ^#([A-Za-z0-9_]+) ]]; then
            current_vntr="${BASH_REMATCH[1]}"
            if [[ "${current_vntr^^}" == "$target_vntr" ]]; then
                found_vntr=1
            else
                found_vntr=0
            fi
            continue
        fi

        # Skip lines outside target VNTR or that are comments
        [[ "$found_vntr" -eq 0 || "$line" =~ ^# ]] && continue

        # Parse key and value
        key=$(echo "$line" | awk '{print $1}')
        rest=$(echo "$line" | awk '{$1=""; sub(/^[ \t]+/, ""); print}')
        rest=$(echo "$rest" | tr -d '\r\n')  # Remove any newlines or carriage returns

        varname="${key^^}"
        eval "${varname}=\"\$rest\""

    done < "$CONFIG_FILE"

    # Validate required variables
    validate_vntr_variables
}

validate_vntr_variables() {
    local required_general_vars=(
        "VNTR_BOUNDARY_SEQUENCE_RIGHT"
        "VNTR_BOUNDARY_SEQUENCE_LEFT"
        "VNTR_ASSEMBLY_SIZE"
        "VNTR_MIN_PRODUCT_SIZE"
        "VNTR_REPEAT_SIZE"
        "VNTR_MOTIFS"
        "VNTR_COLORS"
    )

    local hg38_vars=(
        "VNTR_COORDINATES_HG38_CHR"
        "VNTR_COORDINATES_HG38_START"
        "VNTR_COORDINATES_HG38_END"
    )

    local chm13_vars=(
        "VNTR_COORDINATES_CHM13_CHR"
        "VNTR_COORDINATES_CHM13_START"
        "VNTR_COORDINATES_CHM13_END"
    )

    local missing=0
    local missing_vars=()

    # --- Check general required vars ---
    for var in "${required_general_vars[@]}"; do
        if ! printenv "$var" >/dev/null && [[ -z "${!var+x}" ]]; then
            missing_vars+=("$var")
            missing=1
        fi
    done

    # --- Check coordinate groups (but don't unset yet) ---
    local hg38_missing=0
    for var in "${hg38_vars[@]}"; do
        if ! printenv "$var" >/dev/null && [[ -z "${!var+x}" ]]; then
            missing_vars+=("$var")
            hg38_missing=1
        fi
    done

    local chm13_missing=0
    for var in "${chm13_vars[@]}"; do
        if ! printenv "$var" >/dev/null && [[ -z "${!var+x}" ]]; then
            missing_vars+=("$var")
            chm13_missing=1
        fi
    done

    # --- Require at least one full coordinate group ---
    if [[ "$hg38_missing" -eq 1 && "$chm13_missing" -eq 1 ]]; then
        missing=1
    fi

    # --- Report and exit if anything is missing ---
    if [[ "$missing" -ne 0 ]]; then
        echo "❌ ERROR: The following required variable(s) are missing:"
        for var in "${missing_vars[@]}"; do
            echo "   - $var"
        done
        echo "❌ Exiting due to missing configuration for VNTR '$target_vntr'." >&2
        [[ "${BASH_SOURCE[0]}" != "${0}" ]] && exit 1 || exit 1
    fi

    # --- Now unset incomplete coordinate groups ---
    if [[ "$hg38_missing" -eq 1 ]]; then
        for var in "${hg38_vars[@]}"; do
            unset "$var"
        done
    fi

    if [[ "$chm13_missing" -eq 1 ]]; then
        for var in "${chm13_vars[@]}"; do
            unset "$var"
        done
    fi
}

commit_vntr_variables() {
    local vntr_name="$1"
    local config_file="$CONFIG_FILE"  # Assumes CONFIG_FILE is set in your environment

    if [[ -z "$vntr_name" ]]; then
        echo "❌ Please provide a VNTR name to save."
        exit 1
    fi

    local section_header="#${vntr_name^^}"  # Uppercase section header

    # Check if section header exists (case-insensitive)
    if grep -i -q "^#${vntr_name}\b" "$config_file"; then
        echo "Section $section_header already exists in $config_file. Not appending."
        return 0
    fi

    echo "Appending VNTR variables for $section_header to $config_file..."

    {
        echo ""
        echo "$section_header"
        # Print all VNTR_* variables that are set
        compgen -v | grep '^VNTR_' | while read -r var; do
            # Only print if variable is set and not empty
            if [[ -n "${!var}" ]]; then
                printf "%s %s\n" "$var" "${!var}"
            fi
        done
    } >> "$config_file"

    echo "Variables appended."
}

