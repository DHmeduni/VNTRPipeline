#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh

# --- Argument check ---
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <input_file> <threshold> <separation> <min_frequency> <min_loss> <max_loss>"
    exit 1
fi

input_file=$1
threshold=$2
separation=$3
min_frequency=$4
min_loss=$5
max_loss=$6

# Optional: if find_maxima.py is in the same directory
cd "$(dirname "$0")"

# Activate conda environment
conda activate python-env || { echo "Error: Could not activate conda environment 'python-env'"; exit 1; }

# --- Call the Python function and capture its output ---
maxima=$(python - <<END
import sys
from find_maxima import find_maxima_from_file

# Values passed from bash are injected using string formatting below
input_file = "${input_file}"
threshold = int("${threshold}")
separation = int("${separation}")
min_frequency = int("${min_frequency}")
min_loss = float("${min_loss}")
max_loss = float("${max_loss}")

try:
    maxima = find_maxima_from_file(input_file, threshold, separation, min_frequency, min_loss, max_loss)

    if len(maxima) == 2:
        print(f"{maxima[0]} {maxima[1]}")
    elif len(maxima) == 1:
        print(f"{maxima[0]}")
    else:
        print("Error: No maxima found")
except Exception as e:
    print(f"Error: {str(e)}")
END
)



# Check if maxima were found
if [[ $maxima == Error* ]]; then
    echo "$maxima"
    exit 1
fi

# Parse the maxima into variables
echo "$maxima"