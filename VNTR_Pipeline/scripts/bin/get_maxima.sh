#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh

# Check if correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

input_file=$1
threshold=$2
separation=$3

cd "$(dirname "$0")"

conda activate python-env

# Call the Python function and capture its output
maxima=$(python - <<END
import sys
from find_maxima import find_maxima_from_file

try:
    input_file = "$input_file"
    threshold = int("$threshold")
    separation = int("$separation")
    
    maxima = find_maxima_from_file(input_file, threshold, separation)
    
    if len(maxima) == 2:
        print(f"{maxima[0]} {maxima[1]}")
    elif len(maxima) == 1:
        print(f"{maxima[0]}")
    else:
        print("Error: Could not find one or two maxima.")
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