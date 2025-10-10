#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh

# Check if correct number of arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input.bam> <range>"
    echo "Example: $0 input.bam 100"
    echo "This will filter reads for haplotype1 and for haplotype2"
    echo "based on 2 maxima larger than 2500 bp found in the sequencing data"
    echo "with a minimum seperation of 60bp !"
    echo ""
    echo "If two maxima are not found, then get_maxima throws an error"
    echo "and this script exits with exit(1). I've got to call this script"
    echo "and let it process the error, so that if the error comes whatshap"
    echo "used instead. Still working on this."
    exit 1
fi

input_bam=$1

LOGFILE="$(dirname "$input_bam")/script.log"

# Extract filename without extension
filename=$(basename "$input_bam" .bam)

script=$(readlink -f "${BASH_SOURCE[0]}")
script_path=$(dirname "$script")
cd "$script_path" || exit


threshold=$(( (VNTR_MIN_PRODUCT_SIZE / 10) * 10 ))
separation=$(( (VNTR_REPEAT_SIZE / 10) * 10 ))
range=$(( separation / 2 ))

read length1 length2 < <(./get_maxima.sh "$input_bam" "$threshold" "$separation" "${MIN_FREQUENCY:-20}")

echo "$LENGTH_1" >> "$LOGFILE"
echo "$LENGTH_2" >> "$LOGFILE"
#read -p "Press any key to continue..."

if [[ -z "$LENGTH_1" ]]; then
  if [[ $length1 == "Error:" ]]; then
    echo "$length1 $length2" >&2
    exit 1
  elif [[ -z $length2 ]]; then
    echo "$length1 $length2" >&2
    homozygote_marker=Y
  fi
elif [[ -n "$LENGTH_1" ]]; then
    length1="$LENGTH_1"
    length2="$LENGTH_2"
    unset LENGTH_1
    unset LENGTH_2
    echo "$length1 $length2" >&2
fi
#read -p "Press any key to continue..."

echo "$length1" >> "$LOGFILE"
echo "$length2" >> "$LOGFILE"


echo "Sample: $(basename "$input_bam") Length1: $length1, Length2: $length2" >> "$LOGFILE"



mkdir -m 777 "${input_bam%.*}_length_haplotypes"
cd "${input_bam%.*}_length_haplotypes" || exit

echo "$threshold $range $separation $MIN_FREQUENCY" > "${input_bam%.*}"_length_haplotypes/"${filename}"_length_histogram.txt

cd "${input_bam%.*}_length_haplotypes" || exit

conda activate python-env
log python "$script_path"/find_maxima.py -i "$input_bam" \
  -b 10 -o "${input_bam%.*}"_length_haplotypes/"${filename}"_length_histogram.png \
  -t "$threshold" -s "$separation" \
  --x_max 7000 $( [ -n "$MIN_FREQUENCY" ] && echo "-m $MIN_FREQUENCY" )
  
conda deactivate


# Calculate range boundaries
min1=$((length1 - range))
max1=$((length1 + range))
min2=$((length2 - range))
max2=$((length2 + range))

unset length1
unset length2

conda activate samtools-env

# Filter and create haplotype1 BAM file
log bash -c "samtools view -h \"$input_bam\" \
  | awk -v min=\"$min1\" -v max=\"$max1\" 'substr(\$0,1,1)==\"@\" || (length(\$10) >= min && length(\$10) <= max)' \
  | samtools view -b \
  | samtools sort -o \"${filename}_haplotype1.bam\""

log samtools index "${filename}_haplotype1.bam"

log samtools view "$input_bam" | awk -v min="$min1" -v max="$max1" 'substr($0,1,1)=="@" || (length($10) >= min && length($10) <= max)' | wc -l >> "${input_bam%.*}"_length_haplotypes/"${filename}"_length_histogram.txt

# Filter and create haplotype2 BAM file
if [[ "$homozygote_marker" != "Y" ]]; then
    log samtools view -h "$input_bam" | awk -v min="$min2" -v max="$max2" 'substr($0,1,1)=="@" || (length($10) >= min && length($10) <= max)' | samtools view -b | samtools sort -o "${filename}_haplotype2.bam"
    log samtools index "${filename}_haplotype2.bam"

    log samtools view "$input_bam" | awk -v min="$min2" -v max="$max2" 'substr($0,1,1)=="@" || (length($10) >= min && length($10) <= max)' | wc -l >> "${input_bam%.*}"_length_haplotypes/"${filename}"_length_histogram.txt
fi

echo "Filtering, sorting, and indexing complete." | tee -a "${input_bam%.*}"_length_haplotypes/"${filename}"_length_histogram.txt >> "$LOGFILE"
echo "Haplotype1: ${min1}-${max1} bp" | tee -a "${input_bam%.*}"_length_haplotypes/"${filename}"_length_histogram.txt >> "$LOGFILE"
if [[ "$homozygote_marker" != "Y" ]]; then
echo "Haplotype2: ${min2}-${max2} bp" | tee -a "${input_bam%.*}"_length_haplotypes/"${filename}"_length_histogram.txt >> "$LOGFILE"
fi

cd "$(dirname "$0")" || exit
