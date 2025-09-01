#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh

############################################################
# Help                                                     #
############################################################
Help() { 
    # Display Help
    echo "#################################################"
    echo "VNTR Pipeline Med Uni Wien Institute Medical Genetics"
    echo "This program calls loss-of-function mutations"
    echo "in VNTR regions."
    echo "#################################################"
    echo "The pipeline is currently only built for"
    echo "ACAN and MUC1 VNTRs, other options will be"
    echo "added in future versions."
    echo ""
    echo "Samples must be presented as either bam or fastq."
    echo "Output folder is a subfolder of input folder."
    echo "Test Data from HG001 to HG004 PCR reactions is"
    echo "found in the folder /test_data and should be copied to"
    echo "an external drive if the docker container is not"
    echo "executed in an interactive mode."
    echo "#################################################"
    echo "Syntax: script name [-a|b|c|d|f|h|i|l|m|n|p|r|s|t|u|v|w|V] <image_name> "
    echo "options:"
    echo "-h, --help    Print this Help."
    echo "-c            Commit Variable to cfg file (VNTR_NAME)"
    echo "-i            Input folder path/Input file"
    echo "-o            Output folder path"
    echo "-r            Reference Name hg38_p14|chm13v2.0"
    echo "-p            PCR or WGS [0|1]"
    echo "-v            VNTR (MUC1 or ACAN, more will be added)"
    echo "-q/--quiet (future)"
    echo "--version (future)"
    echo ""
    echo ""
    echo "Example Usage:"
    echo ""
    echo "           {script_name} -p 0 -o output -i /path/to/bam/or/fastq -r reference"
    echo "$(dirname $0)"
}

DELETE_TMP=Y
variant_caller_options=()
output_folder="tmp_$(date +%Y%m%d_%H%M%S)"

# Get the options
while getopts "c:hi:m:o:p:r:qs:v:V-:" option; do
    case $option in
        h) # display Help
            Help
            exit ;;
        c) # commit VNTR variables
            commit_variables="$OPTARG" ;;
        i) # input folder/file (bams, fastqs)
            input_folder="$(readlink -f "$OPTARG")"
            variant_caller_options+=("-i" "$(readlink -f "$OPTARG")") ;;
        m) # demultiplex_barcodes
            variant_caller_options+=("-m" "$OPTARG") ;;
        o) # output folder (absolute)
            output_folder_final="$OPTARG" ;;
        r) # reference name (absolute)
            reference="$OPTARG" 
            variant_caller_options+=("-r" "/references/"$OPTARG".fa") ;;
        s) # sample_sheet (absolute path)
            sample_sheet="$OPTARG" 
            variant_caller_options+=("-s" "/$input_folder/"$OPTARG"");;
        q) # quiet option
            quiet_opt="-q" ;;
        p) # pcr or wgs [pcr 0 | wgs 1]
            pcr_or_wgs="$OPTARG"
            variant_caller_options+=("-p" "$OPTARG");;
        v) # which vntr to look at
            vntr="$OPTARG"
            variant_caller_options+=("-v" "${vntr^^}")
            case "${vntr^^}" in
                MUC1|ACAN) ;;
                *) echo "Invalid Option: $vntr Please see Help" >&2; exit 1 ;;
            esac ;;
        V) # VNTR pipeline
            whole_pipeline=Y ;;
        -) # Parse long option (long option branch)
             case "$(echo $OPTARG | awk -F'=' '{print $1}')" in
                help) #-c long-option
                    Help 
                    exit ;;
                quiet)
                    quiet_opt="-q"
                    ;;
                version)
                    echo "$script_name"
                    exit 1 ;;
                *)
                    echo "Invalid long option: --${OPTARG}" >&2
                    exit 1 ;;
            esac ;;
        \?) # Invalid option
            echo "Error: Invalid option -$OPTARG" >&1
            Help
            exit ;;
    esac
done

script=$(readlink -f "${BASH_SOURCE[0]}")
script_path=$(dirname "$script")

# Input config file
if [[ -z $script_path ]]; then
    script_path="$(dirname $0)"
fi

LIB_PATH="$(dirname $script_path)"/lib
export LIB_PATH
CONFIG_FILE="$LIB_PATH""/vntr_variables.cfg"

if [[ $vntr ]] && [[ -z $commit_variables ]]; then
    source "$script_path"/config_parse.sh
    parse_vntr_config "$vntr" "MUC1" || {
    echo "Failed to parse VNTR config!"
    exit 1
    }
    validate_vntr_variables
elif [[ -n $commit_variables ]]; then
    source "$script_path"/config_parse.sh
    validate_vntr_variables
    commit_vntr_variables "$commit_variables"
else
    source "$script_path"/config_parse.sh
    validate_vntr_variables || {
    echo "VNTR missing!"
    exit 1
    }
fi

# Shift away processed options
shift $((OPTIND - 1))

# --- Parse key=value arguments safely ---
for arg in "$@"; do
  if [[ "$arg" =~ ^[A-Za-z_][A-Za-z0-9_]*=.+$ ]]; then
    key=${arg%%=*}
    val=${arg#*=}   
    # Remove surrounding quotes if present
    if [[ "$val" =~ ^\".*\"$ || "$val" =~ ^\'.*\'$ ]]; then
      val="${val:1:-1}"
    fi
    declare "$key=$val"
  else
    echo "⚠️ Warning: Unrecognized argument '$arg'" >&2
  fi
done

for var in $(compgen -v | grep '^[A-Z_]\+$'); do
    export "$var"
done
for var in $(compgen -v | grep '^VNTR_'); do
  export "$var"
done


# output is path or folder for tmp data
if [[ "$output_folder_final" == */* ]]; then
    # output_folder_final is a path
    output_folder="$output_folder_final/$output_folder/"
else
    # output_folder_final is a folder name
    output_folder="$input_folder/$output_folder/"
fi
echo "$output_folder"
if [ -d "$output_folder" ]; then
    echo "❌ ERROR: Directory '$output_folder' already exists."
    exit 1
else
    mkdir -m 777 -p "$output_folder"
    echo "Directory '$output_folder' (and parents) created."
fi


declare -A options
options["1"]="-awc"
options["0"]="-awlc"


options_list="${options["$pcr_or_wgs"]}"

if [[ $WHATSHAP_FORCE == "Y" ]]; then
    options_list="${options[1]}"
fi

variant_caller_options+=("-o" "$(readlink -f "$output_folder")")
variant_caller_options+=("$options_list")

echo "${variant_caller_options[@]}"

 "$script_path"/variant_caller.sh \
"${variant_caller_options[@]}"


if [[ "$pcr_or_wgs" == "0" ]]; then
    mode="pcr"
elif [[ "$pcr_or_wgs" == "1" ]]; then
    mode="wgs"
else
    echo "somethings not right here pick a pcr or wgs"
    exit 1 
fi



for dir in  "$output_folder"/*haplotypes/; do
	"$script_path"/whole_pipeline.sh -d $dir -m "$mode" -l "${vntr^^}"
done

"$script_path"/variant_caller.sh -u -i "$input_folder" -o "$output_folder" -r /references/"$reference".fa



for dir in "$output_folder"/*haplotypes/; do
    if [[ "$output_folder_final" == */* ]]; then
        # output_folder_final is a path
        output_base="$output_folder_final"
    else
        # output_folder_final is a folder name
        output_base="$input_folder/$output_folder_final"
    fi
    # Compose the final output directory name
    output_dir="$output_base/$(basename "${dir%_*_*}")_results"
	mkdir -m 777 -p "$output_dir"
    find "$dir" -type d -name "assembly_mapping" -exec cp -r {} "$output_dir"  \;
    find "$dir" \( -name "*.xlsx" -o -name "*.png" -o -name "*histogram.txt" \) -type f -exec cp {} "$output_dir"  \;
    find "$dir" -name "best_hit_combined*" -type f -exec cp {} "$output_dir/assembly_mapping"  \;
    find "$dir" -name "*###*" -type f -exec cp {} "$output_base"  \;
done

find "$output_folder" \( -name "script.log" -o -name "version.log" \) -type f -exec cp {} "$output_base"  \;

conda activate python-env
python "$script_path"/merge_xlsx.py "$output_folder" LOF_check.xlsx "$output_base"/LOF_results.xlsx
python "$script_path"/merge_xlsx.py "$output_folder" new_and_lof_seqs.xlsx "$output_base"/New_Repeats.xlsx
conda deactivate

if [[ "$DELETE_TMP" == "Y" ]]; then
    rm -r "$output_folder"
fi

chmod -R 777 "$output_folder"