#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh

trap '[[ -n "$output_base" ]] && chmod -R 777 "$output_base"' EXIT

script=$(readlink -f "${BASH_SOURCE[0]}")
script_path=$(dirname "$script")
script_name=$(basename "$script")

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
    echo "the MUC1 VNTR, other options can be"
    echo "added by the user (please see Github)."
    echo ""
    echo "Samples must be presented as either bam or fastq."
    echo "Output can be a subdirectory of input folder or own path."
    echo "Test Data from HG001 to HG004 PCR reactions is"
    echo "found in the folder /test_data and can be copied to"
    echo "an external drive if desired."
    echo "#################################################"
    echo "Syntax: script name [-c|h|i|o|p|r|v|] <image_name> "
    echo "options:"
    echo "-h, --help    Print this Help."
    echo "-c            Commit Variables to cfg file (VNTR_NAME)"
    echo "-i            Input folder path/Input file"
    echo "-o            Output folder path"
    echo "-r            Reference Name hg38|t2t"
    echo "-p            PCR or WGS [pcr|wgs]"
    echo "-v            VNTR (MUC1 or ACAN, more will be added)"
    echo "-q/--quiet (future)"
    echo "--version (future)"
    echo ""
    echo "#################################################"
    echo "#add for other VNTRs"
    echo "For other VNTRs/motifs, instead of -v, the following additional options have to be set when executing VNTR_pipeline."
    echo "These variables can be committed to the configuration file by adding the option -c {VNTR NAME} and must be added after the above options."
    echo "Note: any of the below options used with the –v {VNTR_NAME} option will override the values saved to the configuration file."
    echo ""
    echo "VNTR_BOUNDARY_SEQUENCE_LEFT   Sequence at beginning of VNTR (will be used to extract the VNTR, one mismatch allowed)"
    echo "VNTR_BOUNDARY_SEQUENCE_RIGHT  Sequence at end of VNTR (will be used to extract the VNTR, one mismatch allowed)"
    echo "VNTR_ASSEMBLY_SIZE            Approximate size of PCR Product (e.g.: 500, 1k, 3k, 10k,), used for PCR data, minimum should be 500"
    echo "VNTR_MIN_PRODUCT_SIZE         Minimum PCR Product Size expected, necessary for filtering"
    echo "VNTR_REPEAT_SIZE              The Number of bases found in one repeat unit of the VNTR"
    echo "VNTR_COORDINATES_T2T          T2T-chm13v2.0 Reference coordinates of the VNTR for filtering in the format [chr]:[start]-[end]"
    echo "VNTR_COORDINATES_HG38         GRCh38.p14 Reference coordinates of the VNTR for filtering in the format [chr]:[start]-[end]"
    echo "VNTR_MOTIFS                   Path to VNTR Motif file, containing Sequencing, Alphanumeric Designation and Color code, tab separated .txt file, example files are available on github"
    echo ""
    echo "#Following options may be added:"
    echo ""
    echo "TMP_DELETE=N                  Retain the temporary files for troubleshooting"
    echo "ALL_FIGURES=Y                 Only produces figures, based on already assembled and trimmed sequences of *best_hit.fasta files, recursively to a depth of one subfolder from a folder that is given as input"
    echo "ALL_FIGURES=A                 Additionally to the analysis of each sample a combined figure is produced"
    echo "NON_CODING=Y                  Analyse VNTRs in non-coding regions, LoF prediction is skipped"
    echo "VNTR_ALL=Y                    Analyse all VNTRs (assemblies or polished sequences the top three most common sequences per haplotype) found by the workflow (pseudogenes/duplications)"
    echo "VNTR_PACBIO=Y                 Process PacBio WGS Data"
    echo "WHATSHAP_FORCE=Y              Force Whatshap haplotyping"
    echo "CONFIG_FILE                   Define another path for the configuration file (CONFIG_FILE=/path/to/file)"
    echo "LENGTH_1                      Define the shorter of two lengths for length based haplotyping (e.g. LENGTH_1=2500), default is automatic"
    echo "LENGTH_2                      Define the longer of two lengths for length based haplotyping (e.g LENGTH_2=3000), default is automatic"
    echo "MIN_FREQUENCY                 Temporarily lower the minimum frequency threshold for length based haplotyping (default MIN_FREQUENCY=20)"
    echo ""
    echo ""
    echo "#################################################"
    echo ""
    echo "Example Usage:"
    echo ""
    echo "           {"$script_name"} -o /path/to/output -i /path/to/bam/or/fastq -p pcr -r hg38 -v MUC1"
    echo "$(dirname $0)"
}

ALL_FIGURES=N
DELETE_TMP=Y
variant_caller_options=()
output_folder="tmp_$(date +%Y%m%d_%H%M%S)"

# Input config file
if [[ -z $script_path ]]; then
    script_path="$(dirname $0)"
fi

LIB_PATH="$(dirname $script_path)"/lib
export LIB_PATH
CONFIG_FILE="$LIB_PATH""/vntr_variables.cfg"

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
            if [[ "${OPTARG,,}" == "hg38" ]]; then
                reference="hg38_p14"
            elif [[ "${OPTARG,,}" == "t2t" ]]; then
                reference="chm13v2.0"
            else
                echo "Reference not known, please check syntax"
                exit 1
            fi
            variant_caller_options+=("-r" "/references/"$reference".fa") ;;
        s) # sample_sheet (absolute path)
            sample_sheet="$OPTARG" 
            variant_caller_options+=("-s" "/$input_folder/"$OPTARG"") ;;
        q) # quiet option
            quiet_opt="-q" ;;
        p) # pcr or wgs [pcr 0 | wgs 1]
            if [[ "${OPTARG,,}" == "pcr" ]]; then
                pcr_or_wgs="0"
            elif [[ "${OPTARG,,}" == "wgs" ]]; then
                pcr_or_wgs="1"
            else
                echo "PCR or WGS, please check syntax"
                exit 1
            fi
            variant_caller_options+=("-p" "$pcr_or_wgs");;
        v) # which vntr to look at
            vntr="${OPTARG^^}"  # Make uppercase
            variant_caller_options+=("-v" "$vntr")
            # Check $CONFIG_FILE for line containing #$vntr
            if ! grep -q "#$vntr" "$CONFIG_FILE"; then
                echo "Error: #$vntr not found in $CONFIG_FILE" >&2
                exit 1
            fi
            ;;
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


if [[ -n $vntr ]] && [[ -z $commit_variables ]]; then
    source "$script_path"/config_parse.sh
    parse_vntr_config "$vntr" || {
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
    echo "VNTR variables missing!"
    exit 1
    }
fi

for var in $(compgen -v | grep '^[A-Z0-9_]\+$'); do
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

if [[ "$output_folder_final" == */* ]]; then
    # output_folder_final is a path
    output_base="$output_folder_final"
else
    # output_folder_final is a folder name
    output_base="$input_folder/$output_folder_final"
fi

create_all_figures() {
    declare "ALL_FIGURES_CREATE=Y"
    export ALL_FIGURES_CREATE
    mkdir -m 777 -p "$output_base/output_TRViz/"
    if [[ "$VNTR_ALL" == "Y" ]]; then
        #find "$1" -name "*result_comi.fasta" -type f -exec cat {} + > "$output_base/output_TRViz/all_figures_result_combined.fasta" -path "$output_folder" -prune -o
        find "$1" -maxdepth 2 -name "*result.fasta" -type f -exec awk '
            /^>/ {
            n++
            if (n > 3) exit
            print $1"_"$2
            next
            }
            { print }
            ' {} \; > "$output_base"/output_TRViz/all_figures_result_combined.fasta
    else    
        find "$1" -maxdepth 2 -name "*best_hit.fasta" -type f -exec cat {} + > "$output_base/output_TRViz/all_figures_best_hit_combined.fasta"
        #-path "$output_folder" -prune -o 
    fi
    "$script_path"/whole_pipeline.sh -d $output_base -m "$mode" -l "${vntr^^}"
    if [[ -d "$output_base/output_TRViz" ]]; then
        mkdir -p -m 777 "$output_base"/all_figures
        mv "$output_base"/output_TRViz/best_trviz_fig.png "$output_base"/all_figures
        if [[ "$ALL_FIGURES" != "Y"  ]]; then
            rm -r "$output_base"/output_TRViz
            rm -r "$output_base"/tmp_motif.txt
            rm -r "$output_base"/motifs_evaluation.xlsx
        else
            mv "$output_base"/output_TRViz/new_and_lof_seqs.xlsx "$output_base"/all_figures
            rm -r "$output_base"/output_TRViz
            rm -r "$output_base"/tmp_motif.txt
            mv "$output_base"/motifs_evaluation.xlsx "$output_base"/all_figures
        fi
    fi
}

if [[ "$ALL_FIGURES" == "Y" ]]; then
    create_all_figures "$input_folder"
    rm -r $output_folder
    chmod -R 777 $(dirname "$output_folder")
    exit 0
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

if [[ "$mode" == "pcr" ]] || [[ "$ASSEMBLY" == "Y" ]]; then
    "$script_path"/variant_caller.sh -u -i "$input_folder" -o "$output_folder" -r /references/"$reference".fa -p "$pcr_or_wgs"
fi



for dir in "$output_folder"/*haplotypes/; do
    # Compose the final output directory name
    output_dir="$output_base/$(basename "${dir%_*_*}")_results"
	mkdir -m 777 -p "$output_dir"
    find "$dir" \( -name "*.xlsx" -o -name "*fig.png" -o -name "*histogram.png" -o -name "*histogram.txt" \) -type f -exec cp {} "$output_dir"  \;
    find "$dir" -name "*###*" -type f -exec cp {} "$output_base"  \;
    find "$dir" -name "*best_hit.fasta" -type f -exec cp {} "$output_dir" \;    
    find "$dir" -name "*result.fasta" -type f -exec cp {} "$output_dir" \;
    if [[ "$mode" == "pcr" ]] || [[ "$ASSEMBLY" == "Y" ]]; then
        mkdir -m 777 -p "$output_dir/assembly_mapping"
        find "$dir" -name "best_hit_combined*" -type f -exec cp {} "$output_dir/assembly_mapping"  \;
        find "$dir" -name "combined_contigs*" -type f -exec cp {} "$output_dir/assembly_mapping"  \;
        find "$dir" -name "result_combined*" -type f -exec cp {} "$output_dir/assembly_mapping"  \;
        find "$dir" -type d -name "assembly_mapping" -exec cp -r {} "$output_dir"  \;
    fi
done

find "$output_folder" \( -name "script.log"  \) -type f -exec cp {} "$output_base"  \;

conda activate python-env
python "$script_path"/merge_xlsx.py "$output_folder" new_and_lof_seqs.xlsx "$output_base"/All_Samples_new_and_lof_seqs.xlsx
python "$script_path"/merge_xlsx.py "$output_folder" seqs_distribution.xlsx "$output_base"/All_Samples_seqs_distribution.xlsx
conda deactivate

if [[ "$ALL_FIGURES" != "N" ]]; then
    create_all_figures "$output_base"
fi

if [[ "$DELETE_TMP" == "Y" ]]; then
    rm -r "$output_folder"
    find "$output_base" -type f -name tmp* -exec rm {} \;
else
    find "$output_base" -type f -name tmp* -exec rm {} \;
fi


chmod -R 777 $(dirname "$output_base")