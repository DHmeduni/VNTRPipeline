#!/bin/bash

# --- Environment Setup ---
source /opt/conda/etc/profile.d/conda.sh

# --- Global Variables ---
script=$(readlink -f "${BASH_SOURCE[0]}")
script_path=$(dirname "$script")
if [[ -z $script_path ]]; then
    script_path="$(dirname $0)"
fi
lib_path="$(dirname $script_path)"/lib
source "$script_path"/config_parse.sh

# --- Trap for Permissions ---
trap '[[ -n "$output_folder" ]] && chmod -R 777 "$output_folder"' EXIT


############################################################
# Help                                                     #
############################################################
Help() { 
    # Display Help
    echo "#################################################"
    echo "This script will basecall (future), map, variant call and"
    echo "phase Nanopore sequencing data."
    echo "#################################################"
    echo "If this program is to be used with VNTR"
    echo "Pipeline, please use the T2T Reference and MUC1"
    echo "or ACAN options for the coordinates."
    echo "coordinates and reference will be automatically found."
    echo "Please use the image vntr_pipeline, latter versions"
    echo "will include basecalling and other functions."
    echo "#################################################"
    echo "Syntax: script name [-a|b|c|d|f|h|i|l|m|n|p|r|s|t|u|v|w|V] <image_name> "
    echo "options:"
    echo "-h, --help    Print this Help."
    echo "-a            Use Minimap2 for alignment"
    echo "-b            Dorado Basecaller (Meth, Sup, Hac, Fast) Default:Meth (future version)"
    echo "-c            Use clair3 variant caller"
    echo "-d            Script Folder (unneeded at moment)"
    echo "-f            Use bcftools (future option)"
    echo "-i            Input folder (Absolute, depends on what is being processed)"
    echo "-m            Demutliplex (1|0) (determines if demultiplex or for loop)"
    echo "-o            Output folder (absolute)"
    echo "-r            Reference Path. (Absolute)"
    echo "-w            Phase with Whatshap"
    echo "-p            PCR or WGS [0|1]"
    echo "-s            SV calling with Sniffles (future option)"
    echo "-v            VNTR (MUC1 or ACAN, more can be added)"
    echo "-V            VNTR pipeline (at this level unneeded)"
    echo "-q/--quiet"
    echo "--version"
    echo ""
    echo ""
    echo "Example Usage:"
    echo ""
    echo "           {script_name} -acmws -o /path/to/output -i /path/to/bam/or/fastq -r /path/to/reference"
}

parse_options() {
    while getopts "ab::cd:fhi:lm:o:p:r:qs:uwv:V-:" option; do
        case $option in
            a) # alignment Minimap2
                minimap=Y ;;
            b) # basecaller
                basecaller="${OPTARG:-Meth}"
                case "${basecaller,,}" in
                    meth|sup|hac|fast) ;;
                    *) echo "Invalid Option: $basecaller Please see Help" >&2; exit 1 ;;
                esac ;;
            c) # variant calling clair3
                clair=Y ;;
            d) # script folder mount
                script_folder="$(readlink -f "$OPTARG")" ;;
            f) # variant calling bcftools
                bcftools=Y ;;
            h) #help flag
                Help
                exit ;;
            i) # input folder/file (bams, fastqs)
                input_folder="$(readlink -f "$OPTARG")" ;;
            l) # haplotype_by_length
                length_haplotype=Y ;;
            m) # demultiplex_barcodes
                demultiplex="$OPTARG" ;;
            o) # output folder (absolute)
                output_folder="$(readlink -f "$OPTARG")" ;;
            p) # pcr of wgs [pcr 0| wgs 1]
                pcr_or_wgs="$OPTARG" ;;
            r) # reference path (absolute)
                reference_path="$OPTARG" ;;
            q) # quiet option
                quiet_opt="-q" ;;
            s) # sample_sheet (absolute path)
                sample_sheet="$OPTARG" ;;
            u) # mapping of the sequences to the assemblies
                assembly_mapping=Y ;;
            w) # phase with whatshap
                whatshap=Y
                whatshap_now=Y;;
            v) # which vntr to look at
                vntr="${OPTARG^^}"  # Make uppercase
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
}

##################################################
###       Functions
##################################################

log() {
    # Log only the command itself with expanded variables (but not its output)
    echo "\$ $*" >> "$LOGFILE"
    # Execute the command, redirect stdout and stderr to the console, but NOT log to the file
    #"$@" 2>&1 > /proc/self/fd/1 2>/proc/self/fd/2
    "$@"
}

version_log() {
    # Log the version of the program used
    local program="$1"  
    shift          

    local version_output
    version_output=$("$program" "$@" 2>&1 | head -n 1)  # Capture both stdout and stderr

    echo "$program $*:" >> "$VERSION_LOGFILE"
    echo "$version_output" >> "$VERSION_LOGFILE"

}


vntr_coordinates() {

    if echo "$reference_path" | grep -qE "(hg38|chm13)"; then
        local ref=$(echo "$reference_path" | grep -oE "(hg38|chm13)" | head -n 1)
        [[ $ref == "chm13" ]] && ref="t2t"
    else
        echo "Neither hg38 nor t2t found in the path." 
        Help
        exit 1
    fi
    
    eval "VNTR_COORDINATES_${ref^^}_CHR=\${VNTR_COORDINATES_${ref^^}%%:*}"
    eval "temp=\"\${VNTR_COORDINATES_${ref^^}#*:}\""
    eval "VNTR_COORDINATES_${ref^^}_START=\${temp%%-*}"
    eval "VNTR_COORDINATES_${ref^^}_END=\${temp#*-}"

    local chr_name="VNTR_COORDINATES_"${ref^^}"_CHR"
    local start_name="VNTR_COORDINATES_"${ref^^}"_START"
    local end_name="VNTR_COORDINATES_"${ref^^}"_END"

    if [[ -z $end_name ]]; then
        echo "Error: no variables coordinates"
        exit 1
    fi

    local chr="${!chr_name}"
    local start="${!start_name}"
    local end="${!end_name}"

     if [[ -z $end ]]; then
        echo "Error: no variables coordinates values"
        exit 1
    fi
    
    if [[ "$pcr_or_wgs" == "1" ]]; then
        local region="${chr}:$(( start > 25000 ? start - 25000 : 0 ))-$(( end + 25000 ))"
    elif [[ "$pcr_or_wgs" == "0" ]]; then
        local region="${chr}:${start}-${end}"
    elif [[ $assembly_mapping != "Y" ]]; then
        echo "Error: Missing required parameter pcr_or_wgs."
        Help
        exit 1
    fi

    echo "$region"
}


process_bam_to_fastq() {
    conda activate samtools-env
    for bam_file in "$1"/*.bam; do
        if [[ -f "$bam_file" ]]; then
            fastq_file="$2"/"$(basename "${bam_file%.bam}")".fastq
            log samtools fastq -T '*' "$bam_file" > "$fastq_file"
            #log samtools fastq "$bam_file" > "$fastq_file"
        fi
    done
    conda deactivate
}

run_minimap2() {
    conda activate minimap-env
    version_log minimap2 --version | head -n 1
    local ref_path="${2:-$reference_path}"  # Use $2 if provided, otherwise default to reference_path
    
    log minimap2 -Yax lr:hq -t 32 -y --MD -d "${ref_path%.bam}.mmi" "$ref_path"

    for fastq_file in "$1"/*.fastq; do
        if [[ -f "$fastq_file" ]]; then
            sam_file="${fastq_file%.fastq}.sam"
            local sam_path="${3:-$sam_file}"
            log minimap2 -Yax lr:hq -t 32 -y --MD \
                -o "$sam_path" \
                "${ref_path%.bam}.mmi" \
                "$fastq_file"
        fi
    done
    conda deactivate
}

run_demultiplexing() {
    conda activate python-env
    version_log python --version | head -n 1
    log python "$(dirname "$0")/demux_split_barcodes.py" -i "$1" -o "$output_folder" -b "$sample_sheet" -p 16 -rs
    conda deactivate
    for i in "$output_folder"/*.fastq; do
        if [ $(stat -c %s "$i") -lt 500000 ]; then rm "$i"; fi
    done
}


process_sam_to_bam() {
    local input_sam="$1"
    local output_prefix="${input_sam%.*}"

    conda activate samtools-env

    version_log samtools --version | head -n 1
    
    log samtools view -bS -@ 32 "$input_sam" | \
    log samtools sort -@ 32 -o "${output_prefix}.tmp.bam" -
    log samtools index -@ 32 "${output_prefix}.tmp.bam"

    local coordinates=$(vntr_coordinates)

    if [[ -n "$coordinates" ]] && [[ "$assembly_mapping" != "Y" ]]; then
        log samtools view -b -@ 32 "${output_prefix}.tmp.bam" "$coordinates" | \
        log samtools sort -@ 32 -o "${output_prefix}.bam" -
        log samtools index -@ 32 "${output_prefix}.bam"
    else
        log mv "${output_prefix}.tmp.bam" "${output_prefix}.bam"
        log mv "${output_prefix}.tmp.bam.bai" "${output_prefix}.bam.bai"
    fi
    
    if [[ -f "${output_prefix}.tmp.bam" ]]; then
        rm "${output_prefix}.tmp.bam" "${output_prefix}.tmp.bam.bai"
    fi
    conda deactivate
}

extract_lengths() {
    export -f log
    #local vntr="${2:-muc1}"
    local script_path="$(dirname "$0")/extract_lengths.sh"

    # Log the command and capture the data and set exid code
    #output="$(log "$script_path" "$1" 2>&1)"
    output="$("$script_path" "$1" 2>&1)"
    exit_code=$?
    echo "Extract lengths condition : "$exit_code"" >> "$LOGFILE"

    # Log the output
    echo "$output" >> "$LOGFILE"
    last_line=$(echo "$output" | tail -n 1)
    echo "Here is the last_line: $last_line" >> "$LOGFILE"
    echo "Here is the word count of last_line $(wc -w <<< "$last_line")" >> "$LOGFILE"
    first_line=$(echo "$output" | head -n 1)
    echo "Here is the last_line: $first_line" >> "$LOGFILE"
    echo "Here is the word count of last_line $(wc -w <<< "$first_line")" >> "$LOGFILE"


    if [[ $(echo "$first_line" | wc -w) -lt 2 ]]; then
        echo "Only one or no maxima" >> "$LOGFILE"
        whatshap_now=Y
        return 1       
    fi
}



run_clair3() {
    local input_bam="$1"
    local output_dir="${input_bam%.*}"
    local ref_path="${2:-$reference_path}"

    conda activate clair3-env

    version_log /opt/bin/run_clair3.sh --version | head -n 1

    if [[ "$pcr_or_wgs" == "1" ]] && [[ "$VNTR_PACBIO" != "Y" ]]; then 
        log /opt/bin/run_clair3.sh \
        --bam_fn="$input_bam" \
        --ref_fn="$ref_path" \
        --threads=32 \
        --platform="ont" \
        --model_path="/opt/rerio/clair3_models/r1041_e82_400bps_sup_v500" \
        --output="$output_dir"_clair3
    elif [[ "$VNTR_PACBIO" == "Y" ]]; then
        log /opt/bin/run_clair3.sh \
        --bam_fn="$input_bam" \
        --ref_fn="$ref_path" \
        --threads=32 \
        --platform="hifi" \
        --model_path="/opt/models/hifi_revio" \
        --output="$output_dir"_clair3
    else
        log /opt/bin/run_clair3.sh \
        --bam_fn="$input_bam" \
        --ref_fn="$ref_path" \
        --threads=32 \
        --platform="ont" \
        --model_path="/opt/rerio/clair3_models/r1041_e82_400bps_sup_v500" \
        --output="$output_dir"_clair3 \
        --var_pct_full=1 \
        --ref_pct_full=1 \
        --enable_variant_calling_at_sequence_head_and_tail \
        --var_pct_phasing=1
    fi

    conda deactivate
}

run_bcftools() {

    local input_bam="$1"
    local ref_path="$2"
    local output_vcf="${input_bam%.*}"_bcf.vcf

    conda activate bcftools-env
    version_log bcftools --version | head -n 1
    log bcftools mpileup -Ov -o "$output_vcf".tmp -f "$2" "$1" 
    log bcftools call -mv -Ov -o "$output_vcf" "$output_vcf".tmp
    rm "$output_vcf".tmp
    conda deactivate
}

run_sniffles() {
    conda activate sniffles-env
    local input_bam="$1"
    local output_vcf="${input_bam%.*}"_sv.vcf
    version_log sniffles --version | head -n 1
    log sniffles --input "$input_bam" --vcf "$output_vcf"

    conda deactivate

}

run_whatshap() {
    local input_bam="$1"
    local output_dir="${input_bam%.*}_whatshap_haplotypes"
    local base_name="$(basename "${input_bam%.*}")"
    
    mkdir -p "$output_dir"
    
    conda activate whatshap-env
    version_log whatshap --version | head -n 1
    log whatshap phase -o "$output_dir/$base_name.bam.phased.vcf" \
        -r "$reference_path" --indels --ignore-read-groups \
        "${input_bam%.*}_clair3/merge_output.vcf.gz" \
        "$input_bam"
    conda deactivate

    conda activate bcftools-env
    version_log bcftools --version | head -n 1
    log bcftools sort "$output_dir/$base_name.bam.phased.vcf" -Oz -o "$output_dir/$base_name.bam.phased.vcf.gz"
    log bcftools index "$output_dir/$base_name.bam.phased.vcf.gz"
    conda deactivate

    conda activate whatshap-env
    version_log whatshap --version | head -n 1
    log whatshap haplotag -o "$output_dir/$base_name.bam.phased.haplotag.bam" --ignore-read-groups \
        --reference "$reference_path" --output-haplotag-list "$output_dir/$base_name.bam.phased.haplotag.tsv" \
        "$output_dir/$base_name.bam.phased.vcf.gz" "$input_bam"
    log whatshap split --output-h1 "$output_dir/$base_name.bam.phased.haplotype1.bam" \
        --output-h2 "$output_dir/$base_name.bam.phased.haplotype2.bam" \
        --discard-unknown-reads "$output_dir/$base_name.bam.phased.haplotag.bam" \
        "$output_dir/$base_name.bam.phased.haplotag.tsv"
    conda deactivate

    conda activate samtools-env
    version_log samtools --version | head -n 1
    log samtools index "$output_dir/$base_name.bam.phased.haplotype1.bam" 
    log samtools index "$output_dir/$base_name.bam.phased.haplotype2.bam" 
    log samtools index "$output_dir/$base_name.bam.phased.haplotag.bam"
    conda deactivate
}

################################################
####     MAIN
################################################

main() {
    parse_options "$@"
    VNTR_PACBIO=N

    # --- Parameter checks ---
    if [[ "$assembly_mapping" != "Y" && "$demultiplex" == "1" ]]; then
        for v in input_folder output_folder reference_path vntr pcr_or_wgs sample_sheet; do
            if [[ -z "${!v}" ]]; then
                echo "Error: Missing required parameter: $v"
                Help
                exit 1
            fi
        done
    fi

    LOGFILE="$output_folder/script.log"
    VERSION_LOGFILE="$output_folder/version.log"

    # --- Helper: process all files matching a glob with a callback ---
    process_files() {
        local pattern="$1"
        local callback="$2"
        shift 2
        shopt -s nullglob
        local files=($pattern)
        shopt -u nullglob
        for f in "${files[@]}"; do
            [[ -f "$f" ]] && "$callback" "$f" "$@"
        done
    }

    # --- Minimap2 block ---
    if [[ $minimap == "Y" ]]; then
        if [[ ! -v demultiplex || $demultiplex == "0" ]]; then
            if compgen -G "$input_folder/*.bam" > /dev/null; then
                process_bam_to_fastq "$input_folder" "$output_folder"
                run_minimap2 "$output_folder"
            elif compgen -G "$input_folder/*.fastq" > /dev/null; then
                cp "$input_folder"/*.fastq "$output_folder"
                run_minimap2 "$output_folder"
            elif compgen -G "$input_folder/*.fastq.gz" > /dev/null; then
                for file in "$input_folder"/*.fastq.gz; do
                    gzip -d -c "$file" > "$output_folder/$(basename "${file%.gz}")"
                done
                run_minimap2 "$output_folder"
            elif [[ -f "$input_folder" && "$input_folder" == *.fastq.gz ]]; then
                zcat "$input_folder" > "$output_folder/$(basename "${input_folder%.fastq.gz}").fastq"
                run_minimap2 "$output_folder"
            elif [[ -f "$input_folder" && "$input_folder" == *.fastq ]]; then
                cp "$input_folder" "$output_folder"
                run_minimap2 "$output_folder"
            elif [[ -f "$input_folder" && "$input_folder" == *.bam ]]; then
                cp "$input_folder" "$output_folder"
                process_bam_to_fastq "$output_folder" "$output_folder"
                rm -f "$output_folder/$(basename "$input_folder")"
                run_minimap2 "$output_folder"
            fi
        elif [[ $demultiplex == "1" ]]; then
            if compgen -G "$input_folder/*.fastq" > /dev/null; then
                run_demultiplexing "$input_folder"
            elif compgen -G "$input_folder/*.bam" > /dev/null; then
                process_bam_to_fastq "$input_folder" "$input_folder"
                run_demultiplexing "$input_folder/$(basename "${input_folder%.bam}").fastq"
            fi
            run_minimap2 "$output_folder"
        fi
    fi

    # --- Process SAM to BAM ---
    process_files "$output_folder/*.sam" process_sam_to_bam

    # --- Length Haplotype Extraction ---
    declare -A result_array
    if [[ "$length_haplotype" == "Y" && "$pcr_or_wgs" == "0" ]]; then
        whatshap_now=N
        for bam in "$output_folder"/*.bam; do
            [[ -f "$bam" ]] || continue
            extract_lengths "$bam"
            exit_code=$?
            if [[ "$exit_code" -eq 1 ]]; then
                echo "Length QC function detected an issue: $whatshap_now"
                result_array["$bam"]=1
            else
                result_array["$bam"]=0
                echo "Here is something happening: $bam"
            fi
        done
    fi

    # --- Clair3 Variant Calling ---
    run_clair3_for_bams() {
        for bam in "$output_folder"/*.bam; do
            [[ -f "$bam" ]] || continue
            "$@" "$bam"
        done
    }
   
    if [[ "$clair" == "Y" ]]; then
        if [[ "$whatshap_now" == "Y" && ( "$length_haplotype" == "Y" || "$pcr_or_wgs" == "0" ) && "$benchmarking" != "1" ]]; then
            for bam in "${!result_array[@]}"; do
                if [[ ${result_array["$bam"]} == "1" && -f "$bam" ]]; then
                    run_clair3 "$bam"
                fi
            done
        elif [[ ( "$whatshap_now" == "Y" && "$pcr_or_wgs" == "1" ) || "$benchmarking" == "1" ]]; then
            for bam in "$output_folder"/*.bam; do
                [[ -f "$bam" ]] || continue
                run_clair3 "$bam"
            done        
        fi
    fi

    # --- WhatsHap Phasing ---
    if [[ "$whatshap" == "Y" ]]; then
        if [[ "$whatshap_now" == "Y" && ( "$length_haplotype" == "Y" || "$pcr_or_wgs" == "0" ) && "$benchmarking" != "1" ]]; then
            for bam in "${!result_array[@]}"; do
                if [[ ${result_array["$bam"]} == "1" && -f "$bam" ]]; then
                    run_whatshap "$bam"
                fi
            done
        elif [[ ( "$whatshap_now" == "Y" && "$pcr_or_wgs" == "1" ) || "$benchmarking" == "1" ]]; then
            for bam in "$output_folder"/*.bam; do
                [[ -f "$bam" ]] || continue
                run_whatshap "$bam"
            done   
        fi
    fi


    # --- Assembly Mapping ---
    if [[ "$assembly_mapping" == "Y" ]]; then
        for dir in "$output_folder"/*haplotypes/; do
            [[ -d "$dir" ]] || continue
            if compgen -G "$dir/output_TRViz" > /dev/null; then
                if [[ "$(basename "$dir")" =~ length_haplotypes$ ]]; then
                    if [ -f "$dir/output_TRViz/best_hit_combined.fasta" ] && \
                        awk '/^>/ {if (l && l != seq) exit 1; l=seq; seq=0; next} {seq+=length} END {exit (l == seq) ? 0 : 1}' "$dir/output_TRViz/best_hit_combined.fasta"; then
                        touch "$dir/###ERROR____$(basename "$dir")_____ERROR###___length_contig____#######"
                    fi
                fi
                mkdir -p -m 777 "$dir/assembly_mapping"
                cat "$dir"/*.fastq > "$dir/assembly_mapping/$(basename "$dir").fastq"
                if [[ -f "$dir/assembly_mapping/$(basename "$dir").fastq" ]] && [[ "$pcr_or_wgs" == 0 ]]; then
                    if [[ -f "$dir/output_TRViz/best_hit_combined.fasta" ]]; then
                        run_minimap2 "$dir/assembly_mapping" \
                        "$dir/output_TRViz/best_hit_combined.fasta" \
                        "$dir/assembly_mapping/$(basename "$dir").sam"
                    elif [[ -f "$dir/output_TRViz/result_combined.fasta" ]]; then
                        run_minimap2 "$dir/assembly_mapping" \
                        "$dir/output_TRViz/result_combined.fasta" \
                        "$dir/assembly_mapping/$(basename "$dir").sam"
                    fi
                fi
                if [[ -f "$dir/assembly_mapping/$(basename "$dir").fastq" ]] && [[ "$pcr_or_wgs" == 1 ]]; then
                    find "$dir" -name "*contigs.fasta" -type f | while read -r f; do
                        fname=$(basename "$f" .fasta)   # strip .fasta extension
                        awk -v prefix="$fname" '
                        /^>/ {print ">" prefix "_" substr($0,2); next}
                        {print}
                        ' "$f" >> "$dir/combined_contigs.fasta"
                    done
                    run_minimap2 "$dir/assembly_mapping" \
                        "$dir/combined_contigs.fasta" \
                        "$dir/assembly_mapping/$(basename "$dir").sam"
                fi

                conda activate samtools-env
                samtools faidx "$dir/output_TRViz/best_hit_combined.fasta" || true
                samtools faidx "$dir/combined_contigs.fasta" || true
                samtools faidx "$dir/output_TRViz/result_combined.fasta" || true
                conda deactivate



                if [[ -f "$dir/assembly_mapping/$(basename "$dir").sam" ]]; then
                    process_sam_to_bam "$dir/assembly_mapping/$(basename "$dir").sam"
                fi

                if [[ "$pcr_or_wgs" == 0 ]] && [[ -f "$dir/output_TRViz/best_hit_combined.fasta" ]]; then
                    run_bcftools "$dir/assembly_mapping/$(basename "$dir").bam" "$dir/output_TRViz/best_hit_combined.fasta"
                    run_sniffles "$dir/assembly_mapping/$(basename "$dir").bam"
                elif [[ "$pcr_or_wgs" == 0 ]] && [[ -f "$dir/output_TRViz/result_combined.fasta" ]]; then
                    run_bcftools "$dir/assembly_mapping/$(basename "$dir").bam" "$dir/output_TRViz/result_combined.fasta"
                    run_sniffles "$dir/assembly_mapping/$(basename "$dir").bam"
                fi
                if [[ "$pcr_or_wgs" == 1 ]]; then
                    run_bcftools "$dir/assembly_mapping/$(basename "$dir").bam" "$dir/combined_contigs.fasta"
                    run_sniffles "$dir/assembly_mapping/$(basename "$dir").bam"
                fi
                
                conda activate bcftools-env
                bcf_error="$(bcftools view --no-header "$dir/assembly_mapping/$(basename "$dir")_bcf.vcf")"
                sniffles_error="$(bcftools view --no-header "$dir/assembly_mapping/$(basename "$dir")_sv.vcf")"
                if [[ -n "$bcf_error" || -n "$sniffles_error" ]]; then
                    touch "$dir/###ERROR____$(basename "$dir")_____ERROR###___vcf_non_zero____#######"
                fi
                conda deactivate
            fi
        done
    fi

    chmod -R 777 "$output_folder"
}


############################################################
# Script Entrypoint
############################################################
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi