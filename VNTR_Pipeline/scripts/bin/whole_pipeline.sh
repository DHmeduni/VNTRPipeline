#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh


#get absolute path were current script and all other side scripts are
script=$(readlink -f $0)
script_path=$(dirname $script)


helpFunction()
{
   echo ""
   echo "Usage: $0 -d work_dir -g genome_size -s start_seq -e end_seq -f motifs -m mode"
   echo -e "\t-d Absolute path of the working directory including *haplotype*.bam files to analyse, in case of WGS or AS, reads must be filtered for the vntr region before"
   echo -e "\t-g Genome size for de novo assembly e.g. '3k'"
   echo -e "\t-s Boundary sequence at the beginning of VNTR, first sequence of motif files stored in script directory"
   echo -e "\t-e Boundary sequence at the end of VNTR, last sequence of motif files stored in script directory"
   echo -e "\t-f Absolute path to .txt file, containaing one motif sequence per line"
   echo -e "\t-l Locus that is to be analyzed. -g, -s, -e and -f are automatically set. At the moment required."
   echo -e "\t-m Either 'pcr' or 'wgs' or for established PCRs 'ACAN' or 'MUC1'. At the moment required."
   exit 1 # Exit script after printing help
}

while getopts "d:g:s:e:m:f:l:" opt
do
   case "$opt" in
      d ) work_dir="$OPTARG" ;;
      g ) genome_size="$OPTARG" ;;
      s ) start_seq="$OPTARG" ;;
      e ) end_seq="$OPTARG" ;;
      f ) motifs="$OPTARG" ;;
      l ) locus="$OPTARG" ;;
      m ) mode="$OPTARG" ;;
	  ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done



genome_size="$VNTR_ASSEMBLY_SIZE"
start_seq="$VNTR_BOUNDARY_SEQUENCE_LEFT"
end_seq="$VNTR_BOUNDARY_SEQUENCE_RIGHT"
motifs=$LIB_PATH/"$VNTR_MOTIFS"


if [[ "$mode" == "wgs" ]]
then
	genome_size="50k"
fi


# Print helpFunction in case parameters are empty

if [ -z "$work_dir" ] || [ -z "$genome_size" ] || [ -z "$start_seq" ]  || [ -z "$end_seq" ]  || [ -z "$motifs" ] || [ -z "$mode" ]; then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

if { [[ "$mode" != "pcr" ]] && [[ "$mode" != "wgs" ]]; }  ||  { [[ "${locus^^}" != "ACAN" ]] && [[ "${locus^^}" != "MUC1" ]]; }; then
	echo "the mode parameter -m must be 'pcr', 'wgs' , 'ACAN' or 'MUC1'!"
	exit 1
fi

# Begin script in case all parameters are correct
echo "$work_dir"
echo "$genome_size"
echo "$start_seq"
echo "$end_seq"
echo "$mode"
echo "$script_path"
echo "$motifs"


#run canu pipeline
if [ "$mode" == "pcr" ]
then
	$script_path/canu_pcr.sh $genome_size $start_seq $end_seq $work_dir
fi

if [ "$mode" == "wgs" ]
then
	$script_path/canu_wgs.sh $genome_size $start_seq $end_seq $work_dir
fi


conda activate python-env
#analyse assemblies and compute combined fasta

Rscript $script_path/count_seqs_of_canu_trimmed.R $work_dir


cat $work_dir/output_TRViz/*best_hit.fasta > $work_dir/output_TRViz/best_hit_combined.fasta


#trviz
echo "$work_dir"

#TRVIz
python "$script_path"/trviz_script.py "$work_dir"/output_TRViz/ "$LIB_PATH"/"$VNTR_MOTIFS"


#Run statistics


Rscript $script_path/analyse_seqs_from_trviz_for_pipeline.R $LIB_PATH/$VNTR_MOTIFS $work_dir
Rscript $script_path/statistics_of_canu_trimmed_reads.R $work_dir
Rscript $script_path/trviz_plot_modified.R $work_dir $LIB_PATH/$VNTR_COLORS



conda deactivate