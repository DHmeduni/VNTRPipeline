#!/bin/bash

genome_size=$1
start_seq=$2
end_seq=$3
work_dir=$4

source /opt/conda/etc/profile.d/conda.sh

all_bam=$(find "$work_dir" -maxdepth 1 -name "*haplotype*.bam")
echo $all_bam
for file in $all_bam
do
	input_base=$(basename $file .bam)
	echo $input_base
	
	#bam2fq
	conda activate samtools-env
	samtools fastq $file > "$work_dir"/$input_base.fastq
	
	
	#correct reads with canu and build consensus
	mkdir -p "$work_dir"/$input_base
	
	canu -d "$work_dir"/"$input_base" -p "$input_base" genomeSize="$genome_size" -nanopore-raw "$work_dir"/"$input_base.fastq" -stopOnLowCoverage=0 -minInputCoverage=2
	conda deactivate
    chmod a+rwx -R "$work_dir"/"$input_base"
	
	#extract seqs from .trimmed file
	start_seq_rev=$(echo $start_seq | tr 'ATCGatcg' 'TAGCtagc' | rev)
	end_seq_rev=$(echo $end_seq | tr 'ATCGatcg' 'TAGCtagc' | rev)
	awk -v ORS= '/^>/ { $0 = (NR==1 ? "" : RS) $0 RS } END { printf RS }1' $work_dir/$input_base/$input_base.contigs.fasta > $work_dir/$input_base/$input_base.contigs_no_new_line.fasta
	# shellcheck disable=SC1087
	forward=$(cat $work_dir/$input_base/$input_base.contigs_no_new_line.fasta | grep -o "$start_seq[^)]*$end_seq")
	# shellcheck disable=SC1087
	reverse=$(cat $work_dir/$input_base/$input_base.contigs_no_new_line.fasta | grep -o "$end_seq_rev[^)]*$start_seq_rev")
 
	counter=1
	for line in $forward
	do
		printf ">$counter\n$line\n" >>  $work_dir/$input_base/$input_base.trimmedReads_for.fasta
		counter=$((counter+1))
	done
	counter=1
	for line in $reverse
	do
		printf ">$counter\n$line\n" >>  $work_dir/$input_base/$input_base.trimmedReads_rev.fasta
		counter=$((counter+1))
	done
	
done

