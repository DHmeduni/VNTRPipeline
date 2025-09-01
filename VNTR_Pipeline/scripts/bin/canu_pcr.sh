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
	#input_dir=$(dirname $file)
	echo $input_base
	
		
	#bam2fq
	conda activate samtools-env
	samtools fastq $file > "$work_dir"/$input_base.fastq
	
	
	#correct reads with canu and build consensus
	
	canu -d "$work_dir"/"$input_base" -p "$input_base" genomeSize="$genome_size" -nanopore-raw "$work_dir"/"$input_base.fastq" -readSamplingCoverage=100 -contigFilter="2 0 1.0 0.5 0"
	conda deactivate
    chmod a+rwx -R "$work_dir"/"$input_base"
	
	#extract seqs from .trimmed file
	start_seq_rev=$(echo $start_seq | tr 'ATCGatcg' 'TAGCtagc' | rev)
	end_seq_rev=$(echo $end_seq | tr 'ATCGatcg' 'TAGCtagc' | rev)
	forward=$(zcat $work_dir/$input_base/$input_base.trimmedReads.fasta.gz | grep -o "$start_seq[^)]*$end_seq")
	reverse=$(zcat $work_dir/$input_base/$input_base.trimmedReads.fasta.gz | grep -o "$end_seq_rev[^)]*$start_seq_rev")

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

