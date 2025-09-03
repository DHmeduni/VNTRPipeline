# VNTRPipeline
VNTR Analysis Pipeline for PCR and WGS datasets.

This pipeline is associated with the research paper: (Link will be added upon paper acceptance)

## Introduction

This Pipeline can be fed with BAM or FASTQ formatted sequencing 
information to detect and genotype
VNTR variants in Amplicon or Native DNA Long-Read Sequencing Data.

## Table of Contents
- [Docker Image](#docker-image)
- [Test Data](#test-data)
- [Usage](#usage)
- [Options](#options)
- [Output](#output)
- [Quality Control](#quality-control)
- [Project Workflow](#project-workflow)
- [License](#license)


## Docker Image
The necessary prerequisites, scripts and variables for the analysis of 
MUC1 and ACAN VNTRs are already precompiled in a docker image.

https://github.com/users/DHmeduni/packages/container/package/vntr_pipeline

```
docker pull ghcr.io/dhmeduni/vntr_pipeline:latest
```

## Test Data

Data to test the workflow (PCR Amplicons Sequencing data of the MUC1 VNTR 
from HG001 through HG004) is inside the docker image (test_data)

The pipeline can be used to test in interactive mode.

```
docker run -it ghcr.io/dhmeduni/vntr_pipeline:latest
```

And the pipeline started as follows:

```
VNTR_pipeline -i /test_data -o analysis -p 0 -v MUC1 -r chm13v2.0
```

## Usage

The container can be executed using data in a mounted drive:

```
docker run -v /directory/to/link:/data ghcr.io/dhmeduni/vntr_pipeline:latest VNTR_pipeline -i /test_data -o /data/analysis -p 0 -v MUC1 -r chm13v2.0
```

## Options

```
-h, --help    Print Help
--version    (future)
q/--quiet    (future)

# required
-i            Input folder path/Input file name (BAM or FASTQ)
-o            Output folder path
-r            hg38 or t2t (reference genome)
-p            pcr or wgs
-v            MUC1 or ACAN, (motif seqs and additional parameters are alredy set)


#For other VNTRs/motifs, instead of -v following additional options have to be set
when executing VNTR_pipeline (MUC1 as an example):

VNTR_BOUNDARY_SEQUENCE_LEFT=AAGGAGACTTCGGCTACCCAGAGAAGTTCAGTGCCCAGCTCTACTGAGAAGAATGCTGTG \
VNTR_BOUNDARY_SEQUENCE_RIGHT=GGCTCCACCGCCCCTCCAGTCCACAATGTCACCTCGGCCTCAGGCTCTGCATCAGGCTCA \
VNTR_ASSEMBLY_SIZE=5k \
VNTR_MIN_PRODUCT_SIZE=2000 \
VNTR_REPEAT_SIZE=60 \
VNTR_COORDINATES_HG38_CHR=chr1 \
VNTR_COORDINATES_HG38_START=155188487 \
VNTR_COORDINATES_HG38_END=155192239 \
VNTR_COORDINATES_CHM13_CHR=chr1 \
VNTR_COORDINATES_CHM13_START=154328103 \
VNTR_COORDINATES_CHM13_END=154330802 \
VNTR_MOTIFS=motifs_muc1_with_char_sheme.txt \
VNTR_COLORS=motifs_muc1_with_char_sheme.txt \

These variables can be committed to the configuration file by adding the option -c {VNTR NAME}.
```
## Output

The output files can be found in the input directory, under the given output directory name (analysis),
in the directory named for the sample.

```
/input_folder/output_folder/Sample_Name_results
```

The following files are then of particular interest:
- [best_trviz_fig.png](best_trviz_fig.pdf)
- [new_and_lof_seq.xlsx](new_and_lof_seq.pdf)


## Quality Control

Assembly Information is found in the directory /assembly_mapping
- best_hit_combined.fasta (Reference)
- Sample_Name.bam (Alignment File)
- Sample_Name_sv.vcf (Sniffles Variant File)
- Sample_Name_bcf.vcf (BCF Variant File)

If either the Sniffles or BCF Variant file are non-zero (contain variants),
then an error message is outputed into the output folder with error description:
```
###ERROR____Sample_Name_____ERROR###___vcf_non_zero____#######
```
There are many reasons for this (poor read quality, poor sequencing results, etc.)
setting DELETE_TMP=N can allow you to troublshoot the problem, the tmp directory contains all files

If the length based separation of haplotypes results in an incorrect separation,
then an error message is outputed into the output folder with error description:
```
###ERROR____Sample_Name_____ERROR###___length_contig____#######
```
There are several reasons for this (poor read quality, poor sequencing results, etc.)
Setting WHATSHAP_FORCE=Y will force the haplotypes to be phased using variant data.
  
## Project Workflow
![Alt text](/VNTRPipeline_workflow.png?raw=true "Project workflow")

! Basecalling and Demultiplexing is performed externally to the VNTR_pipeline

## License

MIT License

Copyright (c) 2025 MedGenMedUniWien

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE





### Code
All code used to generate the figures and analyses are available upon request
