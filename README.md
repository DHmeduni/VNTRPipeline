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


## Installation
The necessary prerequisites, scripts and variables for the analysis of
*MUC1* and *ACAN* VNTRs are already precompiled in a docker image. \
[Docker](https://docs.docker.com/get-started/get-docker/) need to be installed on your system.

```
docker pull ghcr.io/dhmeduni/vntr_pipeline:latest
```

## Test Data

Data to test the workflow (PCR Amplicons Sequencing data of the MUC1 VNTR 
from HG001 through HG004) is inside the docker image ( in the test_data folder)

```
docker run -v /output_directory/to/link:/data ghcr.io/dhmeduni/vntr_pipeline:latest VNTR_pipeline -i /test_data -o /data/analysis -r chm13 -p pcr -v MUC1 
```


## Usage

The container can be executed using data in a mounted drive (example usage):

```
docker run -v /directory/to/link:/data ghcr.io/dhmeduni/vntr_pipeline:latest VNTR_pipeline -i /data -o /data/analysis -r chm13 -p pcr -v MUC1 
```

## Options

```
-h, --help    Print Help

# required
-i            Input folder path or Input file name (.bam or fastq/fastq.gz format)
-o            Output folder path
-r            hg38 or chm13 (reference genome)
-p            pcr or wgs
-v            MUC1 or ACAN, (motif seqs and additional parameters are already set)


#For other VNTRs/motifs, instead of -v, following additional options have to be set when executing VNTR_pipeline.
These variables can be committed to the configuration file by adding the option -c {VNTR NAME}:
VNTR_BOUNDARY_SEQUENCE_LEFT    Sequence at beginning of VNTR (will be used to cut out the VNTR, one mismatch allowed)
VNTR_BOUNDARY_SEQUENCE_RIGHT   Sequence at end of VNTR (will be used to cut out the VNTR, one mismatch allowed)
VNTR_ASSEMBLY_SIZE             Approximate size of VNTR (e.g.: 500, 1k, 3k, 10k,), used for PCR data
VNTR_MIN_PRODUCT_SIZE          Minimum PCR Product Size expected, necessary for filtering
VNTR_REPEAT_SIZE               The Number of bases found in one repeat unit of the VNTR
VNTR_COORDINATES_HG38          T2T-chm13v2.0 Referernce coordinates of the VNTR for filtering
VNTR_COORDINATES_CHM13         GRCh38.p14 Referernce coordinates of the VNTR for filtering
VNTR_MOTIFS                    Path to VNTR Motif file, containing Sequencing, Aplhanumeric Designation and Color code

#Following options may be added:
DELETE_TMP=Y               Allows user to retain the temporary files for troubleshooting
ALL_FIGURES=Y              Only produces Figures, based on already assembled and trimmed sequences in *best_hit.fasta files,
                           recursively to a depth of one subfolder from a folder that is given as input
NON_CODING=Y               Analyse VNTR’s in non-coding regions, LoF prediction is skipped
VNTR_PACBIO=Y              Process PacBio WGS Data
WHATSHAP_FORCE=Y           Force Whatshap haplotyping
VNTR_ALL=Y                 Analyse all VNTR’s (assemblies or polished sequences) found by the workflow (pseudogenes/duplications)
CONFIG_FILE                Define another path for the configuration file (CONFIG_FILE=/path/to/file)
LENGTH_1                   Define the shorter of two lengths for length based haplotyping (e.g. LENGTH_1=2500), default is automatic
LENGTH_2                   Define the longer of two lengths for length based haplotyping (e.g LENGTH_2=3000), default is automatic
MIN_FREQUENCY              Temporarly lower the minimum frequency threshold for length based haplotyping (default MIN_FREQUENCY=20)

```
Pre-saved options vor *ACAN* and *MUC1* can be found here:
- [configuration_file](https://github.com/DHmeduni/VNTRPipeline/blob/f7c60b42f65db005dcbfd38ac3a87cf833541033/VNTR_Pipeline/scripts/lib/vntr_variables.cfg#L1C1-L18C44)

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
