# VNTRPipeline
VNTR Analysis Pipeline for PCR and WGS datasets.

This pipeline is associated with the research paper: (Link will be added upon paper acceptance)

## Introduction

This Pipeline can be fed with BAM or FASTQ formatted sequencing 
information to detect and genotype
VNTR variants in Amplicon or Native DNA Long-Read Sequencing Data.

## Table of Contents
- [Installation](#Installation)
- [Options](#options)
- [Test Data](#test-data)
- [Usage](#usage)
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

## Options

```
-h, --help    Print Help

# required
-i            Input folder path (all existing bam and fastq/fastq.gz files are used for analysis) or specific input file name directly
-o            Output folder path
-r            hg38 or t2t (reference genome)
-p            pcr or wgs
-v            MUC1 or ACAN, (motif seqs and additional parameters are already set)

#For other VNTRs/motifs, instead of -v, the following additional options have to be set when executing VNTR_pipeline.
These variables can be committed to the configuration file by adding the option -c {VNTR NAME} and must be added after the above options.
Note: any of the below options used with the –v {VNTR_NAME} option will override the values for the current run.
VNTR_BOUNDARY_SEQUENCE_LEFT    Sequence at beginning of VNTR (will be used to extract the VNTR, one mismatch allowed)
VNTR_BOUNDARY_SEQUENCE_RIGHT   Sequence at end of VNTR (will be used to extract the VNTR, one mismatch allowed)
VNTR_ASSEMBLY_SIZE             Approximate size of PCR Product (e.g.: 500, 1k, 3k, 10k,), used for PCR data, minimum should be 500
VNTR_MIN_PRODUCT_SIZE          Minimum PCR Product Size expected, necessary for filtering
VNTR_REPEAT_SIZE               The Number of bases found in one repeat unit of the VNTR
VNTR_COORDINATES_T2T           T2T-chm13v2.0 Reference coordinates of the VNTR for filtering in the format [chr]:[start]-[end]
VNTR_COORDINATES_HG38          GRCh38.p14 Reference coordinates of the VNTR for filtering in the format [chr]:[start]-[end]
VNTR_MOTIFS                    Path to VNTR Motif file, containing Sequencing, Aplhanumeric Designation and Color code,
                               tab separated .txt file, example files are available on github

#Following options may be added:
DELETE_TMP=N               Retain the temporary files for troubleshooting
ALL_FIGURES=Y              Only produces figures, based on already assembled and trimmed sequences of *best_hit.fasta files,
                           recursively to a depth of one subfolder from a folder that is given as input
NON_CODING=Y               Analyse VNTRs in non-coding regions, LoF prediction is skipped
VNTR_PACBIO=Y              Process PacBio WGS Data
WHATSHAP_FORCE=Y           Force Whatshap haplotyping
VNTR_ALL=Y                 Analyse all VNTRs (assemblies or polished sequences the top three most common sequences per haplotype)
                           found by the workflow (pseudogenes/duplications)
CONFIG_FILE                Define another path for the configuration file (CONFIG_FILE=/path/to/file)
LENGTH_1                   Define the shorter of two lengths for length based haplotyping (e.g. LENGTH_1=2500), default is automatic
LENGTH_2                   Define the longer of two lengths for length based haplotyping (e.g LENGTH_2=3000), default is automatic
MIN_FREQUENCY              Temporarily lower the minimum frequency threshold for length based haplotyping (default MIN_FREQUENCY=20)

```
Pre-saved options for *ACAN* and *MUC1* can be found here:
- [configuration_file](https://github.com/DHmeduni/VNTRPipeline/blob/f7c60b42f65db005dcbfd38ac3a87cf833541033/VNTR_Pipeline/scripts/lib/vntr_variables.cfg#L1C1-L18C44)

## Test Data

Data to test the workflow (PCR Amplicons Sequencing data of the MUC1 VNTR 
from HG001 through HG004) is inside the docker image (inside test_data)

```
docker run -v /output_directory/to/link:/data ghcr.io/dhmeduni/vntr_pipeline:latest VNTR_pipeline -i /test_data -o /data/analysis -r t2t -p pcr -v MUC1 
```


## Usage

The container can be executed using data in a mounted drive (example usage),

```
docker run --rm -v /directory/to/link:/data ghcr.io/dhmeduni/vntr_pipeline:latest VNTR_pipeline -i /data -o /data/analysis -r t2t -p pcr -v MUC1
```

or in interactive mode:
```
docker run -it -v /directory/to/link:/data ghcr.io/dhmeduni/vntr_pipeline:latest VNTR_pipeline
VNTR_pipeline -i /data -o /data/analysis -r t2t -p pcr -v MUC1
```
For other VNTRs following parameters are necessary (example for CEL WGS data)

```
VNTR_pipeline -i <input_folder containing HG002_CEL_WGS.fastq> -o <output_folder> -r hg38 -p wgs VNTR_BOUNDARY_SEQUENCE_LEFT=TATCTGGCGCTGCCCACAGTGACCGACCAG VNTR_BOUNDARY_SEQUENCE_RIGHT=AAGGAAGCTCAGATGCCTGCAGTCATTAGGTTTTAG VNTR_ASSEMBLY_SIZE=560 VNTR_MIN_PRODUCT_SIZE=200 VNTR_REPEAT_SIZE=33 VNTR_COORDINATES_HG38=chr9:133071168-133071728 VNTR_MOTIFS=<folder_to>/motifs_cel_with_char_scheme.txt (VNTR_ALL=Y)
```

## Output
**Files in \*results Folder (one folder per sample)**
| **_File name_**                        | **_Description_**                                                                                                                                                                                                                                                                                                        |
| -------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| \*_haplotype1_best_hit.fasta           | Assembled/Consensus sequence of haplotype 1. The most common sequence found after error correction in PCR and filtering using boundary sequences and the first hit (closest to boundary sequence left) in WGS for haplotype 1                                                                                            |
| \*_haplotype2_best_hit.fasta           | Assembled/Consensus sequence of haplotype 2. The most common sequence found after error correction in PCR and filtering using boundary sequences and the first hit (closest to boundary sequence left) in WGS for haplotype 2                                                                                            |
| \*_haplotype1_result.fasta             | All unique sequences with counts found after filtering with the boundary sequences in PCR, or all hits found in WGS for haplotype 1                                                                                                                                                                                      |
| \*_haplotype2_result.fasta             | All unique sequences with counts found after filtering with the boundary sequences in PCR, or all hits found in WGS for haplotype 2                                                                                                                                                                                      |
| \*_length_histogram.png                | If PCR a read-length histogram is produced after filtering for the alignment region                                                                                                                                                                                                                                      |
| \*_length_histogram.txt                | If PCR a text file containing variables used and results of length-based haplotype separation (min_product_size, range_for_length_filtering, repeat_size, min_frequency \\ number of reads/entries für halpotype 1  \\ number of reads/entries für halpotype 2 \\ size range of haplotype 1 \\ size range of haplotype 2 |
| best_trviz_fig.png                     | VNTR composition figure                                                                                                                                                                                                                                                                                                  |
| motifs_evaluation.xlsx                 | All found VNTR motifs including their occurrences and ID (assigned letter)                                                                                                                                                                                                                                               |
| new_and_lof_seqs.xlsx                  | All new and loss-of-function motifs and their positions within the VNTR are shown here                                                                                                                                                                                                                                   |
| seq_distribution.xlsx                  | Distribution of unique VNTR sequences after assembly/correction                                                                                                                                                                                                                                                          |
| **_Part of Assembly_Mapping folder:_** |                                                                                                                                                                                                                                                                                                                          |
| \*haplotypes.bam                       | Uncorrected haplotype assigned Reads aligned to the best_hit_combined.fasta file                                                                                                                                                                                                                                         |
| \*bcf.vcf                              | BCF variant file from assembly alignment                                                                                                                                                                                                                                                                                 |
| \*sv.vcf                               | Sniffles2 variant file from assembly alignment                                                                                                                                                                                                                                                                           |
| \*best_hit_combined.fasta              | Reference file for alignment of filtered reads                                                                                                                                                                                                                                                                           |
**Other files/folders**
| all_figures folder containing \*best_trviz_fig.png | Graphic VNTR output of all haplotypes and samples                                        |
| All_Samples_new_and_lof_seqs.xlsx                  | All new and loss-of-function motifs and their positions combined for all samples         |
| All_Samples_seqs_distribution.xlsx                 | Distribution of unique VNTR sequences after assembly/correction combined for all samples |
| script.log                                         | Log file                                                                                 |

The following files are then of particular interest:
- [best_trviz_fig.png](best_trviz_fig.pdf)
- [new_and_lof_seq.xlsx](new_and_lof_seq.pdf)


## Quality Control

The file *seqs_distribution.xlsx* in the *results folder shows the distribution of all corrected and trimmed sequences per haplotype. When analyzing PCR data, the most common sequence should represent a majority of all sequences found (first value compared with all subsequent values). If the first and second values are not significantly different, then there are two sequences with similar frequency present in one haplotype, indicating a problem with haplotype separation or the presence of a mosaic background. When analyzing WGS data, the file indicates the number of assembled and extractd VNTR sequences.

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
