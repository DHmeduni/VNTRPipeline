# VNTRPipeline
VNTR Analysis Pipeline for PCR and WGS datasets


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
docker run -v /directory/to/link:/data ghcr.io/dhmeduni/vntr_pipeline:latest VNTR_pipeline -i /data -o analysis -p 0 -v MUC1 -r chm13v2.0
```

## Options

```
-h, --help    Print Help
-i            Input folder path
-o            Output folder path
-r            Reference Name hg38_p14|chm13v2.0
-p            PCR or WGS [0|1]
-v            MUC1 or ACAN, more VNTRs will be added
-q/--quiet (future)
--version (future)
```

## Output

The output files can be found in the input directory, under the given output directory name (analysis),
in the directory named for the sample.

```
/input_folder/output_folder/Sample_Name_results
```

The following files are then of particular interest:
- best_trviz_fig.png
- new_and_lof_seq.xlsx


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
If the length based separation of haplotypes results in an incorrect separation,
then an error message is outputed into the output folder with error description:
```
###ERROR____Sample_Name_____ERROR###___length_contig____#######
```


  
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
