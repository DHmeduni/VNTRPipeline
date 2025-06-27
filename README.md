# VNTRPipeline
VNTR Analysis Pipeline for PCR and WGS datasets


## Introduction
Please refer to our publication at BioRxiv. https://doi.org

This Pipeline can be fed with BAM or FASTQ formatted sequencing 
information to detect and genotype
VNTR variants in Amplicon or Native DNA Long-Read Sequencing Data.

## Table of Contents
- [Docker Image](#docker-image)
- [Use](#use)
- [Output](#output)
- [Project Workflow](#project-workflow)
- [License](#license)


## Docker Image
The necessary prerequisites, scripts and variables for the analysis of 
MUC1 and ACAN VNTRs is already precompiled in a docker image.

https://github.com/users/DHmeduni/packages/container/package/vntr_pipeline

```
docker pull ghcr.io/dhmeduni/vntr_pipeline:latest
```

## Use

Data to test the workflow (PCR Amplicons Sequencing data of the MUC1 VNTR 
from HG001 through HG004) is inside the docker image (test_data)

The pipeline can be used to test in interactive mode.

```
docker run --it ghcr.io/dhmeduni/vntr_pipeline:latest
```

And the pipeline started as follows:

```
/scripts/VNTR_pipeline.sh -i /test_data -o analysis -p 0 -v MUC1 -r chm13v2.0
```

The container can be executed using data in a mounted drive:

```
docker run -v /directory/to/link:/data ghcr.io/dhmeduni/vntr_pipeline:latest /scripts/VNTR_pipeline.sh -i /data -o analysis -p 0 -v MUC1 -r chm13v2.0
```

## Output

The output files can be found in the input directory, under the given output directory name (analysis),
in the directory containing the file name that was analyzed ending in haplotypes, then in the
subdirectory output_TRViz.

```
/input_folder/output_folder/sample_haplotypes/output_TRViz
```

The following files are then of interest:
- best_trviz_fig_test.png
- new_andlof_seq.xlsx

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
All code used to generate the figures and analyses in this manuscript
are available in this repository, in the `scripts` repo
