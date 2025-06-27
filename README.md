# VNTRPipeline
VNTR Analysis Pipeline for PCR and WGS datasets


## Introduction
Please refer to our publication at BioRxiv. https://doi.org/10.1101/2023.09.08.556789

This Pipeline can be fed with BAM or FASTQ formatted sequencing 
information to detect and genotype
VNTR variantions in Amplicon or Native DNA Long-Read Sequencing Platforms

## Table of Contents
- [Docker Image](#docker-image)
- [Project Workflow](#project-workflow)
- [License](#license)


## Docker Image
The necessary prerequisites, scripts and variables for the analysis of 
MUC1 and ACAN VNTRs is already precompiled in a docker image.

https://github.com/users/DHmeduni/packages/container/package/vntr_pipeline

```
docker pull ghcr.io/medgenmeduniwien/vntr_pipeline:latest
```


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
