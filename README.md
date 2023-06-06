# pipebox

<!-- SOFTWARE pipebox -->
<!-- VERSION 1.0.0 -->

Pipeline in a Docker container

- [Background](#background)
  - [Motivation](#motivation)
  - [Preliminaries](#preliminaries)
- [Demo](#demo)
  - [Build](#build)
  - [Run](#run)
  - [Results](#results)
- [How it works](#how-it-works)
  - [The `Dockerfile`](#the-dockerfile)
  - [The pipeline script](#the-pipeline-script)
  - [The post-processing script](#the-post-processing-script)
- [For developers](#for-developers)

<!--

# Build and ship: 
# export GHCR_PAT="YOUR_TOKEN_HERE"
# echo $GHCR_PAT | docker login ghcr.io -u USERNAME --password-stdin
docker build build --no-cache --tag pipebox:latest --tag ghcr.io/colossal-compsci/pipebox:latest .
docker push ghcr.io/colossal-compsci/pipebox:latest

-->

## Background

### Motivation

This repository demonstrates how to build an entire data analysis pipeline that runs inside a Docker container. 

Docker is a virtualization technology that can be used to bundle an application and all its dependencies in a virtual container that can be distributed and deployed to run reproducibly on any Windows, Linux or MacOS operating system. Docker is commonly used in computational biology to bundle a tool and all its dependencies such that it can essentially be used as an executable (see the [StaPH-B Docker](https://github.com/StaPH-B/docker-builds) and [slimbioinfo](https://github.com/stephenturner/slimbioinfo) projects for examples). For instance, installing something as simple as `bwa` on a modern ARM Macbook is challenging because of the way M1 Macs compile C code. However, a Docker container can be used where `bwa` can be essentially replaced with `docker run bwa ...`.

However, while containers are typically used to deploy microservices, containerization can also be used to deploy an entire data analysis pipeline. I illustrate this in pipebox using a fairly trivial data analysis pipeline involving read mapping and variant calling.

Inputs: 

1. Paired end reads
1. A reference genome
1. An output base filename

What the pipeline does:

1. If the reference genome isn't indexed, it'll index it with both `samtools faidx` and `bwa index`.
1. Very simple QC on the sequencing reads using `seqtk`.
1. Use `bwa mem` to map reads to the indexed reference genome.
1. Very simple QC on the alignments with `samtools stats`.
1. Use `bcftools mpileup` and `bcftools call` to call variants.
1. Very simple QC on the alignments with `samtools stats`.
1. Very simple QC on the variant calls with `bcftools stats`.
1. Creates a simple single page PDF "report" with a few plots from the QC steps.

All of these data processing steps happen _inside the container_. You provide a volume mount to the container so that you can mount the local filesystem to a location inside the running container. Data processing happens inside the container. Output files are also written to a location that is mounted to the container (e.g., the present working directory). This pipeline will run on any system where Docker is available. 

### Preliminaries

This repo assumes a minimal working understanding of containerization, including the difference between an image and a container, and how they're each created. If you're familiar with all these concepts, you can skip to the [Demo](#demo) section below.

A Docker _image_ is like a blueprint or a snapshot of what a container will look like when it runs. Images are _created from code_ using a Dockerfile, which is a set of instructions used to create the image. This is important, because with Dockerfiles we can _define infrastructure in code_ -- code which can be easily version controlled and worked on collaboratively in a platform like GitHub. When we have our instructions for building an image written out in a Dockerfile, we can create an image using `docker build`.

Once we have an image, we can "rehydrate" that image and bring to life a running _container_ using `docker run`. When a container is initialized, it spins up, does whatever it needs to do, and is destroyed. A container could serve a long-running process, such as a webservice or database. Here, we're spinning up a container that runs a script _inside that container_. When the script runs to completion, the container is destroyed. From a single image we can use `docker run` to spin up one running container or one thousand running containers, all starting from the identical image. 

![Docker images versus containers.](img/dockerfile-image-container.png)

Another important concept is a _volume mount_. A running container is a virtualized operating system with its own software and importantly, its own _filesystem_ inside the container. If we have a plain Alpine linux container and run something like `docker run alpine ls`, this will start the alpine container and run `ls` in the default working directory, which is `/`, _inside the container_. That is, you'll see the usual system folders like `/bin`, `/dev`, `/etc`, and others. If you were to run `docker run alpine ls /home` you'll see _nothing_ - because the container has no users other than root! 

This concept is important because if you're running tools inside the container, you must provide a way to have files on the host system (i.e., your computer or VM) to be accessible to the running container. To do this we use the `--volume` or more commonly the short `-v` flag, which allows us to mount local directories to the container. For instance, `docker run -v /home/turner/mydata:/data` would mount `/home/turner/mydata` from my laptop to `/data` inside the running container. If I have a script that expects to see data in `/data` inside the container, this would work well. 

Here's a common trick that allows you to mount contents of the present working directory on your machine to the running container using shell substitution with `$(pwd)`. Using `-v $(pwd):(pwd)` mounts the present working directory on the host system to a path of the same name inside the container, and `-w $(pwd)` _sets_ the working directory inside the container to the same directory specified by the volume mount. 

Further suggested reading:

- **StaPH-B docker user guide**: https://staphb.org/docker-builds/. This is an excellent guide to using Docker in bioinformatics written by a group of bioinformaticians working in public health labs. Start with the chapters in the upper-right. Specifically useful are the [running containers chapter](https://staphb.org/docker-builds/run_containers/) which provides additional detail on volume mounts, and [developing containers chapter](https://staphb.org/docker-builds/make_containers/) which discusses creating Dockerfiles, building images, and initializing containers. 
- **"Ten simple rules for writing Dockerfiles for reproducible data science."** _PLoS computational biology_ 16.11 (2020): e1008316. DOI:[10.1371/journal.pcbi.1008316](https://doi.org/10.1371/journal.pcbi.1008316). This paper provides a great overview on containerization with a specific focus on the importance of Dockerfiles for reproducible data science and bioinformatics workflows.
- **"pracpac: Practical R Packaging with Docker."** arXiv:2303.07876 (2023). DOI:[10.48550/arXiv.2303.07876](https://doi.org/10.48550/arXiv.2303.07876). Shameless self-promo: I wrote this paper and the software it describes for automating the building and deployment of pipelines-as-Docker-containers with specific attention to pipelines which require a custom-built R package as part of the pipeline. The software is implemented as an R package and it's R-centric (Shiny, MLOps with tidymodels, etc.), but isn't limited solely to R and R-related tools, as demonstrated in the package vignettes.


## Demo

### Build

First, get this repository and build the pipebox image.

```sh
git clone git@github.com:colossal-compsci/pipebox.git
cd pipebox
docker build --tag pipebox .
```

Alternatively, [create a personal access token](https://github.com/settings/tokens) with read permissions for the GitHub Container Registry, login (once), then [pull the container directly from the GitHub Container Registry](https://ghcr.io/colossal-compsci/pipebox).

```sh
## Run once
# export GHCR_PAT="YOUR_TOKEN_HERE"
# echo $GHCR_PAT | docker login ghcr.io -u USERNAME --password-stdin

# Pull the image
docker pull ghcr.io/colossal-compsci/pipebox

# Give it a short name so you can run it more easily
docker tag ghcr.io/colossal-compsci/pipebox pipebox
```

### Run

Running the container with no arguments prints a minimal help message (we use `--rm` to destroy the container when it exits). The pipebox container is running the [pipebox.sh](src/pipebox.sh) script as its `ENTRYPOINT`, which runs the entire pipeline.

```sh
docker run --rm pipebox
```

```
Usage: pipebox.sh <read1> <read2> <reffa> <outbase>
```

Let's run it on the [test data in this repository](testdata). The `-v $(pwd):(pwd)` mounts the present working directory on the host system to a path of the same name inside the container, while the `-w $(pwd)` _sets_ the working directory inside the container to the same directory specified by the volume mount. These two flags make it easy to make data on the host accessible to the running container.

```sh
cd testdata
docker run --rm -v $(pwd):$(pwd) -w $(pwd) pipebox \
    SRR507778_1.fastq.gz \
    SRR507778_2.fastq.gz \
    yeastref.fa.gz \
    results/SRR507778
```

### Results

You should see all the files in the `testdata/results` folder:

- [`SRR507778.seqtkfqchk.tsv`](testdata/results/SRR507778.seqtkfqchk.tsv): Results from running `seqtk fqchk` after interleaving the paired FASTQ files.
- [`SRR507778.bam`](testdata/results/SRR507778.bam): Sorted alignment.
- [`SRR507778.samtoolsstats.tsv`](testdata/results/SRR507778.samtoolsstats.tsv): Results from running `samtools stats` on this alignment.
- [`SRR507778.vcf.gz`](testdata/results/SRR507778.vcf.gz): Results from variant calling against the reference genome.
- [`SRR507778.bcfstats.tsv`](testdata/results/SRR507778.bcfstats.tsv): Results from running `bcftools stats` on these variant calls.
- [`SRR507778.metrics.pdf`](testdata/results/SRR507778.metrics.pdf): "Report" compiling a few outputs from each of the tools above (PDF format).
- [`SRR507778.metrics.png`](testdata/results/SRR507778.metrics.png): "Report" compiling a few outputs from each of the tools above (PNG format).

![Report compiling a few outputs from each of the tools above.](testdata/results/SRR507778.metrics.png)

