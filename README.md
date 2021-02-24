# Workflow template

This repository contains a Nextflow workflow template and associated Docker
container build. The workflow also supports using conda environments as an
alternative software isolation method to Docker.

## Quickstart

### Building the container

> This step is not necessary if you intend to run the workflow using
> conda environments.

The Docker container image can be built with the following command:

```bash
CONTAINER_TAG=ontresearch/workflow-template
docker build \
    -t ${CONTAINER_TAG} -f Dockerfile \
    --build-arg BASEIMAGE=ontresearch/base-workflow-image:v0.1.0 \
    .
```

The `BASEIMAGE` argument here can be changed to use an alternative image.

### Running the workflow

The template includes a simple workflow that outputs a file with the lengths
of sequences contained in a .fastq.gz file.

**Running the workflow with Docker containers**

To run the workflow using Docker containers supply the `-profile standard`
argument to `nextflow run`:

```
OUTPUT=workflow-template
nextflow run main.nf \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --fastq test_data/reads.fq.gz \
    --out_dir ${OUTPUT}
```

The output of the pipeline will be found in `./workflow-template` for the above
example. This directory contains the nextflow working directories alongside
the two primary outputs of the pipeline.

**Using conda environments**

To run the workflow backed by conda environments, simply provide the
`-profile conda` argument to `nextflow run`.

```
# run the pipeline with the test data
OUTPUT=workflow-template
nextflow run main.nf \
    -w ${OUTPUT}/workspace \
    -profile conda \
    --fastq test_data/reads.fq.gz \
    --out_dir ${OUTPUT}
```

This will create a conda environment with all required software within the
workspace directory. When running multiple analyses on distinct datasets
it may not be desirable to have Nextflow create a conda environment for each
analysis. To avoid the situation editing the file `nextflow.config` will
be necessary. Search for the term `cacheDir` and set this to a directory
where you wish the conda environment to be placed.
