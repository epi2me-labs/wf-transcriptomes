# Workflow template

This repository contains a nextflow workflow template and associated
container build.

## Quickstart

```bash
# build the container
CONTAINER_TAG=template-workflow
docker build \
    -t ${CONTAINER_TAG} -f Dockerfile \
    --build-arg BASEIMAGE=epi2melabs/base-workflow-image:latest \
    .
```

The `BASEIMAGE` argument here can be changed to use an alternative image.


The template includes a simple workflow that outputs a file with the lengths
of sequences contained in a .fastq.gz file.
```
# run the pipeline with the test data
OUTPUT=template-workflow
nextflow run workflow.nf \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --reads test_data/reads.fq.gz \
    --out_dir ${OUTPUT}
```

The output of the pipeline will be found in `./template-workflow` for the above
example. This directory contains the nextflow working directories alongside
the two primary outputs of the pipeline.
