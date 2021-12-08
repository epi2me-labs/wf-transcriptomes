#!/usr/bin/env bash

# A couple of simple tests with different combinations of CLI options

if [[ "$#" -ne 1 ]]; then
    echo "Please supply path to out_dir"
    exit 1
fi

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd $SCRIPT_DIR;

singledir="test_data/fastq"
multisampledir="test_data/demultiplexed_fastq"

echo $singledir
echo $multisampledir
results=()

# Test1 single sample
OUTPUT=$1/test1;
nextflow run . --fastq $singledir --ref_genome test_data/SIRV_150601a.fasta \
--ref_annotation test_data/SIRV_isofroms.gtf -profile conda --out_dir ${OUTPUT} -w ${OUTPUT}/workspace -resume;

results+=("test1: $?")

# Test2 multiple samples demultiplexed
OUTPUT=$1/test2;
nextflow run . --fastq $multisampledir --ref_genome test_data/SIRV_150601a.fasta \
--ref_annotation test_data/SIRV_isofroms.gtf -profile conda --out_dir ${OUTPUT} -w ${OUTPUT}/workspace \
--sample_sheet test_data/sample_sheet -resume;

results+=("test2: $?")

# Test3 single sample. No reference annotation
OUTPUT=$1/test3;
nextflow run . --fastq $singledir --ref_genome test_data/SIRV_150601a.fasta \
 -profile conda --out_dir ${OUTPUT} -w ${OUTPUT}/workspace -resume;

results+=("test3: $?")

# Test4 single sample. force split_bam to make multiple alignment bundles
OUTPUT=$1/test4;
nextflow run . --fastq $singledir  --ref_genome test_data/SIRV_150601a.fasta \
--ref_annotation test_data/SIRV_isofroms.gtf -profile conda --out_dir ${OUTPUT} -w ${OUTPUT}/workspace \
--bundle_min_reads 5 -resume;

results+=("test4: $?")

echo "Exit status codes for each test"
for value in "${results[@]}"; do
     echo "${value}"
done