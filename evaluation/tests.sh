#!/usr/bin/env bash

# A few simple tests with different combinations of CLI options
# Run from within an appropriate active conda environment

if [[ "$#" -lt 1 ]]; then
    echo "usage: tests.sh <outdir> [nextflow.config]"
    exit 1
fi

if [[ "$#" -eq 1 ]]; then
    config=''
fi

if [[ "$#" -eq 2 ]]; then
  config="-c $2";
fi

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd $SCRIPT_DIR/../;

singledir="test_data/fastq"
multisampledir="test_data/demultiplexed_fastq"

# This is for when using SIRV dataset with non-canonical spice junctions
#"--minimap2_opts '-uf --splice-flank=no'"
results=()

# Reference based tests
OUTPUT=$1/reference_single_dir;
nextflow run . --fastq $singledir $config --ref_genome test_data/SIRV_150601a.fasta --minimap2_opts '-uf --splice-flank=no' \
--ref_annotation test_data/SIRV_isoforms.gtf -profile local --out_dir ${OUTPUT} -w ${OUTPUT}/workspace -resume;
r=$?
results+=("$(basename $OUTPUT): $r")

OUTPUT=$1/multiple_samples;
nextflow run . --fastq $multisampledir $config --ref_genome test_data/SIRV_150601a.fasta --minimap2_opts '-uf --splice-flank=no'\
--ref_annotation test_data/SIRV_isoforms.gtf -profile local --out_dir ${OUTPUT} -w ${OUTPUT}/workspace \
--sample_sheet test_data/sample_sheet.csv -resume;
r=$?
results+=("$(basename $OUTPUT): $r")

OUTPUT=$1/reference_no_ref_annotation;
nextflow run . --fastq $singledir $config --ref_genome test_data/SIRV_150601a.fasta --minimap2_opts '-uf --splice-flank=no'\
 -profile local --out_dir ${OUTPUT} -w ${OUTPUT}/workspace -resume;
r=$?
results+=("$(basename $OUTPUT): $r")

# Force split_bam to make multiple alignment bundles
OUTPUT=$1/reference_frce_split_bam;
nextflow run . --fastq $singledir  $config --ref_genome test_data/SIRV_150601a.fasta --minimap2_opts '-uf --splice-flank=no'\
--ref_annotation test_data/SIRV_isoforms.gtf -profile local --out_dir ${OUTPUT} -w ${OUTPUT}/workspace \
--bundle_min_reads 5 -resume;
r=$?
results+=("$(basename $OUTPUT): $r")

echo "Exit status codes for each test"
for value in "${results[@]}"; do
     echo "${value}"
done
