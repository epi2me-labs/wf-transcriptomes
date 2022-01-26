#!/usr/bin/env bash

# Usage: ./run_evaluation_dmel.sh pathto/outputdir

# See the isONcorrect paper https://www.nature.com/articles/s41467-020-20340-8 where this dataset is described


if [[ "$#" -lt 1 ]]; then
    echo "usage: run_evaluation_dmel.sh <outdir> [nextflow.config]"
    exit 1
fi

if [[ "$#" -eq 1 ]]; then
    config=''
fi

if [[ "$#" -eq 2 ]]; then
  config="-c $2";
fi

OUTDIR=$1;

FASTQ_URL="http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/005/ERR3588905/ERR3588905_1.fastq.gz"
REF_URL="http://ftp.ensembl.org/pub/release-99/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa.gz"
GFF_URL="http://ftp.ensembl.org/pub/release-99/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.28.99.gff3.gz"

DATA_DIR="$OUTDIR/data"
READS_DIR="$DATA_DIR/reads"
FASTQ="$READS_DIR/ERR3588905_1.fastq.gz"
REF="$DATA_DIR/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa"
GFF="$DATA_DIR/Drosophila_melanogaster.BDGP6.28.99.gff3"

mkdir -p $READS_DIR

if [ ! -f $REF ];
then (echo "downloading reference genome"; cd $DATA_DIR; curl -L -C - -O $REF_URL); gzip -d ${REF}.gz
fi

if [ ! -f $GFF ];
then
	(echo "downloading reference annotation"; cd $DATA_DIR; curl -L -C - -O $GFF_URL); gzip -d ${GFF}.gz
fi

if [ ! -f $FASTQ ];
then (echo "downloading reads"; cd $READS_DIR; curl -L -C - -O $FASTQ_URL); gzip -d ${FASTQ}.gz
fi


OUT_REF="$OUTDIR/ref"
OUT_DENOVO="$OUTDIR/denovo"


nextflow run ../ --fastq $READS_DIR  $config \
--ref_genome $REF --ref_annotation $GFF -profile local --out_dir $OUT_REF --minimap2_opts '-uf --splice-flank=no' \
-w $OUT_REF/workspace -resume;

echo "Doing de novo evaluation"
nextflow run ../ --fastq $READS_DIR  $config --denovo -profile local --out_dir $OUT_DENOVO \
-w $OUT_DENOVO/workspace -resume;
