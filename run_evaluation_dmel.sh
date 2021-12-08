#!/bin/bash

# Usage: ./run_evaluation_dmel.sh pathto/outputdir

# See the isONcorrect paper https://www.nature.com/articles/s41467-020-20340-8 where this dataset is described

OUTDIR=$1;

FASTQ_URL="http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/005/ERR3588905/ERR3588905_1.fastq.gz"
REF_URL="http://ftp.ensembl.org/pub/release-99/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa.gz"
GFF_URL="http://ftp.ensembl.org/pub/release-99/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.28.99.gff3.gz"

RESULTS_DIR="$OUTDIR/results"
DATA_DIR="$OUTDIR/data"
READS_DIR="$DATA_DIR/reads"
FASTQ="$READS_DIR/ERR3588905_1.fastq"
REF="$DATA_DIR/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa"
GFF="$DATA_DIR/Drosophila_melanogaster.BDGP6.28.99.gff3"

echo $READS_DIR;

rm -fr $OUT_DIR/results
mkdir -p $OUTDIR/data

if [ ! -f $REF ];
then (cd $DATA_DIR; curl -L -C - -O $REF_URL); gzip -d ${REF}.gz
fi

if [ ! -f $GFF ]
then
	(cd $DATA_DIR; curl -L -C - -O $GFF_URL); gzip -d ${GFF}.gz
fi

if [ ! -f $FASTQ ];
then (cd $READS_DIR; curl -L -C - -O $FASTQ_URL); gzip -d ${FASTQ}.gz
fi

nextflow run ../wf-isoforms --fastq $READS_DIR  \
--reference_genome $REF --annotation $GFF -profile conda --out_dir $OUTDIR \
-w $OUTDIR/workspace -resume
