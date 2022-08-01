#!/bin/bash
set -e

git clone https://github.com/Oshlack/JAFFA.git &&
cd JAFFA
git checkout 24b1c3b

cp ../subworkflows/JAFFAL/install_linux64.sh .
./install_linux64.sh

# JAFFA uses tools.groovy to locate binaries. We modify it as we only need a small subset or they are already
# included in our env

# Tools to be compiled from src/
bin=$(realpath tools/bin)
echo "//Tools built locally" >> tools.groovy
declare -a tools=("reformat" "extract_seq_from_fasta" "make_simple_read_table" "process_transcriptome_align_table" "make_3_gene_fusion_table")
for b in "${tools[@]}"; do
  echo "$b=\"$bin/$b\"" >> tools.groovy
done
echo "minimap2=\"minimap2\"" >> tools.groovy
