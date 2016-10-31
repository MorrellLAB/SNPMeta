#!/usr/bin/env bash

set -e
set -u
set -o pipefail

#   Fill in values for the variables according to your needs.

#   Path to the executable
#   This will be something like blastn, or tblastx...
BLAST=
#   The blast database we will use
#   You can point it to a local copy of an NCBI database
DB=
#   The directory containing the SNPs to blast
TO_BLAST=
#   The number of hits to return total for each query
NUM_HITS=
#   The maximum e-value to return
EVALUE=
#   cd into the directory, so we can run a for-loop easily
cd ${TO_BLAST}

#   For each SNP...
for SNP in *
do
#   replace the BACK .fasta with .xml
#   This is a construct only present in recent versions of Bash
    OUT_BLAST=${SNP/%.fasta/.xml}
    #   The next few lines actually perform the search
    #   If you want to BLAST against NCBI's database, add the '-remote' switch
    $BLAST -query $SNP\
    -db $DB\
    -out ${OUT_BLAST}\
    -outfmt 5\
    -max_target_seqs ${NUM_HITS}\
    -evalue ${EVALUE}
done
