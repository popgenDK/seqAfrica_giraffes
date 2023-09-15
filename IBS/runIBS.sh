#!/bin/bash

REF=$1 #reference genome
SITES=$2 # sites file
CHR=$3 #list of chr
OUT=$4 ##Final results dir
BAMS='metadata/67allInd_bamList.txt'

angsd -sites $SITES -rf $CHR -GL 2 -minMapQ 30 -minQ 30 -out $OUT -nthreads 20 \
-doGlf 2 -doMajorMinor 1 -skipTriallelic 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 -bam $BAMS \
-doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1  2>&1 | tee -a $OUT.angsd.log
