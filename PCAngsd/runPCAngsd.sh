#!/bin/bash

# https://github.com/Rosemeis/pcangsd
# http://www.popgen.dk/software/index.php/PCAngsd

PC=$1 #PCAngsd path
SITES=$2 #sites file
BAMNAME=$3 #bamname
DIRNAME=$4 #dirname
resultsdir=$5 #results dir
BEAGLE=$6 #beagle file out
out=$7 #out
log=$8 #log
BAMS='metadata/54unrelatedInd_bamList.txt'

# make beagle file input for PCAngsd
angsd -sites $SITES -GL 2 -minMapQ 30 -minQ 30 -out $BEAGLE -nthreads 4 \
-doGlf 2 -doMajorMinor 1 -skipTriallelic 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 \
-bam $BAMS 2>&1 | tee -a $DIRNAME/logs/$BAMNAME.log

#-------
#Plot PCA
$PC --beagle $BEAGLE --threads 8 -o $out 2>&1 > $log #PC3
$PC -e12 --beagle $BEAGLE --threads 8 -o $out.e12 2>&1 > $log.e12.log #npops-1
