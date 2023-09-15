#!/bin/bash

#Variables
EVALADMIX=$1 #evalAdmix path
bgl=$2 #beagle file
resdir=$3 #results dir
outdir=$4 #output directory
outname=$5 #output basename

for k in `ls -1 $resdir | grep qopt` #qopt k paths
do
    qfile="$resdir/$k"
    ffilegz="$resdir/`basename $qfile .qopt_conv`.fopt_conv.gz"
    ffile="`dirname $ffilegz`/`basename $ffilegz .gz`.txt"
    adx=`basename $ffile | sed 's/NGSadmix_minMaf005.//' | sed 's/.fopt_conv.txt//'`
    out=$outdir/${outname}_K${adx}
    gunzip -c $ffilegz > $ffile
    $EVALADMIX -beagle $bgl -qname $qfile -fname $ffile -P 20 -o $out.corres 2>&1 | tee -a $out.log
done
