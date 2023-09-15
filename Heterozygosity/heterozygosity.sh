#!/bin/bash

REF=$1 #reference genome
SITES=$2 #sites file
CHR=$3 #list of chr
bamdir=$4 #bam directory

#run angsd with filters
cat $samples | while read -r LINE
do
angsd -i $bamdir/${LINE}.Warthog.bam -sites $SITES \
-rf $CHR -anc $REF -dosaf 1 -gl 2 -ref $REF -C 50 \
-minMapQ 30 -minQ 30 -nThreads 1 -out ${LINE}

#realSFS
$realSFS ${LINE}.saf.idx -P 10 > ${LINE}.ml
done

#Concatenate ML data
for i in `ls -1 *.ml`
do
sampleID=`basename $i .ml`
cat $i | awk -v sampleID=$sampleID '{print sampleID, $0}' >> $outml
done
