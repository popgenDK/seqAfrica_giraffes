#!/bin/bash
FILE=$1 # Path to bed file
K=$2 # Number of ancestral populations
N=$3 # Number of runs (different seeds)
OUT=$4 # Output prefix
THREADS=$5 # Number of threads to use
TMP=`basename $FILE`
DEF=${TMP::-4}

for i in `seq 1 $N`
do
        admixture $FILE $K -j$THREADS -s $i > $OUT.K$K.s$i.log
        mv $DEF.$K.Q $OUT.K$K.s$i.Q
        mv $DEF.$K.P $OUT.K$K.s$i.P
        gzip $OUT.K$K.s$i.P
done
