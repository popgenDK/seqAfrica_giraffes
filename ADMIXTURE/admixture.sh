#!/bin/bash
file=$1 # PLINK prefix filename
num=$2  # Number of ADMIXTURE runs per K
threads=$3 # Number of threads
out=$4 # Output directory
K=$5 # Number of ancestral populations
base=$(basename $file)

# File to store loglikelihoods
touch ${out}/${base}.K$K.loglikes

# Run ADMIXTURE for num times
for s in `seq 1 $num`
do
		admixture ${file}.bed $K -s $s -j${threads} > ${out}/${base}.K${K}.s${s}.log
		mv ${base}.${K}.Q ${out}/${base}.K${K}.s${s}.Q
		mv ${base}.${K}.P ${out}/${base}.K${K}.s${s}.P
		grep ^Loglikelihood ${out}/${base}.K${K}.s${s}.log | cut -f2 -d" " >> ${out}/${base}.K$K.loglikes
done

# Find best solution
best=$(nl ${out}/${base}.K$K.loglikes | sort -k2 -n -r | head -n 1 | cut -f1 | tr -d ' ')
cp ${out}/${base}.K${K}.s${best}.Q ${out}/${base}.K${K}.s${best}.best.Q
cp ${out}/${base}.K${K}.s${best}.P ${out}/${base}.K${K}.s${best}.best.P
