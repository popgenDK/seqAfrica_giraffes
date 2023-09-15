#!/bin/bash

file=$1 #input beagle file
nfile=$2 #file name 
num=$3  #maximum number of iterations
P=$4 #number of threads/cores used
K=$5 #number of populations
NGSA=$6 - #NGSadmix path

out=$out/$K #output directory

mkdir -p $out

for f in `seq $num`
do
    echo -n -e $f"\t"
    echo $file
    echo $K
    #Run NGSadmix
    $NGSA -minMaf 0.05 -likes $file -seed $f -K $K -P $P -printInfo 1 -o $out/$nfile.$K.$f
    grep "like=" $out/$nfile.$K.$f.log | cut -f2 -d " " | cut -f2 -d "=" >> $out/$nfile.$K.likes
    #Convergence criteria - #r = 10* no. of millions of SNPs
    CONV=`Rscript -e "r<-read.table('$out/$nfile.$K.likes');r<-r[order(-r[,1]),];cat(sum(r[1]-r<148),'\n')"` 

    #Check convergence
    if [ $CONV -gt 2 ]  # 3 files meeting convergence criteria
    then
	cp $out/$nfile.$K.$f.qopt $out/$nfile.$K.$f.qopt_conv
	cp $out/$nfile.$K.$f.fopt.gz $out/$nfile.$K.$f.fopt_conv.gz
	cp $out/$nfile.$K.$f.log $out/$nfile.$K.$f.log_conv
	break
    fi

done
