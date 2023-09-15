#!/bin/bash
file=$1 #plink filename without .bed
num=$2  #number of iterations
P=$3 #number of threads
out=$4 #output directory (no /)
K=$5 #number of populations
bfile=`basename $file`
echo num = $num
echo p = $P
echo out = $out
echo K = $K
if [ -f $file.bed ]
then
   echo file = $file
else
   echo File $file.bed does not exist.
fi

touch $out/$bfile.$K.likes.tmp
rm $out/$bfile.$K.likes.tmp

for f in `seq 1 $num`
        do
    echo -n -e $f"\t" >> $out/$bfile.$K.likes.tmp
    admixture $file.bed $K -s $f -j$P > $out/$bfile.$K.log_$f
    mv $bfile.$K.Q $out/$bfile.$K.Q_$f
    mv $bfile.$K.P $out/$bfile.$K.P_$f
    grep ^Loglikelihood $out/$bfile.$K.log_$f | cut -f2 -d" " >> $out/$bfile.$K.likes.tmp
    CONV=`Rscript -e "r<-read.table('$out/$bfile.$K.likes.tmp');r<-r[order(-r[,2]),];cat(sum(r[1,2]-r[,2]<5),'\n')"`

    if [ $CONV -gt 2 ]
    then
        cp $out/$bfile.$K.Q_$f $out/$bfile.$K.Q_conv
        cp $out/$bfile.$K.P_$f $out/$bfile.$K.P_conv
        cp $out/$bfile.$K.log_$f $out/$bfile.$K.log_conv
        break
        fi
done
cat $out/$bfile.$K.likes.tmp | sort -k2 -n -r > $out/$bfile.$K.likes
