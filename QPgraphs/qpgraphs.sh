#!/bin/bash

# Filters for polymorphic sites in ingroup and non-polymorphic sites in outgroup
# Allows for a --geno filter for within the ingroup (max missingness)

# Software prerequisites
# PLINK

#########
Input variables - EXAMPLES
inputplink="/path/to/"
inputinds="/path/to/individuals.txt" # individual names have to be in first column
outgroupkey="OJoh" # Should match sample names of outgroup (only one outgroup for now)
genofilter="--geno 0.1" #decimal e.g. 0.1 means <10% missingness in ingroup
outputplink="/path/to/`basename $inputplink`.outfile" #Final output plink filename

###################################
#Make tmp and logs directories
#outputdir="`dirname $outputplink`"
mkdir -p $outputdir/tmp
mkdir -p $outputdir/polyFilter_logs

echo ${outputdir}

###################################
#Filter for individuals, including outgroups.
cd ${outputdir}
#inds="${outputdir}/tmp/`basename $inputinds .txt`.indfilters.txt"
cut -f1 ${inputinds} | awk '{print $1,$1}' OFS="\t" > ${inds}

#Filter for individuals we want (i.e. MAF>0.0001 and geno 0)
#outplinkInd="$outputdir/tmp/`basename $inputplink`.indfilters.tmp1"

#PLINK FILE OUT 1

plink \
    --make-bed \
    --bfile ${inputplink} \
    --keep ${inds} \
    --maf 0.000000000000000001 \
    --out ${outplinkInd} \
    --allow-extra-chr \
    --set-missing-var-ids @:#:\$1:\$2 \
    $genofilter

mv $outputdir/tmp/*.log  $outputdir/polyFilter_logs 

###################################
#Keep only polymorphic sites in species i.e. MAF>0.
#outplinkpolyIngroup="$outplinkInd.ingroupfilters.tmp2a"
#ingroupIDs="$outputdir/tmp/`basename $inds .txt`.noOutgroups.txt"

#Subset for ingroup individuals only
cat $inds | grep -v $outgroupkey > $ingroupIDs

plink \
    --make-just-bim \
    --bfile $outplinkInd \
    --keep $ingroupIDs \
    --maf 0.000000000000000001 \
    --allow-extra-chr \
    --out $outplinkpolyIngroup \
    --set-missing-var-ids @:#:\$1:\$2

mv $outputdir/tmp/*.log  $outputdir/polyFilter_logs 

################################
####Extract sites that are polymorphic in ingroup
cut -f2 $outplinkpolyIngroup.bim > $inputInPlink.txt

#PLINK FILE OUT 2

plink \
    --make-bed \
    --bfile $outplinkInd \
    --extract $inputInPlink.txt \
    --allow-extra-chr \
    --out $outplinkpolyIngroup2 \
    --set-missing-var-ids @:#:\$1:\$2

mv $outputdir/tmp/*.log  $outputdir/polyFilter_logs 

#Separately filter for non-polymorphic sites in outgroup - MAF=0 (max-maf <0.0000000001).
#outplinkIndNonpolyout="$outplinkInd.outgroupfilters.tmp2b"
#outgroupIDs="$outputdir/tmp/`basename $inds .txt`.Outgroups.txt"

#Subset for outgroup
cat $inds | grep $outgroupkey > $outgroupIDs

plink \
    --make-just-bim \
    --bfile $outplinkpolyIngroup2 \
    --keep $outgroupIDs \
    --max-maf 0.000000000000000001 \
    --allow-extra-chr \
    --out $outplinkIndNonpolyout \
    --set-missing-var-ids @:#:\$1:\$2

mv $outputdir/tmp/*.log  $outputdir/polyFilter_logs 

###################################
#Generate new plink file containing filtered sites from original plink file 
#inputOutPlink="$outputdir/tmp/`basename $outputplink`.Input.tmp3"

#Generate new bim file containing intersection of filtered files 
cut -f2 $outplinkIndNonpolyout.bim > $inputOutPlink.txt

#PLINK FILE OUT 3
plink \
    --make-bed \
    --bfile ${inputplink} \
    --extract $inputOutPlink.txt \
    --out $outputplink \
    --allow-extra-chr
 
#rm -rf $outputdir/tmp
######FIN
