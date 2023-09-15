#!/usr/bin/env bash

# module load plink/1.9.0

ind_plink(){
    sam=$1
    if [ -d ${sam}_plink ]; then
        /bin/rm -r ${sam}_plink
    fi
    mkdir ${sam}_plink
    grep -w ${sam} ${target_ind} > ${sam}_plink/${sam}.list
    plink --keep ${sam}_plink/${sam}.list --allow-extra-chr --bfile filtered.maf5.geno5 --exclude filtered.maf5.geno5.hwe.OHE_ge0.5.list --keep-allele-order --make-bed --homozyg --homozyg-window-het 3 --homozyg-window-missing 20 --out ${sam}_plink/${sam}_maf5_geno5_ge0.5 &>${sam}_plink/${sam}_maf5_geno5_ge0.5.plog
}

export target_ind=$1 #/path/to/target/individual/list/xxx.fam
export -f ind_plink

awk '{print $1}' ${target_ind} | parallel -j 3 "ind_plink {}" :::: -


sam=$(head -n 1 /path/to/target/individual/list/xxx.fam)
head -n 1 ${sam}_plink/${sam}_maf5_geno5_ge0.5.hom > all.hom; cat */*.hom | grep -v "DENSITY" >>all.hom
