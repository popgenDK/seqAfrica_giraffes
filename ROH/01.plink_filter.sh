#!/usr/bin/env bash

# module load plink/1.9.0

pre_bfile=$1 #/path/to/preprocessed/bfile
target_ind=$2 #/path/to/target/individual/list/xxx.fam

plink --allow-extra-chr --bfile ${pre_bfile} --keep ${target_ind} --keep-allele-order --chr-set 16 --maf 0.05 --geno 0.05 --het --hardy --ibc --missing --make-bed --set-missing-var-ids "@:#:\$1:\$2" --out filtered.maf5.geno5

awk '$7>=0.5' filtered.maf5.geno5.hwe |awk '{print $2}' |sed '1d' > filtered.maf5.geno5.hwe.OHE_ge0.5.list
