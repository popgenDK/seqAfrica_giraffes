#!/usr/bin/env bash

crd=$(pwd -P)

plot_roh(){
    sam=$1
    Rscript ROH.v0.1.1.R --homfile ${sam}_maf5_geno5_ge0.5.hom -p ${sam}_maf5_geno5_ge0.5.bed -s ${sam}
}

export -f plot_roh
awk '{print $1}' /path/to/target/individual/list/xxx.fam | parallel -j 3 "plot_roh {}" :::: -
