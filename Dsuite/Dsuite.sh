#!/bin/bash

# Define the paths to Dsuite and Bcftools
Dsuite="/path/to/Dsuite"
Bcftools="/path/to/bcftools"

# Replace these paths with your actual input VCF and output file names
input_vcf="/path/to/input.vcf.gz"
noSNPs=12345  # Replace 12345 with the actual number of SNPs in the vcf

# Run Bcftools and Dsuite Dtrios
$Bcftools view "$input_vcf" | $Dsuite Dtrios -l "$noSNPs" stdin SETS.txt -t tree.nwk

# Run Dsuite Fbranch
$Dsuite Fbranch tree.nwk SETS.txt > Fbranch_out.txt

# Replace these paths with your actual python script and arguments
python_script="/path/to/python3" 
dtools_script="/path/to/dtools.py"

# Run the Python script
$python_script "$dtools_script" Fbranch_out.txt tree.nwk  --tree-label-size  8
