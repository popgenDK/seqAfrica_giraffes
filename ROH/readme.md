# ROH (Runs of Homozygosity)

## Introduction

The pipeline for ROH analysis of individuals in giraffes project.

## Requirement

- GNU parallel
- plink == 1.9.0
- R >= 4.0
- R library:
  - snpStats
  - windowscanr
  - tidyverse
  - RColorBrewer
  - argparse
  - ggplot2
  - data.table
  - ggthemes

## Steps

Before running the pipeline, you may want to preprocess your SNP data, like quality control or data filtering.

Please generate the input files for this pipeline using PLINK (ver 1.9.0, with `--make-bed`). Output should be files like `/path/to/preprocessed/bfile.bed` (and other suffixes).

Please provide a list of individuals which you're interested in. It should be in `.fam` (PLINK) format, like `/path/to/target/individual/list/xxx.fam`. Make sure the **FID and IID should be the same**.

Then please run ROH by following steps:

**(NOTE: please check the plink parameters in `.sh` scripts, or modify plotting parameters in R scripts if they don't fit your data well.)**

```bash
# step 01: filter out invalid sites
bash 01.plink_filter.sh /path/to/preprocessed/bfile /path/to/target/individual/list/xxx.fam
```

```bash
# step 02: ROH for each individual
bash 02.individual_plink.sh /path/to/target/individual/list/xxx.fam
```

```bash
# step 03: draw ROH region plot across each chromosomes of each individual
bash 03.ROH_plots.sh
```

```bash
# step 04: draw ROH proportion plot of all individuals
bash 04.ROH_proportion.sh
```

