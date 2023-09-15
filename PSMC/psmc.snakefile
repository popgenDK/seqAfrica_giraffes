
import os
from pathlib import Path
# https://github.com/lh3/psmc

def add_ext(path, *args, sep=".", t=Path):
    p = str(path)
    for x in args:
        p = p + sep + x
    return t(p)

CURR_DIR = os.getcwd()

PSMC_PLOT = os.path.join(CURR_DIR, "plot_psmc.R")

VCFUTILS=config["VCFUTILS"] ## "/path/to/vcfutils.pl"
BCFTOOLS=config["BCFTOOLS"] ## "/path/to/bcftools"
PSMC_DIR=config["PSMC_DIR"] ## "/path/to/psmc"
FQ2PSMCFA=os.path.join(PSMC_DIR,"utils/fq2psmcfa")
PSMC=os.path.join(PSMC_DIR,"psmc")

FA=config["REF"] ## "/path/to/reference.fa"
MIN_MQ=config["MIN_MQ"]
MIN_BQ=config["MIN_BQ"] ## 30
ALLELE_SUP = config["ALLELE_SUPPORT"]
BED=config["GOODBED"]

BAMS = config["SAMPLES"]
DEPTHS = config["DEPTHS"]

OUTMAIN = Path(config["OUTDIR"])
OD_VCF = OUTMAIN / "vcf" / "{sample}"
OD_PSMCIN = OUTMAIN / "PSMC_IN" / "{sample}"
OD_PSMCOUT = OUTMAIN / "PSMC_OUT" / "{sample}"
OD_PLOT = OUTMAIN / "plot" / "res"

wildcard_constraints:
    sample = "|".join(BAMS.keys()),

rule all:
    input:
        expand(add_ext(OD_PSMCOUT, "psmc", t=str), sample=BAMS.keys()),
        add_ext(OD_PLOT, "png"),
        OUTMAIN / "heterozygosity.txt",


rule gen_bcftools_genome_wide:
    input:
        lambda wildcards: BAMS[wildcards.sample]
    output:
        temp(add_ext(OD_VCF, "bcf.gz"))
    threads: 2
    shell:
        "{BCFTOOLS} mpileup --threads {threads} --full-BAQ -T {BED} -Q {MIN_BQ}  -q {MIN_MQ} "
        " -O u --fasta-ref {FA} --per-sample-mF "
        "-a FORMAT/AD,FORMAT/DP {input} | "
        " {BCFTOOLS} call -Ob -o {output} --threads {threads} -c"

rule filter_bcf:
    input:
        add_ext(OD_VCF, "bcf.gz")
    output:
        add_ext(OD_VCF, "di_mono_allelic", "bcf.gz")
    params:
        mindepth=lambda wildcards: DEPTHS[wildcards.sample][0],  # 30/3
        maxdepth=lambda wildcards: DEPTHS[wildcards.sample][1],  # 30*2
    threads: 2
    shell: """
    {BCFTOOLS} view --threads {threads} -i 'strlen(REF)==1 & (strlen(ALT)==1 || ALT=".") &  FMT/DP>={params.mindepth} & FMT/DP<={params.maxdepth}' -M 2 -Ou {input} | \
    {BCFTOOLS} view --threads {threads} -i '(GT=="het" & FMT/AD[*:0]>={ALLELE_SUP} & FMT/AD[*:1]>={ALLELE_SUP} ) || GT=="hom"' -Ob -o {output}
    """

rule gen_fq_mindepthx:
    input:
        add_ext(OD_VCF, "di_mono_allelic", "bcf.gz")
    output:
        add_ext(OD_PSMCIN, "fq.gz")
    shell:
        "{BCFTOOLS} view  {input} | {VCFUTILS} vcf2fq | gzip > {output}"

rule gen_psmcfa:
    input:
        add_ext(OD_PSMCIN, "fq.gz")
    output:
        add_ext(OD_PSMCIN, "psmcfa")
    shell:
        """
        {FQ2PSMCFA} -q 20 {input} > {output}
        """

## https://github.com/lh3/psmc
rule run_psmc:
    input:
        add_ext(OD_PSMCIN, "psmcfa")
    output:
        add_ext(OD_PSMCOUT, "psmc")
    shell:
        """{PSMC} -N25 -t15 -r5 -p "4+25*2+4+6" -o {output} {input}"""

#Rscript ../plot_psmc.R psmc.files groups.txt 1e-8 9.5 test.png
rule psmc_plot_prep:
    input:
        expand(add_ext(OD_PSMCOUT, "psmc",t=str), sample=BAMS.keys())
    output:
        add_ext(OD_PLOT, "psmcs"),
        add_ext(OD_PLOT, "groups"),
    run:
        import os
        with open(output[0], 'w') as fh, open(output[1], 'w') as fh2:
            for psmc in input:
                print(psmc, file=fh)
                print(os.path.basename(psmc)[:7], file=fh2)

rule psmc_plot:
    input:
        rules.psmc_plot_prep.output
    output:
        add_ext(OD_PLOT, "png"),
    params:
        mut = 1e-8,
        gen = 7.5,
    shell:
        "Rscript {PSMC_PLOT} {input} {params.mut} {params.gen} {output}"
## BOOTSTRAPS will appear here, magically

rule AC_count:
    input:
        add_ext(OD_VCF, "di_mono_allelic", "bcf.gz")
    output:
        add_ext(OD_VCF, "di_mono_allelic", "AC")
    params:
        awk = '{a[$1]++} END {for (allele in a){print allele, a[allele]}}'
    shell:
        "{BCFTOOLS} query -f '%INFO/AC1\n' {input} | awk '{params.awk}' > {output}"

rule hetero:
    input:
        add_ext(OD_VCF, "di_mono_allelic", "AC")
    output:
        add_ext(OD_VCF, "di_mono_allelic", "het")
    run:
        total = 0
        with open(input[0], 'r') as fh:
            for line in fh:
                (ac, counts) = line.rstrip().split()
                ac = int(ac)
                counts = int(counts)
                total += counts
                if ac == 1:
                    het = counts
        with open(output[0], 'w') as fh:
            print(wildcards.sample, round(het / float(total), 6), file=fh)

rule merge_hetero:
    input:
        expand(add_ext(OD_VCF, "di_mono_allelic", "het", t=str), sample=BAMS.keys()),
    output:
        OUTMAIN / "heterozygosity.txt"
    shell:
        "cat {input} > {output}"


## parser to collect data for all samples.
