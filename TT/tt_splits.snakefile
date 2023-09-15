BCFTOOLS="/usr/bin/bcftools"
BGZIP="/usr/bin/bgzip"
TABIX="/usr/bin/tabix"
RSCRIPT="Rscript"

SCRIPTSDIR=config["scriptsdir"]
TT=os.path.join(SCRIPTSDIR, "tt.R")

OUTMAIN=config["outmain"]

rule all_tt:
    input:
        expand(os.path.join(OUTMAIN, "estimates", "tt", "{pair}.tt.{f}"),
               pair=config["do_tt"].keys(),
               f=["params.res", "split.years"]
        ),

rule select_polarized_sites:
    input:
        bcf = config["bcf"]
    output:
        sites = os.path.join(OUTMAIN, "sites", "ancestral_supported_sites.txt.bgz")
    params:
        outgroupsamples = ",".join(config["outgroups"])
    shell:
        """{BCFTOOLS} view -s {params.outgroupsamples} {input.bcf} | {BCFTOOLS} view  -e 'GT=="1/1" || GT=="0/1" || GT=="./."'  | {BCFTOOLS} query -f '%CHROM\t%POS\n' | {BGZIP} -c > {output.sites}"""

rule index_polarized_sites:
    input:
        sites = os.path.join(OUTMAIN, "sites", "ancestral_supported_sites.txt.bgz")
    output:
        index = os.path.join(OUTMAIN, "sites", "ancestral_supported_sites.txt.bgz.tbi")
    shell:
        """{TABIX} -b 2 -e 2 {input.sites}"""

rule get_genomewide_2dsfs:
    input:
        bcf = config["bcf"],
        sites = os.path.join(OUTMAIN, "sites", "ancestral_supported_sites.txt.bgz"),
        sites_index = os.path.join(OUTMAIN, "sites", "ancestral_supported_sites.txt.bgz.tbi"),
    params:
        s1 = lambda wildcards: config["do_tt"][wildcards.pair][0],
        s2 = lambda wildcards: config["do_tt"][wildcards.pair][1],
    output:
        sfs = os.path.join(OUTMAIN, "2dsfs", "{pair}_unfolded.2dsfs")
    shell:
        """{BCFTOOLS} view -s {params.s1},{params.s2} -T {input.sites} {input.bcf} | {BCFTOOLS} view  -e 'GT=="./."'  | {BCFTOOLS} query -f "[%GT ]\n" | awk '{{a[$1"_"$2]++}} END{{for(g in a){{print g, a[g]}}}}' > {output.sfs}"""

rule get_tt_estimates:
    input:
        sfs = os.path.join(OUTMAIN, "2dsfs", "{pair}_unfolded.2dsfs")
    output:
        os.path.join(OUTMAIN, "estimates", "tt", "{pair}.tt.params.res"),
        os.path.join(OUTMAIN, "estimates", "tt", "{pair}.tt.split.years"),
    params:
        outprefix = os.path.join(OUTMAIN, "estimates", "tt", "{pair}"),
        g = config["g"],
        mu = config["mu"]
    shell:
        """
        {RSCRIPT} {TT} {input.sfs} {params.outprefix} {params.mu} {params.g}
"""
