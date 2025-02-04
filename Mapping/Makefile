.PHONEY = adapters

# Needed to support `source`, since default shell (sh) doesn't and for advanced globs
SHELL := /bin/bash -O extglob -c

# Prefix that ensures that commands run with activated modules
ACTIVATED := source install/activate.sh && nice -19 ionice -c 3

# Override by running `make BATCH=batch_1`
BATCH := batch_1
# Find YAML files for this batch. Requires extglob to match ${BATCH}.mapping.yaml and ${BATCH}.mapping.blah.yaml
YAMLFILES := ${BATCH}.mapping.*yaml

ifeq ($(strip ${BATCH}), batch_4)
PERFECT_SAMPLE := GCamCaP20387
#else ifeq ($(strip ${BATCH}), batch_2)
#PERFECT_SAMPLE :=
else
PERFECT_SAMPLE := TODO
endif


help:
	@echo 'Usage: make <adapters/mapping/filtering> BATCH=$${batch}'
	@echo "  e.g. make mapping BATCH=batch_1"

adapters:
	# Support the case where output is a symlink to a removed folder
	mkdir -p $(shell readlink -f ${BATCH}.adapters)
	${ACTIVATED} python3 install/scripts/identify_adapters.py $(YAMLFILES) --output ${BATCH}.adapters
	${ACTIVATED} python3 install/scripts/classify_adapters.py ${BATCH}.adapters | tee ${BATCH}.adapters.tsv

mapping:
	# Support the case where output is a symlink to a removed folder
	mkdir -p $(shell readlink -f ${BATCH}.raw_bams)
	mkdir -p ${BATCH}.raw_bams/logs
	${ACTIVATED} paleomix bam run \
		--destination ${BATCH}.raw_bams/ \
		--temp-root ${BATCH}.raw_bams/temp \
		--log-file ${BATCH}.raw_bams/logs/${HOSTNAME}_$(shell date +"%Y%m%d_%H%M%S").txt \
		$(YAMLFILES)

filtering:
	${ACTIVATED} snakemake --snakefile install/scripts/Snakefile.filter_bams --config batch=${BATCH} -j 32 -p



#############################################################################
## 2. Preliminary QC of FASTQ files

${BATCH}.qc.prior.yaml: ${YAMLFILES}
	# Generate configuration file
	${ACTIVATED} python3 install/qc_workflow/scripts/make_config.py \
		--paleomix_yaml ${YAMLFILES} \
		--output_yaml "$@" \
		--outputdir_snk "${BATCH}.qc"

${BATCH}.qc.main.yaml: ${YAMLFILES}
	# Generate configuration file
	${ACTIVATED} python3 install/qc_workflow/scripts/make_config.py \
		--paleomix_yaml ${YAMLFILES} \
		--bams ${BATCH}/*.bam \
		--close_ref_name "RothschildsGiraffe" \
		--distant_ref_name "Okapi" \
		--perfect_bam_name "${PERFECT_SAMPLE}" \
		--output_yaml "$@" \
		--outputdir_snk "${BATCH}.qc"
	#	 --exclude A B C


QC_PIPELINE := ${ACTIVATED} snakemake --use-conda -p -c ${CORES} \
	--conda-prefix ./install/qc_workflow.conda/ \
	--snakefile ./install/qc_workflow/Snakefile \
	--rerun-incomplete --keep-going

qc_1_fastq: ${BATCH}.qc.prior.yaml ./${BATCH}.qc/
	${QC_PIPELINE} --configfile "${BATCH}.qc.prior.yaml" -- RUN_multiqc_pre

qc_2_mapstat: ${BATCH}.qc.main.yaml ./${BATCH}.qc/
	${QC_PIPELINE} --configfile "${BATCH}.qc.main.yaml" -- RUN_json

qc_3_all: ${BATCH}.qc.main.yaml ./${BATCH}.qc/
	${QC_PIPELINE} --configfile "${BATCH}.qc.main.yaml" -- RUN_all
