#!/bin/bash
WORKDIR=$PWD

SMKDIR=$PWD
export SMKDIR

CLUSTRCONF="conf/cluster.json"
CONFIGFILE="conf/config.yaml"

# This command use --profile argument

snakemake --jobs 500 -k --profile slurm_profile -s Snakefile --verbose --cluster-config "${CLUSTRCONF}" --configfile "${CONFIGFILE}" --stats "${WORKDIR}"/logs/snakemake_stats.json all |& tee -a "${WORKDIR}"/logs/run_snakemake_"$TSTAMP".log

#snakemake --jobs 500 -k --profile slurm_profile -s Snakefile --verbose --cluster-config "${CLUSTRCONF}" --configfile "${CONFIGFILE}" --stats "${WORKDIR}"/logs/snakemake_stats.json all |& tee -a "${WORKDIR}"/logs/run_snakemake_"$TSTAMP".log

# This command DON'T use --profile argument
#snakemake --jobs 500 -k --latency-wait 120 --max-jobs-per-second 2 --restart-times 0 --cluster-config "${CLUSTRCONF}" --configfile "${CONFIGFILE}" --jobname "{jobid}.{cluster.name}" all

## END ##
