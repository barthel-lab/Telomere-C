#!/bin/bash
WORKDIR="~/Telomere-C"

SMKDIR="~/Telomere-C"
export SMKDIR

CLUSTRCONF="conf/cluster.json"
CONFIGFILE="conf/config.yaml"
#snakemake --jobs 500 -k --latency-wait 120 --max-jobs-per-second 2 --restart-times 0 --cluster-config cluster.json --configfile config.yaml --jobname "{jobid}.{cluster.name}" --drmaa " -S /bin/bash -j {cluster.j} -M {cluster.M} -m {cluster.m} -q {cluster.queu} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem}gb -e ${WORKDIR}/{cluster.stderr} -o ${WORKDIR}/{cluster.stdout}" all
snakemake --profile sumner -s Snakefile --verbose --cluster-config "${CLUSTRCONF}" --configfile "${CONFIGFILE}" --stats "${WORKDIR}"/logs/snakemake_stats.json all |& tee -a "${WORKDIR}"/logs/run_snakemake_"$TSTAMP".log
## END ##
