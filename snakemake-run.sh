#!/bin/bash
WORKDIR="/projects/barthf/Telomere-C/Telo-C"
snakemake --jobs 500 -k --latency-wait 120 --max-jobs-per-second 2 --restart-times 0 --cluster-config cluster.json --configfile config.yaml --jobname "{jobid}.{cluster.name}" --drmaa " -S /bin/bash -j {cluster.j} -M {cluster.M} -m {cluster.m} -q {cluster.queu} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime} -l mem={cluster.mem}gb -e ${WORKDIR}/{cluster.stderr} -o ${WORKDIR}/{cluster.stdout}" all
## END ##
