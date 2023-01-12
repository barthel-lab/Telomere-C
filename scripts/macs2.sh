#!/bin/bash
# usage ./macs2.sh <lib name> <path/to/treated.MQ30.bam> <path/to/input.MQ30.bam>
Sname=$1
treat=$2
input=$3

macs2 callpeak \
-t $treat \
-c $input \
--outdir results/align/macs2/${Sname}/ \
-g hs \
-n ${Sname} \
-B \
-q 0.01

