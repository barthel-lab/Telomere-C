#!/bin/bash
WD='results/align/UmiDeDup'
scripts/SES/SES_norm.sh ${WD}/VA13_capture_TCCAACGC.realn.mdup.MQ30.bam ${WD}/VA13_input_CTTGGTAT.realn.mdup.MQ30.bam VA13
scripts/SES/SES_norm.sh ${WD}/IMR90_caputre_CGTTAGAA.realn.mdup.MQ30.bam ${WD}/IMR90_input_TACCGAGG.realn.mdup.MQ30.bam IMR90
scripts/SES/SES_norm.sh ${WD}/WI38_capture_GATTCTGC.realn.mdup.MQ30.bam ${WD}/WI38_input_AGCCTCAT.realn.mdup.MQ30.bam WI38
scripts/SES/SES_norm.sh ${WD}/U2OS_capture_GTCGGAGC.realn.mdup.MQ30.bam ${WD}/U2OS_input_ACTAAGAT.realn.mdup.MQ30.bam U2OS
scripts/SES/SES_norm.sh ${WD}/HeLa_capture_TTACAGGA.realn.mdup.MQ30.bam ${WD}/HeLa_input_CCGTGAAG.realn.mdup.MQ30.bam HeLa
