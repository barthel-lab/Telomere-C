#!/bin/bash
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8g
#SBATCH --time 1-00:00:00
WD='results/align/UmiDeDup'
bam1=($(ls results/align/UmiDeDup/|grep input.*MQ30.bam$))
bam2=($(ls results/align/UmiDeDup/|grep cap.*MQ30.bam$))

for i in {0..5};
do outName=$(echo ${bam1[$i]}|sed 's/.input_[A-Z].*$//g');
plotEnrichment -b ${WD}/${bam1[$i]} ${WD}/${bam2[$i]} \
--BED data/repeatMasker.ITS.rmLowMap.bed \
data/TPE_candidate_gene.coor.all.bed \
--labels input TeloC \
--regionLabels "ITS" "TPE_genes" \
--ignoreDuplicates \
--blackListFileName data/T2T.excluderanges.noTelo.bed \
--variableScales \
--outRawCounts ${outName}.rawCounts.txt \
-o ${outName}.enrichment.png \
--numberOfProcessors 10;
done

