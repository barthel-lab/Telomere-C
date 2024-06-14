#!/bin/bash
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8g
#SBATCH --time 1-00:00:00

# Star from SES normalized bw data, which save our time

# Create compare bw
WD="results/align/SES/*"
ips=($(ls results/align/SES/*/|grep cap.*MQ30.SES.bw$))
input=($(ls results/align/SES/*/|grep input.*MQ30.SES.bw$))

for i in {0..5}
do echo "Start bigwigCompare..."; 
name=$(echo ${input[$i]}|sed 's/.input_[A-Z].*$//g')
bigwigCompare -b1 ${WD}/${ips[$i]} \
-b2 ${WD}/${input[$i]} \
--numberOfProcessors 10 \
--binSize 100 \
--verbose \
--outFileFormat bigwig \
--outFileName $(echo "${name}.log2ratio.bw");
done

new_bw=$(ls|grep log2ratio.bw)

# Generate intermdeiate matrix file from bw
computeMatrix scale-regions  -S $(echo $new_bw) \
-R data/repeatMasker.ITS.rmLowMap.bed \
-b 3000 -a 3000 \
--missingDataAsZero \
-o ITS.matrix.mat.gz

# Plot heatmap and profile
plotHeatmap -m ITS.matrix.mat.gz \
--startLabel "Start" \
--endLabel "End" \
-out ITS.TSS.heatmap.png

