#!/bin/bash
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8g
#SBATCH --time 1-00:00:00

# Usage: regions_vs_signal.sh [.run_signal.bw] [.run_peaks.bed] 

# Inputs
reads=$1
peak=$2
regions='/home/ychen/yc_data/annotation/CHM13v2.promoter_regions.CDS.bed'

# Outputs
mkdir -p results/align/promoter_analysis
outdir='results/align/promoter_analysis/'
region_out=$(echo ${outdir}$(basename ${peak} | sed 's/-capture.realn.mdup.MQ30.run_peaks/.hitted_promoters/g'))
matrix=$(echo ${outdir}$(basename ${reads} | sed 's/-capture.realn.mdup.MQ30.run_signal.bw/.matrix.txt/g'))

# Find peak hitted promoters
bedtools intersect \
-a ${peak} \
-b ${regions} -wb | \
awk 'BEGIN{OFS="\t"} ; {print $7,$8,$9,$10,$11,$12}' | sort | uniq \
> ${region_out}

# Make matrix
row_name_reads=$(basename ${reads} | sed 's/-capture.realn.mdup.MQ30.run_signal.bw/-Telomere-C/g')
row_name_hitpromoter=$(basename ${region_out})
cell=$(basename ${reads} | sed 's/-capture.realn.mdup.MQ30.run_signal.bw//g')
cell_factor=$(basename ${reads} | sed 's/.realn.mdup.MQ30.run_signal.bw//g')

echo -e "name\ttype\tfile\tfactor\tcell" > ${matrix}
echo -e "${row_name_reads}\treads\t${reads}\t${cell_factor}\t${cell}" >> ${matrix}
echo -e "${row_name_hitpromoter}\tregions\t${region_out}\tpromoters\t${cell}" >> ${matrix}

# Draw lineplot
ext=3000
bin=600
step=300

rgt-viz lineplot ${matrix} -o ${outdir}lineplot -t lineplot_${cell} -col reads -c regions -e ${ext} -center midpoint -add_region_number -bs ${bin} -ss ${step}
