#!/bin/bash

# Usage SES_norm.sh <TeloC.bam> <input.bam> <lib name>

window="data/chm13v2.0.w1000.bed"
ip_bam=$1
input_bam=$2
lib_name=$3
bin=100
# Intermediate file
ip_bed=$(basename ${ip_bam}|sed 's/.bam/.count.bed/g')
input_bed=$(basename ${input_bam}|sed 's/.bam/.count.bed/g')

# Output file
out=$(basename ${input_bam}|sed 's/.bam/.SES.bw/g')
out_ip=$(basename ${ip_bam}|sed 's/.bam/.SES.bw/g')
out_dir="results/align/SES"

#  Checking if window file exist
if [ ! -f "${window}" ]; then
    echo "${window} no found"
    echo "creating..."
    bedtools makewindows -g /labs/barthel/references/CHM13v2/chm13v2.0.size.genome -w 1000 > ${window}
fi

# Checking read counts for capture
if [ ! -f  "${ip_bed}" ]; then
    echo "Count reads in each window for capture"
    samtools view -b -F 1024 ${ip_bam} | \
    bedtools intersect -a ${window} -b stdin -c > ${ip_bed}
fi

# Checking read counts for input
if [ ! -f "${input_bed}" ]; then
    echo "Count reads in each window for input"
    samtools view -b -F 1024 ${input_bam} | \
    bedtools intersect -a ${window} -b stdin -c > ${input_bed}
fi

echo "Calculating scaling factor"
echo "Be patient..."
Rscript scripts/SES/SES_alpha.r ${ip_bed} ${input_bed}

# Apply scale factor to correct input signal
scale_factor=$(cat scale_factor.txt)
echo "Scaling factor is $(cat scale_factor.txt)"

echo "Normalyzing input data, and export in bigwig"
echo "Be patient..."
mkdir -p ${out_dir}/${lib_name}
bamCoverage \
        -b ${input_bam}  \
        -o ${out_dir}/${lib_name}/${out} \
        -of bigwig \
        --binSize ${bin} \
        --numberOfProcessors 10 \
        --ignoreDuplicates \
        --scaleFactor ${scale_factor} \
        --normalizeUsing None

bamCoverage \
        -b ${ip_bam}  \
        -o ${out_dir}/${lib_name}/${out_ip} \
        -of bigwig \
        --binSize ${bin} \
        --numberOfProcessors 10 \
        --ignoreDuplicates \
        --scaleFactor 1 \
        --normalizeUsing None

mv scale_factor.txt ${out_dir}/${lib_name}
mv Signal_extraction_scaling.tiff ${out_dir}/${lib_name}
rm ${ip_bed} ${input_bed}
echo 'done'
