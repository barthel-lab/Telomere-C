#!/bin/bash
WD='results/align/UmiDeDup'
cells=(BJ IMR90 WI38 U2OS HeLa VA13)
scale=(0.6085872 0.7008441 0.7201158 0.7081338 0.7002526 0.7027948)
mkdir -p results/align/dTelo
OUT='results/align/dTelo/'

for i in {0..5};
do input=$(ls $WD/${cells[$i]}*|grep MQ30|grep input|grep -v bai);
capture=$(ls $WD/${cells[$i]}*|grep MQ30|grep -v input|grep -v bai);
bamCoverage \
        -b ${input}   \
        -o ${OUT}/${cells[$i]}_input.SES.bin.100K.bed  \
        -of bedgraph \
        --binSize 100000 \
        --numberOfProcessors 10 \
        --ignoreDuplicates \
        --scaleFactor ${scale[$i]} \
        --normalizeUsing None;

bamCoverage \
        -b ${capture}   \
        -o ${OUT}/${cells[$i]}_capture.SES.bin.100K.bed  \
        -of bedgraph \
        --binSize 100000 \
        --numberOfProcessors 10 \
        --ignoreDuplicates \
        --scaleFactor 1 \
        --normalizeUsing None;

bedtools intersect \
-a /home/ychen/yc_data/CHM13v2_ref_by_chr/chm13.v2.100K.bed \
-b ${OUT}/${cells[$i]}_input.SES.bin.100K.bed -wa -wb|awk '{print $1 "\t" $2 "\t" $3 "\t" $7}' > ${OUT}/${cells[$i]}_input.SES.bin.100000.window.bed;

bedtools intersect \
-a /home/ychen/yc_data/CHM13v2_ref_by_chr/chm13.v2.100K.bed \
-b ${OUT}/${cells[$i]}_capture.SES.bin.100K.bed -wa -wb|awk '{print $1 "\t" $2 "\t" $3 "\t" $7}' > ${OUT}/${cells[$i]}_capture.SES.bin.100000.window.bed;

done

rm *.SES.bin.100K.bed
