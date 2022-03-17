# This is for AluI cutted lib
# Inputs
# ${Sname}.realn.mdup.MQ30.bam
# AluI_RE.bed
# Usage: make_fragment.sh [lib name(str)] 
# fragment size default is 15K

Sname=$1

RE_bed='data/AluI_RE.bed'
IN="results/align/macs2/${Sname}/${Sname}.realn.mdup.MQ30.bam"
OUT_DIR="results/align/vFragments/${Sname}/"
OUT_stat="${OUT_DIR}/${Sname}.Telo-RE.stat"
fSize=15000

mkdir -p ${OUT_DIR}

# L-end: reverse reads, -f 16, end with AG & intersect with AluI RE sites 
samtools view -b -f 16 -r 'Telo-unknown' -r 'unknown-Telo' ${IN}|bedtools intersect -a stdin -b ${RE_bed} -wa|samtools view|awk '$10 ~ /AG$/ {print $1}'|sort|uniq > "${OUT_DIR}/${Sname}.realn.mdup.MQ30.L-end.readNames"
samtools view -b -N ${OUT_DIR}/${Sname}.realn.mdup.MQ30.L-end.readNames -o ${OUT_DIR}/${Sname}.realn.mdup.MQ30.L-end.bam ${IN}

# R-end: forward reads, -F 16, start with CT & intersect with AluI RE sites
samtools view -b -F 16 -r 'Telo-unknown' -r 'unknown-Telo' ${IN}|bedtools intersect -a stdin -b ${RE_bed} -wa|samtools view|awk '$10 ~ /^CT/ {print $1}'|sort|uniq > "${OUT_DIR}/${Sname}.realn.mdup.MQ30.R-end.readNames"
samtools view -b -N ${OUT_DIR}/${Sname}.realn.mdup.MQ30.R-end.readNames -o ${OUT_DIR}/${Sname}.realn.mdup.MQ30.R-end.bam ${IN}

# Extend L-end forward read by max fragment size and write into bed file
samtools view -b -F 1040 ${OUT_DIR}/${Sname}.realn.mdup.MQ30.L-end.bam|bedtools bamtobed -i stdin|awk -v size=${fSize} '{print $1 "\t" $2 "\t" $2+size}' > "${OUT_DIR}/${Sname}.realn.mdup.MQ30.L-end.Fext.bed"

# Get R-end reverse reads and write as bed file
samtools view -b -f 16 -F 1024 ${OUT_DIR}/${Sname}.realn.mdup.MQ30.R-end.bam|bedtools bamtobed -i stdin|awk '{print $1 "\t" $2 "\t" $3}' > "${OUT_DIR}/${Sname}.realn.mdup.MQ30.R-end.Rev.bed"

# Fragment assembling. start: L-end start coordinate, end: R-end end coordinate
bedtools intersect -a  ${OUT_DIR}/${Sname}.realn.mdup.MQ30.L-end.Fext.bed -b  ${OUT_DIR}/${Sname}.realn.mdup.MQ30.R-end.Rev.bed -wa -wb|awk '$6-$2>300 {print $1 "\t" $2 "\t" $6}' > "${OUT_DIR}/${Sname}.realn.mdup.MQ30.fragments.bed"


# Generate metrics
MQ30=0
MQ30_TelTel=0
MQ30_TelRE=0
MQ30_RERE=0
MQ30_other=0
frag=0


MQ30=$(samtools view ${IN} |wc -l)
MQ30_TelTel=$(samtools view -r 'Telo-Telo' ${IN} |wc -l)
MQ30_TelTel_cent=$(echo "scale=2;$MQ30_TelTel/$MQ30"|bc)

MQ30_TelRE=$((\
$(samtools view ${OUT_DIR}/${Sname}.realn.mdup.MQ30.L-end.bam|wc -l)+\
$(samtools view ${OUT_DIR}/${Sname}.realn.mdup.MQ30.R-end.bam|wc -l)\
))
MQ30_TelRE_cent=$(echo "scale=2;$MQ30_TelRE/$MQ30"|bc)

samtools view -b -F 16 ${IN}|bedtools intersect -a stdin -b ${RE_bed} -wa|samtools view|awk '$10 ~ /^CT/ {print $1}'|sort|uniq > "${OUT_DIR}/RE-X.readNames"
samtools view -b -f 16 -N ${OUT_DIR}/RE-X.readNames ${IN}|bedtools intersect -a stdin -b ${RE_bed} -wa|samtools view|awk '$10 ~ /AG$/ {print $1}'|sort|uniq > "${OUT_DIR}/RE-RE.readNames"
MQ30_RERE=$(samtools view -N ${OUT_DIR}/RE-RE.readNames ${IN}|wc -l)
MQ30_RERE_cent=$(echo "scale=2;$MQ30_RERE/$MQ30"|bc)

MQ30_other=$((\
$MQ30 - $MQ30_TelTel - $MQ30_TelRE - $MQ30_RERE\
))
MQ30_other_cent=$(echo "scale=2;$MQ30_other/$MQ30"|bc)
frag=$(wc -l <  ${OUT_DIR}/${Sname}.realn.mdup.MQ30.fragments.bed)


echo $"\
Total ${Sname}.realn.mdup.MQ30.bam reads: $MQ30
Tel + Tel: $MQ30_TelTel $MQ30_TelTel_cent
Tel + RE: $MQ30_TelRE $MQ30_TelRE_cent
RE + RE: $MQ30_RERE $MQ30_RERE_cent
Ohter: $MQ30_other $MQ30_other_cent
Fragments: $frag" \
> ${OUT_stat}


