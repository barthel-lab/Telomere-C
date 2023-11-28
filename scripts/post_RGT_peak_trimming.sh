#!/bin/bash
# Patch the bug from RGT peakcaller by removing chrM and corrdinate that start > end
# Usage:
# peak_trimming.sh <untrimed_peak.bed>

IN=$1
OUT=$(echo $IN | sed s/.bed/.trim.bed/g)

awk '{if ($2<$3 && $1 != "chrM") print $0}' $IN > $OUT


