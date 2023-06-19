# Usage: python RGT_peakCaller.py capture.bam input.bam
# Reference: https://reg-gen.readthedocs.io/en/latest/rgt/tutorial-peak-calling.html
from rgt.Util import GenomeData
from rgt.helper import get_chrom_sizes_as_genomicregionset
from rgt.THOR.get_extension_size import get_extension_size
from rgt.CoverageSet import CoverageSet
import sys
import os
# 0. Configure
# I/O 
bamfile = sys.argv[1]
bamfile_input = sys.argv[2]
outdir = 'result/align/RGT_peakCall/'
outtxt = outdir + os.path.basename(bamfile).replace('.bam','.run_peaks.rs20.bed')
outbw =  outdir + os.path.basename(bamfile).replace('.bam','.run_signal.rs20.bw')

# Peak calling filters
minRd = 20
pvalue = 0.05

# You must complete the Configuration of Genomic Data in your first time of running
# Please check: https://reg-gen.readthedocs.io/en/latest/rgt/setup_data.html
g = GenomeData('CHM13v2')
regionset = get_chrom_sizes_as_genomicregionset(g.get_chromosome_sizes())

# 1. Preprocessing Steps
# Compute the fragment size of the reads
ext, _ = get_extension_size(bamfile, start=0, end=300, stepsize=5)

# Create an instance of CoverageSet
cov = CoverageSet('IP coverage', regionset)
cov_input = CoverageSet('input-dna coverage', regionset)

cov.coverage_from_bam(bam_file=bamfile, extension_size=ext)
cov_input.coverage_from_bam(bam_file=bamfile_input, extension_size=0)

# Normalize against the GC-content followed by subtracting the input-DNA data from the IP data
cov.norm_gc_content(cov_input.coverage, g.get_genome(), g.get_chromosome_sizes())
cov.subtract(cov_input)

# 2. Peak Calling (Binomial distribution)
from numpy import sum, mean
s = sum(cov.overall_cov)
p = mean(cov.overall_cov[cov.overall_cov > 0]) / s

# Ignore bins with low number of reads
from scipy.stats import binom_test
def filter_peaks(c, empirical_s, empirical_p, pvalue_theshold=pvalue, min_reads=minRd):
    p = binom_test((c, s-c), p=empirical_p) if c > min_reads else 1
    return True if p < pvalue_theshold else False

from rgt.GenomicRegionSet import GenomicRegionSet
res = GenomicRegionSet('identified_peaks')

# Create an instance of GenomicRegionSet. Use e+1 as ending coordinate, such that the considered bin overlap in one position
from rgt.GenomicRegion import GenomicRegion
for i, c in enumerate(cov.overall_cov):
    if filter_peaks(c, s, p):
        chrom, s, e = cov.index2coordinates(i, regionset)
        res.add(GenomicRegion(chrom, s, e+1))

# Merge consecutive bin to one single peak
res.merge()

# Output the identified peaks and computed ChIP-seq signal
res.write(outtxt)
cov.write_bigwig(outbw, g.get_chromosome_sizes())
