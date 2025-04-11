# Usage: python RGT_peakCaller.py capture.bam input.bam
# Reference: https://reg-gen.readthedocs.io/en/latest/rgt/tutorial-peak-calling.html
from rgt.Util import GenomeData
from rgt.helper import get_chrom_sizes_as_genomicregionset
from rgt.THOR.get_extension_size import get_extension_size
from rgt.THOR.DualCoverageSet import DualCoverageSet
from rgt.CoverageSet import CoverageSet
import rgt.THOR
import sys
import os
from collections import defaultdict
from scipy.stats import binomtest
from rgt.THOR.tracker import Tracker

# 0. Configure
# I/O 
bamfile = sys.argv[1]
bamfile_input = sys.argv[2]
outdir = 'results/align/RGT_peakCall/'
outtxt = outdir + os.path.basename(bamfile).replace('.bam','.run_peaks.bed')
outbw =  outdir + os.path.basename(bamfile).replace('.bam','.run_signal.bw')
#os.mkdir(outdir)

# You must complete the Configuration of Genomic Data in your first time of running
# Please check: https://reg-gen.readthedocs.io/en/latest/rgt/setup_data.html
g = GenomeData('CHM13v2')
regionset = get_chrom_sizes_as_genomicregionset(g.get_chromosome_sizes())

# 1. Preprocessing Steps
# Compute the fragment size of the reads
ext, _ = get_extension_size(bamfile, start=0, end=300, stepsize=5)

chrom_sizes_dict={}
with open(g.get_chromosome_sizes()) as f:
  for line in f:
    line = line.strip()
    line = line.split('\t')
    c, e = line[0], int(line[1])
    chrom_sizes_dict[c] = e

dims = [len(bamfile), len(bamfile_input)]
tracker_instance = Tracker(p=0, bamfiles=['file1','file2'], genome=g.get_genome(), chrom_sizes=g.get_chromosome_sizes(), dims=dims, inputs=None, options='', version='')
# Create an instance of CoverageSet
cov = DualCoverageSet(name='coverage',region=regionset,genome_path=g.get_genome(),binsize=500, \
                      stepsize=300, rmdup=True,file_1=bamfile, ext_1=ext,file_2=bamfile, ext_2=ext, \
                      input_1=bamfile_input, ext_input_1=0, input_factor_1=None, input_2=bamfile_input, ext_input_2=0, input_factor_2=None, \
                      chrom_sizes=g.get_chromosome_sizes(), verbose=True, norm_strategy=1, no_gc_content=False, deadzones=None, \
                      factor_input_1=None, factor_input_2=None, chrom_sizes_dict=chrom_sizes_dict, debug=None, tracker=tracker_instance)
#cov = CoverageSet('IP coverage', regionset)
#cov_input = CoverageSet('input-dna coverage', regionset)

cov.coverage_from_bam(bam_file=bamfile, extension_size=ext, binsiz=500, stepsize=300)
cov_input.coverage_from_bam(bam_file=bamfile_input, extension_size=0, binsize=500, stepsize=300)

# Normalize against the GC-content followed by subtracting the input-DNA data from the IP data
cov.norm_gc_content(cov_input.coverage, g.get_genome(), g.get_chromosome_sizes())
cov.subtract(cov_input)

# 2. Peak Calling (Binomial distribution)
# This section is modified from https://reg-gen.readthedocs.io/en/latest/rgt/tutorial-peak-calling.html
# The working assumption is that the reads falling into a bin follow a Binomial distribution
from numpy import sum, mean
#s is the total number of mapped reads (after GC normalization, and subtracting the input-DNA data from the IP data)
s = sum(cov.overall_cov)
#p is the expected probability of a read landing in a given bin (ignoring bins with zero counts)
p = mean(cov.overall_cov[cov.overall_cov > 0]) / s
#pvalue = 1e-32
pvalue = 1e-24
print(f"overall s: {s}")
print(f"overall p: {p}")
#print(f"overall p value: {pvalue}")
#Given s and p, we can apply the binomial test to determine whether a particular bin has a higher read count than expected by chance
#We set a fixed pvalue threshold of 1e-24
#Using the code below, we determine what bin coverage (cov_threshold) is necessary in order to have a pvalue lower than 1e-24
i=0
while True:
  i+=1
  pval = binomtest(i, s, p).pvalue
  if pval < pvalue:
    cov_threshold = i
    break

print(f"new required Rd: {cov_threshold}")

#The above code works in a genome-wide manner, where we use the overall s and p values from the entire genome
#However, we can calculate chromosome specific bin coverage thresholds
#For each chromosome, we store the chromosome specific number of mapped reads (s), expected probability of a read landing in a given bin (p) and number of bins with nonzero counts
chrom_to_s = defaultdict(int)
chrom_to_ps = defaultdict(int)
chrom_to_numposcov = defaultdict(int)
for i, c in enumerate(cov.overall_cov):
  chrom, start, end = cov.index2coordinates(i, regionset)
  chrom_to_s[chrom] += c
  if c > 0:
    chrom_to_numposcov[chrom] += 1
    chrom_to_ps[chrom] += c

#We calculate the chromosome specific probabilies (p)
chrom_to_p = {}
for chrom in chrom_to_s:
  if chrom_to_numposcov[chrom] > 0:
    chrom_to_p[chrom] = chrom_to_ps[chrom]/chrom_to_numposcov[chrom]/chrom_to_s[chrom]
  else:
    print(chrom)
    chrom_to_p[chrom] = 0

print(chrom_to_s)
print(chrom_to_p)

#we calculate the chromosome specific bin coverage/read depth, using the fixed pvalue threshold of 1e-24
#We can also set a minRd, if we want to ensure a minimum bin coverage (that may be higher than what is necessary to satisfy the pvalue threshold)
#Setting minRd to 0 ensures we only care about the pvalue threshold
minRd = 0
chrom_to_rd = {}
for chrom in chrom_to_s:
  s = chrom_to_s[chrom]
  p = chrom_to_p[chrom]
  if s==0:
    continue
  i = 0
  while True:
    i += 1
    if i > s:
      break
    pval = binomtest(i, s, p).pvalue
    if pval < pvalue:
      chrom_to_rd[chrom] = max(i,minRd)
      #chrom_to_rd[chrom] = i
      break

print(chrom_to_rd)

# Ignore bins with low number of reads
#from scipy.stats import binom_test
from scipy.stats import binomtest
def filter_peaks(c, empirical_s, empirical_p, pvalue_theshold=pvalue, min_reads=minRd):
    #p = binom_test((c, s-c), p=empirical_p) if c > min_reads else 1
    p = binomtest(c, empirical_s, empirical_p).pvalue if c > min_reads else 1
    return True if p < pvalue_theshold else False

#This functions returns true as long as the number of reads mapping to a bin is larger than the bin coverage threshold for the corresponding chromosome
def filter_peaks_cov(c, chrom, chrom_to_rd):
  if chrom not in chrom_to_rd:
    return False
  rd = chrom_to_rd[chrom]
  if c >= rd:
    return True
  else:
    return False

from rgt.GenomicRegionSet import GenomicRegionSet
res = GenomicRegionSet('identified_peaks')

# Create an instance of GenomicRegionSet. Use e+1 as ending coordinate, such that the considered bin overlap in one position
from rgt.GenomicRegion import GenomicRegion
for i, c in enumerate(cov.overall_cov):
    chrom, start, end = cov.index2coordinates(i, regionset)
    s = chrom_to_s[chrom]
    if s==0:
      continue
    p = chrom_to_p[chrom]
    #if filter_peaks(c, s, p):
    if filter_peaks_cov(c, chrom, chrom_to_rd):
        #chrom, s, e = cov.index2coordinates(i, regionset)
        res.add(GenomicRegion(chrom, start, end+1))

# Merge consecutive bin to one single peak
res.merge()

# Output the identified peaks and computed ChIP-seq signal
res.write(outtxt)
cov.write_bigwig(outbw, g.get_chromosome_sizes())
