library(tidyverse)
library(ggplot2)

setwd("/Volumes/Helix-Common/USERS/barthf/Telomere-C/Telo-C/")

## Bedtools was used to bin 100Kb windows across the genome
## Coverage was counted using the sample bam file in each of these bins

## Version 1
#dat <- read.delim("sandbox/b37.100Kb.windows.counts.bedg", header = FALSE, as.is = TRUE)
#colnames(dat) <- c("chr", "start", "stop", "counts_tel")
#dat <- dat %>% filter(chr %in% c(as.character(1:22),"X"))

## Version 3
dat <- read.delim("results/align/bedtools/A2780-GT20.counts.gc.bed", header = TRUE, as.is = TRUE)
colnames(dat) <- c("chr", "start", "stop", "counts_tel", "gc")
dat <- dat %>% filter(chr %in% c(as.character(1:22),"X"))

## GC-normalize
## 1. Create GC bins
dat$gc_grp <- cut(dat$gc, breaks = c(0, 0.3, 0.35, 0.40, 0.45, 0.50, 0.55, 0.65))
## 2. Check if GC bins are normally distributed
ggplot(dat, aes(x = log2(counts_tel), color = gc_grp)) + geom_density()
## 3. Compute Mean/SD log2 counts per GC-bin and
## 4. Compute Z-scores
## Dropping super low GC bins
dat <- dat %>% 
  group_by(gc_grp) %>% 
  mutate(mean_log2_counts_tel = mean(log2(counts_tel+1)), sd_log2_counts_tel = sd(log2(counts_tel+1))) %>% 
  ungroup() %>%
  filter(gc >= 0.30) %>%
  mutate(z = (log2(counts_tel+1) - mean_log2_counts_tel) / sd_log2_counts_tel)

## Plot the coverage per bin for each chromosome
ggplot(dat, aes(x = start, y = z)) + 
  geom_point(size = 0.01) +
  facet_wrap(~chr, scales = "free_x")# + coord_cartesian(ylim = c(0,5000))

ggplot(dat, aes(x = start, y = counts_tel)) + 
  geom_point(size = 0.01) +
  facet_wrap(~chr, scales = "free_x")# + coord_cartesian(ylim = c(0,5000))

## RNA
dat2 <- read.delim("results/align/bedtools/SRR8616019.counts.rna.gc.gencode.bed", header = FALSE, as.is = TRUE)
colnames(dat2) <- c("chr", "start", "stop", "counts_tx", "gc", "prop_tx")
dat2 <- dat2 %>% filter(chr %in% c(as.character(1:22),"X"))

ggplot(dat2, aes(x = start, y = (counts_tx/((prop_tx+0.0001) * 100000)))) + geom_line() + facet_wrap(~chr, scales = "free_x") #+ coord_cartesian(ylim = c(0,500))


colnames(dat) <- c("chr", "start", "stop", "counts_tel", "gc")
colnames(dat2) <- c("chr", "start", "stop", "counts_tx", "gc", "prop_tx")
dat3 <- dat %>% left_join(dat2)

ggplot(dat3, aes(x = counts_tel, y = 0.1*(counts_tx/(prop_tx+0.0001) * 100000)))) + geom_point() + coord_cartesian()
cor.test(dat$counts_tel, 0.1*(dat3$counts_tx/((dat3$prop_tx+0.0001) * 100000)), method = "s")

fisher.test()

dat3 <- dat3 %>% filter(complete.cases(counts_tel, counts_tx)) %>% mutate(exp = counts_tx/((prop_tx+0.0001) * 100000))

## END ##
