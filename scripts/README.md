# Script: LinkerQC.R

## Reference
[4C-seq from beginning to end: A detailed protocol for sample preparation and data analysis (Krijger, P. H. L. et al, 2020, *Methods*)](https://pubmed.ncbi.nlm.nih.gov/31351925/)

## Description
In the pipeline, the fastq files are subgrouped by its 4C PCR primers combination (i.e. FtFb-FtFb, FtFb-unknown, unknown-unknown...). The `*clipreads.metrics.txt` files contain information of reads and clipped reads in each subgroup. This information is useful to help evaluate 4C-seq library. 

- Reads with two linkers reflect the **successful linker ligation** at the both ends
- Reads with one linker reflect **partial ligation** OR **incompleted fragmentation** at an one end
- Reads with no linker reflect **failed linker ligation**
- Clipped reads, are mostly reverse-complement 4C PCR primers at the 3' end of fragments due to short inserts, reflect **self-ligation**

## Specification

Inputs: All `*.clipreads.metrics.txt` files

Output: A .linkerQC.txt table with 3 columns: Two_linker, One_linker, No_linker describing bellow information.

|Specification            |Description|
|:---                     |:---       |
|Number of examined reads |Number of sequenced reads in each subgroup|
|Percent of examined reads|Percentage of subgroup reads in total sequenced reads|
|Number of clipped reads  |Number of reads containing reverse-complement 4C PCR primers at the 3' end were clipped in each subgroup|
|Percent of clipped reads |Percentage of clipped reads in each subgroup reads|
|Number of examined bases |Number of sequenced bases in each subgroup|
|Percent of examined bases|Percentage of subgroup bases in total sequenced bases|
|Number of clipped bases  |Number of bases containing reverse-complement 4C PCR primers at the 3’ end were clipped in each subgroup|
|Percent of clipped bases |Percentage of clipped bases in each subgroup bases|


