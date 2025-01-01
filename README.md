# Telomere-C ver 1.0

# Introduction

This pipeline processes Telomere-C sequencing data barcoded with unique molecular identifiers (UMIs) and produces output files, including Telomere-C peak coordinates and normalized read counts for downstream analysis. It also generates quality control (QC) metrics using GATK and FastQC, which are compatible with MultiQC to provide a comprehensive QC report.

The workflow is managed by Snakemake to ensure reproducible and scalable data analyses. It also prevents overwriting and allows resumption from interrupted jobs.

## Key Workflow Steps

**1. UMI Extraction**. Raw sequencing data is processed to extract UMI information within the reads. UMIs are critical for identifying and removing duplicate reads.

**2. Adapter Marking**. Illumina adapter sequences are marked to enable their removal during alignment, preserving the integrity of the sequence data.

**3. BWA Alignment**. Processed reads are aligned to a reference genome using BWA-MEM. The pipeline merges the aligned reads with metadata from the pre-alignment steps.

**4. UMI Grouping and Deduplication**. Reads are grouped by UMI and deduplicated to eliminate PCR artifacts, ensuring that only unique reads are retained for analysis.

**5. Reads Normalization**. Aligned reads are normalized using the Signal Extraction Scaling (SES) method and adjusted for GC-content using functions provided by the [Regulatory Genomics Toolbox](https://reg-gen.readthedocs.io/en/latest/).

**6. Peak Calling**. For peak calling, a binomial distribution is applied to identify regions of significant enrichment. Candidate peaks are filtered based on p-value and dynamic coverage threshold.

# Graphic workflow of Telomere-C
![Workflow](dag.svg)
---

# Prerequisites

## Conda/Miniconda
[Installing Miniconda
](https://docs.anaconda.com/miniconda/install/)

## Snakemake v7.25.0
[Full install instructions for Snakemake](https://snakemake.readthedocs.io/en/v7.25.0/snakefiles/best_practices.html).

**Build and activate environment**
```
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake=7.25.0 python=3.11
conda activate snakemake
```

## RGT - Regulatory Genomics Toolbox
[Full install instruction of RGT](https://reg-gen.readthedocs.io/en/latest/rgt/installation.html).

**Quick installation**
```
pip install cython numpy scipy
pip install RGT
```
Please refer [Configuration of Genomic Data](https://reg-gen.readthedocs.io/en/latest/rgt/setup_data.html) to configure Genomic Data in the home directory. python = 3.11 is recommended.

## Other packages
Telomere-C required bellowing packages and suggested versions:

|packge_name|current version|document|
|:---   |:---            |:---|
|cutadapt|3.4|https://cutadapt.readthedocs.io/en/stable/|
|gatk4|4.2.2.0|https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4|
|bwa|0.7.17|https://github.com/lh3/bwa|
|samtools|1.13|http://www.htslib.org/|
|telseq|0.0.2|https://github.com/zd1/telseq|
|bedtools|2.30.0|https://bedtools.readthedocs.io/en/latest/|
|macs2|2.2.7.1|https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html|
|deepTools|3.5.4|https://test-argparse-readoc.readthedocs.io/en/latest/index.html|
|fastQC|0.12.0|https://www.bioinformatics.babraham.ac.uk/projects/fastqc/|
|umi_tools|1.1.1|https://umi-tools.readthedocs.io/en/latest/QUICK_START.html|

## Indexes
- BWA indexes
```
bwa index <refernce.fasta>
```

- fasta indexes

- reference genome dictionary
```
gatk CreateSequenceDictionary -R "<refernce.fasta>" -O "<refernce.dict>"
```

If your refence fasta files are compressed (i.e. .gz or .bzip), use gatk tools to get normalized and uncompressed refence fasta file.
```
gatk NormalizeFasta -I <ref.fasta.gz> -O <ref.norm.fastq>
```

# Configures

## Reference genome
Modify the path to reference fasta in the `Snakefile` line 3 
```
ref_fasta="/home/references/CHM13v2/chm13v2.0.fasta"
```
Note: the index files should be in the same directory of reference fasta

## Configure `slurm_profile` using Cookiecutter (Only for the first time of running)
We use [cookiecutter](https://github.com/cookiecutter/cookiecutter) to configure slurm for Telomere-C pipeline.

Installation
```
pip install cookiecutter
```

Example setting for Slurm:
```
cookiecutter cookiecutter_slurm_tempate
```

And follow the setting below:
```
profile_name [slurm]: slurm_profile
sbatch_defaults []: --nodes=1 --cpus-per-task=8  --time=72:00:00 --mem=12G
cluster_config []:
Select advanced_argument_conversion:
1 - no
2 - yes
Choose from 1, 2 [1]: 2
cluster_name []: 
```

## Configure `scripts/primary_peak_call.py`

In lines 20-22, make sure the name of Genome data exists in `~/rgtdata` and is [configured](https://reg-gen.readthedocs.io/en/latest/rgt/setup_data.html) 
```
20 # You must complete the Configuration of Genomic Data in your first time of running
21 # Please check: https://reg-gen.readthedocs.io/en/latest/rgt/setup_data.html
22 g = GenomeData('CHM13v2')
                   ^^^^^^^
```

## Configure `snakemake-run.sh`
In the line 13, check if the path of `--conda-prefix` is correct. You can use command `conda env list` to find parent directory of current conda enviromnet. 

In this case, `/tgen_labs/barthel/software` is the `--conda-prefix`
```
conda env list
telomereC.py3.1      * /tgen_labs/barthel/software/miniforge3/envs/telomereC.py3.1
                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^
```

## Prepare the `fastqList.txt` File
Create a tab-separated file named fastqList.txt with the following format:
- Sample Name: Must include either `-input` or `-capture`.
- PATH1: Provide the full paths to the Fastq files for read1.
- PATH2: Provide the full paths to the Fastq files for read2.

Example:
```
Cell-input   /home/fastq/Cell-input.R1.fastq   /home/fastq/Cell-input.R2.fastq  
Cell-capture   /home/fastq/Cell-capture.R1.fastq   /home/fastq/Cell-capture.R2.fastq  
```

# Quick Start
Execute the following command to initiate a dry-run:
```
sh snakemake-run.sh  
```

Modify the Script for a full run by removing the `-n` option in snakemake-run.sh, line 12, then execute the command agagin
```
sh snakemake-run.sh  
```

# Output files
All Output files were stored in the subdirectory of `results/alig/`

## Main outputs
deduplicated BAM: 
- `results/align/UmiDeDup/<Sampe Name>-input.realn.mdup.MQ30.bam`
- `results/align/UmiDeDup/<Sampe Name>-capture.realn.mdup.MQ30.bam`

unnormalized BigWig:
- `results/align/bamCoverage/<Sampe Name>-input.realn.mdup.MQ30.norm.100bp.bigwig`
- `results/align/bamCoverage/<Sampe Name>-capture.realn.mdup.MQ30.norm.100bp.bigwig`

normalized BigWig files: 
- `results/align/RGT_peakCall/<Sampe Name>-capture.realn.mdup.MQ30.run_signal.bw`

called peaks: 
- `results/align/RGT_peakCall/<Sampe Name>-capture.realn.mdup.MQ30.run_peaks.merge.bed`
 
# FAQ

**Q: No peaks are called in `results/align/RGT_peakCall/`.**

A: This issue occurs when one of the input files during read normalization, either `input.bam` or `capture.bam`, is missing due to a timeout.

Verify that the following files exist:
- `results/align/UmiDeDup/<Sample Name>-input.realn.mdup.MQ30.bam`
- `results/align/UmiDeDup/<Sample Name>-capture.realn.mdup.MQ30.bam`

Remove the `results/align/RGT_peakCall/` directory, and re-run the pipeline using the command:
```
sh snakemake-run.sh
```

