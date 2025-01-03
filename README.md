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
pip install python numpy scipy
pip install RGT
```
Please refer to [Configuration of Genomic Data](https://reg-gen.readthedocs.io/en/latest/rgt/setup_data.html) to configure Genomic Data in the home directory. python = 3.11 is recommended.

Example setting of for the non-default reference genome.

`~/rgtdata/data.config.user`
```
[CHM13v2]
genome: /tgen_labs/barthel/references/CHM13v2/chm13v2.0.fasta
chromosome_sizes: /tgen_labs/barthel/references/CHM13v2/chm13v2.0.size.genome
gene_regions: CHM13v2/CHM13v2.gene_regions.bed
annotation: CHM13v2/CHM13v2.annotation.gtf
gene_alias: hg38/alias_human.txt
```
- genome: The fasta file of the reference genome.
- chromosome_sizes: The size of the chromosome in the reference genome, can be acquired from the fasta index.
- gene_regions: Corrdinate of genes, can be acuqired by acquired by converting the gft file. e.g. `gtf2bed` from `GFFUtils`.
- annotation: Annotation of reference genome in gft format.
- gene_alias: Optional for our pipeline. You could assign it from existing gene_alias.

Files in `~/rgtdata/CHM13v2/`
```
CHM13v2.annotation.gtf
CHM13v2.gene_regions.bed
```

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
For Slurm users, we use [cookiecutter](https://github.com/cookiecutter/cookiecutter) to configure Slurm for the Telomere-C pipeline.

Installation
```
pip install cookiecutter
```
Remove the `slurm_profile` and build a new one.
```
rm -rf slurm_profile
cookiecutter gh:Snakemake-Profiles/slurm
```
Follow the setting below:
```
 [1/17] profile_name (slurm): slurm_profile
 [2/17] Select use_singularity
 1 - False
 2 - True
 Choose from [1/2] (1): 
 [3/17] Select use_conda
 1 - False
 2 - True
 Choose from [1/2] (1): 2
 [4/17] jobs (500): 
 [5/17] restart_times (0): 3
 [6/17] max_status_checks_per_second (10): 
 [7/17] max_jobs_per_second (10): 
 [8/17] latency_wait (5): 
 [9/17] Select print_shell_commands
 1 - False
 2 - True
 Choose from [1/2] (1): 
 [10/17] sbatch_defaults (): --nodes=1 --cpus-per-task=8  --time=72:00:00 --mem=12G
 [11/17] cluster_sidecar_help (Use cluster sidecar. NB! Requires snakemake >= 7.0! Enter to continue...): 
 [12/17] Select cluster_sidecar
 1 - yes
 2 - no
 Choose from [1/2] (1): 
 [13/17] cluster_name (): 
 [14/17] cluster_jobname (%r_%w): 
 [15/17] cluster_logpath (logs/slurm/%r/%j): 
 [16/17] cluster_config_help (The use of cluster-config is discouraged. Rather, set snakemake CLI options in the 
profile configuration file (see snakemake documentation on best practices). Enter to continue...): 
 [17/17] cluster_config (): 
```

## Configure `scripts/primary_peak_call.py`

In lines 20-22, make sure the name of Genome data exists in `~/rgtdata` and is [configured](https://reg-gen.readthedocs.io/en/latest/rgt/setup_data.html) 
```
20 # You must complete the Configuration of Genomic Data in your first time of running
21 # Please check: https://reg-gen.readthedocs.io/en/latest/rgt/setup_data.html
22 g = GenomeData('CHM13v2')
 ^^^^^^^
```

For non-slurm user, uncomment line 16 and comment the line 13
```
15 # This command DON'T use --profile argument
16 snakemake --jobs 500 -k --latency-wait 120 --max-jobs-per-second 2 --restart-times 0 --cluster-config "${CLUSTRCONF}" --configfile "${CONFIGFILE}" --jobname "{jobid}.{cluster.name}" all
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

Or remove specific intermediate data. For example:
```
rm results/align/RGT_peakCall/intermediate/<Sample Name>*.token
rm results/align/RGT_peakCall/.<Sample Name>*.token
```
