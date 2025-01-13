# Telomere-C ver 1.0b (Modified for ChIP-seq, single end)

# Introduction

This pipeline processes Telomere-C sequencing data barcoded with unique molecular identifiers (UMIs) and produces output files, including Telomere-C peak coordinates and normalized read counts for downstream analysis. It also generates quality control (QC) metrics using GATK and FastQC, which are compatible with MultiQC to provide a comprehensive QC report.

The workflow is managed by Snakemake to ensure reproducible and scalable data analyses. It also prevents overwriting and allows resumption from interrupted jobs.

## Key Workflow Steps

**1. Adapter Marking**. Illumina adapter sequences are marked to enable their removal during alignment, preserving the integrity of the sequence data.

**2. BWA Alignment**. Processed reads are aligned to a reference genome using `BWA-MEM`. The pipeline merges the aligned reads with metadata from the pre-alignment steps.

**3. Deduplication**. Duplicated reads are marked by GATK's `MarkDuplicates` to eliminate PCR artifacts, ensuring that only unique reads are retained for analysis.

**4. Reads Normalization**. Aligned reads are normalized using deepTools's `bamCoverage` with RPKM method, bin = 100 bp. 

**5. Peak Calling**. For peak calling, `mcas2 callpeak` is applied to identify regions of significant peaks with cutoff `-q 0.01`.

# Graphic workflow of Telomere-C
![Workflow](dag.rule.pdf)
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

## Other packages
Telomere-C required bellowing packages and suggested versions:

|packge_name|current version|document|
|:---   |:---            |:---|
|cutadapt|3.4|https://cutadapt.readthedocs.io/en/stable/|
|gatk4|4.2.2.0|https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4|
|bwa|0.7.17|https://github.com/lh3/bwa|
|samtools|1.13|http://www.htslib.org/|
|bedtools|2.30.0|https://bedtools.readthedocs.io/en/latest/|
|macs2|2.2.7.1|https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html|
|deepTools|3.5.4|https://test-argparse-readoc.readthedocs.io/en/latest/index.html|
|fastQC|0.12.0|https://www.bioinformatics.babraham.ac.uk/projects/fastqc/|

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

For non-slurm user, uncomment line 16 and comment the line 13
```
15 # This command DON'T use --profile argument
16 snakemake --jobs 500 -k --latency-wait 120 --max-jobs-per-second 2 --restart-times 0 --cluster-config "${CLUSTRCONF}" --configfile "${CONFIGFILE}" --jobname "{jobid}.{cluster.name}" all
```

## Prepare the `fastqList.txt` File
Create a tab-separated file named fastqList.txt with the following format:
- Sample Name: Sample name
- PATH1: Provide the full paths to the Fastq files for read1.

Example:
```
Cell-rep1   /home/fastq/Cell-rep1.fastq
Cell-rep2   /home/fastq/Cell-rep2.fastq
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
All Output files were stored in the subdirectory of `results/align/`

## Main outputs
deduplicated BAM: 
- `results/align/markduplicates/<Sampe Name>-rep1.realn.mdup.MQ30.bam`
- `results/align/markduplicates/<Sampe Name>-rep2.realn.mdup.MQ30.bam`

normalized BigWig:
- `results/align/bamCoverage/<Sampe Name>-rep1.realn.mdup.MQ30.norm.100bp.bigwig`
- `results/align/bamCoverage/<Sampe Name>-rep2.realn.mdup.MQ30.norm.100bp.bigwig`

called peaks: 
- `results/align/macs2/<Sampe Name>/<Sampe Name>-rep1_peaks.narrowPeak`
- `results/align/macs2/<Sampe Name>/<Sampe Name>-rep2_peaks.narrowPeak`