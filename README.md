# Telomere-C
Telomere interaction analysis
---

# To start:

1. [Build environment](https://github.com/fpbarthel/Telomere-C/tree/deploy_on_dback#1-build-enviroment)
2. [Check dependent packages](https://github.com/fpbarthel/Telomere-C/tree/deploy_on_dback#2-check-dependent-packages)
3. [Configure Telomere-C](https://github.com/fpbarthel/Telomere-C/tree/deploy_on_dback#3-configure-telomere-c)
4. [Configure Snakemake --profile to run on Slurm (Optional)](https://github.com/fpbarthel/Telomere-C/tree/deploy_on_dback#4-configure-snakemake---profile-optional)

5. Enter compute node

        srun --pty --time=72:00:00 /bin/bash

6. Activate condaStart snakemake environment

        conda activate snakemake

7. run the job

        sh snakemake-run.sh

# 1 Build enviroment

## 1.1 Install Conda
To installed Anaconda, see [Conda install instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

The current version is ***conda 4.10.3***

## 1.2 Install Snakemake
To use snakemake, you have install a Conda-based Python3 distribution, here recommend Mambaforge. The current versions are ***snakemake 6.8.0*** and ***mamba 0.15.3***

*This example assume you're running on a Linux x86_64 machine.* 

For other optionals, see [full install instruction](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). 

### Download and install Mambaforge

    wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
    bash Mambaforge-Linux-x86_64.sh

### Build and activate enviroment

    mamba create -c conda-forge -c bioconda -n snakemake snakemake
    conda activate snakemake

You can check snakemake version by typing
    
    snakemake --version

# 2 Check dependent packages

Telomere-C required bellowing packages:

|packge_name|current version|document|
|:---   |:---            |:---|
|cutadapt|3.4|https://cutadapt.readthedocs.io/en/stable/|
|gatk4|4.2.2.0|https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4|
|bwa|0.7.17|https://github.com/lh3/bwa|
|samtools|1.13|http://www.htslib.org/|
|telseq|0.0.2|https://github.com/zd1/telseq|
|bedtools|2.30.0|https://bedtools.readthedocs.io/en/latest/|
|macs2|2.2.7.1|https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html|
|sra-tools|2.8.0|https://github.com/ncbi/sra-tools/wiki|
|r-base| 4.0.5|https://www.r-project.org/|
|r-tidyverse|1.3.1|https://www.tidyverse.org/|

These dependencies will be automatically downloaded and installin the first run. Environment settings were described in `/Telomere-C/envs/` .yaml files.

## 2.1 Configure SRA-tools 

For full instruction, see https://github.com/ncbi/sra-tools/wiki/05.-Toolkit-Configuration for instructions.

To enter config page, type:

    vdb-config -i

, then check on options: 

- **'Enable Remote Access'**
- **'Enable Local Caching'** 

Then, set SRA tools' **Default Import Path**. The default value should be

`/home/<user>/ncbi/public` or 

`~/ncbi/public`

# 3 Configure Telomere-C

## 3.1 Download Telomere-C

    git clone https://github.com/fpbarthel/Telomere-C.git

## 3.2 Prepare input files

0. Go to working directory

        cd Telomere-C

1. Copy or create a soft link pair-end **sequencing fastq files** to `data/fastq`. For example:

        ln -s <your sequence_R1.fastq> data/fastq/
        ln -s <your sequence_R2.fastq> data/fastq/

2. Copy or create a soft link **reference fastq files**, **index files**, and **dictionary file** to `data/ref`.For example:

        ln -s <refernce.fasta> data/ref/
        ln -s <refernce.fasta.btw> data/ref/
        ln -s <refernce.fasta.pac> data/ref/
        ln -s <refernce.fasta.ann> data/ref/
        ln -s <refernce.fasta.amb> data/ref/
        ln -s <refernce.fasta.sa> data/ref/
        ln -s <refernce.dict> data/ref/


3. Create **`ncbi/public`** soft link to the working directory 

        ln -s ~/ncbi /PATH/TO/Telomere-C

4. Copy or create soft a link **genome annotation file** to `data/`.

        ln -s <annotaion.BED/GFF/VCF> data/

**Note1:** To generate index files:

    bwa index <refernce.fasta>

See [bwa](https://github.com/lh3/bwa#getting-started) for the detailed instruction.

**Note2:** To generate dictionary files:

    gatk CreateSequenceDictionary -R "<refernce.fasta>" -O "<refernce.dict>"

**Note3:** If your refence fasta files are compressed (i.e. .gz or .bzip), use gatk tools to get normalized and uncompressed refence fasta file.

    gatk NormalizeFasta -I <ref.fasta.gz> -O <ref.norm.fastq>
    
See [gatk CreateSequenceDictionary](https://gatk.broadinstitute.org/hc/en-us/articles/360036729911-CreateSequenceDictionary-Picard-) for the detailed instruction.

## 3.3 Configuration

To start configure Telomere-C, use text editor to open Snakefile.

    cd Telomere-C
    vim Snakefile

**#==| Start of Configure |==#**

### Sample name

Give an any string as your sample name. This will be the prefix of your output files.

`Sname = "<string>"`

### Sequence fastq files

Assign path of your sequence .fastq files, both read1 and read2. Recommend to put your sequence file in `data/fastq/` director.

`input_read1= "<data/fastq/R1.fastq.gz>"`

`input_read2= "<data/fastq/R2.fastq.gz>"`

### Reference fastq genome

Assign a path of reference genome .fasta file. Recommend to put the file in `data/ref/` director.

`ref_fasta= "<data/ref/genome_ref.fasta>"`

### sra refence genome

`sra_ref_fasta="<data/ref/sra_genome_ref.fasta>"`

### sra refence genome annoation 

`sra_anno = "<sra_anno.bed>"`

### Bin file

Assign a path of bin bed file for counting coverage per bin. The default file is in the `data/`. 

`bin = "<data/bin.windows.bed>"`

`sra_bin="<data/srq_bin.windows.bed>"`

To get bin file, first, generate the genome file describing size of chromosome. It should tab delimited and structured as follows: 
        
        <chromName><TAB><chromSize>

        For example, Human (hg19):
        chr1	249250621
        chr2	243199373
        ...
        chr18_gl000207_random	4262

typeing:

        cat <data/ref/refernce.dict>|awk '{print $2 "\t" $3}'|sed 's/SN\:chr//g'|sed 's/LN\://g'|tail -n +2> <refernce.genome>

Second generate bin bed file:

        bedtools makewindows -g <refernce.genome> -w <window_size>|cat <bin.windows.bed>.bed

### Optimal  Group ID

The default setting is

`opt_rg= "data/optimalrg.txt"`

### Metadata 

Metadata will be written into Sam file.

The READ_GROUP_NAME will be assigned automatically according to telomere adapter combination.

#### Sequecning platform

The platform type (e.g. illumina, solid) to insert into the read group header

`platform = "<string>"`

#### Sequecning platform unit number 

The platform unit (often run_barcode.lane) to insert into the read group header

`unit = "<string>"`

#### Library name

The library name to place into the LB attribute in the read group header

`lib = "<string>"`

#### Sequencing date

Date the run was produced, to insert into the read group header

`date = "<string>"`

#### Sequencing center

The sequencing center from which the data originated

`center = "<string>"`

### SRA ID

Give an SRA ID for downloading SRA file, a RNA-seq data.

`SRAID = "<SRA ID>"`

### Telomere-C adapter

A fasta file describing 5' anchored telomere-C primers sequence. This file is required for adapter trimming. The default file is:

`ADPT = "data/telomerec.cutadapt.fasta"`

### Telomeric clip reads

A fasta files describing telomere-C PCR primers sequence and reverse-complement sequence on the 3' end of fragments due to short inserts. The default file is:

`CLIP = "data/telomerec.clipreads.fasta"`

**#==| END of Configure|==#**

# 4 Configure Snakemake --profile (Optional)

*The snakemake-run.sh provides two option of command lines for running Telomere-C with or without using `--profile` option. Simply add "#" to the head of unwanted command line.*

Adapting Snakemake to a particular environment (i.e. cluster systems like Slurm, SGE, Grid middleware, cloud computing) can entail many flags and options. Therefore, since Snakemake 4.1, it is possible to specify a configuration profile to be used to obtain default options. See [detailed doc](https://snakemake.readthedocs.io/en/stable/executing/cli.html?highlight=--profile#profiles) 

To get profile template, see [Snakemake-profiles project](https://github.com/snakemake-profiles/doc)

Example setting for Slurm:

    $ cookiecutter slurm
    profile_name [slurm]: slurm_profile
    sbatch_defaults []: --nodes=1 --cpus-per-task=2  --time=72:00:00 --mem=12G
    cluster_config []:
    Select advanced_argument_conversion:
    1 - no
    2 - yes
    Choose from 1, 2 [1]: 2
    cluster_name []: 

## Output files

Output files were stored in the subdirector of `results/alig/`

### Main results

|Results|Files|
|:---|:---|
|BWA-MEM aligned reads and index|`bwa/<sample_name>/<sample_name>.<adapter_combination>.aln.bam`|
||`bwa/<sample_name>/<sample_name>.<adapter_combination>.aln.bai`|
|A single bam file by combining adapter_combination read groups with mark of PCR duplication. AND its index as well as matrix|`markduplicates/<sample_name>.realn.mdup.bam`|
||`markduplicates/<sample_name>.realn.mdup.bai`|
||`markduplicates/<sample_name>.metrics.txt`|
|Peaks calling results using MACS2|`macs2/<sample_name>/<sample_name>_peaks.narrowPeak`|
||`macs2/<sample_name>/<sample_name>_peaks.xls`|
||`macs2/<sample_name>/<sample_name>_control_lambda.bdg`|
||`macs2/<sample_name>/<sample_name>_treat_pileup.bdg`|
||`macs2/<sample_name>/<sample_name>_model.r`|
||`macs2/<sample_name>/<sample_name>_summits.bed`|
|Quantification of telomere sequences using TelSeq|`telseq/<sample_name>.telseq.txt`|
|Covert count read coverage per bin from Mapq 30 filtered bam file AND GC content annotated bed file|`bedtools/A2780-GT20.counts.bed`|
||`bedtools/AA2780-GT20.counts.gc.bed`|
|Convert SRR RNA-seq file to BAM and create index|`samdump/<SRR_ID>.bam`|
||`samdump/<SRR_ID>.bam.bai`|
|Count read coverage per bin for RNAseq, annotation of GC contents, and Gencode GTF annotation|`bedtools/<SRR_ID>.counts.rna.bed`|
||`bedtools/<SRR_ID>.counts.rna.gc.bed`|
||`bedtools/<SRR_ID>.counts.rna.gc.gencode.bed`|

### Interndeiate files
|Results|Files|
|:---|:---|
|5' anchored adapter cut reads|`cutadapt/<sample_name>/<sample_name>.<adapter_combination>.fastq.gz`|
|Unaligment bam file converted from fastq|`ubam/<sample_name>/<sample_name>.<adapter_combination>.unaligned.bam`|
|adapter marked Unaligment bam file and its matrix|`markadapters/<sample_name>/<sample_name>.<adapter_combination>.markadapters.bam`|
||`markadapters/<sample_name>/<sample_name>.<adapter_combination>.markadapters.metrics.txt`|
|3' adapter clipped reads and its matrix|`clipreads/<sample_name>/<sample_name>.<adapter_combination>.clipreads.bam`|
||`clipreads/<sample_name>/<sample_name>.<adapter_combination>.clipreads.metrics.txt`|
|Mapq 30 filtered bam file and its index|`macs2/<sample_name>/<sample_name>.realn.mdup.MQ30.bam`|
||`macs2/<sample_name>/<sample_name>.realn.mdup.MQ30.bam.bai`|

