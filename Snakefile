##"results/align/cutadapt/A2780-GT20/A2780-GT20.unknown-unknown.1.fastq.gz"        

## 5' 
#cutadapt -e 0.07 --no-indels -g file:data/ref/telomerec.cutadapt.fasta -G file:data/ref/telomerec.cutadapt.fasta -o sandbox/{name1}-{name2}.1.fastq.gz -p sandbox/{name1}-{name2}.2.fastq.gz data/fastq/A2780-0_5M_GT20-05417_CAATTAAC-CGAGATAT_S1_R1_001.fastq.gz data/fastq/A2780-0_5M_GT20-05417_CAATTAAC-CGAGATAT_S1_R2_001.fastq.gz
#cutadapt -e 0.16 -g file:data/ref/telomerec.cutadapt.fasta -G file:data/ref/telomerec.cutadapt.fasta -o sandbox2/{name1}-{name2}.1.fastq.gz -p sandbox2/{name1}-{name2}.2.fastq.gz data/fastq/A2780-0_5M_GT20-05417_CAATTAAC-CGAGATAT_S1_R1_001.fastq.gz data/fastq/A2780-0_5M_GT20-05417_CAATTAAC-CGAGATAT_S1_R2_001.fastq.gz

## 3'
#cutadapt -e 0.07 --no-indels -a file:data/ref/telomerec.cutadapt2.fasta -A file:data/ref/telomerec.cutadapt2.fasta -o sandbox2/{name1}-{name2}.1.fastq.gz -p sandbox2/{name1}-{name2}.2.fastq.gz data/fastq/A2780-0_5M_GT20-05417_CAATTAAC-CGAGATAT_S1_R1_001.fastq.gz data/fastq/A2780-0_5M_GT20-05417_CAATTAAC-CGAGATAT_S1_R2_001.fastq.gz


## Define adapter combinations
import itertools
adapters = ['FtFb','FtRb','RtFb','RtRb','unknown']
adapters_comb = list(itertools.product(adapters, repeat=2))
adapters_comb_str = ["{}-{}".format(r1,r2) for (r1,r2) in adapters_comb]

## 3' adapters are assumed using -g and -G for read1 and read2, respectively
## The adapters are anchored (indicated by the ^ in the fasta file)
## This means that the adapters are fixed to the 3' end

##TODO list
## 1. change to relatvive path and test the script if it is working
## 2. pull 'file path' from a userConfig file
## 3. pull patameters from userConfig


rule cutadapt:
### TODO pell value from config file
    input:
        R1 = "data/fastq/A2780-0_5M_GT20-05417_CAATTAAC-CGAGATAT_S1_R1_001.fastq.gz",
        R2 = "data/fastq/A2780-0_5M_GT20-05417_CAATTAAC-CGAGATAT_S1_R2_001.fastq.gz",
        adapters = "data/ref/telomerec.cutadapt.fasta"
    output:
        R1 = expand("results/align/cutadapt/{{aliquot_barcode}}/{{aliquot_barcode}}.{adapt}.1.fastq.gz", adapt = adapters_comb_str),
        R2 = expand("results/align/cutadapt/{{aliquot_barcode}}/{{aliquot_barcode}}.{adapt}.2.fastq.gz", adapt = adapters_comb_str)
    params:
        o = lambda wildcards: "results/align/cutadapt/{aliquot_barcode}/{aliquot_barcode}.{{name1}}-{{name2}}.1.fastq.gz".format(aliquot_barcode = wildcards.aliquot_barcode),
        p = lambda wildcards: "results/align/cutadapt/{aliquot_barcode}/{aliquot_barcode}.{{name1}}-{{name2}}.2.fastq.gz".format(aliquot_barcode = wildcards.aliquot_barcode)
    log:
        "logs/align/cutadapt/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/align/cutadapt/{aliquot_barcode}.txt"
    message:
        "Splitting FASTQ file by telomere bridge linker adapter permutations\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:"""cutadapt \
            -e 0.16 \
            -g file:{input.adapters} \
            -G file:{input.adapters} \
            -o {params.o} \
            -p {params.p} \
            {input.R1} {input.R2} \
            > {log} 2>&1"""

rule fq2ubam:
    input:
        R1 = "results/align/cutadapt/{aliquot_barcode}/{aliquot_barcode}.{adapt}.1.fastq.gz",
        R2 = "results/align/cutadapt/{aliquot_barcode}/{aliquot_barcode}.{adapt}.2.fastq.gz"
    output:
        "results/align/ubam/{aliquot_barcode}/{aliquot_barcode}.{adapt}.unaligned.bam"
    params:
### TODO pull value from congfig file
        RGID = "{adapt}.G5B3N.1",
        RGPL = "ILLUMINA",
        RGPU = "G5B3NM02838.1",
        RGLB = "BAR65723",
        RGDT = "20200321",
        RGSM = "A2780-GT20",
        RGCN = "JAX"
    log:
        "logs/align/fq2ubam/{aliquot_barcode}.{adapt}.log"
    benchmark:
        "benchmarks/align/fq2ubam/{aliquot_barcode}.{adapt}.txt"
    message:
        "Converting FASTQ file to uBAM format\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Index permutation: {wildcards.adapt}"
    shell:"""gatk --java-options -Xmx6g FastqToSam \
            --FASTQ={input.R1} \
            --FASTQ2={input.R2} \
            --OUTPUT={output} \
            --READ_GROUP_NAME=\"{params.RGID}\" \
            --PLATFORM_UNIT=\"{params.RGPU}\" \
            --SAMPLE_NAME=\"{params.RGSM}\" \
            --PLATFORM=\"{params.RGPL}\" \
            --LIBRARY_NAME=\"{params.RGLB}\" \
            --SEQUENCING_CENTER=\"{params.RGCN}\" \
            --SORT_ORDER=queryname \
            --TMP_DIR=/fastscratch/barthf/temp \
            > {log} 2>&1"""

rule markadapters:
    input:
        "results/align/ubam/{aliquot_barcode}/{aliquot_barcode}.{adapt}.unaligned.bam"
    output:
        bam = "results/align/markadapters/{aliquot_barcode}/{aliquot_barcode}.{adapt}.markadapters.bam",
        metric = "results/align/markadapters/{aliquot_barcode}/{aliquot_barcode}.{adapt}.markadapters.metrics.txt"
    log: 
        "logs/align/markadapters/{aliquot_barcode}.{adapt}.log"
    benchmark:
        "benchmarks/align/markadapters/{aliquot_barcode}.{adapt}.txt"
    message:
        "Adding XT tags. This marks Illumina Adapters and allows them to be removed in later steps\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Index permutation: {wildcards.adapt}"
    shell:
        "gatk --java-options -Xmx6g MarkIlluminaAdapters \
            --INPUT={input} \
            --OUTPUT={output.bam} \
            --METRICS={output.metric} \
            --TMP_DIR=/fastscratch/barthf/temp \
            > {log} 2>&1"

rule clipreads:
    input:
        ubam = "results/align/markadapters/{aliquot_barcode}/{aliquot_barcode}.{adapt}.markadapters.bam",
### TODO pull value from configfile
        clipseq = "data/ref/telomerec.clipreads.fasta"
    output:
        ubam = "results/align/clipreads/{aliquot_barcode}/{aliquot_barcode}.{adapt}.clipreads.bam",
        stats = "results/align/clipreads/{aliquot_barcode}/{aliquot_barcode}.{adapt}.clipreads.metrics.txt"
    log: 
        "logs/align/markadapters/{aliquot_barcode}.{adapt}.log"
    benchmark:
        "benchmarks/align/markadapters/{aliquot_barcode}.{adapt}.txt"
    message:
        "Clip reads. This clips the mostly reverse-complement 4C PCR primers on the 3' end of fragments due to short inserts\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Index permutation: {wildcards.adapt}"
    shell:
        "gatk --java-options -Xmx6g ClipReads \
            -I {input.ubam} \
            -O {output.ubam} \
            --clip-sequences-file {input.clipseq} \
            --output-statistics {output.stats} \
            > {log} 2>&1"

rule samtofastq_bwa_mergebamalignment:
    input:
        bam = "results/align/clipreads/{aliquot_barcode}/{aliquot_barcode}.{adapt}.clipreads.bam",
### TODO pull value from configfile
        ref = "data/ref/human_g1k_v37_decoy.fasta"
    output:
        bam = "results/align/bwa/{aliquot_barcode}/{aliquot_barcode}.{adapt}.aln.bam",
        bai = "results/align/bwa/{aliquot_barcode}/{aliquot_barcode}.{adapt}.aln.bai"
    log: 
        "logs/align/samtofastq_bwa_mergebamalignment/{aliquot_barcode}.{adapt}.log"
    threads: 12
    benchmark:
        "benchmarks/align/revertsam/{aliquot_barcode}.{adapt}.txt"
    message:
        "BAM to FASTQ --> BWA-MEM --> Merge BAM Alignment.\n"
        "The first step converts the reverted BAM to an interleaved FASTQ, removing Illumina "
        "adapters. The output is then piped to BWA-MEM and aligned. Aligned reads are merged "
        "with the original pre-aligned BAM to preserve original metadata, including read groups.\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Index permutation: {wildcards.adapt}"
    shell:
        "gatk --java-options '-Dsamjdk.buffer_size=131072 -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx128m' SamToFastq \
            --INPUT={input.bam} \
            --FASTQ=/dev/stdout \
            --CLIPPING_ATTRIBUTE=XT \
            --CLIPPING_ACTION=2 \
            --INTERLEAVE=true \
            --NON_PF=true \
            --TMP_DIR=/fastscratch/barthf/temp | \
         bwa mem -M -t {threads} -p {input.ref} /dev/stdin | \
         gatk --java-options '-Dsamjdk.buffer_size=131072 -Dsamjdk.use_async_io=true -Dsamjdk.compression_level=1 -XX:+UseStringCache -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx5000m' MergeBamAlignment \
            --ALIGNED_BAM=/dev/stdin \
            --UNMAPPED_BAM={input.bam} \
            --OUTPUT={output.bam} \
            --REFERENCE_SEQUENCE={input.ref} \
            --CREATE_INDEX=true \
            --ADD_MATE_CIGAR=true \
            --CLIP_ADAPTERS=false \
            --CLIP_OVERLAPPING_READS=true \
            --INCLUDE_SECONDARY_ALIGNMENTS=true \
            --MAX_INSERTIONS_OR_DELETIONS=-1 \
            --PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
            --ATTRIBUTES_TO_RETAIN=XS \
            --TMP_DIR=/fastscratch/barthf/temp \
            > {log} 2>&1"

rule markduplicates:
    input:
        expand("results/align/bwa/A2780-GT20/A2780-GT20.{adapt}.aln.bam", adapt = adapters_comb_str)
    output:
        bam = "results/align/markduplicates/{aliquot_barcode}.realn.mdup.bam",
        bai = "results/align/markduplicates/{aliquot_barcode}.realn.mdup.bai",
        metrics = "results/align/markduplicates/{aliquot_barcode}.metrics.txt"
    params:
        max_records = 6000000
    log:
        "logs/align/markduplicates/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/align/markduplicates/{aliquot_barcode}.txt"
    message:
        "Readgroup-specific BAM files are combined into a single BAM. "
        "Potential PCR duplicates are marked.\n"
        "Sample: {wildcards.aliquot_barcode}"
    run:
        multi_input = " ".join(["--INPUT=" + s for s in input])
        shell("gatk --java-options -Xmx6g MarkDuplicates \
            {multi_input} \
            --OUTPUT={output.bam} \
            --METRICS_FILE={output.metrics} \
            --CREATE_INDEX=true \
            --TMP_DIR=/fastscratch/barthf/temp \
            --MAX_RECORDS_IN_RAM={params.max_records} \
            > {log} 2>&1")


rule callpeaks:
    input:
        bam = "results/align/markduplicates/{aliquot_barcode}.realn.mdup.bam",
        rg = "data/ref/optimalrg.txt"
    output:
        filtered_bam = "results/align/macs2/{aliquot_barcode}/{aliquot_barcode}.realn.mdup.MQ30.bam",
        peaks = "results/align/macs2/{aliquot_barcode}/{aliquot_barcode}_peaks.xls"
    params:
        outdir = "results/align/macs2/{aliquot_barcode}/"
    log:
        "logs/align/macs2/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/align/macs2/{aliquot_barcode}.txt"
    message:
        "Calling peaks using MACS2.\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:"""
        samtools view -R {input.rg} -b -q 30 {input.bam} > {output.filtered_bam}; \
         samtools index {output.filtered_bam}; \
         macs2 callpeak \
            -t {output.filtered_bam} \
            --outdir {params.outdir} \
            -f BAM \
            -g hs \
            -n {wildcards.aliquot_barcode} \
            -B \
            -q 0.01 \
            > {log} 2>&1"""

rule telseq:
    input:
        "results/align/markduplicates/{aliquot_barcode}.realn.mdup.bam"
    output:
        "results/align/telseq/{aliquot_barcode}.telseq.txt"
    log:
        "logs/align/telseq/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/align/telseq/{aliquot_barcode}.txt"
    message:
        "Quantification of telomere sequences using TelSeq\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:"""
        telseq -r 151 -k 7 -o {output} {input} \
            > {log} 2>&1"""

rule bedtools_count:
    input:
        bam = "results/align/macs2/{aliquot_barcode}/{aliquot_barcode}.realn.mdup.MQ30.bam",
        windows = "data/ref/b37.100Kb.windows.bed"
    output:
        "results/align/bedtools/{aliquot_barcode}.counts.bed"
    log:
        "logs/align/bedtools/{aliquot_barcode}.counts.log"
    benchmark:
        "benchmarks/align/bedtools/{aliquot_barcode}.counts.txt"
    message:
        "Count read coverage per bin\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:"""
        bedtools intersect -a {input.windows} \
                   -b {input.bam} \
                   -c -sorted \
            > {output} \
            2> {log}"""

rule bedtools_gc:
    input:
        bed = "results/align/bedtools/{aliquot_barcode}.counts.bed",
### TODO pull value from configfile
        fa = "data/ref/human_g1k_v37_decoy.fasta" 
    output:
        "results/align/bedtools/{aliquot_barcode}.counts.gc.bed"
    log:
        "logs/align/bedtools/{aliquot_barcode}.counts.gc.log"
    benchmark:
        "benchmarks/align/bedtools/{aliquot_barcode}.counts.gc.txt"
    message:
        "Annotate bins by GC content\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:"""
        bedtools nuc -fi {input.fa} \
             -bed {input.bed} \
            | cut -f 1-4,6 \
            > {output} \
            2> {log}"""

## For this step to work, one needs to install the SRA toolkit (sra-tools on conda)
## We also need to setup the SRA environment and load it with keyfiles
## See https://ncbi.github.io/sra-tools/install_config.html for instructions
## The output base path needs to match the directory in SRA config (vdb-config -i)
rule prefetch:
    output:
### TODO change to relative path
        "~/ncbi/public/sra/{sraid}.sra"
    params:
        sraid = "SRR8616019"
    log:
        "logs/align/prefetch/{sraid}.log"
    message:
        "Downloading SRR file\n"
        "SRA ID: {wildcards.sraid}"
    shell:"""
        prefetch {params.sraid} \
            > {log} 2>&1"""


## SAM-dump (eg. convert SRR into SAM) downloaded file
rule samdump:
    input:
### TODO change to relative path
        srr = "~/ncbi/public/sra/{sraid}.sra"
    output:
        sam = temp("results/align/samdump/{sraid}.sam"),
        bam = "results/align/samdump/{sraid}.bam",
        bai = "results/align/samdump/{sraid}.bam.bai"
    log:
        "logs/align/samdump/{sraid}.log"
    message:
        "Convert SRR file to BAM and create index\n"
        "SRA ID: {wildcards.sraid}"
    shell:"""
        sam-dump {input.srr} \
            > {output.sam} \
            2> {log};
        samtools view -S -b {output.sam} > {output.bam} 2> {log};
        samtools index {output.bam} > {log} 2>&1;"""

## See https://www.biostars.org/p/92744/ for a short overview of steps
## To get reads per bin
rule bedtools_count_rna:
    input:
        bam = "results/align/samdump/{sraid}.bam",
### TODO pull value from configfile
        windows = "data/ref/b37.100Kb.windows.bed"
    output:
        "results/align/bedtools/{sraid}.counts.rna.bed"
    log:
        "logs/align/bedtools/{sraid}.counts.rna.log"
    message:
        "Count read coverage per bin for RNAseq\n"
        "Sample: {wildcards.sraid}"
    shell:"""
        bedtools intersect -a {input.windows} \
                   -b {input.bam} \
                   -c -sorted \
            > {output} \
            2> {log}"""

rule bedtools_gc_rna:
    input:
        bed = "results/align/bedtools/{sraid}.counts.rna.bed",
### TODO pull value from configfile
        fa = "data/ref/human_g1k_v37_decoy.fasta"
    output:
        "results/align/bedtools/{sraid}.counts.rna.gc.bed"
    log:
        "logs/align/bedtools/{sraid}.counts.rna.gc.log"
    message:
        "Annotate bins by GC content\n"
        "Sample: {wildcards.sraid}"
    shell:"""
        bedtools nuc -fi {input.fa} \
             -bed {input.bed} \
            | cut -f 1-4,6 \
            > {output} \
            2> {log}"""

rule bedtools_gencode_rna:
    input:
        bed = "results/align/bedtools/{sraid}.counts.rna.gc.bed",
### TODO pull value from configfile
        gtf = "data/ref/gencode.v19.flattened.captured.sorted.bed"
    output:
        "results/align/bedtools/{sraid}.counts.rna.gc.gencode.bed"
    log:
        "logs/align/bedtools/{sraid}.counts.rna.gc.gencode.log"
    message:
        "Annotate bins by overlap with Gencode GTF\n"
        "Sample: {wildcards.sraid}"
    shell:"""
        bedtools annotate -i {input.bed} \
            -files {input.gtf} \
            > {output} \
            2> {log}"""
### TODO can we generate these by a script
rule all:
    input: "results/align/macs2/A2780-GT20/A2780-GT20_peaks.xls", "results/align/telseq/A2780-GT20.telseq.txt", "results/align/bedtools/A2780-GT20.counts.gc.bed", "results/align/samdump/SRR8616019.bam", "results/align/bedtools/SRR8616019.counts.rna.gc.bed", "results/align/bedtools/SRR8616019.counts.rna.gc.gencode.bed"

# ## END ##
