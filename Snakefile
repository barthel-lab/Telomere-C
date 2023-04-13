# Configure
# Reference genome (Make sure coresponding index files are in the same director)
ref_fasta="/labs/barthel/references/CHM13v2/chm13v2.0.fasta"

# Define sample names
import pandas as pd
fastqls = pd.read_csv("fastqList.txt", sep='\t', header=None, names=["name","R1","R2"])
fastqls.index = fastqls['name']
Sname = pd.Series(fastqls['name'])

# Start
rule ExtractUmis:
    input:
        R1 = lambda wildcards: fastqls.loc[wildcards.aliquot_barcode][1],
        R2 = lambda wildcards: fastqls.loc[wildcards.aliquot_barcode][2]
    output:
        f1 = "results/align/ExtractUmis/{aliquot_barcode}/{aliquot_barcode}_primary_processed.1.fastq.gz",
        f2 = "results/align/ExtractUmis/{aliquot_barcode}/{aliquot_barcode}_primary_processed.2.fastq.gz"
    log:
        "logs/align/ExtractUmis/{aliquot_barcode}.pre.ExtractUmis.R1.log"
    message:
        "Extracting Umi information from fastq"
    shell:"""
             module add umi_tools/1.1.1
             umi_tools extract -I {input.R1} --bc-pattern="^(?P<umi_1>.{{5}})(?P<discard_1>.{{2}}).*" \
             --read2-in={input.R2} \
             --stdout={output.f1} --read2-out={output.f2} \
             --extract-method=regex > {log} 2>&1
          """

rule RemoveEmptyReads:
    input:
        f1 = "results/align/ExtractUmis/{aliquot_barcode}/{aliquot_barcode}_primary_processed.1.fastq.gz",
        f2 = "results/align/ExtractUmis/{aliquot_barcode}/{aliquot_barcode}_primary_processed.2.fastq.gz"
    output:
        R1 = "results/align/ExtractUmis/{aliquot_barcode}/{aliquot_barcode}_processed.1.fastq.gz",
        R2 = "results/align/ExtractUmis/{aliquot_barcode}/{aliquot_barcode}_processed.2.fastq.gz"
    log:
        "logs/align/ExtractUmis/{aliquot_barcode}.ExtractUmis.R1.log"
    shell:"""
             cutadapt -m 1 -o {output.R1} -p {output.R2} {input.f1} {input.f2};
          """

rule fq2ubam:
    input:
        R1 = "results/align/ExtractUmis/{aliquot_barcode}/{aliquot_barcode}_processed.1.fastq.gz",
        R2 = "results/align/ExtractUmis/{aliquot_barcode}/{aliquot_barcode}_processed.2.fastq.gz"
    output:
        "results/align/ubam/{aliquot_barcode}/{aliquot_barcode}.unaligned.bam"
    params:
        RGCN = 'TGen',
        RGSM = lambda wildcards: wildcards.aliquot_barcode
    log:
        "logs/align/fq2ubam/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/align/fq2ubam/{aliquot_barcode}.txt"
    message:
        "Converting FASTQ file to uBAM format\n"
        "Sample: {wildcards.aliquot_barcode}\n"
    conda:
        "envs/gatk4.yaml"
    shell:"""gatk --java-options -Xmx6g FastqToSam \
            --FASTQ {input.R1} \
            --FASTQ2 {input.R2} \
            --OUTPUT {output} \
            --SEQUENCING_CENTER \"{params.RGCN}\" \
            --SAMPLE_NAME \"{params.RGSM}\" \
            --SORT_ORDER queryname \
            --TMP_DIR Temp \
            > {log} 2>&1"""

rule markadapters:
    input:
        "results/align/ubam/{aliquot_barcode}/{aliquot_barcode}.unaligned.bam"
    output:
        bam = "results/align/markadapters/{aliquot_barcode}/{aliquot_barcode}.markadapters.bam",
        metric = "results/align/markadapters/{aliquot_barcode}/{aliquot_barcode}.markadapters.metrics.txt"
    log: 
        "logs/align/markadapters/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/align/markadapters/{aliquot_barcode}.txt"
    message:
        "Adding XT tags. This marks Illumina Adapters and allows them to be removed in later steps\n"
        "Sample: {wildcards.aliquot_barcode}\n"
    conda:
        "envs/gatk4.yaml"
    shell:
        """gatk --java-options -Xmx6g MarkIlluminaAdapters \
            --INPUT {input} \
            --OUTPUT {output.bam} \
            --METRICS {output.metric} \
            --TMP_DIR Temp \
            > {log} 2>&1"""

# QC
rule fastqc:
    input:
        "results/align/ubam/{aliquot_barcode}/{aliquot_barcode}.unaligned.bam"
    output:
        "results/align/fastqc_preclip/{aliquot_barcode}/{aliquot_barcode}.unaligned_fastqc.html"
    params:
        dir = "results/align/fastqc_preclip/{aliquot_barcode}"
    log:
        "logs/align/fastqc_preclip/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/align/fastqc_preclip/{aliquot_barcode}.txt"
    message:
        "Running FASTQC (pre-clipping)\n"
        "Sample: {wildcards.aliquot_barcode}\n"
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc \
            --extract \
            -o {params.dir} \
            -f bam \
            {input} \
            > {log} 2>&1"

rule samtofastq_bwa_mergebamalignment:
    input:
        bam = "results/align/markadapters/{aliquot_barcode}/{aliquot_barcode}.markadapters.bam",
        ref = ref_fasta
    output:
        bam = "results/align/bwa/{aliquot_barcode}/{aliquot_barcode}.aln.bam",
        bai = "results/align/bwa/{aliquot_barcode}/{aliquot_barcode}.aln.bai"
    log: 
        "logs/align/samtofastq_bwa_mergebamalignment/{aliquot_barcode}.log"
    threads: 12
    resources:
         mem_mb=128728
    benchmark:
        "benchmarks/align/revertsam/{aliquot_barcode}.txt"
    message:
        "BAM to FASTQ --> BWA-MEM --> Merge BAM Alignment.\n"
        "The first step converts the reverted BAM to an interleaved FASTQ, removing Illumina "
        "adapters. The output is then piped to BWA-MEM and aligned. Aligned reads are merged "
        "with the original pre-aligned BAM to preserve original metadata, including read groups.\n"
        "Sample: {wildcards.aliquot_barcode}\n"
    conda:
        "envs/bwa.yaml"
    shell:
        """gatk --java-options '-Dsamjdk.buffer_size=131072 -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx128m' SamToFastq \
            --INPUT {input.bam} \
            --FASTQ /dev/stdout \
            --CLIPPING_ATTRIBUTE XT \
            --CLIPPING_ACTION 2 \
            --INTERLEAVE true \
            --NON_PF true \
            --TMP_DIR Temp | \
         bwa mem -M -t {threads} -p {input.ref} /dev/stdin | \
         gatk --java-options '-Dsamjdk.buffer_size=131072 -Dsamjdk.use_async_io=true -Dsamjdk.compression_level=1 -XX:+UseStringCache -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx5000m' MergeBamAlignment \
            --ALIGNED_BAM /dev/stdin \
            --UNMAPPED_BAM {input.bam} \
            --OUTPUT {output.bam} \
            --REFERENCE_SEQUENCE {input.ref} \
            --CREATE_INDEX true \
            --ADD_MATE_CIGAR true \
            --CLIP_ADAPTERS false \
            --CLIP_OVERLAPPING_READS true \
            --INCLUDE_SECONDARY_ALIGNMENTS true \
            --MAX_INSERTIONS_OR_DELETIONS -1 \
            --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
            --ATTRIBUTES_TO_RETAIN XS \
            --TMP_DIR Temp \
            > {log} 2>&1"""

rule UmiReadGroup:
    input:
       "results/align/bwa/{aliquot_barcode}/{aliquot_barcode}.aln.bam"
    output:
        bam = "results/align/UmiReadGroup/{aliquot_barcode}_mapped_grouped.bam",
        groups = "results/align/UmiReadGroup/{aliquot_barcode}_mapped_grouped.txt",
        sort = "results/align/UmiReadGroup/{aliquot_barcode}_sorted_mapped_grouped.bam"
    log:
        "logs/align/UmiReadGroup/{aliquot_barcode}_UMIs.mapped_grouped.log"
    shell:
        """ module add umi_tools/1.1.1
            umi_tools group -I {input} --paired --group-out={output.groups} \
            --output-bam -S {output.bam} --umi-group-tag RX > {log} 2>&1
            samtools sort {output.bam} -o {output.sort}
            samtools index {output.sort}"""

rule UmiDeDup:
    input:
        "results/align/UmiReadGroup/{aliquot_barcode}_sorted_mapped_grouped.bam"        
    output:
        dedup = "results/align/UmiDeDup/{aliquot_barcode}.realn.mdup.bam",
        log = "results/align/UmiDeDup/{aliquot_barcode}.UmiDeDup.log"
    log:
        "logs/align/UmiDeDup/{aliquot_barcode}.UmiDeDup.log"
    params:
        extra=r"'{aliquot_barcode}'"
    shell:
        """ module add umi_tools/1.1.1
            umi_tools dedup -I {input} \
            --log={output.log} \
            --paired -S {output.dedup} --output-stats={params.extra} > {log} 2>&1"""

# QC
rule alignmetrics:
    input:
        bam = "results/align/UmiDeDup/{aliquot_barcode}.realn.mdup.bam",
        ref = ref_fasta 
    output:
        "results/align/alignmetrics/{aliquot_barcode}.AlignMetrics.txt"
    log:
        "logs/align/alignmetrics/{aliquot_barcode}.AlignMetrics.log"
    benchmark:
        "benchmarks/align/alignmetrics/{aliquot_barcode}.AlignMetrics.txt"
    message:
        "Computing Alignment Summary Metrics\n"
        "Sample: {wildcards.aliquot_barcode}"
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk --java-options -Xmx6g CollectAlignmentSummaryMetrics \
            -R {input.ref} \
            -I {input.bam} \
            -O {output} \
            --METRIC_ACCUMULATION_LEVEL READ_GROUP \
            > {log} 2>&1"
# QC
rule multiplemetrics:
    input:
        bam = "results/align/UmiDeDup/{aliquot_barcode}.realn.mdup.bam",
        ref = ref_fasta
    output:
        "results/align/multiplemetrics/{aliquot_barcode}.alignment_summary_metrics"
    log:
        "logs/align/multiplemetrics/{aliquot_barcode}.MultipleMetrics.log"
    benchmark:
        "benchmarks/align/multiplemetrics/{aliquot_barcode}.MultipleMetrics.txt"
    message:
        "Computing Multiple Metrics\n"
        "Sample: {wildcards.aliquot_barcode}"
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk --java-options -Xmx6g CollectMultipleMetrics \
            -R {input.ref} \
            -I {input.bam} \
            -O results/align/multiplemetrics/{wildcards.aliquot_barcode} \
            --METRIC_ACCUMULATION_LEVEL READ_GROUP \
            > {log} 2>&1"
# QC
rule collectinsertsizemetrics:
    input:
        "results/align/UmiDeDup/{aliquot_barcode}.realn.mdup.bam"
    output:
        txt = "results/align/insertmetrics/{aliquot_barcode}.insertmetrics.txt",
        pdf = "results/align/insertmetrics/{aliquot_barcode}.insertmetrics.pdf"
    log:
        "logs/align/insertmetrics/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/align/insertmetrics/{aliquot_barcode}.txt"
    message:
        "Collect Insert Size Metrics\n"
        "Sample: {wildcards.aliquot_barcode}"
    conda:
        "envs/gatk4.yaml"
    shell:"""
        gatk CollectInsertSizeMetrics \
            -I {input} \
            -O {output.txt} \
            -H {output.pdf} \
            -M 0.5 \
            --METRIC_ACCUMULATION_LEVEL ALL_READS \
            --METRIC_ACCUMULATION_LEVEL SAMPLE \
            --METRIC_ACCUMULATION_LEVEL LIBRARY \
            --METRIC_ACCUMULATION_LEVEL READ_GROUP \
            2> {log}"""

rule telseq:
    input:
        "results/align/UmiDeDup/{aliquot_barcode}.realn.mdup.bam"
    output:
        "results/align/telseq/{aliquot_barcode}.telseq.txt"
    log:
        "logs/align/telseq/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/align/telseq/{aliquot_barcode}.txt"
    message:
        "Quantification of telomere sequences using TelSeq\n"
        "Sample: {wildcards.aliquot_barcode}"
    conda:
        "envs/telseq.yaml"
    shell:"""
        telseq -r 151 -k 7 -o {output} {input} \
            > {log} 2>&1"""

rule MQ30_n_blacklist_filter:
    input:
        bam = "results/align/UmiDeDup/{aliquot_barcode}.realn.mdup.bam",
        blacklist = "/labs/barthel/references/CHM13v2/T2T.excluderanges.noTelo.bed"
    output:
        filtered_bam = "results/align/UmiDeDup/{aliquot_barcode}.realn.mdup.MQ30.bam",
        filtered_bai ="results/align/UmiDeDup/{aliquot_barcode}.realn.mdup.MQ30.bam.bai"
    log:
        "logs/align/MQ30/{aliquot_barcode}.log"
    message:
        "apply Q30 filter to bam files"
    shell:"""
        bedtools intersect -v -abam {input.bam} -b {input.blacklist} | \
        samtools view -b -q 30 > {output.filtered_bam}; \
        samtools index {output.filtered_bam}"""

# Analysis method
# So far, I do macs2 by another script manually. ./scripts/macs2.sh <treated> <input>
# I keep this rule to make sure the pipeline goes well
rule callpeaks:
    input:
        bam = "results/align/UmiDeDup/{aliquot_barcode}.realn.mdup.MQ30.bam",
        input = "results/align/UmiDeDup/{aliquot_barcode}.realn.mdup.MQ30.bam"
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
    conda:
        "envs/macs2.yaml"
    shell:"""
        # macs2 callpeak \
        #    -t {output.filtered_bam} \
        #    --outdir {params.outdir} \
        #    -g hs \
        #    -n {wildcards.aliquot_barcode} \
        #    -B \
        #    -q 0.01 \
        #    > {log} 2>&1
        echo 'do nothing;"""

# Analysis method
rule bamCoverge:
    input:
        "results/align/UmiDeDup/{aliquot_barcode}.realn.mdup.MQ30.bam"
    output:
        "results/align/bamCoverage/{aliquot_barcode}.realn.mdup.MQ30.norm.100bp.bigwig"
    params:
        bin = 100,
        normalization = "RPKM"
    log:
        "logs/align/bamCoverage/{aliquot_barcode}.log"
    message:
        "Normalize reads by RPKM and bin"
    shell:"""bamCoverage \
        -b {input} \
        -o {output} \
        -of bigwig \
        --binSize {params.bin} \
        --numberOfProcessors 10 \
        --ignoreDuplicates \
        --scaleFactor 1 \
        --normalizeUsing {params.normalization} \
        > {log} 2>&1"""

# QC
rule wgsmetrics:
    input:
        bam = "results/align/UmiDeDup/{aliquot_barcode}.realn.mdup.MQ30.bam",
        ref = ref_fasta 
    output:
        "results/align/wgsmetrics/{aliquot_barcode}.WgsMetrics.txt"
    log:
        "logs/align/wgsmetrics/{aliquot_barcode}.WgsMetrics.log"
    benchmark:
        "benchmarks/align/wgsmetrics/{aliquot_barcode}.WgsMetrics.txt"
    message:
        "Computing WGS Metrics\n"
        "Sample: {wildcards.aliquot_barcode}"
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk --java-options -Xmx6g CollectWgsMetrics \
            -R {input.ref} \
            -I {input.bam} \
            -O {output} \
            --USE_FAST_ALGORITHM false \
            > {log} 2>&1"

# Define QC output files of workflows
rule all:
   input:
        expand("results/align/macs2/{name}/{name}_peaks.xls",name=Sname),
        expand("results/align/telseq/{name}.telseq.txt",name=Sname),
        expand("results/align/bamCoverage/{aliquot_barcode}.realn.mdup.MQ30.norm.100bp.bigwig",aliquot_barcode=Sname),
        expand("results/align/fastqc_preclip/{aliquot_barcode}/{aliquot_barcode}.unaligned_fastqc.html",aliquot_barcode=Sname),
        expand("results/align/wgsmetrics/{aliquot_barcode}.WgsMetrics.txt",aliquot_barcode=Sname),
        expand("results/align/alignmetrics/{aliquot_barcode}.AlignMetrics.txt",aliquot_barcode=Sname),
        expand("results/align/multiplemetrics/{aliquot_barcode}.alignment_summary_metrics",aliquot_barcode=Sname),
        expand("results/align/insertmetrics/{aliquot_barcode}.insertmetrics.txt",aliquot_barcode=Sname),
        expand("results/align/insertmetrics/{aliquot_barcode}.insertmetrics.pdf",aliquot_barcode=Sname)


## END ##
