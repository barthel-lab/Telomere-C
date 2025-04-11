rule ExtractUmis:
    input:
        R1 = get_TelomereC_r1,
        R2 = get_TelomereC_r2
    output:
        f1 = "results/ExtractUmis/{sample}/{sample}-{enrich}_primary_processed.1.fastq.gz",
        f2 = "results/ExtractUmis/{sample}/{sample}-{enrich}_primary_processed.2.fastq.gz"
    message:
        "Extracting Umi information from fastq"
    conda:
        "telomereC.py3.1"
    shell:"""
        umi_tools extract -I {input.R1} --bc-pattern="^(?P<umi_1>.{{5}})(?P<discard_1>.{{2}}).*" \
            --read2-in={input.R2} \
            --stdout={output.f1} --read2-out={output.f2} \
            --extract-method=regex
        """

rule RemoveEmptyReads:
    input:
        f1 = "results/ExtractUmis/{sample}/{sample}-{enrich}_primary_processed.1.fastq.gz",
        f2 = "results/ExtractUmis/{sample}/{sample}-{enrich}_primary_processed.2.fastq.gz"
    output:
        R1 = "results/ExtractUmis/{sample}/{sample}-{enrich}_processed.1.fastq.gz",
        R2 = "results/ExtractUmis/{sample}/{sample}-{enrich}_processed.2.fastq.gz"
    conda:
        "telomereC.py3.1"
    shell:"""
        cutadapt -m 1 -o {output.R1} -p {output.R2} {input.f1} {input.f2};
        """

rule fq2ubam:
    input:
        R1 = "results/ExtractUmis/{sample}/{sample}-{enrich}_processed.1.fastq.gz",
        R2 = "results/ExtractUmis/{sample}/{sample}-{enrich}_processed.2.fastq.gz"
    output:
        "results/ubam/{sample}/{sample}-{enrich}.unaligned.bam"
    params:
        RGCN = 'TGen',
        RGSM = lambda wildcards: wildcards.sample
    message:
        "Converting FASTQ file to uBAM format\n"
    conda:
        "telomereC.py3.1"
    shell:"""
        gatk --java-options -Xmx6g FastqToSam \
            --FASTQ {input.R1} \
            --FASTQ2 {input.R2} \
            --OUTPUT {output} \
            --SEQUENCING_CENTER \"{params.RGCN}\" \
            --SAMPLE_NAME \"{params.RGSM}\" \
            --SORT_ORDER queryname \
            --TMP_DIR Temp
        """

rule markadapters:
    input:
        "results/ubam/{sample}/{sample}-{enrich}.unaligned.bam"
    output:
        bam = "results/markadapters/{sample}/{sample}-{enrich}.markadapters.bam",
        metric = "results/markadapters/{sample}/{sample}-{enrich}.markadapters.metrics.txt"
    message:
        "Adding XT tags. This marks Illumina Adapters and allows them to be removed in later steps\n"
    conda:
        "telomereC.py3.1"
    shell:"""
        gatk --java-options -Xmx6g MarkIlluminaAdapters \
            --INPUT {input} \
            --OUTPUT {output.bam} \
            --METRICS {output.metric} \
            --TMP_DIR Temp
        """

# QC
rule fastqc:
    input:
        "results/ubam/{sample}/{sample}-{enrich}.unaligned.bam"
    output:
        "results/fastqc_preclip/{sample}/{sample}-{enrich}.unaligned_fastqc.html"
    params:
        dir = "results/fastqc_preclip/{sample}"
    message:
        "Running FASTQC (pre-clipping)\n"
    conda:
        "telomereC.py3.1"
    shell:"""
        fastqc \
            --extract \
            -o {params.dir} \
            -f bam \
            {input}
        """

rule samtofastq_bwa_mergebamalignment:
    input:
        bam = "results/markadapters/{sample}/{sample}-{enrich}.markadapters.bam",
        ref = ref_fasta
    output:
        bam = "results/bwa/{sample}/{sample}-{enrich}.aln.bam",
        bai = "results/bwa/{sample}/{sample}-{enrich}.aln.bai"
    threads: 12
    resources:
         mem_mb=250000
    message:
        "BAM to FASTQ --> BWA-MEM --> Merge BAM Alignment.\n"
        "The first step converts the reverted BAM to an interleaved FASTQ, removing Illumina "
        "adapters. The output is then piped to BWA-MEM and aligned. Aligned reads are merged "
        "with the original pre-aligned BAM to preserve original metadata, including read groups.\n"
        "Sample: {wildcards.sample}\n"
    conda:
        "telomereC.py3.1"    
    shell:"""
        gatk --java-options '-Dsamjdk.buffer_size=131072 -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx128m' SamToFastq \
            --INPUT {input.bam} \
            --FASTQ /dev/stdout \
            --CLIPPING_ATTRIBUTE XT \
            --CLIPPING_ACTION 2 \
            --INTERLEAVE true \
            --NON_PF true \
            --TMP_DIR Temp | \
        bwa mem -M -t {threads} -p {input.ref} /dev/stdin | \
        gatk --java-options '-Dsamjdk.buffer_size=131072 -Dsamjdk.use_async_io=true -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx5000m' MergeBamAlignment \
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
            --TMP_DIR Temp
        """

rule UmiReadGroup:
    input:
        bam = "results/bwa/{sample}/{sample}-{enrich}.aln.bam"
    output:
        bam = "results/UmiReadGroup/{sample}-{enrich}_mapped_grouped.bam",
        groups = "results/UmiReadGroup/{sample}-{enrich}_mapped_grouped.txt",
        sort = "results/UmiReadGroup/{sample}-{enrich}_sorted_mapped_grouped.bam"
    conda:
        "telomereC.py3.1"
    shell:"""
        umi_tools group -I {input.bam} --paired --group-out={output.groups} \
            --output-bam -S {output.bam} --umi-group-tag RX;
            samtools sort {output.bam} -o {output.sort};
            samtools index {output.sort};
        """

rule UmiDeDup:
    input:
        "results/UmiReadGroup/{sample}-{enrich}_mapped_grouped.bam"        
    output:
        dedup = "results/UmiDeDup/{sample}-{enrich}.realn.mdup.bam",
        log = "results/UmiDeDup/{sample}-{enrich}.UmiDeDup.log"
    params:
        extra=r"'{sample}'"
    conda:
        "telomereC.py3.1"  
    shell:"""
        umi_tools dedup -I {input} \
            --log={output.log} \
            --paired -S {output.dedup} --output-stats={params.extra}
        """

# QC
rule alignmetrics:
    input:
        bam = "results/UmiDeDup/{sample}-{enrich}.realn.mdup.bam",
        ref = ref_fasta 
    output:
        "results/alignmetrics/{sample}-{enrich}.AlignMetrics.txt"
    message:
        "Computing Alignment Summary Metrics\n"
        "Sample: {wildcards.sample}"
    conda:
        "telomereC.py3.1"  
    shell:"""
        gatk --java-options -Xmx6g CollectAlignmentSummaryMetrics \
            -R {input.ref} \
            -I {input.bam} \
            -O {output} \
            --METRIC_ACCUMULATION_LEVEL READ_GROUP
        """
# QC
rule multiplemetrics:
    input:
        bam = "results/UmiDeDup/{sample}-{enrich}.realn.mdup.bam",
        ref = ref_fasta
    output:
        "results/multiplemetrics/{sample}-{enrich}.alignment_summary_metrics"
    message:
        "Computing Multiple Metrics\n"
        "Sample: {wildcards.sample}"
    conda:
        "telomereC.py3.1" 
    shell:"""
        gatk --java-options -Xmx6g CollectMultipleMetrics \
            -R {input.ref} \
            -I {input.bam} \
            -O results/align/multiplemetrics/{wildcards.sample} \
            --METRIC_ACCUMULATION_LEVEL READ_GROUP
        """
# QC
rule collectinsertsizemetrics:
    input:
        "results/UmiDeDup/{sample}-{enrich}.realn.mdup.bam"
    output:
        txt = "results/insertmetrics/{sample}-{enrich}.insertmetrics.txt",
        pdf = "results/insertmetrics/{sample}-{enrich}.insertmetrics.pdf"
    message:
        "Collect Insert Size Metrics\n"
    conda:
        "telomereC.py3.1" 
    shell:"""
        gatk CollectInsertSizeMetrics \
            -I {input} \
            -O {output.txt} \
            -H {output.pdf} \
            -M 0.5 \
            --METRIC_ACCUMULATION_LEVEL ALL_READS \
            --METRIC_ACCUMULATION_LEVEL SAMPLE \
            --METRIC_ACCUMULATION_LEVEL LIBRARY \
            --METRIC_ACCUMULATION_LEVEL READ_GROUP
        """

rule telseq:
    input:
        "results/UmiDeDup/{sample}-{enrich}.realn.mdup.bam"
    output:
        "results/telseq/{sample}-{enrich}.telseq.txt"
    message:
        "Quantification of telomere sequences using TelSeq\n"
        "Sample: {wildcards.sample}"
    conda:
        "telomereC.py3.1" 
    shell:"""
        telseq -r 151 -k 7 -o {output} {input}
        """

rule MQ30_n_blacklist_filter:
    input:
        bam = "results/UmiDeDup/{sample}-{enrich}.realn.mdup.bam",
        blacklist = blacklist
    output:
        filtered_bam = "results/UmiDeDup/{sample}-{enrich}.realn.mdup.MQ30.bam",
        filtered_bai ="results/UmiDeDup/{sample}-{enrich}.realn.mdup.MQ30.bam.bai"
    message:
        "apply Q30 filter to bam files"
    conda:
        "telomereC.py3.1" 
    shell:"""
        bedtools intersect -v -abam {input.bam} -b {input.blacklist} | \
        samtools view -b -q 30 > {output.filtered_bam}; \
        samtools index {output.filtered_bam}
        """
# Analysis method
rule bamCoverge:
    input:
        "results/UmiDeDup/{sample}-{enrich}.realn.mdup.MQ30.bam"
    output:
        "results/bamCoverage/{sample}-{enrich}.realn.mdup.MQ30.norm.100bp.bigwig"
    params:
        bin = 100,
        normalization = "RPKM"
    message:
        "Normalize reads by RPKM and bin"
    conda:
        "telomereC.py3.1" 
    shell:"""
        bamCoverage \
            -b {input} \
            -o {output} \
            -of bigwig \
            --binSize {params.bin} \
            --numberOfProcessors 10 \
            --ignoreDuplicates \
            --scaleFactor 1 \
            --normalizeUsing {params.normalization}
        """

# QC
rule wgsmetrics:
    input:
        bam = "results/UmiDeDup/{sample}-{enrich}.realn.mdup.MQ30.bam",
        ref = ref_fasta 
    output:
        "results/wgsmetrics/{sample}-{enrich}.WgsMetrics.txt"
    message:
        "Computing WGS Metrics\n"
        "Sample: {wildcards.sample}"
    conda:
        "telomereC.py3.1" 
    shell:"""
        gatk --java-options -Xmx6g CollectWgsMetrics \
            -R {input.ref} \
            -I {input.bam} \
            -O {output} \
            --USE_FAST_ALGORITHM false
        """

rule callpeaks1:
    input:
        expand("results/UmiDeDup/{{sample}}-{enrich}.realn.mdup.MQ30.bam",sample=samples, enrich=['input','capture'])
    output:
        bw = "results/RGT_peakCall/{sample}/{sample}.normal.run_signal.bw",
        bed = temp("results/RGT_peakCall/{sample}/{sample}.normal.run_peaks.bed")
    params:
        bam1 = "results/UmiDeDup/{{sample}}-capture.realn.mdup.MQ30.bam",
        bam2 = "results/UmiDeDup/{{sample}}-input.realn.mdup.MQ30.bam",
        peak_caller = "scripts/primary_peak_call.py"
    threads: 4
    resources:
         mem_mb=204800
    message:
        "Calling peaks using RGT_peakCall.\n"
    conda:
        "telomereC.py3.1" 
    shell:"""
    python {params.peak_caller} {params.bam1} {params.bam2}
    """

rule callpeaks2::
    input:
        "results/RGT_peakCall/{sample}/{sample}.normal.run_peaks.bed"
    output:
        trim = temp("results/RGT_peakCall/{sample}/{sample}.normal.run_peaks.trim.bed"),
        merge = "results/RGT_peakCall/{sample}/{sample}.normal.run_peaks.merge.bed"
    message:
        "Merge called peaks if they are close\n"
    params:
        distance = 100
    conda:
        "telomereC.py3.1" 
shell:"""
    awk '{{if ($2<$3 && $1 != "chrM") print $0}}' {input} | sort -k1,1 -k2,2n > {output.trim}
    bedtools merge -d 100 -i  {output.trim} | awk '{{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}}' > {params.merge}
    """

    
