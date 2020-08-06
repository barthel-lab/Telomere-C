rule fq2ubam:
    input:
        R1 = "data/fastq/A2780-0_5M_GT20-05417_CAATTAAC-CGAGATAT_S1_R1_001.fastq.gz",
        R2 = "data/fastq/A2780-0_5M_GT20-05417_CAATTAAC-CGAGATAT_S1_R2_001.fastq.gz"
    output:
        "results/align/ubam/{aliquot_barcode}.{readgroup}.unaligned.bam"
    params:
        RGID = "G5B3N.1",
        RGPL = "ILLUMINA",
        RGPU = "G5B3NM02838.1",
        RGLB = "BAR65723",
        RGDT = "20200321",
        RGSM = "A2780-GT20",
        RGCN = "JAX",
        mem = 6
    threads:
        1
    log:
        "logs/align/fq2ubam/{aliquot_barcode}.{readgroup}.log"
    benchmark:
        "benchmarks/align/fq2ubam/{aliquot_barcode}.{readgroup}.txt"
    message:
        "Converting FASTQ file to uBAM format\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "gatk --java-options -Xmx{params.mem}g FastqToSam \
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
            --TMP_DIR={config[tempdir]} \
            > {log} 2>&1"

rule markadapters:
    input:
        "results/align/ubam/A2780-GT20.G5B3N.1.unaligned.bam"
    output:
        bam = "results/align/markadapters/{aliquot_barcode}.{readgroup}.markadapters.bam",
        metric = "results/align/markadapters/{aliquot_barcode}.{readgroup}.markadapters.metrics.txt"
    threads:
        1
    params:
        mem = 6
    log: 
        "logs/align/markadapters/{aliquot_barcode}.{readgroup}.log"
    benchmark:
        "benchmarks/align/markadapters/{aliquot_barcode}.{readgroup}.txt"
    message:
        "Adding XT tags. This marks Illumina Adapters and allows them to be removed in later steps\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "gatk --java-options -Xmx{params.mem}g MarkIlluminaAdapters \
            --INPUT={input} \
            --OUTPUT={output.bam} \
            --METRICS={output.metric} \
            --TMP_DIR={config[tempdir]} \
            > {log} 2>&1"

rule clipreads:
    input:
        ubam = "results/align/markadapters/{aliquot_barcode}.{readgroup}.markadapters.bam",
        clipseq = "data/ref/telomerec.fasta"
    output:
        ubam = "results/align/clipreads/{aliquot_barcode}.{readgroup}.clipreads.bam",
        stats = "results/align/clipreads/{aliquot_barcode}.{readgroup}.clipreads.metrics.txt"
    params:
        mem = 6
    threads:
        1
    log: 
        "logs/align/markadapters/{aliquot_barcode}.{readgroup}.log"
    benchmark:
        "benchmarks/align/markadapters/{aliquot_barcode}.{readgroup}.txt"
    message:
        "Clip reads. This clips the 4C PCR primers.\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "gatk --java-options -Xmx{params.mem}g ClipReads \
            -I {input.ubam} \
            -O {output.ubam} \
            --clip-sequences-file {input.clipseq} \
            --output-statistics {output.stats} \
            > {log} 2>&1"

rule samtofastq_bwa_mergebamalignment:
    input:
        bam = "results/align/clipreads/{aliquot_barcode}.{readgroup}.clipreads.bam"
    output:
        bam = "results/align/bwa/{aliquot_barcode}.{readgroup}.aln.bam",
        bai = "results/align/bwa/{aliquot_barcode}.{readgroup}.aln.bai"
    threads:
        1
    params:
        mem = 6
    log: 
        "logs/align/samtofastq_bwa_mergebamalignment/{aliquot_barcode}.{readgroup}.log"
    benchmark:
        "benchmarks/align/revertsam/{aliquot_barcode}.{readgroup}.txt"
    message:
        "BAM to FASTQ --> BWA-MEM --> Merge BAM Alignment.\n"
        "The first step converts the reverted BAM to an interleaved FASTQ, removing Illumina "
        "adapters. The output is then piped to BWA-MEM and aligned. Aligned reads are merged "
        "with the original pre-aligned BAM to preserve original metadata, including read groups.\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "module load bwa/0.7.9a; \
         gatk --java-options {config[samtofastq_java_opt]} SamToFastq \
            --INPUT={input.bam} \
            --FASTQ=/dev/stdout \
            --CLIPPING_ATTRIBUTE=XT \
            --CLIPPING_ACTION=2 \
            --INTERLEAVE=true \
            --NON_PF=true \
            --TMP_DIR={config[tempdir]} | \
         bwa mem -M -t {threads} -p {config[reference_fasta]} /dev/stdin | \
         gatk --java-options {config[mergebamalignment_java_opt]} MergeBamAlignment \
            --ALIGNED_BAM=/dev/stdin \
            --UNMAPPED_BAM={input.bam} \
            --OUTPUT={output.bam} \
            --REFERENCE_SEQUENCE={config[reference_fasta]} \
            --CREATE_INDEX=true \
            --ADD_MATE_CIGAR=true \
            --CLIP_ADAPTERS=false \
            --CLIP_OVERLAPPING_READS=true \
            --INCLUDE_SECONDARY_ALIGNMENTS=true \
            --MAX_INSERTIONS_OR_DELETIONS=-1 \
            --PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
            --ATTRIBUTES_TO_RETAIN=XS \
            --TMP_DIR={config[tempdir]} \
            > {log} 2>&1"

rule markduplicates:
    input:
        "results/align/bwa/A2780-GT20.G5B3N.1.aln.bam"
    output:
        bam = temp("results/align/markduplicates/{aliquot_barcode}.realn.mdup.bam"),
        bai = temp("results/align/markduplicates/{aliquot_barcode}.realn.mdup.bai"),
        metrics = "results/align/markduplicates/{aliquot_barcode}.metrics.txt"
    params:
        max_records = 6000000,
        mem = 6
    threads:
        1
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
        shell("gatk --java-options -Xmx{params.mem}g MarkDuplicates \
            {multi_input} \
            --OUTPUT={output.bam} \
            --METRICS_FILE={output.metrics} \
            --CREATE_INDEX=true \
            --TMP_DIR={config[tempdir]} \
            --MAX_RECORDS_IN_RAM={params.max_records} \
            > {log} 2>&1")

rule callpeaks:
    input:
        "results/align/markduplicates/{aliquot_barcode}.realn.mdup.bam"
    output:
        peaks = "results/align/macs2/{aliquot_barcode}/{aliquot_barcode}_peaks.xls"
    params:
        outdir = "results/align/macs2/{aliquot_barcode}/",
        mem = 6
    threads:
        1
    log:
        "logs/align/macs2/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/align/macs2/{aliquot_barcode}.txt"
    message:
        "Calling peaks using MACS2.\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "module load MACS; \
         macs2 callpeak \
            -t {input} \
            --outdir {params.outdir} \
            -f BAM \
            -g hs \
            -n {wildcards.aliquot_barcode} \
            -B \
            -q 0.01 \
            > {log} 2>&1"   

rule all:
    input: "results/align/macs2/A2780-GT20/A2780-GT20_peaks.xls"

## END ##