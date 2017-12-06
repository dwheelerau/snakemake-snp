from os.path import join
import glob
from datetime import datetime, date, time
configfile: "./config.yaml"

# globals
now = datetime.now().strftime('%Y-%m-%d-%H-%M')
date_string = "{}".format(now)

# from config.yaml file
THREADS = config['THREADS']
GENOMEFTP = config['GENOMEFTP']

# other params
DIRS = ['bams', 'genome', 'logs', 'clean_reads', 'vcf']
SAMPLE_PATH = config['sample_path']
SAMPLE_TAIL = config['sample_tail']
SAMPLES, = glob_wildcards(SAMPLE_PATH)
PATTERN_R1 = config['pat_r1']
PATTERN_R2 = config['pat_r2']

rule all:
    input:
        DIRS,
        "genome/genome.fa",
        "genome/genome.1.bt2",
        expand("logs/{sample}.aln_report.txt", sample=SAMPLES),
        expand("bams/{sample}.sorted.md.bam", sample=SAMPLES),
        expand("vcf/{sample}.vcf.gz.tbi", sample=SAMPLES),
        "vcf/merged.vcf.gz"
        #expand("bams/{sample}.sorted.bam", sample=SAMPLES),
        #"logs/qc_report.txt",
        #expand("clean_reads/{sample}_R1.cln.fastq.gz", sample=SAMPLES),
        #expand("data/{sample}{tail}", sample=SAMPLES, tail=SAMPLE_TAIL),
        #expand("bams/{sample}.bam", sample=SAMPLES),

rule clean:
    shell:
        """
        rm -f genome/genome.*
        rm -f logs/*
        rm -f bams/*
        rm -f clean_reads/*gz
        rm -f vcf/*
        """

rule setup:
    output: DIRS
    shell:
        "mkdir -p "+' '.join(DIRS)

rule get_genome:
    output:
        "genome/genome.fa"
    params:
        genome=GENOMEFTP
    log:
        "logs/genome_download.log"
    shell:
        """
        wget {params.genome} -O {output}.gz -o {log}
        gunzip {output}.gz
        """

rule genome_index:
    input:
        "genome/genome.fa"
    output:
        "genome/genome.1.bt2"
    log:
        "logs/bowtie2-build.log"
    shell:
        "bowtie2-build {input} genome/genome > {log}"

rule sample_qc:
    input:
        r1 = join("data", PATTERN_R1),
        r2 = join("data", PATTERN_R2),
    params:
        minlen=config['minlen'],
        qtrim=config['qtrim'],
        trimq=config['trimq'],
        ktrim=config['ktrim'],
        kwin=config['kwin'],
        mink=config['mink'],
        hdist=config['hdist'],
        adapt=config['adaptors']
    output:
        r1_out = "clean_reads/{sample}_R1.cln.fastq.gz",
        r2_out = "clean_reads/{sample}_R2.cln.fastq.gz",
    log:
        "logs/read_qc.log"
    shell:
        "bbduk.sh in1={input.r1} in2={input.r2} out1={output.r1_out} out2={output.r2_out} "
        "minlen={params.minlen} qtrim={params.qtrim} trimq={params.trimq} "
        "ktrim={params.ktrim} k={params.kwin} mink={params.mink} "
        "ref={params.adapt} hdist={params.hdist} 2>&1 | tee -a {log}"

#rule qc_report:
#    input:
#        "clean_reads/{sample}_R2.cln.fastq.gz"
#    output:
#        "logs/qc_report.txt"
#    shell:
#        "python scripts/report_qc.py logs/read_qc.log | tee {output}"

rule bowtie_aln:
    input:
        r1 = "clean_reads/{sample}_R1.cln.fastq.gz",
        r2 = "clean_reads/{sample}_R2.cln.fastq.gz",
        genome="genome/genome.1.bt2",
    output:
        "bams/{sample}.bam"
    params:
        thr=THREADS,
        mates=config['mates'],
        mode=config['alnmode'],
        sensitivity=config['sensitivity']
    log:
        "logs/mapping.log"
    shell:
        "bowtie2 --rg 'ID:{wildcards.sample}\tSM:{wildcards.sample}' "
        "-p {params.thr} -x genome/genome --met-file {log} {params.mode} "
        "{params.mates} {params.sensitivity} "
        "-1 {input.r1} -2 {input.r2} | samtools view -bS - > {output}"

rule sort_bam:
    input:
        "bams/{sample}.bam"
    output:
        "bams/{sample}.sorted.bam"
    shell:
        "samtools sort -o {output} {input}"

rule aln_report:
    input:
        "bams/{sample}.sorted.bam"
    output:
        aln="logs/{sample}.aln_report.txt",
    shell:
        """
        samtools flagstat {input} > {output.aln}
        echo alignment report for sample {output.aln}
        grep 'mapped' {output.aln}
        python scripts/report_qc.py logs/read_qc.log | tee logs/qc_report.txt
        """

rule mark_duplicates:
    input:
        "bams/{sample}.sorted.bam"
    output:
        bam="bams/{sample}.sorted.md.bam",
        info="bams/{sample}.markduplicates.txt"
    log:
        "logs/markduplicates.log"
    shell:
        "picard MarkDuplicates I={input} "
        "O={output.bam} M={output.info} VALIDATION_STRINGENCY=LENIENT 2> {log}"

rule call_snp:
    input:
        "bams/{sample}.sorted.md.bam"
    output:
        "vcf/{sample}.vcf"
    params:
        ploidy=config['ploidy']
    shell:
        "freebayes -f genome/genome.fa -b {input} -p {params.ploidy} > {output}"

rule compress_vcf:
    input:
        "vcf/{sample}.vcf"
    output:
        "vcf/{sample}.vcf.gz"
    shell:
        "bgzip {input}"

rule index_vcf:
    input:
        "vcf/{sample}.vcf.gz"
    output:
        "vcf/{sample}.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"

rule merge_vcf:
    input:
        expand("vcf/{sample}.vcf.gz.tbi", sample=SAMPLES)
    output:
        "vcf/merged.vcf.gz"
    shell:
        "vcf-merge vcf/*gz > {output}"
