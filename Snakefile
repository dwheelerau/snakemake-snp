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
DIRS = ['bams', 'genome', 'logs']

rule all:
    input:
        DIRS,
        "genome/genome.fa",
        "genome/genome.1.bt2"

rule clean:
    shell:
        """
        rm -f genome/genome.*
        rm -f logs/*
        rm -f bams/*bam
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
