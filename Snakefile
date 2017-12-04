from os.path import join
import glob
from datetime import datetime, date, time
configfile: "./config.yaml"

# globals
now = datetime.now().strftime('%Y-%m-%d-%H-%M')
date_string = "{}".format(now)

# from config.yaml file
THREADS = config['THREADS']
MINOVLEN = config['MINOVLEN']
MAXDIFF = config['MAXDIFF']
MAXEE = config['MAXEE']
MINLEN = config['MINLEN']
MAXLEN = config['MAXLEN']
MAXNS = config['MAXNS']
MINQ = config['MINQ'] # new
CLUSTERMODE = config['CLUSTERMODE']
CLUSTERID = config['CLUSTERID']
REF_DB = config['REF_DB']
TAX_REF = config['TAX_REF']
SAM_DEPTH = config['SAM_DEPTH']
MAP_FILE = config['MAP_FILE']
TAX_METH = config['TAX_METH']

rule all:
    input:
        "sample.otu_table.txt",
        "sample.rename.otu_table.txt",
        "fastqc/",
        "tax_summary_proportion/",
        "sample.rename.otu_table_summary.txt",
        "all.otus.tre",
        "coreout/"

rule cluster_otus:
    input:
        "scripts/run_otu.sh"
    output:
        "sample.otu_table.txt", "all.otus.fasta"
    log:
        stdout=expand("run_otu.{date_string}.stdout",date_string=date_string),
        stderr=expand("run_otu.{date_string}.stderr",date_string=date_string)
    params:
        threads=THREADS,
        minovlen=MINOVLEN,
        maxdiff=MAXDIFF,
        maxee=MAXEE,
        minlen=MINLEN,
        maxlen=MAXLEN,
        maxns=MAXNS,
        minq=MINQ,
        clustermode=CLUSTERMODE,
        clusterid=CLUSTERID,
        ref_db=REF_DB
    shell:
        "scripts/run_otu.sh "
        "-t {params.threads} -m {params.minovlen} -d {params.maxdiff} "
        "-e {params.maxee} -i {params.minlen} -l {params.maxlen} "
        "-n {params.maxns} -q {params.minq} -c '{params.clustermode}' -r {params.clusterid} "
        "-b '{params.ref_db}' > >(tee -a {log.stdout}) 2> >(tee -a {log.stderr} >&2)"

