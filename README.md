# Snakemake pipeline for snp calling using freebayes

## requires  
sudo apt-get install vcftools  

## clone workflow into working directory
git clone https://github.com/dwheelerau/snakemake-snp
cd snakemake-snp

## edit config and workflow as needed
vim config.yaml

## install dependencies into isolated environment
conda env create -n snp --file environment.yml

## activate environment
source activate snp

## execute workflow
snakemake all --cores 24
