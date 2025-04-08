#!/bin/bash

eval "$(conda shell.bash hook)"

if [ ! -d Example ]; then
        mkdir Example
fi

wget -c https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR270/031/SRR27010831/SRR27010831_1.fastq.gz -O Example/HG001.fastq.gz

wget -c https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O Example/genome.fasta.gz

gunzip Example/genome.fasta.gz

if [ ! -f Example/chr_list.txt ]; then
        for x in {1..22..1}; do echo $x; done | cat - <(echo "X") > Example/chr_list.txt
fi

bash get_config.sh -d f -r $PWD/Example/genome.fasta -m 24 -s 24 -p ONT -i $PWD/Example -o $PWD/Example -l $PWD/Example/chr_list.txt

conda activate snakemake

python main.py -c $PWD/Example/snake_config.yml

bash run_snakemake.sh -w $PWD/Example -t 24
