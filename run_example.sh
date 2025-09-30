#!/bin/bash

eval "$(conda shell.bash hook)"

if [ ! -d Example ]; then
        mkdir Example
fi

wget -c https://zenodo.org/records/17065234/files/simulated_data_1.fastq.gz -O Example/simulated_data_1.fastq.gz

wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz -O Example/genome.fasta.gz

gunzip Example/genome.fasta.gz

grep ">" Example/genome.fasta | sed 's/>//g;s/ .*$//g' | paste -d "\t" - <(for x in {1..22..1}; do echo $x; done | cat - <(echo -e "X\nY\nMT")) > header_pairs.txt

conda activate anomaly

seqkit replace -p '^(\S+)(.+?)$' -r '{kv}' -k header_pairs.txt Example/genome.fasta > temp && mv temp Example/genome.fasta

conda activate base

if [ ! -f Example/chr_list.txt ]; then
        for x in {1..22..1}; do echo $x; done | cat - <(echo -e "X\nY") > Example/chr_list.txt
fi

bash get_config.sh -d f -r $PWD/Example/genome.fasta -m 24 -s 24 -p ONT -i $PWD/Example -o $PWD/Example -l $PWD/Example/chr_list.txt

conda activate anomaly

python main.py -c $PWD/Example/snake_config.yml

bash run_snakemake.sh -w $PWD/Example -t 24
