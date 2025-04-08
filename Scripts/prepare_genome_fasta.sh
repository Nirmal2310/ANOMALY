#!/bin/bash

input_fasta=$1

path=$(conda info --base)

cwd=$PWD

script_dir=$(dirname "$(readlink -f "$0")")

data_dir=$(echo "$script_dir" | sed 's/Scripts/DATA/g')

fasta_name=$(basename "$input_fasta")

fasta_concat=$(echo $fasta_name | sed 's/.fa/_concatenated.fa/g;s/.fasta/_concatenated.fa/g;s/.fna/_concatenated.fa/g')

cp "$input_fasta" "$data_dir"/"$fasta_name"

source $path/bin/activate blast

cd "$data_dir"

sed "s/>.*$/>MT/g" "$fasta_name" > temp && mv temp "$fasta_name"

seqkit concat -w 0 -t dna "$fasta_name" "$fasta_name" > "$fasta_concat"

makeblastdb -in "$fasta_concat" -parse_seqids -blastdb_version 5 -title "$fasta_concat" -dbtype nucl -out "$fasta_concat"

source ~/.bashrc

source $path/bin/activate base

grep  -qF "export MT_DATA=\"$data_dir/$fasta_concat\"" ~/.bashrc || echo "export MT_DATA=\"$data_dir/$fasta_concat\"" >> ~/.bashrc && source ~/.bashrc

cd $cwd
