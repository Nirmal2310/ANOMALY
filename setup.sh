#!/bin/bash

if which conda >/dev/null; then
        
        echo "Conda Exist"

else
        source ~/.bashrc
        
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh \
        && chmod +x miniconda.sh && bash miniconda.sh -b -p miniconda
        
        base_dir=$(echo $PWD)
        
        export PATH=$base_dir/miniconda/bin:$PATH
        
        source ~/.bashrc
        
        echo -e "$base_dir/miniconda/etc/profile.d/conda.sh" >> ~/.profile
        
        conda init bash

fi

path=$(conda info --base)

base_dir=$PWD

script_dir=$(dirname "$(readlink -f "$0")")

data_dir=$(echo "$script_dir" | sed 's/$/\/DATA/g')

## Installing Required Tools for the Pipeline

if { conda env list |  grep -w "snakemake"; } > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name snakemake --file Envs/snakemake.txt

fi

if { conda env list |  grep -w "minimap2"; } > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name minimap2 --file Envs/minimap2.txt

fi

if { conda env list |  grep -w "sniffles"; } > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name sniffles --file Envs/sniffles.txt

fi

if { conda env list |  grep -w "blast"; } > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name blast --file Envs/blast.txt

fi

if { conda env list |  grep -w "bedtools"; } > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name bedtools --file Envs/bedtools.txt

fi

if { conda env list |  grep -w "renv"; } > /dev/null 2>&1; then

        echo "Environment Exist"

else

        conda create --name renv --file Envs/renv.txt

fi

## Creating the BLASTn Index For Mitochondrial Genome

if [ ! -d "$data_dir" ]; then

 	mkdir "$data_dir"

	wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz -O "$data_dir"/Homo_sapiens.GRCh38.MT.fa.gz

	gunzip "$data_dir/"Homo_sapiens.GRCh38.MT.fa.gz

 	source $path/bin/activate blast

  	seqkit concat -w 0 -t DNA "$data_dir"/Homo_sapiens.GRCh38.MT.fa "$data_dir"/Homo_sapiens.GRCh38.MT.fa > "$data_dir"/Homo_sapiens.GRCh38.MT_concatenated.fa

	makeblastdb -in "$data_dir"/Homo_sapiens.GRCh38.MT_concatenated.fa -parse_seqids -blastdb_version 5 -title Homo_sapiens.GRCh38.MT_concatenated.fa -dbtype nucl -out "$data_dir"/Homo_sapiens.GRCh38.MT_concatenated.fa

	source ~/.bashrc

	cd $base_dir

	source $path/bin/activate base
	
fi

grep -qF "export MT_DATA=\"$data_dir/Homo_sapiens.GRCh38.MT_concatenated.fa\"" ~/.bashrc || echo "export MT_DATA=\"$data_dir/Homo_sapiens.GRCh38.MT_concatenated.fa\"" >> ~/.bashrc && source ~/.bashrc
