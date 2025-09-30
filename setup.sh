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

if { conda env list |  grep -w "anomaly"; } > /dev/null 2>&1; then

        conda list -n anomaly --explicit > _current_env.txt

        if diff -q $script_dir/Envs/anomaly.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name anomaly --file $script_dir/Envs/anomaly.txt -y && rm -r _current_env.txt
                conda list -n anomaly --explicit > $script_dir/Envs/anomaly.txt
        fi

else

        conda create --name anomaly --file $script_dir/Envs/anomaly.txt
        conda list -n anomaly --explicit > $script_dir/Envs/anomaly.txt

fi

if { conda env list |  grep -w "snakemake"; } > /dev/null 2>&1; then

        conda list -n snakemake --explicit > _current_env.txt

        if diff -q $script_dir/Envs/snakemake.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name snakemake --file $script_dir/Envs/snakemake.txt -y && rm -r _current_env.txt
                conda list -n snakemake --explicit > $script_dir/Envs/snakemake.txt
        fi

else

        conda create --name snakemake --file $script_dir/Envs/snakemake.txt
        conda list -n snakemake --explicit > $script_dir/Envs/snakemake.txt

fi

## Creating the BLASTn Index For Mitochondrial Genome

if [ ! -d "$data_dir" ]; then

 	mkdir "$data_dir"

	wget https://www.ebi.ac.uk/ena/browser/api/fasta/CP068254.1 -O "$data_dir"/Homo_sapiens.CHM13.MT.fasta

 	source $path/bin/activate blast

  	seqkit concat -w 0 -t DNA "$data_dir"/Homo_sapiens.CHM13.MT.fasta "$data_dir"/Homo_sapiens.CHM13.MT.fasta > "$data_dir"/Homo_sapiens.CHM13.MT_concatenated.fasta

	makeblastdb -in "$data_dir"/Homo_sapiens.CHM13.MT_concatenated.fasta -parse_seqids -blastdb_version 5 -title Homo_sapiens.CHM13.MT_concatenated.fasta -dbtype nucl -out "$data_dir"/Homo_sapiens.CHM13.MT_concatenated.fasta

	source ~/.bashrc

	cd $base_dir

	source $path/bin/activate base
	
fi

grep -qF "export MT_DATA=\"$data_dir/Homo_sapiens.CHM13.MT_concatenated.fasta\"" ~/.bashrc || echo "export MT_DATA=\"$data_dir/Homo_sapiens.CHM13.MT_concatenated.fasta\"" >> ~/.bashrc && source ~/.bashrc
