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

if { conda env list |  grep -w "minimap2"; } > /dev/null 2>&1; then

        conda list -n minimap2 --explicit > _current_env.txt

        if diff -q $script_dir/Envs/minimap2.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name minimap2 --file $script_dir/Envs/minimap2.txt -y && rm -r _current_env.txt
                conda list -n minimap2 --explicit > $script_dir/Envs/minimap2.txt
        fi

else

        conda create --name minimap2 --file $script_dir/Envs/minimap2.txt
        conda list -n minimap2 --explicit > $script_dir/Envs/minimap2.txt

fi

if { conda env list |  grep -w "sniffles"; } > /dev/null 2>&1; then

        conda list -n sniffles --explicit > _current_env.txt

        if diff -q $script_dir/Envs/sniffles.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name sniffles --file $script_dir/Envs/sniffles.txt -y && rm -r _current_env.txt
                conda list -n sniffles --explicit > $script_dir/Envs/sniffles.txt
        fi

else

        conda create --name sniffles --file $script_dir/Envs/sniffles.txt
        conda list -n sniffles --explicit > $script_dir/Envs/sniffles.txt

fi

if { conda env list |  grep -w "blast"; } > /dev/null 2>&1; then

        conda list -n blast --explicit > _current_env.txt

        if diff -q $script_dir/Envs/blast.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name blast --file $script_dir/Envs/blast.txt -y && rm -r _current_env.txt
                conda list -n blast --explicit > $script_dir/Envs/blast.txt
        fi

else

        conda create --name blast --file $script_dir/Envs/blast.txt
        conda list -n blast --explicit > $script_dir/Envs/blast.txt

fi

if { conda env list |  grep -w "bedtools"; } > /dev/null 2>&1; then

        conda list -n bedtools --explicit > _current_env.txt

        if diff -q $script_dir/Envs/bedtools.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name bedtools --file $script_dir/Envs/bedtools.txt -y && rm -r _current_env.txt
                conda list -n bedtools --explicit > $script_dir/Envs/bedtools.txt
        fi

else

        conda create --name bedtools --file $script_dir/Envs/bedtools.txt
        conda list -n bedtools --explicit > $script_dir/Envs/bedtools.txt

fi

if { conda env list |  grep -w "renv"; } > /dev/null 2>&1; then

        conda list -n renv --explicit > _current_env.txt

        if diff -q $script_dir/Envs/renv.txt _current_env.txt > /dev/null; then
                echo "Environment exists and up to date." && rm -r _current_env.txt
        else
                conda create --name renv --file $script_dir/Envs/renv.txt -y && rm -r _current_env.txt
                conda list -n renv --explicit > $script_dir/Envs/renv.txt
        fi

else

        conda create --name renv --file $script_dir/Envs/renv.txt
        conda list -n renv --explicit > $script_dir/Envs/renv.txt

fi

## Creating the BLASTn Index For Mitochondrial Genome

if [ ! -d "$data_dir" ]; then

 	mkdir "$data_dir"

	wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz -O "$data_dir"/Homo_sapiens.GRCh38.MT.fasta.gz

	gunzip "$data_dir/"Homo_sapiens.GRCh38.MT.fasta.gz

 	source $path/bin/activate blast

  	seqkit concat -w 0 -t DNA "$data_dir"/Homo_sapiens.GRCh38.MT.fasta "$data_dir"/Homo_sapiens.GRCh38.MT.fasta > "$data_dir"/Homo_sapiens.GRCh38.MT_concatenated.fasta

	makeblastdb -in "$data_dir"/Homo_sapiens.GRCh38.MT_concatenated.fasta -parse_seqids -blastdb_version 5 -title Homo_sapiens.GRCh38.MT_concatenated.fasta -dbtype nucl -out "$data_dir"/Homo_sapiens.GRCh38.MT_concatenated.fasta

	source ~/.bashrc

	cd $base_dir

	source $path/bin/activate base
	
fi

grep -qF "export MT_DATA=\"$data_dir/Homo_sapiens.GRCh38.MT_concatenated.fasta\"" ~/.bashrc || echo "export MT_DATA=\"$data_dir/Homo_sapiens.GRCh38.MT_concatenated.fasta\"" >> ~/.bashrc && source ~/.bashrc
