#!/bin/bash

bam_file=$1

out_nuclear_bam=$2

out_mt_bam=$3

threads=$4

chr_list=$5

path=$(conda info --base)

source $path/bin/activate anomaly

MT_HEADER=$(cat $MT_DATA | grep ">" | sed 's/>//g')

sambamba view -t $threads -F "supplementary and ref_name!='$MT_HEADER'" $bam_file | \
    rg -j $threads -P 'SA:Z:(?:[^;]*;)*MT,' | \
	cat <(samtools view -H $bam_file) - | \
	samtools view -@ $threads -bS | samtools sort -@ $threads -o $out_nuclear_bam && \
	samtools index -@ $threads $out_nuclear_bam

sambamba view -t $threads -F "supplementary" -h $bam_file MT | samtools sort -@ $threads -o $out_mt_bam

source $path/bin/activate base