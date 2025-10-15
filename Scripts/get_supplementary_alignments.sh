#!/bin/bash

bam_file=$1

out_bam=$2

threads=$3

chr_list=$4

path=$(conda info --base)

source $path/bin/activate anomaly

MT_HEADER=$(cat $MT_DATA | grep ">" | sed 's/>//g')

sambamba view -t $threads -F "ref_name!='$MT_HEADER'" $bam_file | \
    rg -j $threads -P 'SA:Z:(?:[^;]*;)*MT,' | \
	cat <(samtools view -H $bam_file) - | \
	samtools view -@ $threads -bS | samtools sort -@ $threads -o $out_bam && \
	samtools index -@ $threads $out_bam

source $path/bin/activate base