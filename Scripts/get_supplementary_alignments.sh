#!/bin/bash

bam_file=$1

out_bam=$2

threads=$3

chr_list=$4

path=$(conda info --base)

source $path/bin/activate anomaly

sambamba view -t $threads $bam_file | \
    rg -j $threads -P 'SA:Z:(?:[^;]*;)*MT,' | \
    awk -v ref_list="$chr_list" '
    BEGIN {
    	while ((getline line < ref_list) > 0) {
     		keys[line] = 1
		}
 		close(ref_list)
  	}
  	$3 in keys {
   		print
    }' | cat <(samtools view -H $bam_file) - | \
	samtools view -@ $threads -bS | samtools sort -@ $threads -o $out_bam && \
	samtools index -@ $threads $out_bam

source $path/bin/activate base
