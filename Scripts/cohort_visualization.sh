#!/bin/bash

helpFunction()
{
   echo "Usage: cohort_visualization.sh -s <output svg file> -p <output png file> -r <reference fasta> -v <valid chromosome list>"
   echo -e "\t-s <str> Output SVG file name."
   echo -e "\t-p <str> Output PNG file name."
   echo -e "\t-r <str> Reference Nuclear FASTA file."
   echo -e "\t-v <str> list containing valid chromosome fasta headers"
   exit 1 # Exit script after printing help
}

while getopts "s:p:r:v:" opt
do
        case "$opt" in
        s )
                svg="$OPTARG"
                ;;
        p )
                png="$OPTARG"
                ;;
        r )
                ref_genome="$OPTARG"
                ;;
        v )
                ref_headers="$OPTARG"
                ;;
        ? ) helpFunction ;;
        esac
done

if [ -f "$MT_DATA" ]; then
    
    mt_length=$(cat "$MT_DATA" | awk '/^>/{if (l!="") print; l=0; next}{l+=length($0)}END{print l/2}')

else
    
    echo "No MT_DATA variable found in bashrc. Kindly check and modify the bashrc."

fi

if [ -z "$output" ] && [ -z "$svg" ] && [ -z "$png" ] && [ -z "$ref_genome" ] && [ -z "$ref_headers" ]; then
    
    echo -e "Please Provide the result directory path, SVG file name, PNG file name,\nReference Nuclear Genome FASTA and, List containing valid chromosomal headers."
    
    helpFunction
fi

script_path=$(dirname "$(readlink -f "$0")")

ref_genome_full=$(realpath "$ref_genome")

result_dir_full=$(dirname "$(readlink -f "$svg")")

ref_headers_full=$(realpath "$ref_headers")

ref_genome_index=$(echo "$ref_genome_full" | sed 's/$/.fai/g')

Rscript "$script_path"/circos_numts_all.R -o "$result_dir_full" -s "$svg" -p "$png" -r "$ref_headers_full" -i "$ref_genome_index" -l "$mt_length"