#!/bin/bash

helpFunction()
{
     echo "Usage: bash get_insertion_calls.sh -i <input.gz> -o <output.txt> -r <ref_headers.txt> -t <threads>"
     echo -e "\t-i <file> Input VCF file."
     echo -e "\t-o <file> Output text file."
     echo -e "\t-r <file> Reference headers file."
     echo -e "\t-t <threads> Number of threads."
     exit 1
}

while getopts "i:o:r:t:" opt
do
     case "$opt" in
     i )
          input="$OPTARG"
          ;;
     o )
          output="$OPTARG"
          ;;
     r )
          chr_list="$OPTARG"
          ;;
     t )
          threads="$OPTARG"
          ;;
     ? )  helpFunction ;;
     esac
done

if [ ! -f "$input" ] && [ ! -f "$output" ] && [ ! -f "$chr_list" ]
     then
     echo "Please provide the required arguments"
     helpFunction
fi

CONDA_BASE=$(conda info --base)

source "${CONDA_BASE}/bin/activate" anomaly

tmp_insertion_file=$(mktemp)

tmp_support_file=$(mktemp)

process_region() {
     local chr="$1"
     local start="$2"
     local end="$3"
     local support="$4"
     local input="$5"

     awk -v chr="$chr" -v start="$start" -v end="$end" -v support="$support" '
        BEGIN{FS="\t";OFS="\t"}
        $1==chr && $2>=start && $3<=end && $6==support { print $0 }
    ' "$input"
}

zcat "$input" | \
    grep -v "#" | \
    grep "INS" | \
    awk 'BEGIN { FS="\t"; OFS="\t" }
         { print $1, $2, $5, $8 }' | \
    awk 'BEGIN { FS="\t"; OFS="\t" }
         { 
           split($4, a, ";"); 
           print $1, $2, $2+1, $3, a[3], a[5] 
         }' | \
    sed 's/SVLEN=//g;s/SUPPORT=//g' | \
     awk -F "\t" '{ if($4 != "<INS>") print }' | \
     grep -Ff $chr_list - > "$tmp_insertion_file"

awk 'BEGIN{FS="\t";OFS="\t"}{
     print $1,$2,$3,$4,$5,$6
}' "$tmp_insertion_file" | sortBed | bedtools merge -i stdin -d 5 -c 6 -o max > "$tmp_support_file"

export -f process_region

while IFS=$'\t' read -r chr start end support; do
     echo "process_region $chr $start $end $support $tmp_insertion_file"
done < "$tmp_support_file" | parallel -j $threads {} | \
awk 'BEGIN{FS="\t";OFS="\t"}{if($4!="<INS>") print $1, $2, $4, $5}' > "$output"

rm -r "$tmp_insertion_file" "$tmp_support_file"