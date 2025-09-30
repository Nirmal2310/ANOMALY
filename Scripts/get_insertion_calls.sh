#!/bin/bash
# Usage: ./script.sh <input.gz> <output.txt> <ref_headers.txt>

# Check for the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input.gz> <output.txt> <ref_headers.txt>"
    exit 1
fi

CONDA_BASE=$(conda info --base)

source "${CONDA_BASE}/bin/activate" anomaly

input="$1"

output="$2"

chr_list="$3"

tmp_insertion_file=$(mktemp)

tmp_support_file=$(mktemp)


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
    awk -v ref_list="$chr_list" '
    BEGIN {
          while ((getline line < ref_list) > 0) {
               keys[line] = 1
          }
          close(ref_list)
    }
    $1 in keys {
     print }' | awk -F "\t" '{ if($3 != "<INS>") print }' > "$tmp_insertion_file"

awk 'BEGIN{FS="\t";OFS="\t"}{
     print $1,$2,$3,$4,$5,$6
}' "$tmp_insertion_file" | sortBed | bedtools merge -i stdin -d 5 -c 6 -o max > "$tmp_support_file"

while IFS=$'\t' read -r chr start end support; do
     awk -v chr="$chr" -v start="$start" -v end="$end" -v support="$support" '
     BEGIN{FS="\t";OFS="\t"}
     $1==chr && $2>=start && $3<=end && $6==support {
          print $0
     }' "$tmp_insertion_file"
done < "$tmp_support_file" | awk 'BEGIN{FS="\t";OFS="\t"}{if($4!="<INS>") print $1, $2, $4, $5}' > "$output"

rm -r "$tmp_insertion_file" "$tmp_support_file"
