#!/bin/bash

# Script to process supplementary alignments from BAM files
# Usage: ./get_potential_numts_from_sa.sh <input_bam> <intermediate_output> <final_output>

set -euo pipefail

# Input validation
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <input_bam> <intermediate_output> <final_output> <threads> <supporting_reads> <len_cutoff>"
    exit 1
fi

# Input parameters
BAM_FILE="$1"
OUT_FILE="$2"
TARGET_FILE="$3"
THREADS="$4"
SUPPORTING_READS="$5"
LEN_CUTOFF="$6"

# Conda environment setup
CONDA_BASE=$(conda info --base)

source "${CONDA_BASE}/bin/activate" anomaly

# Function to process CIGAR strings
process_cigar() {
    awk '
    BEGIN { FS = "\t"; OFS = "\t" }
    {
        has_d = 0
        for (i = 1; i <= NF; i++) {
            if ($i ~ /D/) {
                has_d = 1
                break
            }
        }
        if (has_d == 0)
            print $0, "0D"
        else
            print $0
    }'
}

# Function to extract match and deletion values
extract_values() {
    awk '
    BEGIN { FS = "\t"; OFS = "\t" }
    {
        split($0, a, "\t")
        match_cigar = ""
        del_cigar = "0"
        for (i = 1; i <= length(a); i++) {
            if (a[i] ~ /M/)
                match_cigar = a[i]
            else if (a[i] ~ /D/)
                del_cigar = a[i]
        }
        sub("M", "", match_cigar)
        sub("D", "", del_cigar)
        print $2, $3 - 1, $3, $1, $4, $5, $5 + match_cigar + del_cigar, match_cigar + del_cigar
    }'
}

# Main processing pipeline

num_line=$(sambamba view "$BAM_FILE" | wc -l)

if [ "$num_line" -gt 2 ]; then
    sambamba view -t "$THREADS" "$BAM_FILE" \
        | awk 'BEGIN {FS=OFS="\t"} {
            for (i=1; i<=NF; i++) {
                if ($i ~ /SA:Z:/) print $1, $3, $4, $i
            }
        }' \
        | sed 's/SA:Z://g' \
        | sed 's/;/\t/g' \
        | awk 'BEGIN {FS=OFS="\t"} {
            split($0,a,"\t")
            for(i=4;i<=length(a);i++) {
                if(a[i]~/MT/) {
                    print $1, $2, $3, a[i]
                }
            }}' \
        | sed 's/,/\t/g' \
        | awk '{gsub(/([A-Z])/, "&\t", $7); print}' \
        | sed 's/ /\t/g' \
        | process_cigar \
        | extract_values \
        | awk -v len_cutoff="$LEN_CUTOFF" 'BEGIN{FS=OFS="\t"}{if($8>=len_cutoff) print $1,$2,$3,$4,$5,$6,$7}' > "$OUT_FILE"

    awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$4}' "$OUT_FILE" | sort -k1 -n | uniq | sortBed \
    | bedtools merge -d 100 -i stdin -c 4 -o count \
    | awk -v read_cutoff="$SUPPORTING_READS" 'BEGIN {FS=OFS="\t"} {if ($4 >= read_cutoff) print $1, $2, $3}' \
    > "$TARGET_FILE"

    echo "Processing complete. Results written to $TARGET_FILE"
else
    echo "MT BAM File is empty."
    touch "$OUT_FILE" "$TARGET_FILE"
fi