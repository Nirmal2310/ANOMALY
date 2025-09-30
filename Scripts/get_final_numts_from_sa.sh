#!/bin/bash

# Script to process genomic coordinates and find overlapping regions
# Usage: ./script.sh <input_file> <output_file> <potential_numts>

set -euo pipefail

# Input validation
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_file> <output_file> <potential_numts>"
    exit 1
fi

# Input parameters
input_file="$1"

output_file="$2"

potential_numts="$3"

tmp_file=$(mktemp)

# Ensure the temporary file is cleaned up on the script exit
trap 'rm -f "$tmp_file"' EXIT

# Set up conda environment
CONDA_BASE=$(conda info --base)

source "${CONDA_BASE}/bin/activate" anomaly

# Function to process coordinates and find overlapping regions
process_coordinates() {
    local chr="$1"
    local start="$2"
    local end="$3"
    
    awk -v chr="$chr" -v start="$start" -v end="$end" '
        BEGIN { FS = "\t"; OFS = "\t" }
        {
            if ($1 == chr && $2 >= start && $3 <= end)
                print $5, $6, $7, chr, start
        }
    ' "$input_file" \
    | awk '
        BEGIN { FS = OFS = "\t" }
        {
            if ($2 > $3)
                print $1, $3, $2, $4, $5
            else
                print $0
        }
    ' \
    | sortBed \
    | bedtools merge -i stdin -d 50 -c 4,5 -o distinct,distinct \
    | awk '
        BEGIN { FS = "\t"; OFS = "\t" }
        { print $4, $5, $1, $2 + 1, $3 - 1 }
    '
}

# Function to process output based on the number of lines
process_output() {
    local line_count="$1"
    
    if [ "$line_count" -gt 1 ]; then
        awk '
            BEGIN { FS = "\t"; OFS = "\t" }
            {
                if (max == "" || $4 > max) max = $4
                if (min == "" || $5 < min) min = $5
            }
            END { print $1, $2, $2+1, $3, min, max, "With Control" }
        ' "$tmp_file"
    else
        awk 'BEGIN{FS="\t";OFS="\t"}{print $1, $2, $2+1, $3, $4, $5, "Without Control"}' "$tmp_file"
    fi
}

# Main processing loop

num_line=$(cat "$input_file" | wc -l)

if [ "$num_line" -gt 0 ]; then
    while read -r line; do
        # Extract fields from the current line
        IFS=$'\t' read -r chr start end <<< "$line"
        
        # Process coordinates and save to a temporary file
        process_coordinates "$chr" "$start" "$end" > "$tmp_file"
        
        # Get line count and process output accordingly
        line_count=$(wc -l < "$tmp_file")
        process_output "$line_count"
        
    done < "$potential_numts" > "$output_file"

    echo "Processing complete. Results are written to $output_file"
else
    echo "Input File is Empty."
    touch "$output_file"
fi
