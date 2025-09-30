#!/bin/bash

# Set strict error handling
set -euo pipefail

# Constants

WORD_SIZE=7

# Function to display usage information
usage() {
    cat << EOF
Usage: $(basename "$0") <input_fasta> <output_file> <threads> <chr_list>

Process FASTA file against MT database using BLAST and filter results.

Arguments:
    input_fasta  - Input FASTA file
    output_file  - Output file for filtered BLAST results
    threads      - Number of threads to use for BLAST
    chr_list     - List containing chromosome reference headers

Requirements:
    - BLAST+ must be installed
    - MT_DATA environment variable must be set in ~/.bashrc
    - List containing chromosome reference headers
EOF
    exit 1
}

# Function to validate inputs
validate_inputs() {
    local input=$1
    local output=$2
    local threads=$3
    local chr_list=$4
    local final_evalue_cutoff=$5

    # Check if the input file exists
    
    if [ ! -f "$input" ]; then
        echo "Error: Input FASTA file does not exist: $input"
        exit 1
    fi

    if [ ! -f "$chr_list" ]; then
        echo "Error: Input Chromosome List file does not exist: $chr_list"
        exit 1
    fi
    
    # Validate thread count
    if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
        echo "Error: Thread count must be a positive integer"
        exit 1
    fi
}

# Function to get MT database path


MT_DATA=$(grep MT_DATA ~/.bashrc | tail -n 1 | sed 's/export MT_DATA="//;s/"//g')


# Function to run BLAST and process results
run_blast_pipeline() {
    local input=$1
    local output=$2
    local threads=$3
    local chr_list=$4
    local final_evalue_cutoff=$5

    # Run BLAST and filter results
    blastn -task blastn \
           -word_size "$WORD_SIZE" \
           -db "$MT_DATA" \
           -query "$input" \
           -num_threads "$threads" -max_target_seqs 6 -max_hsps 6 \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue" | \
    awk -v evalue_cutoff="$final_evalue_cutoff" '
    BEGIN {
        FS=OFS="\t"
    }
    $11 < evalue_cutoff' | \
    awk '
    BEGIN {
        FS=OFS="\t"
    }
    {
        split($1, a, "_")
        print a[1],a[2]-1,a[2],$2,$9,$10,a[4],$4
    }' | \
    awk '
    BEGIN {
        FS=OFS="\t"
    }
    {
        if($5>$6) {
            print $1,$2,$3,$4,$6,$5,$7,$8
        } else {
            print $0
        }
    }' | awk -v ref_list="$chr_list" '
    BEGIN {
        # Read the list line by line.
        while ((getline line < ref_list) > 0) {
            keys[line] = 1
        }
        close(ref_list)
    }
    $1 in keys {
    print
    }' > "$output"
}

# Main script execution
main() {
    # Check for the correct number of arguments
    if [ $# -ne 5 ]; then
        usage
    fi

    local input_fasta=$1
    local out_data=$2
    local threads=$3
    local chr_list=$4
    local final_evalue_cutoff=$5

    # Validate inputs
    validate_inputs "$input_fasta" "$out_data" "$threads" "$chr_list" "$final_evalue_cutoff"

    # Run BLAST pipeline
    run_blast_pipeline "$input_fasta" "$out_data" "$threads" "$chr_list" "$final_evalue_cutoff"
}

# Execute main function with all arguments
main "$@"
