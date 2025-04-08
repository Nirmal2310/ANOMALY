#!/bin/bash

# Script to process and filter BED format genomic data
# Usage: ./numt_concat.sh <input_bed> <output_bed>

set -euo pipefail

# Input validation
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_bed> <output_bed>"
    exit 1
fi

input="$1"

output="$2"

tmp_file=$(mktemp)

if [ -f "$MT_DATA" ]; then
    mt_length=$(cat "$MT_DATA" | awk '/^>/{if (l!="") print; l=0; next}{l+=length($0)}END{print l/2}')
else
    echo "No MT_DATA variable found in bashrc. Kindly modify the bashrc."
fi


get_freq() {
    sortBed -i "$1" \
    | bedtools merge -i stdin -d 10 -c 2 -o count
}

process_freq_6() {
    local in="$1"
    local chr="$2"
    local start="$3"
    local end="$4"
    local mt_length="$5"

    # Filter by chromosome and position
    awk -v chr="$chr" -v start="$start" -v end="$end" '
    BEGIN {FS=OFS="\t"}
    $1==chr && $2>=start && $3<=end {
        print $0
    }' "$in" | \
    awk -v mt_len="$mt_length" '
    BEGIN{FS=OFS="\t"}
    {
        if($5>mt_len && $6>mt_len) {
            print $1,$2,$3,$4,$5-mt_len,$6-mt_len,$7,$8
        } else {
            print $0
        }
    }' | sort -k8 -n | uniq | \
    awk 'BEGIN{FS=OFS="\t"} {
        lines[NR] = $0
        if($2 < chr_min || NR==1) chr_min = $2
        if($3 < chr_max || NR==1) chr_max = $3
        if($5 < mt_min || NR==1) mt_min = $5
        if($6 > mt_max || NR==1) mt_max = $6
        if($7 > query_max || NR==1) query_max = $7
        if($8 > ref_max || NR==1) ref_max = $8
    } END{
        if(NR==3) {
            print $1,chr_min,chr_max,$4,mt_min,mt_max,query_max,ref_max
        } else if(NR==6) {
            print lines[NR]
        }
    }'
}

process_freq_4() {
    local in="$1"
    local chr="$2"
    local start="$3"
    local end="$4"
    local mt_length="$5"

    # Filter by chromosome and position
    awk -v chr="$chr" -v start="$start" -v end="$end" '
    BEGIN {FS=OFS="\t"}
    $1==chr && $2>=start && $3<=end {
        print $0
    }' "$in" | \
    # Adjust values > 16569
    awk -v mt_len="$mt_length" '
    BEGIN {FS=OFS="\t"}
    {
        if ($5>mt_len && $6>mt_len) {
            print $1,$2,$3,$4,$5-mt_len,$6-mt_len,$7,$8
        } else {
            print $0
        }
    }' | sort -k8 -n | uniq | \
    # Calculate min/max values
    awk 'BEGIN {FS=OFS="\t"} {
        lines[NR] = $0
        if($2 < chr_min || NR==1) chr_min = $2
        if($3 < chr_max || NR==1) chr_max = $3
        if($5 < mt_min || NR==1) mt_min = $5
        if($6 > mt_max || NR==1) mt_max = $6
        if($7 > query_max || NR==1) query_max = $7
    }
    END {
        if(NR==2) {
            print $1,chr_min,chr_max,$4,mt_min,mt_max,query_max,mt_max-mt_min
        } else if(NR==4) {
            print lines[NR]
        }
    }'
}

process_freq_3() {
    local in="$1"
    local chr="$2"
    local start="$3"
    local end="$4"
    
    awk -v chr="$chr" -v start="$start" -v end="$end" '
    BEGIN {FS=OFS="\t"}
    $1==chr && $2>=start && $3<=end {
        print $0
    }' "$in" | sort -k8 -n -r | head -n 1
}

process_other_freq() {
    local in="$1"
    local chr="$2"
    local start="$3"
    local end="$4"
    local mt_length="$5"

    awk -v chr="$chr" -v start="$start" -v end="$end" '
    BEGIN {FS=OFS="\t"}
    $1==chr && $2>=start && $3<=end {
        print $0
    }' "$in" | \
    awk -v mt_len="$mt_length" '
    BEGIN {FS=OFS="\t"} {
        if($5>mt_len && $6>mt_len) {
            print $1,$2,$3,$4,$5-mt_len,$6-mt_len,$7,$8
        } else {
            print $0
        }
    }' | sort -k8 -n | uniq | \
    awk 'BEGIN{FS=OFS="\t"} {
        lines[NR] = $0
        if($2 < chr_min || NR==1) chr_min = $2
        if($3 < chr_max || NR==1) chr_max = $3
        if($5 < mt_min || NR==1) mt_min = $5
        if($6 > mt_max || NR==1) mt_max = $6
        if($7 > query_max || NR==1) query_max = $7
        if($8 > ref_max || NR==1) ref_max = $8
    } END{
        if(NR==1) {
            print $1,chr_min,chr_max,$4,mt_min,mt_max,query_max,ref_max
        } else if(NR==2) {
            print lines[NR]
        }
    }'
}

concat_numt() {
    while IFS=$'\t' read -r chr start end freq; do
        case "$freq" in
            6)
                process_freq_6 "$1" "$chr" "$start" "$end" "$3"
                ;;
            
            4)
                process_freq_4 "$1" "$chr" "$start" "$end" "$3"
                ;;
            3)
                process_freq_3 "$1" "$chr" "$start" "$end" "$3"
                ;;
            2)
                process_other_freq "$1" "$chr" "$start" "$end" "$3"
                ;;
        esac
    done < "$2"
}

if [ -s "$input" ]; then
    # process bed file
    
    get_freq "$input" > "$tmp_file"
    
    concat_numt "$input" "$tmp_file" "$mt_length" | awk 'BEGIN{FS="\t";OFS="\t"} {
        if(($8*100)/$7 >= 70) print $1,$2,$3,$4,$5,$6
    }' | awk -v mt_len="$mt_length" 'BEGIN{FS=OFS="\t"} {
        if($6>mt_len) {
            print $1, $2, $3, $4, $5, $6-mt_len, "With Control"
        }
            else { print $0, "Without Control"
            }
    }' | awk 'BEGIN{FS=OFS="\t"} {
    	if($5>$6) {
    	   print $1, $2, $3, $4, $6, $5, $7
    	} else { print $1, $2, $3, $4, $5, $6, $7}
    	}' > "$output"

    rm -r "$tmp_file"

else
    echo "Input file is empty."
    touch "$output"
    exit 0
fi
