#!/bin/bash

# Script to process and filter BED format genomic data
# Usage: ./numt_concat.sh <input_bed> <output_bed>

set -euo pipefail

# Input validation
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_bed> <output_bed> <coverage_cutoff>"
    exit 1
fi

input="$1"

output="$2"

coverage_cutoff="$3"

tmp_file=$(mktemp)

tmp_file_new=$(mktemp)
    
concat_numt_loc=$(mktemp)

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
    local coverage_cutoff="$6"
    tmp_file_2=$(mktemp)
    passed_calls=$(mktemp)

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
    }' | sort -k8 -n | uniq > "$tmp_file_2"

    awk -F "\t" -v coverage_cutoff="$coverage_cutoff" '{if(($8*100/$7)>=coverage_cutoff) print $0}' "$tmp_file_2" > "$passed_calls"

    if [ -s "$passed_calls" ]; then
       lines=$(cat "$passed_calls" | wc -l)
        if [ "$lines" -gt 1 ]; then
            sort -k8 -n -r "$passed_calls" | head -n 1
            rm -r "$passed_calls"
        else
            cat "$passed_calls"
            rm -r "$passed_calls"
        fi
    else 
        process_concatenated_numts "$tmp_file_2" "$chr" "$start" "$end" "$mt_length"
    fi
}

process_freq_4() {
    local in="$1"
    local chr="$2"
    local start="$3"
    local end="$4"
    local mt_length="$5"
    local coverage_cutoff="$6"
    tmp_file_2=$(mktemp)
    passed_calls=$(mktemp)

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
    }' | sort -k8 -n | uniq > "$tmp_file_2"

    awk -F "\t" -v coverage_cutoff="$coverage_cutoff" '{if(($8*100/$7)>=coverage_cutoff) print $0}' "$tmp_file_2" > "$passed_calls"

    if [ -s "$passed_calls" ]; then
        lines=$(cat "$passed_calls" | wc -l)
        if [ "$lines" -gt 1 ]; then
            sort -k8 -n -r "$passed_calls" | head -n 1
            rm -r "$passed_calls"
        else
            cat "$passed_calls"
            rm -r "$passed_calls"
        fi
    else
        process_concatenated_numts "$tmp_file_2" "$chr" "$start" "$end" "$mt_length"
    fi
}

process_concatenated_numts() {
    local in="$1"
    local chr="$2"
    local start="$3"
    local end="$4"
    local mt_length="$5"
    mt_tmp_file=$(mktemp)

    awk 'BEGIN{FS="\t";OFS="\t"}{print $4, $5, $6}' "$in" | sortBed | bedtools merge -i stdin -d 100 > "$mt_tmp_file"

    line_count=$(cat "$mt_tmp_file" | wc -l)

    if [ "$line_count" -eq 1 ]; then 
        mt_start=$(awk -F "\t" '{print $2}' "$mt_tmp_file")
        mt_end=$(awk -F "\t" '{print $3}' "$mt_tmp_file")
        awk -v mt_start="$mt_start" -v mt_end="$mt_end" '
        BEGIN {FS=OFS="\t"}
        {
            lines[NR] = $0
            if($2 < chr_min || NR==1) chr_min = $2
            if($3 < chr_max || NR==1) chr_max = $3
            if($7 > query_max || NR==1) query_max = $7
            if($8 > ref_max || NR==1) ref_max = $8
        } END {
            print $1, chr_min, chr_max,$4, mt_start, mt_end, query_max, ref_max
        }' "$in"
    else 
        while IFS=$'\t' read -r mt mt_start mt_end; do
            awk -v mt_start="$mt_start" -v mt_end="$mt_end" '
            BEGIN {FS=OFS="\t"}
            {
                lines[NR] = $0
                if($2 < chr_min || NR==1) chr_min = $2
                if($3 < chr_max || NR==1) chr_max = $3
                if($7 > query_max || NR==1) query_max = $7
                ref_sum += $8
            } END {
                print $1, chr_min, chr_max, $4, mt_start, mt_end, query_max, ref_sum
            }' "$in"
        done < "$mt_tmp_file"
    
    fi
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

process_freq_2() {
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

process_other_freq() {
    local in="$1"
    local chr="$2"
    local start="$3"
    local end="$4"
    local mt_length="$5"
    local coverage_cutoff="$6"
    tmp_file_2=$(mktemp)
    passed_calls=$(mktemp)

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
    }' | sort -k8 -n | uniq > "$tmp_file_2"

    awk -F "\t" -v coverage_cutoff="$coverage_cutoff" '{if(($8*100/$7)>=coverage_cutoff) print $0}' "$tmp_file_2" > "$passed_calls"

    if [ -s "$passed_calls" ]; then
        lines=$(cat "$passed_calls" | wc -l)
        if [ "$lines" -gt 1 ]; then
            sort -k8 -n -r "$passed_calls" | head -n 1
            rm -r "$passed_calls"
        else
            cat "$passed_calls"
            rm -r "$passed_calls"
        fi
    else 
        process_concatenated_numts "$tmp_file_2" "$chr" "$start" "$end" "$mt_length"
    fi
}


concat_numt() {
    while IFS=$'\t' read -r chr start end freq; do
        case "$freq" in
            6)
                process_freq_6 "$1" "$chr" "$start" "$end" "$3" "$4"
                ;;
            
            4)
                process_freq_4 "$1" "$chr" "$start" "$end" "$3" "$4"
                ;;
            3)
                process_freq_3 "$1" "$chr" "$start" "$end" "$3"
                ;;
            2)
                process_freq_2 "$1" "$chr" "$start" "$end" "$3"
                ;;
            *)
                if [ "$freq" -gt 3 ] && [ "$freq" -ne 4 ] && [ "$freq" -ne 6 ]; then
                    process_other_freq "$1" "$chr" "$start" "$end" "$3" "$4"
                fi
                ;;
        esac
    done < "$2"
}

if [ -s "$input" ]; then
    # process bed file
    
    get_freq "$input" > "$tmp_file"
    
    concat_numt "$input" "$tmp_file" "$mt_length" "$coverage_cutoff" \
    | awk -v coverage_cutoff="$coverage_cutoff" 'BEGIN{FS="\t";OFS="\t"} {
        if(($8*100)/$7 >= coverage_cutoff) print $1,$2,$3,$4,$5,$6
    }' | awk -v mt_len="$mt_length" 'BEGIN{FS=OFS="\t"} {
        if($6>mt_len) {
            print $1, $2, $3, $4, $5, $6-mt_len, "Control"
        }
        else {
            print $0, "Without Control"
        }
    }' | awk 'BEGIN{FS=OFS="\t"} {
    	if($5>$6) {
    	   print $1, $2, $3, $4, $6, $5, $7
    	} else {
            print $1, $2, $3, $4, $5, $6, $7}
    	}' > "$tmp_file_new"

    cat "$tmp_file_new" | awk -F "\t" '{print $1"\t"$2"\t"$3}' | sort | uniq -cd | awk '{print $2"\t"$3"\t"$4}' > "$concat_numt_loc"

    if [ -s "$concat_numt_loc" ]; then
        while IFS=$'\t' read -r chr start end; do
            awk -v chr="$chr" -v start="$start" -v end="$end" 'BEGIN{FS=OFS="\t"} {
                if ($1==chr && $2==start && $3==end) {
                    print $1, $2, $3, $4, $5, $6, "Concatenated NUMT"
                } else {
                    print $0
                }
            }' "$tmp_file_new"
        done < "$concat_numt_loc" > "$output"
    else
        cat "$tmp_file_new" > "$output"
    fi
    
    rm -r "$tmp_file" "$tmp_file_new" "$concat_numt_loc"

else
    echo "Input file is empty."
    touch "$output"
    exit 0
fi
