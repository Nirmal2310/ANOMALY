#!/bin/bash

helpFunction()
{
   echo "Usage: remove_overlap.sh -i sample_concatenated_NUMTs.txt -s sample_chrM_SA_final.txt"
   echo -e "\t-i <str> Absolute Path of Final NUMTs Calls from Sniffles Insertion Calls."
   echo -e "\t-s <str> Absolute Path of Final NUMTs Calls from Supplementary Alignments."
   echo -e "\t-o <str> Absolute Path of Final NUMT Output."
   exit 1 # Exit script after printing help
}

while getopts "i:s:o:" opt
do
        case "$opt" in
        i )
                ins_numt="$OPTARG"
                ;;
        s )
                sa_numt="$OPTARG"
                ;;
        o )
                final_out="$OPTARG"
                ;;
        ? ) helpFunction ;;
        esac
done

path=$(conda info --base)

source "$path"/bin/activate anomaly

if [ -s "$ins_numt" ] || [ -s "$sa_numt" ]; then
    cat "$ins_numt" "$sa_numt" | sortBed \
    | awk 'BEGIN{FS="\t";OFS="\t"} {
       if($7=="With Control") {
          print $0 | "bedtools merge -i stdin -d 100 -c 2,3,4,5,6,7 -o max,min,distinct,max,min,distinct"
       } else if($7=="Without Control") {
          print $0 | "bedtools merge -i stdin -d 100 -c 2,3,4,5,6,7 -o max,min,distinct,min,max,distinct"
       } else if($7=="Concatenated NUMT") {
          print $1,$2,$3,$2,$3,$4,$5,$6,$7
       }}' \
    | awk 'BEGIN{FS=OFS="\t"} {
        if($9=="Without Control") {
            print $1,$4+1,$6,$7,$8,$8-$7
        } else if ($9=="Control") {
            print $1,$4+1,$6,$7,$8,(16569+$7)-$8
        } else {
            print $1,$4+1,$6,$7,$8,$8-$7
        }
    }' | sort -k1 -n > "$final_out"
elif [ -s "$ins_numt" ]; then
    
    sort -k1 -n "$ins_numt" > "$final_out"

elif [ -s "$sa_numt" ]; then
    
    sort -k1 -n "$sa_numt" > "$final_out"

else
    
    echo "No NUMTs Found."
    
    touch "$final_out"

fi