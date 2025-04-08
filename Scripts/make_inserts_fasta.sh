#!/bin/bash

insertion_calls=$1

inserts_fasta=$2

while read p q r s
do

	echo -e ">${p}_${q}_len_${s}\n$r"

done < "$insertion_calls" > $inserts_fasta
