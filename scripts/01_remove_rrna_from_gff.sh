#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.gtf output.gtf"
    exit 1
fi

input_file=$1
output_file=$2

# Use awk to filter out lines with 'rRNA' biotype and write the rest to the output file
awk 'BEGIN {FS=OFS="\t"} !($0 ~ /gene_biotype "rRNA"/ || $0 ~ /transcript_biotype "rRNA"/ || $0 ~ /gbkey "rRNA"/)' $input_file > $output_file

echo "Filtered GTF file has been saved to $output_file"