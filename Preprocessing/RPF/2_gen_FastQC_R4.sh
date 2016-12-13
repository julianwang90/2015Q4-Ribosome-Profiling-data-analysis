#!/bin/bash

# 2016-01-07
# Generate FastQC reports


cd ~/Downloads/FastQC


./fastqc ~/Documents/Combined/RPF/R4_C1_cut_long.fastq.gz 
echo "FastQC report for Index 9 complete"
./fastqc ~/Documents/Combined/RPF/R4_UV1_cut_long.fastq.gz 
echo "FastQC report for Index 10 complete"
./fastqc ~/Documents/Combined/RPF/R4_C4_cut_long.fastq.gz 
echo "FastQC report for Index 11 complete"
./fastqc ~/Documents/Combined/RPF/R4_UV4_cut_long.fastq.gz
echo "FastQC report for Index 12 complete"

echo "FastQC report processing for Rep4 samples completed"


