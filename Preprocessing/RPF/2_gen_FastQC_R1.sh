#!/bin/bash

# 2016-01-07
# Generate FastQC reports


cd ~/Downloads/FastQC


./fastqc ~/Documents/Combined/RPF/R1_C1_cut_long.fastq.gz
echo "FastQC report for Index 1 complete"
./fastqc ~/Documents/Combined/RPF/R1_UV1_cut_long.fastq.gz 
echo "FastQC report for Index 2 complete"
./fastqc ~/Documents/Combined/RPF/R1_C4_cut_long.fastq.gz 
echo "FastQC report for Index 3 complete"
./fastqc ~/Documents/Combined/RPF/R1_UV4_cut_long.fastq.gz 
echo "FastQC report for Index 4 complete"


echo "FastQC report processing for Rep1 samples completed"


