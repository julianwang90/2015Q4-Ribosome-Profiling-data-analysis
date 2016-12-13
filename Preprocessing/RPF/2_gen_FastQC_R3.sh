#!/bin/bash

# 2016-01-07
# Generate FastQC reports


cd ~/Downloads/FastQC



./fastqc ~/Documents/Combined/RPF/R3_C1_cut_long.fastq.gz
echo "FastQC report for Index 5 complete"
./fastqc ~/Documents/Combined/RPF/R3_UV1_cut_long.fastq.gz 
echo "FastQC report for Index 6 complete"
./fastqc ~/Documents/Combined/RPF/R3_C4_cut_long.fastq.gz 
echo "FastQC report for Index 7 complete"
./fastqc ~/Documents/Combined/RPF/R3_UV4_cut_long.fastq.gz 
echo "FastQC report for Index 8 complete"



echo "FastQC report processing for Rep3 samples completed"


