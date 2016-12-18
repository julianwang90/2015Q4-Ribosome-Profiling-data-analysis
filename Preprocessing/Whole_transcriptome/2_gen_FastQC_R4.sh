#!/bin/bash

# 2016-01-07
# Generate FastQC reports for WT files


cd ~/Downloads/FastQC


./fastqc ~/Documents/Combined/WT/2\ Cut\ using\ index/WT_R4_C1_cut.fastq.gz 
echo "FastQC report for Index 15 complete"
./fastqc ~/Documents/Combined/WT/2\ Cut\ using\ index/WT_R4_UV1_cut.fastq.gz 
echo "FastQC report for Index 16 complete"
./fastqc ~/Documents/Combined/WT/2\ Cut\ using\ index/WT_R4_C4_cut.fastq.gz 
echo "FastQC report for Index 18 complete"
./fastqc ~/Documents/Combined/WT/2\ Cut\ using\ index/WT_R4_UV4_cut.fastq.gz
echo "FastQC report for Index 19 complete"

echo "FastQC report processing for Rep4 samples completed"


