#!/bin/bash

# 2016-01-07
# Generate FastQC reports for WT files


cd ~/Downloads/FastQC



./fastqc ~/Documents/Combined/WT/2\ Cut\ using\ index/WT_R3_C1_cut.fastq.gz
echo "FastQC report for Index 7 complete"
./fastqc ~/Documents/Combined/WT/2\ Cut\ using\ index/WT_R3_UV1_cut.fastq.gz 
echo "FastQC report for Index 12 complete"
./fastqc ~/Documents/Combined/WT/2\ Cut\ using\ index/WT_R3_C4_cut.fastq.gz 
echo "FastQC report for Index 13 complete"
./fastqc ~/Documents/Combined/WT/2\ Cut\ using\ index/WT_R3_UV4_cut.fastq.gz 
echo "FastQC report for Index 14 complete"



echo "FastQC report processing for Rep3 samples completed"


