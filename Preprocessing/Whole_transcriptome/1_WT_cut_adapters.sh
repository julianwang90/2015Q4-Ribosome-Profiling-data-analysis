#!/bin/bash

###################################################################################

#2. WT cut adapters using whole primer
#Using cutadapt, filter out reads smaller than 2nt (bowtie limit), and export cut stats to txt file

###################################################################################

export PATH="~/.local/bin:$PATH"


PrimerStart="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
PrimerEnd="ATCTCGTATGCCGTCTTCTGCTTG"

Index2="CGATGT"
Index4="TGACCA"
Index5="ACAGTG"
Index6="GCCAAT"

Index7="CAGATC"
Index12="CTTGTA"
Index13="AGTCAA"
Index14="AGTTCC"

Index15="ATGTCA"
Index16="CCGTCC"
Index18="GTCCGC"
Index19="GTGAAA"


mkdir cut
# TruSeq Stranded mRNA LT SetA adapters
cutadapt -a $PrimerStart$Index2$PrimerEnd -m 2 -o cut/WT_R1_C1_cut_long.fastq.gz WT_R1_C1.fastq.gz > WT_R1_C1_cut_long_stats.txt
echo "WT Adapter trimming & processing for Index 2 complete"
cutadapt -a $PrimerStart$Index4$PrimerEnd -m 2 -o cut/WT_R1_UV1_cut_long.fastq.gz WT_R1_UV1.fastq.gz > WT_R1_UV1_cut_long_stats.txt
echo "WT Adapter trimming & processing for Index 4 complete"
cutadapt -a $PrimerStart$Index5$PrimerEnd -m 2 -o cut/WT_R1_C4_cut_long.fastq.gz WT_R1_C4.fastq.gz > WT_R1_C4_cut_long_stats.txt
echo "WT Adapter trimming & processing for Index 5 complete"
cutadapt -a $PrimerStart$Index6$PrimerEnd -m 2 -o cut/WT_R1_UV4_cut_long.fastq.gz WT_R1_UV4.fastq.gz > WT_R1_UV4_cut_long_stats.txt
echo "WT Adapter trimming & processing for Index 6 complete"


cutadapt -a $PrimerStart$Index7$PrimerEnd -m 2 -o cut/WT_R3_C1_cut_long.fastq.gz WT_R3_C1.fastq.gz > WT_R3_C1_cut_long_stats.txt
echo "WT Adapter trimming & processing for Index 7 complete"
cutadapt -a $PrimerStart$Index12$PrimerEnd -m 2 -o cut/WT_R3_UV1_cut_long.fastq.gz WT_R3_UV1.fastq.gz > WT_R3_UV1_cut_long_stats.txt
echo "WT Adapter trimming & processing for Index 12 complete"
cutadapt -a $PrimerStart$Index13$PrimerEnd -m 2 -o cut/WT_R3_C4_cut_long.fastq.gz WT_R3_C4.fastq.gz > WT_R3_C4_cut_long_stats.txt
echo "WT Adapter trimming & processing for Index 13 complete"
cutadapt -a $PrimerStart$Index14$PrimerEnd -m 2 -o cut/WT_R3_UV4_cut_long.fastq.gz WT_R3_UV4.fastq.gz > WT_R3_UV4_cut_long_stats.txt
echo "WT Adapter trimming & processing for Index 14 complete"


cutadapt -a $PrimerStart$Index15$PrimerEnd -m 2 -o cut/WT_R4_C1_cut_long.fastq.gz WT_R4_C1.fastq.gz > WT_R4_C1_cut_long_stats.txt
echo "WT Adapter trimming & processing for Index 15 complete"
cutadapt -a $PrimerStart$Index16$PrimerEnd -m 2 -o cut/WT_R4_UV1_cut_long.fastq.gz WT_R4_UV1.fastq.gz > WT_R4_UV1_cut_long_stats.txt
echo "WT Adapter trimming & processing for Index 16 complete"
cutadapt -a $PrimerStart$Index18$PrimerEnd -m 2 -o cut/WT_R4_C4_cut_long.fastq.gz WT_R4_C4.fastq.gz > WT_R4_C4_cut_long_stats.txt
echo "WT Adapter trimming & processing for Index 18 complete"
cutadapt -a $PrimerStart$Index19$PrimerEnd -m 2 -o cut/WT_R4_UV4_cut_long.fastq.gz WT_R4_UV4.fastq.gz > WT_R4_UV4_cut_long_stats.txt
echo "WT Adapter trimming & processing for Index 19 complete"

echo "WT Adapter trimming & processing complete"






