#!/bin/bash

# 2016-01-07
# Batch process all RPF reads using cutadapt to remove the entire PCR index primers
# Using default error rate of 10%

export PATH="~/.local/bin:$PATH"

cd ~/Documents/Combined/RPF/0\ Combined

#Reverse compliment for Index 1-12

PrimerStart="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC"
PrimerEnd="ATCTCGTATGCCGTCTTCTGCTTG"

Index1="ATCACG"
Index2="CGATGT"
Index3="TTAGGC"
Index4="TGACCA"

Index5="ACAGTG"
Index6="GCCAAT"
Index7="CAGATC"
Index8="ACTTGA"

Index9="GATCAG"
Index10="TAGCTT"
Index11="GGCTAC"
Index12="CTTGTA"


cutadapt -a $PrimerStart$Index1$PrimerEnd -o R1_C1_cut_long.fastq.gz R1_C1.fastq.gz > R1_C1_cut_stats.txt
echo "Rev.Comp. Index 1 complete"
cutadapt -a $PrimerStart$Index2$PrimerEnd -o R1_UV1_cut_long.fastq.gz R1_UV1.fastq.gz > R1_UV1_cut_stats.txt
echo "Rev.Comp. Index 2 complete"
cutadapt -a $PrimerStart$Index3$PrimerEnd -o R1_C4_cut_long.fastq.gz R1_C4.fastq.gz > R1_C4_cut_stats.txt
echo "Rev.Comp. Index 3 complete"
cutadapt -a $PrimerStart$Index4$PrimerEnd -o R1_UV4_cut_long.fastq.gz R1_UV4.fastq.gz > R1_UV4_cut_stats.txt
echo "Rev.Comp. Index 4 complete"


cutadapt -a $PrimerStart$Index5$PrimerEnd -o R3_C1_cut_long.fastq.gz R3_C1.fastq.gz > R3_C1_cut_stats.txt
echo "Rev.Comp. Index 5 complete"
cutadapt -a $PrimerStart$Index6$PrimerEnd -o R3_UV1_cut_long.fastq.gz R3_UV1.fastq.gz > R3_UV1_cut_stats.txt
echo "Rev.Comp. Index 6 complete"
cutadapt -a $PrimerStart$Index7$PrimerEnd -o R3_C4_cut_long.fastq.gz R3_C4.fastq.gz > R3_C4_cut_stats.txt
echo "Rev.Comp. Index 7 complete"
cutadapt -a $PrimerStart$Index8$PrimerEnd -o R3_UV4_cut_long.fastq.gz R3_UV4.fastq.gz > R3_UV4_cut_stats.txt
echo "Rev.Comp. Index 8 complete"


cutadapt -a $PrimerStart$Index9$PrimerEnd -o R4_C1_cut_long.fastq.gz R4_C1.fastq.gz > R4_C1_cut_stats.txt
echo "Rev.Comp. Index 9 complete"
cutadapt -a $PrimerStart$Index10$PrimerEnd -o R4_UV1_cut_long.fastq.gz R4_UV1.fastq.gz > R4_UV1_cut_stats.txt
echo "Rev.Comp. Index 10 complete"
cutadapt -a $PrimerStart$Index11$PrimerEnd -o R4_C4_cut_long.fastq.gz R4_C4.fastq.gz > R4_C4_cut_stats.txt
echo "Rev.Comp. Index 11 complete"
cutadapt -a $PrimerStart$Index12$PrimerEnd -o R4_UV4_cut_long.fastq.gz R4_UV4.fastq.gz > R4_UV4_cut_stats.txt
echo "Rev.Comp. Index 12 complete"

echo "Reverse complement adapter processing completed"


