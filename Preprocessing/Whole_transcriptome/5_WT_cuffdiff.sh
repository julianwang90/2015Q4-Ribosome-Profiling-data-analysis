#!/bin/bash

##################################################################################
#5. Cuffmerge and Cuffdiff

# Differential expression analysis for the assembled transcripts

#JW 2016-01-18
##################################################################################


#export PATH="~/Downloads/cufflinks-2.2.1:$PATH"
#export PATH="~/Downloads/samtools-1.3:$PATH"


#Add this to ~/.bashrc and then reload
#PATH=$PATH:/home/julian/Downloads/cufflinks-2.2.1 


genes="/home/julian/Documents/Reference_files/hg38/Genes/genes.gtf"
genome="/home/julian/Documents/Reference_files/hg38/Chromosomes"
	
# Merge assembly files into single gtf

cuffmerge -p 8 -g $genes -s $genome 5_WT_cuffmerge_assembly_list.txt

cuffdiff -p 8 -o cuffdiff_out -b $genome -L C1,UV1 -u merged_asm/merged.gtf \
./WT_R1_C1_cl_out/transcripts.gtf,./WT_R3_C1_cl_out/transcripts.gtf,./WT_R4_C1_cl_out/transcripts.gtf \
./WT_R1_UV1_cl_out/transcripts.gtf,./WT_R3_UV1_cl_out/transcripts.gtf,./WT_R4_UV1_cl_out/transcripts.gtf



