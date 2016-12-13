#!/bin/bash

##################################################################################
# Tophat Alignment
#
# Alignment and with tophat, using input from tRNA and rRNA removed fastq files.
#
##################################################################################

# Add required programs to PATH
export PATH="~/Downloads/bowtie2-2.2.6:$PATH"
export PATH="~/Downloads/tophat:$PATH"
export PATH="~/Downloads/samtools-1.3:$PATH"

# Set sequencing reads location
cd "/mnt/databank_doc/MRC_Tox_Unit/2014Q4 Ribosome Profiling Data/Combined/RPF/4 TopHat"

# Name/address of reference files
genes="/home/julian/Documents/Reference_files/hg38/Genes/genes.gtf"
bt2idx="/home/julian/Documents/Reference_files/hg38/Bowtie2Index/genome"

for d in *.fastq.gz; do
	echo "Processing $d"
		
	# Regex to shorten sample name
	ident=$(echo $d | sed -E "s/^2_([A-Z]+.*)_no_tRNA_rRNA.fastq.gz/\1/")
	
	# Generate output file names 
	th_out=$ident"_tophat_out"
	mkdir $th_out
	
	#-------- Sequence alignment with TopHat
	tophat --GTF $genes -p 8 -g 1 -o $th_out $bt2idx $d

	echo "Processing of $ident complete"

done

