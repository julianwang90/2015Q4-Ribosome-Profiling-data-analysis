#!/bin/bash

##################################################################################
#4. Tophat & Cufflinks
#	Running tophat with fr-firststrand library type (Illumina stranded mRNA library)

##################################################################################

export PATH="~/Downloads/bowtie2-2.2.6:$PATH"
export PATH="~/Downloads/ngsutils-0.5.7/bin:$PATH"
export PATH="~/Downloads/tophat:$PATH"
export PATH="~/Downloads/cufflinks-2.2.1:$PATH"
export PATH="~/Downloads/samtools-1.3:$PATH"


cd rRNA_tRNA_removed

genes="/home/julian/Documents/Reference_files/hg38/Genes/genes.gtf"
bt2idx="/home/julian/Documents/Reference_files/hg38/Bowtie2Index/genome"
transcriptome_idx="/home/julian/Documents/Reference_files/hg38/Transcriptome_idx/known"


# Build transcriptome index
#tophat -p 8 -G $genes --library-type fr-firststrand --transcriptome-index=$transcriptome_idx $bt2idx


for d in *.fastq.gz; do
	echo "Processing $d"
	
	# Regex to shorten sample name
	ident=$(echo $d | sed -E "s/^([A-Z]+.*)_no_rRNA_tRNA.fastq.gz/\1/")
	
	# Generate output file names 
	th_out=$ident"_th-frfs_out"
	mkdir $th_out
	
	#-------- Sequence alignment with TopHat
	tophat -p 8 -g 1 -o $th_out --library-type fr-firststrand --transcriptome-index=$transcriptome_idx $bt2idx $d
	
	echo "Tophat processing of $ident complete"

done

for d in *.fastq.gz; do
	echo "Processing $d"
	
	# Regex to shorten sample name
	ident=$(echo $d | sed -E "s/^([A-Z]+.*)_no_rRNA_tRNA.fastq.gz/\1/")
	
	# Generate output file names 
	th_out=$ident"_th-frfs_out"
	
	cl_out=$ident"_cl-frfs_out"
	cufflinks -p 8 --library-type fr-firststrand -G $genes -o $cl_out $th_out/accepted_hits.bam

	#-------- Create index files for the BAM files
#	samtools index $th_out/accepted_hits.bam
	
	
	echo "Cufflinks processing of $ident complete"

done


