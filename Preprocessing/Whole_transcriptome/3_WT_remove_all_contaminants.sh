#!/bin/bash

##################################################################################
#3. Remove Contaminants for WTs

#Input files from step 2b) Trimmed using whole index primers and processed to remove reads smaller than 2nt (bowtie limit)
#Create bowtie indicies using the abundant contaminants fa files from iGenome and the tRNA database (http://gtrnadb.ucsc.edu)
#Align and remove all contaminants using bowtie in single step. (Prior testing shows only ~0.1% match to tRNAs) 

#JW 2016-01-14 Files have been moved since this script has been run.
#JW 2016-01-15 contaminant removal % is very high using hg38 tRNAs db. Not sure if it's caused by high rRNA levels or tRNA. 
#	Requires further testing, with removal done in 2 separate steps
##################################################################################

export PATH="~/Downloads/bowtie2-2.2.6:$PATH"
export PATH="~/Downloads/ngsutils-0.5.7/bin:$PATH"


cd "/mnt/databank_doc/MRC_Tox_Unit/2014Q4 Ribosome Profiling Data/Combined/WT/2b Cut using complete index primer"


# Name/address of the tRNA fa file from the tRNA database, used to build the bowtie index
tRNAs_fa="hg38-tRNAs.fa"

cd cut_filtered

# Create bowtie indicies
mkdir t_rRNAs_bt2
# Mount appropriate network drives first
hg38humrib="/mnt/databank_doc/MRC_Tox_Unit/Homo_sapiens/UCSC/hg38/Sequence/AbundantSequences/humRibosomal.fa"
hg385s="/mnt/databank_doc/MRC_Tox_Unit/Homo_sapiens/UCSC/hg38/Sequence/AbundantSequences/hum5SrDNA.fa"

bowtie2-build $hg38humrib,$hg385s,$tRNAs_fa t_rRNAs_bt2/t_rRNAs






# Align and remove contaminants

mkdir tRNA_rRNA_removed
for d in *.fastq.gz; do
	echo "Processing $d"

	# Regex to shorten sample name
	out=$(echo $d | sed -E "s/(^[A-Z]+.*)_cut_long.fastq.gz/\1/")
	
	# Generate output file names
	no_tRNA_rRNA="tRNA_rRNA_removed/"$out"_no_tRNA_rRNA.fastq.gz"
	contaminants="tRNA_rRNA_removed/"$out"_contaminants.sam.gz"

	# Filter for rRNA
	bowtie2 -p 8 --un-gz $no_tRNA_rRNA -x t_rRNAs_bt2/t_rRNAs -U <(zcat $d) | gzip > $contaminants
	echo "Processing for $out complete"


done

