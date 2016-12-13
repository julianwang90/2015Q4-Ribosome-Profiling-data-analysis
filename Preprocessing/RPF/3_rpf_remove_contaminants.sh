#!/bin/bash

##################################################################################
# Remove Contaminants for RPFs

#Input files from step 2b) Cut using complete index primers
#Filter out reads smaller than 2nt (bowtie limit, any reads smaller than that are ignored)
#Create bowtie indicies using the abundant contaminants fa files from iGenome and the tRNA database (http://gtrnadb.ucsc.edu)
#Align and remove all contaminants using bowtie in 2 steps, first removing tRNA reads then rRNA reads. 

##################################################################################

# Add bowtie2 to PATH
export PATH="~/Downloads/bowtie2-2.2.6:$PATH"
export PATH="~/Downloads/ngsutils-0.5.7/bin:$PATH"

# Set sequencing reads location
cd "/mnt/databank_doc/MRC_Tox_Unit/2014Q4 Ribosome Profiling Data/Combined/RPF/3 Remove Contaminants"

# Name/address of the reference sequences used to build the bowtie index
tRNAs_fa="hg38-tRNAs.fa"
hg38humrib="/mnt/databank_doc/MRC_Tox_Unit/Homo_sapiens/UCSC/hg38/Sequence/AbundantSequences/humRibosomal.fa"
hg385s="/mnt/databank_doc/MRC_Tox_Unit/Homo_sapiens/UCSC/hg38/Sequence/AbundantSequences/hum5SrDNA.fa"

# Filter short reads
mkdir cut_filtered
for d in *.fastq.gz; do
	echo "Processing $d"

	# Regex to shorten sample name
	out=$(echo $d | sed -E "s/(^[A-Z]+.*)_cut_long.fastq.gz/\1/")
	
	# Generate output file names
	out_final="cut_filtered/"$out"_cut_filtered.fastq.gz"
	
	# Filter out reads smaller than 2nt (bowtie limit)
	fastqutils filter -illumina -size 2 $d | gzip -c > $out_final

	echo "Processing of $out complete"

done

cd cut_filtered


# Create bowtie indicies
mkdir tRNAs_bt2
# Build tRNA bowtie index
bowtie2-build $tRNAs_fa tRNAs_bt2/tRNAs
# Build rRNA bowtie index
mkdir rRNAs_bt2
bowtie2-build $hg38humrib,$hg385s rRNAs_bt2/humRib



mkdir tRNA_rRNA_removed
# Align and remove contaminants
for d in *.fastq.gz; do
	echo "Processing $d"

	# Regex to shorten sample name
	out=$(echo $d | sed -E "s/(^[A-Z]+.*)_cut_filtered.fastq.gz/\1/")
	
	# Generate output file names
	no_rRNA="tRNA_rRNA_removed/1_"$out"_no_rRNA.fastq.gz"
	rRNA="tRNA_rRNA_removed/1_"$out"_rRNA.sam.gz"
	no_rRNA_tRNA="tRNA_rRNA_removed/2_"$out"_no_rRNA_tRNA.fastq.gz"
	no_rRNA_yes_tRNA="tRNA_rRNA_removed/2_"$out"_no_rRNA_yes_tRNA.sam.gz"

	# Filter for rRNA
	bowtie2 -p 8 --un-gz $no_rRNA -x rRNAs_bt2/humRib -U <(zcat $d) | gzip -c > $rRNA
	echo "Processing for $out rRNA complete"
	
	# Filter for tRNA
	bowtie2 -p 8 --un-gz $no_rRNA_tRNA -x tRNAs_bt2/tRNAs -U <(zcat $no_rRNA) | gzip -c > $no_rRNA_yes_tRNA
	echo "Processing for $out tRNA complete"

	echo "Processing of $out complete"

done

