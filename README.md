# 2015Q4-Ribosome-Profiling-data-analysis
2016-01-13 Julian Wang 

R/BASH scripts for analysis of 2015Q4 Ribosome Profiling project

Analysis of illumina NextSeq-500 sequencing data, formed of Ribosome Protected Fragments (RPF) samples
and Whole Transcriptome (WT) samples, with a total of 12 samples for each (from 3 biological replicates and 4 conditions).

## Preprocessing
Preprocessing of sequencing reads was performed using custom BASH scripts, and by following the Tuxedo protocol (doi:10.1038/nprot.2012.016).

#### Processing raw sequencing reads
Raw sequencing reads for each samples were concatenated to form a single file per condition.

#### Trim adapters
Adapters were trimmed using the program cutadapt (Martin, 2011), using the appropriate reference sequences from the illumina Adapter Sequences Nov 2015 guide. Following trimming, all samples were run through FastQC (Andrews, 2016) to verify the resulting reads quality.

#### Remove contaminated reads using bowtie
Removal of contaminating rRNA and tRNA sequencing reads was performed using bowtie2 (v2.2.5), with indicies built using the reference data for abundant human ribosomal sequences from the illumina iGenome *homo sapiens* hg38 UCSC dataset, and from the UCSC tRNA hg38 database (http://gtrnadb.ucsc.edu/).

#### Map reads to reference genome using tophat
Sequence alignment to the reference genome was performed using tophat (v2.1.0), with both the reference gene annotation (GTF) file and the bowtie2 indicies for the genome being provided by the illumina iGenome *homo sapiens* hg38 UCSC dataset.

## Analysis
#### Translatome analysis with Babel
Analysis performed using Babel (https://cran.r-project.org/web/packages/babel/index.html)

#### Transcriptome analysis with DESeq2
Analysis performed using DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
