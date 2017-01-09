#!/bin/bash

######################################################################
# HTSeq_count bash script
#
# Counting reads using HTSeq_count
#
######################################################################

# Set sequencing reads location
cd "/mnt/databank_doc/MRC_Tox_Unit/2015Q4 Ribosome Profiling Data/Combined/RPF/4 Tophat"

# Name/address of reference files
genes="/mnt/databank_doc/MRC_Tox_Unit/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"

# Rep 1
python -m HTSeq.scripts.count -s no -f bam R1_C1_tophat_out/accepted_hits.bam $genes > R1_C1_gc.txt
python -m HTSeq.scripts.count -s no -f bam R1_UV1_tophat_out/accepted_hits.bam $genes > R1_UV1_gc.txt
python -m HTSeq.scripts.count -s no -f bam R1_C4_tophat_out/accepted_hits.bam $genes > R1_C4_gc.txt
python -m HTSeq.scripts.count -s no -f bam R1_UV4_tophat_out/accepted_hits.bam $genes > R1_UV4_gc.txt
echo "HTSeq-count for Rep 1 complete"

# Rep 3
python -m HTSeq.scripts.count -s no -f bam R3_C1_tophat_out/accepted_hits.bam $genes > R3_C1_gc.txt
python -m HTSeq.scripts.count -s no -f bam R3_UV1_tophat_out/accepted_hits.bam $genes > R3_UV1_gc.txt
python -m HTSeq.scripts.count -s no -f bam R3_C4_tophat_out/accepted_hits.bam $genes > R3_C4_gc.txt
python -m HTSeq.scripts.count -s no -f bam R3_UV4_tophat_out/accepted_hits.bam $genes > R3_UV4_gc.txt
echo "HTSeq-count for Rep 3 complete"

# Rep 4
python -m HTSeq.scripts.count -s no -f bam R4_C1_tophat_out/accepted_hits.bam $genes > R4_C1_gc.txt
python -m HTSeq.scripts.count -s no -f bam R4_UV1_tophat_out/accepted_hits.bam $genes > R4_UV1_gc.txt
python -m HTSeq.scripts.count -s no -f bam R4_C4_tophat_out/accepted_hits.bam $genes > R4_C4_gc.txt
python -m HTSeq.scripts.count -s no -f bam R4_UV4_tophat_out/accepted_hits.bam $genes > R4_UV4_gc.txt
echo "HTSeq-count for Rep 4 complete"

