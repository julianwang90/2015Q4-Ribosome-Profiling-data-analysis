#!/bin/bash

# HTSeq_count bash script

genes="/home/julian/Documents/Reference_files/hg38/Genes/genes.gtf"

cd "/mnt/databank_doc/MRC_Tox_Unit/2015Q4 Ribosome Profiling Data/Combined/WT/4c Tophat (fr-firststrand)"

# Rep 1
python -m HTSeq.scripts.count -s reverse -f bam WT_R1_C1_th-frfs_out/accepted_hits.bam $genes > WT_R1_C1_gc.txt
python -m HTSeq.scripts.count -s reverse -f bam WT_R1_UV1_th-frfs_out/accepted_hits.bam $genes > WT_R1_UV1_gc.txt
python -m HTSeq.scripts.count -s reverse -f bam WT_R1_C4_th-frfs_out/accepted_hits.bam $genes > WT_R1_C4_gc.txt
python -m HTSeq.scripts.count -s reverse -f bam WT_R1_UV4_th-frfs_out/accepted_hits.bam $genes > WT_R1_UV4_gc.txt

echo "HTSeq-count for Rep 1 complete"

# Rep 3
python -m HTSeq.scripts.count -s reverse -f bam WT_R3_C1_th-frfs_out/accepted_hits.bam $genes > WT_R3_C1_gc.txt
python -m HTSeq.scripts.count -s reverse -f bam WT_R3_UV1_th-frfs_out/accepted_hits.bam $genes > WT_R3_UV1_gc.txt
python -m HTSeq.scripts.count -s reverse -f bam WT_R3_C4_th-frfs_out/accepted_hits.bam $genes > WT_R3_C4_gc.txt
python -m HTSeq.scripts.count -s reverse -f bam WT_R3_UV4_th-frfs_out/accepted_hits.bam $genes > WT_R3_UV4_gc.txt

echo "HTSeq-count for Rep 3 complete"

# Rep 4
python -m HTSeq.scripts.count -s reverse -f bam WT_R4_C1_th-frfs_out/accepted_hits.bam $genes > WT_R4_C1_gc.txt
python -m HTSeq.scripts.count -s reverse -f bam WT_R4_UV1_th-frfs_out/accepted_hits.bam $genes > WT_R4_UV1_gc.txt
python -m HTSeq.scripts.count -s reverse -f bam WT_R4_C4_th-frfs_out/accepted_hits.bam $genes > WT_R4_C4_gc.txt
python -m HTSeq.scripts.count -s reverse -f bam WT_R4_UV4_th-frfs_out/accepted_hits.bam $genes > WT_R4_UV4_gc.txt

echo "HTSeq-count for Rep 4 complete"

