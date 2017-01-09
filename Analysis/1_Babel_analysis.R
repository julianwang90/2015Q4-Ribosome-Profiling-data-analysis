######################################################################
# Translatome Analysis using babel
# using gene counts calculated by HTSeq-count
# 
######################################################################

library(babel)

# Set project ID and working directory
proj_id <- "Babel"
setwd("/mnt/databank_doc/MRC_Tox_Unit/2015Q4 Ribosome Profiling Data/Combined/R_Analysis/Babel/")

# Import subroutines
source("/mnt/databank_doc/MRC_Tox_Unit/2015Q4 Ribosome Profiling Data/Combined/R_Analysis/0_subroutine_functions.R")


# Import in data
rpfdirloc <- "/mnt/databank_doc/MRC_Tox_Unit/2015Q4 Ribosome Profiling Data/Combined/RPF/4 Tophat/"
wtdirloc <- "/mnt/databank_doc/MRC_Tox_Unit/2015Q4 Ribosome Profiling Data/Combined/WT/4c Tophat (fr-firststrand)/"
importReads(rpfdirloc, wtdirloc)
rm(rpfdirloc,wtdirloc)

######################################################################
# Run babel
######################################################################
options(mc.cores = 8)
set.seed(12345)
test.babel <- babel(wt, rpf, 
                    group = test.group, 
                    nreps = 1e+07,
                    min.rna = 10)


# use biomaRt to obtain gene name descriptions
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
# configure mRNA_FDR and FDR settings
fdr <- 0.1
mrna_fdr <- 0.1
mrna_logfc <- 1.5

# configure xlsx write parameters
library(xlsx)
xlsx_output <- paste(Sys.Date(), proj_id, "summary.xlsx")

#---------  Counting numbers of unusual counts between conditions
between.babel <- test.babel$between

###--------- 1hr time point
c1vsuv1 <- between.babel[[2]]

        # 1hr Transcriptome - Obtain gene descriptions
        print(paste("1hr Transcriptome:", nrow(c1vsuv1[complete.cases(c1vsuv1) & c1vsuv1$mRNA_FDR<mrna_fdr,]), "DEGs detected in total"))
        transcriptome1hr <- c1vsuv1[complete.cases(c1vsuv1) & c1vsuv1$mRNA_FDR<mrna_fdr,]
        Gene_description <- getBM(attributes= c('hgnc_symbol','description'), filters = 'hgnc_symbol', values = transcriptome1hr$Gene, mart = ensembl)
        matched <- match(transcriptome1hr$Gene, Gene_description$hgnc_symbol)
        transcriptome1hr <- cbind(transcriptome1hr, Gene_description[matched,2])
        rm(matched,Gene_description)
        
                # 1hr Transcriptome - Sort into up/down lists, sort, and write to file
                # Down list
                transcriptome1hr_down <- transcriptome1hr[transcriptome1hr$mRNA_logFC > 0,]
                print(paste("1hr Transcriptome:", nrow(transcriptome1hr_down), "DEGs down-regulated"))
                transcriptome1hr_down <- transcriptome1hr_down[order(transcriptome1hr_down$mRNA_FDR),]
                write.xlsx2(transcriptome1hr_down, xlsx_output, sheetName="transcriptome1hr_down", col.names=TRUE, row.names=FALSE, append=FALSE)
                
                # Up list
                transcriptome1hr_up <- transcriptome1hr[transcriptome1hr$mRNA_logFC < 0,]
                print(paste("1hr Transcriptome:", nrow(transcriptome1hr_up), "DEGs up-regulated"))
                transcriptome1hr_up <- transcriptome1hr_up[order(transcriptome1hr_up$mRNA_FDR),]
                write.xlsx2(transcriptome1hr_up, xlsx_output, sheetName="transcriptome1hr_up", col.names=TRUE, row.names=FALSE, append=TRUE)
           
                     
        # 1hr Translatome - Obtain gene descriptions
        print(paste("1hr Translatome:", nrow(c1vsuv1[complete.cases(c1vsuv1) & c1vsuv1$FDR<fdr,]), "DEGs detected in total"))
        translatome1hr <- c1vsuv1[complete.cases(c1vsuv1) & c1vsuv1$FDR<fdr,]
        Gene_description <- getBM(attributes= c('hgnc_symbol','description'), filters = 'hgnc_symbol', values = translatome1hr$Gene, mart = ensembl)
        matched <- match(translatome1hr$Gene, Gene_description$hgnc_symbol)
        translatome1hr <- cbind(translatome1hr, Gene_description[matched,2])
        rm(matched,Gene_description)
                
                # 1hr Translatome - Sort into up/down lists, sort, and write to file
                # Down list
                translatome1hr_down <- translatome1hr[translatome1hr$Direction == "1",]
                print(paste("1hr Translatome:", nrow(translatome1hr_down), "DEGs down-regulated"))
                translatome1hr_down <- translatome1hr_down[order(translatome1hr_down$FDR),]
                write.xlsx2(translatome1hr_down, xlsx_output, sheetName="translatome1hr_down", col.names=TRUE, row.names=FALSE, append=TRUE)
                
                # Up list
                translatome1hr_up <- translatome1hr[translatome1hr$Direction == "-1",]
                print(paste("1hr Translatome:", nrow(translatome1hr_up), "DEGs up-regulated"))
                translatome1hr_up <- translatome1hr_up[order(translatome1hr_up$FDR),]
                write.xlsx2(translatome1hr_up, xlsx_output, sheetName="translatome1hr_up", col.names=TRUE, row.names=FALSE, append=TRUE)               
                

###--------- 4hr time point
c4vsuv4 <- between.babel[[5]]

        # 4hr Transcriptome - Obtain gene descriptions
        print(paste("4hr Transcriptome:", nrow(c4vsuv4[complete.cases(c4vsuv4) & c4vsuv4$mRNA_FDR<mrna_fdr & abs(c4vsuv4$mRNA_logFC)>mrna_logfc,]), "DEGs detected in total"))
        transcriptome4hr <- c4vsuv4[complete.cases(c4vsuv4) & c4vsuv4$mRNA_FDR<mrna_fdr & abs(c4vsuv4$mRNA_logFC)>mrna_logfc,]
        Gene_description <- getBM(attributes= c('hgnc_symbol','description'), filters = 'hgnc_symbol', values = transcriptome4hr$Gene, mart = ensembl)
        matched <- match(transcriptome4hr$Gene, Gene_description$hgnc_symbol)
        transcriptome4hr <- cbind(transcriptome4hr, Gene_description[matched,2])
        rm(matched,Gene_description)
        
                # 4hr Transcriptome - Sort into up/down lists, sort, and write to file
                # Down list
                transcriptome4hr_down <- transcriptome4hr[transcriptome4hr$mRNA_logFC > 0,]
                print(paste("4hr Transcriptome:", nrow(transcriptome4hr_down), "DEGs down-regulated"))
                transcriptome4hr_down <- transcriptome4hr_down[order(transcriptome4hr_down$mRNA_FDR),]
                write.xlsx2(transcriptome4hr_down, xlsx_output, sheetName="transcriptome4hr_down", col.names=TRUE, row.names=FALSE, append=TRUE) 
                # Up list
                transcriptome4hr_up <- transcriptome4hr[transcriptome4hr$mRNA_logFC < 0,]
                print(paste("4hr Transcriptome:", nrow(transcriptome4hr_up), "DEGs up-regulated"))
                transcriptome4hr_up <- transcriptome4hr_up[order(transcriptome4hr_up$mRNA_FDR),]
                write.xlsx2(transcriptome4hr_up, xlsx_output, sheetName="transcriptome4hr_up", col.names=TRUE, row.names=FALSE, append=TRUE)             
                
        
        # 4hr Translatome - Obtain gene descriptions
        print(paste("4hr Translatome:", nrow(c4vsuv4[complete.cases(c4vsuv4) & c4vsuv4$FDR<fdr,]), "DEGs detected in total"))
        translatome4hr <- c4vsuv4[complete.cases(c4vsuv4) & c4vsuv4$FDR<fdr,]
        Gene_description <- getBM(attributes= c('hgnc_symbol','description'), filters = 'hgnc_symbol', values = translatome4hr$Gene, mart = ensembl)
        matched <- match(translatome4hr$Gene, Gene_description$hgnc_symbol)
        translatome4hr <- cbind(translatome4hr, Gene_description[matched,2])
        rm(matched,Gene_description)
        
                # 4hr Translatome - Sort into up/down lists, sort, and write to file
                # Down list
                translatome4hr_down <- translatome4hr[translatome4hr$Direction == "1",]
                print(paste("4hr Translatome:", nrow(translatome4hr_down), "DEGs down-regulated"))
                translatome4hr_down <- translatome4hr_down[order(translatome4hr_down$FDR),]
                write.xlsx2(translatome4hr_down, xlsx_output, sheetName="translatome4hr_down", col.names=TRUE, row.names=FALSE, append=TRUE) 
                # Up list
                translatome4hr_up <- translatome4hr[translatome4hr$Direction == "-1",]
                print(paste("4hr Translatome:", nrow(translatome4hr_up), "DEGs up-regulated"))
                translatome4hr_up <- translatome4hr_up[order(translatome4hr_up$FDR),]
                write.xlsx2(translatome4hr_up, xlsx_output, sheetName="translatome4hr_up", col.names=TRUE, row.names=FALSE, append=TRUE)               
                rm(xlsx_output)
        
                
# Annotate generated up/down gene lists
        # Transcriptome
        annotateExisting(transcriptome1hr_up, output.name="Annotated DEGs", sheet.name="Transcriptome1U")
        annotateExisting(transcriptome1hr_down, output.name="Annotated DEGs", sheet.name="Transcriptome1D", append=TRUE)
        annotateExisting(transcriptome4hr_up, output.name="Annotated DEGs", sheet.name="Transcriptome4U", append=TRUE)
        annotateExisting(transcriptome4hr_down, output.name="Annotated DEGs", sheet.name="Transcriptome4D", append=TRUE)
        
        # Translatome
        annotateExisting(translatome1hr_up, output.name="Annotated DEGs", sheet.name="Translatome1U", append=TRUE)
        annotateExisting(translatome1hr_down, output.name="Annotated DEGs", sheet.name="Translatome1D", append=TRUE)
        annotateExisting(translatome4hr_up, output.name="Annotated DEGs", sheet.name="Translatome4U", append=TRUE)
        annotateExisting(translatome4hr_down, output.name="Annotated DEGs", sheet.name="Translatome4D", append=TRUE)


runstats <- paste("Run on",Sys.Date(),"using RPF dataset. Ran using babel with 1e+07 permutations
Gene counting done using HTSeq-count
RPF: used --stranded=no
WT: used --stranded=reverse")



######################################################################
# Plotting heatmaps of significant genes in translatome
######################################################################

library(d3heatmap)
library(gplots)
library(htmlwidgets)

###--------- 4hr time point
trans4sorted <- translatome4hr[order(translatome4hr$`P-value`),]
trans4SH <- as.data.frame(trans4sorted[,2])
rownames(trans4SH) <- trans4sorted$Gene
colnames(trans4SH) <- "4hr_log2FC"

# Retrieve matched genes from 1hr time point
trans4sig_in_res1 <- match(row.names(trans4SH),c1vsuv1$Gene)
trans4SH[,2] <- c1vsuv1[trans4sig_in_res1,2]
colnames(trans4SH) <- c("4hr_log2FC","1hr_log2FC")

# Reorder the columns so 1hr time point is first
trans4SH <- trans4SH[c(2,1)]

# Select top 20 genes
trans4SH <- trans4SH[1:20,]

d3heatmap(data.matrix(trans4SH), scale = "none", dendrogram='none',col=redgreen(75))

# Saving to HTML file
hm <- d3heatmap(data.matrix(trans4SH), scale = "none", dendrogram='none',col=redgreen(75))
saveWidget(hm, paste0(Sys.Date()," Heatmap Output (v2) - 4hr", ".html"))


###--------- 1hr time point
trans1sorted <- translatome1hr[order(translatome1hr$`P-value`),]
trans1SH <- as.data.frame(trans1sorted[,2])
rownames(trans1SH) <- trans1sorted$Gene
colnames(trans1SH) <- "1hr_log2FC"

# Retrieve matched genes from 1hr time point
trans1sig_in_res4 <- match(row.names(trans1SH),c4vsuv4$Gene)
trans1SH[,2] <- c4vsuv4[trans1sig_in_res4,2]
colnames(trans1SH) <- c("1hr_log2FC","4hr_log2FC")

# Select top 30
trans1SH <- trans1SH[1:20,]

d3heatmap(data.matrix(trans1SH), scale = "none", dendrogram='none',col=redgreen(75))

# Saving to HTML file
hm <- d3heatmap(data.matrix(trans1SH), scale = "none", dendrogram='none',col=redgreen(75))
saveWidget(hm, paste0(Sys.Date()," Heatmap Output (v2) - 1hr", ".html"))





# Save workspace
save.image(file = paste(Sys.Date(), proj_id, "analysis output (RPF -s no, WT -s reverse).RData"))