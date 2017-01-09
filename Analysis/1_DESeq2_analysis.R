######################################################################
# Transcriptome analysis using DESeq2
#
######################################################################

# Set project ID and working directory
proj_id <- "DESeq2 Transcriptome Analysis"
setwd("/mnt/databank_doc/MRC_Tox_Unit/2015Q4 Ribosome Profiling Data/Combined/R_Analysis/DESeq2/")

# Import subroutines
source("/mnt/databank_doc/MRC_Tox_Unit/2015Q4 Ribosome Profiling Data/Combined/R_Analysis/0_subroutine_functions.R")

library("DESeq2")
library("BiocParallel")
register(MulticoreParam(8))

# Import in data & split into 1 & 4hr data frames
rpfdirloc <- "/mnt/databank_doc/MRC_Tox_Unit/2015Q4 Ribosome Profiling Data/Combined/RPF/4b Tophat (tRNA not removed)/"
wtdirloc <- "/mnt/databank_doc/MRC_Tox_Unit/2015Q4 Ribosome Profiling Data/Combined/WT/4c Tophat (fr-firststrand)/"
importReads(rpfdirloc, wtdirloc)
rm(rpfdirloc,wtdirloc)
        rpf1 <- rpf[,c(1,3,5,7,9,11)]
        rpf4 <- rpf[,c(2,4,6,8,10,12)]
        wt1 <- wt[,c(1,3,5,7,9,11)]
        wt4 <- wt[,c(2,4,6,8,10,12)]


# Set colData for wt's
colData1 <- matrix(c(rep(c("control","uv"),3),c(rep("R1",2),rep("R3",2),rep("R4",2))), nrow=6, ncol=2)
row.names(colData1) <- colnames(wt1)
colnames(colData1) <- c("condition","repeat")

colData4 <- colData1
row.names(colData4) <- colnames(wt4)

dds1 <- DESeqDataSetFromMatrix(countData = wt1,
                               colData = as.data.frame(colData1),
                               design = ~ condition)

dds4 <- DESeqDataSetFromMatrix(countData = wt4,
                              colData = as.data.frame(colData4),
                              design = ~ condition)

# Pre-filtering reads
dds1 <- dds1[ rowSums(counts(dds1)) > 1, ]
dds4 <- dds4[ rowSums(counts(dds4)) > 1, ]

# Run & order by p-adj
dds1 <- DESeq(dds1)
res1 <- results(dds1, alpha=0.05)
res1Ordered <- res1[order(res1$padj),]
summary(res1)

dds4 <- DESeq(dds4)
res4 <- results(dds4, alpha=0.05)
res4Ordered <- res4[order(res4$padj),]
summary(res4)
    
# Check numbers of DEGs given p-value
padj <- 0.05
sum(res1$padj < padj, na.rm=TRUE)
sum(res4$padj < padj, na.rm=TRUE)

# MA Plot
plotMA(res1, alpha=0.05, main="DESeq2", ylim=c(-2,2))
plotMA(res4, alpha=0.05, main="DESeq2", ylim=c(-2,2))

# PCA plot 
rld1 <- rlog(dds1)
rld4 <- rlog(dds4)
plotPCA(rld1, intgroup=c("condition","repeat."))
plotPCA(rld4, intgroup=c("condition","repeat."))

# Dispersion plots
plotDispEsts(dds1, main="Dispersion of 1hr counts")
plotDispEsts(dds4, main="Dispersion of 4hr counts")

######################################################################
# Write DESeq2 output to Excel
######################################################################

library(xlsx)

res1sig <- na.omit(as.data.frame(res1)[res1$padj < 0.05,])
res1sig$gene <- row.names(res1sig)
res1up <- res1sig[res1sig$log2FoldChange > 0,]
res1down <- res1sig[res1sig$log2FoldChange < 0,]

res4sig <- na.omit(as.data.frame(res4)[res4$padj < 0.05,])
res4sig$gene <- row.names(res4sig)
res4up <- res4sig[res4sig$log2FoldChange > 0,]
res4down <- res4sig[res4sig$log2FoldChange < 0,]

# Write lists to file
write.xlsx(res1up[order(res1up$log2FoldChange, decreasing = TRUE),], file=paste(Sys.Date(), "DESeq2 Sig Output List.xlsx"), sheetName = paste0("1U (", nrow(res1up)," total)"), append=FALSE)
write.xlsx(res1down[order(res1down$log2FoldChange, decreasing = TRUE),], file=paste(Sys.Date(), "DESeq2 Sig Output List.xlsx"), sheetName = paste0("1D (", nrow(res1down)," total)"), append=TRUE)
write.xlsx(res4up[order(res4up$log2FoldChange, decreasing = TRUE),], file=paste(Sys.Date(), "DESeq2 Sig Output List.xlsx"), sheetName = paste0("4U (", nrow(res4up)," total)"), append=TRUE)
write.xlsx(res4down[order(res4down$log2FoldChange, decreasing = TRUE),], file=paste(Sys.Date(), "DESeq2 Sig Output List.xlsx"), sheetName = paste0("4D (", nrow(res4down)," total)"), append=TRUE)

######################################################################
# Plot heatmaps of the most significant DEGs
######################################################################

library(d3heatmap)
library(gplots)
library(htmlwidgets)

#---------- 1hr time point 
# Sort by significance (padj)
res1sigsorted <- res1sig[order(res1sig[,6]),]
res1sigSH <- as.data.frame(res1sigsorted[,2])
rownames(res1sigSH) <- rownames(res1sigsorted)
colnames(res1sigSH) <- "1hr_log2FC"

# Retrieve matched genes from 4hr time point
res1sig_in_res4 <- match(row.names(res1sigSH),row.names(res4))
res1sigSH[,2] <- res4[res1sig_in_res4,2]
colnames(res1sigSH) <- c("1hr_log2FC","4hr_log2FC")

# Select top 30
res1sigSH <- res1sigSH[1:30,]

d3heatmap(data.matrix(res1sigSH), scale = "none", dendrogram='none',col=redgreen(75))

# Saving to HTML file
hm <- d3heatmap(data.matrix(res1sigSH), scale = "none", dendrogram='none',col=redgreen(75))
saveWidget(hm, paste0(Sys.Date()," Heatmap Output (v2) - 1hr", ".html"))


#---------- 4hr time point 
# Sort by significance (padj)
res4sigsorted <- res4sig[order(res4sig[,6]),]
res4sigSH <- as.data.frame(res4sigsorted[,2])
rownames(res4sigSH) <- rownames(res4sigsorted)
colnames(res4sigSH) <- "4hr_log2FC"

# Retrieve matched genes from 1hr time point
res4sig_in_res1 <- match(row.names(res4sigSH),row.names(res1))
res4sigSH[,2] <- res1[res4sig_in_res1,2]
colnames(res4sigSH) <- c("4hr_log2FC","1hr_log2FC")

# Select top 30
res4sigSH <- res4sigSH[1:30,]

# Re-order columns to place 1hr TP first
res4sigSH <- res4sigSH[c(2,1)]

d3heatmap(data.matrix(res4sigSH), scale = "none", dendrogram='none',col=redgreen(75))

# Saving to HTML file
hm <- d3heatmap(data.matrix(res4sigSH), scale = "none", dendrogram='none',col=redgreen(75))
saveWidget(hm, paste0(Sys.Date()," Heatmap Output (v2) - 4hr", ".html"))




# Save workspace
save.image(file = paste(Sys.Date(), proj_id, "output.RData"))
