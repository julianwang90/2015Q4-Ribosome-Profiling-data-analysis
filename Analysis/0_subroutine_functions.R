######################################################################
# Subroutine functions for R based analysis
#
######################################################################

importReads <- function(rpfdirloc, wtdirloc, suffix="_gc.txt"){
        #--------------- Read in RPFs
        # Read in first sample as expression("log"[2]*"(RPF counts)") and add subsequent samples in as additional columns to the data frame
        rpf <- read.table(paste0(rpfdirloc,"R1_C1",suffix), row.names=1, col.names=c("gene","R1_C1"))
        rpf_r1_c4 <- read.table(paste0(rpfdirloc,"R1_C4",suffix))
        rpf_r1_uv1 <- read.table(paste0(rpfdirloc,"R1_UV1",suffix))
        rpf_r1_uv4 <- read.table(paste0(rpfdirloc,"R1_UV4",suffix))
        
        rpf_r3_c1 <- read.table(paste0(rpfdirloc,"R3_C1",suffix))
        rpf_r3_c4 <- read.table(paste0(rpfdirloc,"R3_C4",suffix))
        rpf_r3_uv1 <- read.table(paste0(rpfdirloc,"R3_UV1",suffix))
        rpf_r3_uv4 <- read.table(paste0(rpfdirloc,"R3_UV4",suffix))
        
        rpf_r4_c1 <- read.table(paste0(rpfdirloc,"R4_C1",suffix))
        rpf_r4_c4 <- read.table(paste0(rpfdirloc,"R4_C4",suffix))
        rpf_r4_uv1 <- read.table(paste0(rpfdirloc,"R4_UV1",suffix))
        rpf_r4_uv4 <- read.table(paste0(rpfdirloc,"R4_UV4",suffix))
        
        # Merge all data frames into one
        rpf$R1_C4 <- rpf_r1_c4[,2]
        rpf$R1_UV1 <- rpf_r1_uv1[,2]
        rpf$R1_UV4 <- rpf_r1_uv4[,2]
        
        rpf$R3_C1 <- rpf_r3_c1[,2]
        rpf$R3_C4 <- rpf_r3_c4[,2]
        rpf$R3_UV1 <- rpf_r3_uv1[,2]
        rpf$R3_UV4 <- rpf_r3_uv4[,2]
        
        rpf$R4_C1 <- rpf_r4_c1[,2]
        rpf$R4_C4 <- rpf_r4_c4[,2]
        rpf$R4_UV1 <- rpf_r4_uv1[,2]
        rpf$R4_UV4 <- rpf_r4_uv4[,2]
        
        # Remove last 5 rows (these contain HTSeq-count stats)
        rpf <<- head(rpf, -5)
        
        # Remove intermediate variables
        rm(list=ls(pattern="rpf_"))
        
        #--------------- Read in WTs
        # Read in first sample as expression("log"[2]*"(WT counts)") and add subsequent samples in as additional columns to the data frame
        wt <- read.table(paste0(wtdirloc,"WT_R1_C1",suffix), row.names=1, col.names=c("gene","R1_C1"))
        wt_r1_c4 <- read.table(paste0(wtdirloc,"WT_R1_C4",suffix))
        wt_r1_uv1 <- read.table(paste0(wtdirloc,"WT_R1_UV1",suffix))
        wt_r1_uv4 <- read.table(paste0(wtdirloc,"WT_R1_UV4",suffix))
        
        wt_r3_c1 <- read.table(paste0(wtdirloc,"WT_R3_C1",suffix))
        wt_r3_c4 <- read.table(paste0(wtdirloc,"WT_R3_C4",suffix))
        wt_r3_uv1 <- read.table(paste0(wtdirloc,"WT_R3_UV1",suffix))
        wt_r3_uv4 <- read.table(paste0(wtdirloc,"WT_R3_UV4",suffix))
        
        wt_r4_c1 <- read.table(paste0(wtdirloc,"WT_R4_C1",suffix))
        wt_r4_c4 <- read.table(paste0(wtdirloc,"WT_R4_C4",suffix))
        wt_r4_uv1 <- read.table(paste0(wtdirloc,"WT_R4_UV1",suffix))
        wt_r4_uv4 <- read.table(paste0(wtdirloc,"WT_R4_UV4",suffix))
        
        # Merge all data frames into one
        wt$R1_C4 <- wt_r1_c4[,2]
        wt$R1_UV1 <- wt_r1_uv1[,2]
        wt$R1_UV4 <- wt_r1_uv4[,2]
        
        wt$R3_C1 <- wt_r3_c1[,2]
        wt$R3_C4 <- wt_r3_c4[,2]
        wt$R3_UV1 <- wt_r3_uv1[,2]
        wt$R3_UV4 <- wt_r3_uv4[,2]
        
        wt$R4_C1 <- wt_r4_c1[,2]
        wt$R4_C4 <- wt_r4_c4[,2]
        wt$R4_UV1 <- wt_r4_uv1[,2]
        wt$R4_UV4 <- wt_r4_uv4[,2]
        
        # Remove last 5 rows (these contain HTSeq-count stats)
        wt <<- head(wt, -5)
        
        # Remove intermediate variables
        rm(list=ls(pattern="wt_"))
        rm(rpfdirloc,wtdirloc,suffix)
        
        # Group by condition
        test.group <<- c("C1","C4","UV1","UV4","C1","C4","UV1","UV4","C1","C4","UV1","UV4")
}

outputGeneList <- function(genes, data, output.name="Output", sheet.name="Sheet", col.names=TRUE, row.names=FALSE, append=FALSE){
        library(xlsx)
        library(biomaRt)
        
        if (output.name=="Output" | sheet.name=="Sheet"){
                cat("WARNING: output.name and/or sheet.name not specified.\nUsing default value(s) instead.\n\n")
        }
        
        xlsx_output <- paste0(Sys.Date(), " ",output.name, ".xlsx")
        ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
        
        genelist <- na.omit(data[genes,])
        Gene_description <- getBM(attributes= c('hgnc_symbol','description'), filters = 'hgnc_symbol', values = genelist$Gene, mart = ensembl)
        matched <- match(genelist$Gene, Gene_description$hgnc_symbol)
        genelist <- cbind(genelist, Gene_description[matched,2])
        rm(matched,Gene_description)
        write.xlsx2(genelist, xlsx_output, sheetName=sheet.name, col.names=col.names, row.names=row.names, append=append)
        
        cat("Matched",nrow(genelist),"Genes out of",length(genes),"from input list.\n",
            "Results successfully written to", xlsx_output)
        
}

plotloglog <- function(wt, rpf, pch=16,cex=0.4){
        # Plot log counts of all samples
        plot.new()
        par(mfrow=c(3,4)) 
        plot(log2(wt[,1]), log2(rpf[,1]), xlab=expression("log"[2]*"(WT counts)"), ylab=expression("log"[2]*"(RPF counts)"), main="R1 C1", pch=pch,cex=cex)
        plot(log2(wt[,3]), log2(rpf[,3]), xlab=expression("log"[2]*"(WT counts)"), ylab=expression("log"[2]*"(RPF counts)"), main="R1 UV1", pch=pch,cex=cex)
        plot(log2(wt[,2]), log2(rpf[,2]), xlab=expression("log"[2]*"(WT counts)"), ylab=expression("log"[2]*"(RPF counts)"), main="R1 C4", pch=pch,cex=cex)
        plot(log2(wt[,4]), log2(rpf[,4]), xlab=expression("log"[2]*"(WT counts)"), ylab=expression("log"[2]*"(RPF counts)"), main="R1 UV4", pch=pch,cex=cex)
        
        plot(log2(wt[,5]), log2(rpf[,5]), xlab=expression("log"[2]*"(WT counts)"), ylab=expression("log"[2]*"(RPF counts)"), main="R3 C1", pch=pch,cex=cex)
        plot(log2(wt[,7]), log2(rpf[,7]), xlab=expression("log"[2]*"(WT counts)"), ylab=expression("log"[2]*"(RPF counts)"), main="R3 UV1", pch=pch,cex=cex)
        plot(log2(wt[,6]), log2(rpf[,6]), xlab=expression("log"[2]*"(WT counts)"), ylab=expression("log"[2]*"(RPF counts)"), main="R3 C4", pch=pch,cex=cex)
        plot(log2(wt[,8]), log2(rpf[,8]), xlab=expression("log"[2]*"(WT counts)"), ylab=expression("log"[2]*"(RPF counts)"), main="R3 UV4", pch=pch,cex=cex)
        
        plot(log2(wt[,9]), log2(rpf[,9]), xlab=expression("log"[2]*"(WT counts)"), ylab=expression("log"[2]*"(RPF counts)"), main="R4 C1", pch=pch,cex=cex)
        plot(log2(wt[,11]), log2(rpf[,11]), xlab=expression("log"[2]*"(WT counts)"), ylab=expression("log"[2]*"(RPF counts)"), main="R4 UV1", pch=pch,cex=cex)
        plot(log2(wt[,10]), log2(rpf[,10]), xlab=expression("log"[2]*"(WT counts)"), ylab=expression("log"[2]*"(RPF counts)"), main="R4 C4", pch=pch,cex=cex)
        plot(log2(wt[,12]), log2(rpf[,12]), xlab=expression("log"[2]*"(WT counts)"), ylab=expression("log"[2]*"(RPF counts)"), main="R4 UV4", pch=pch,cex=cex)

}

plot3dvar <- function(wt, rpf){
        #--------- Checking inter-repeat variability
        library(scatterplot3d)
        
        plot.new()
        par(mfrow=c(2,4)) 
        # RPF C1
        scatterplot3d(log2(rpf[,1]),   # x axis
                      log2(rpf[,2]),     # y axis
                      log2(rpf[,3]),    # z axis
                      pch=".",
                      main="RPF C1",
                      xlab=expression("log"[2]*"(RPF C1 Repeat 1)"),
                      ylab=expression("log"[2]*"(RPF C1 Repeat 2)"),
                      zlab=expression("log"[2]*"(RPF C1 Repeat 3)")
        )
        # RPF UV1
        scatterplot3d(log2(rpf[,7]),   # x axis
                      log2(rpf[,8]),     # y axis
                      log2(rpf[,9]),    # z axis
                      pch=".",
                      main="RPF UV1",
                      xlab=expression("log"[2]*"(RPF UV1 Repeat 1)"),
                      ylab=expression("log"[2]*"(RPF UV1 Repeat 2)"),
                      zlab=expression("log"[2]*"(RPF UV1 Repeat 3)")
        )
        
        # RPF C4
        scatterplot3d(log2(rpf[,4]),   # x axis
                      log2(rpf[,5]),     # y axis
                      log2(rpf[,6]),    # z axis
                      pch=".",
                      main="RPF C4",
                      xlab=expression("log"[2]*"(RPF C4 Repeat 1)"),
                      ylab=expression("log"[2]*"(RPF C4 Repeat 2)"),
                      zlab=expression("log"[2]*"(RPF C4 Repeat 3)")
        )
        
        # RPF UV4
        scatterplot3d(log2(rpf[,10]),   # x axis
                      log2(rpf[,11]),     # y axis
                      log2(rpf[,12]),    # z axis
                      pch=".",
                      main="RPF UV4",
                      xlab=expression("log"[2]*"(RPF UV4 Repeat 1)"),
                      ylab=expression("log"[2]*"(RPF UV4 Repeat 2)"),
                      zlab=expression("log"[2]*"(RPF UV4 Repeat 3)")
        )
        
        ##### WT
        # WT C1
        scatterplot3d(log2(wt[,1]),   # x axis
                      log2(wt[,2]),     # y axis
                      log2(wt[,3]),    # z axis
                      pch=".",
                      main="WT C1",
                      xlab=expression("log"[2]*"(WT C1 Repeat 1)"),
                      ylab=expression("log"[2]*"(WT C1 Repeat 2)"),
                      zlab=expression("log"[2]*"(WT C1 Repeat 3)")
        )
        # WT UV1
        scatterplot3d(log2(wt[,7]),   # x axis
                      log2(wt[,8]),     # y axis
                      log2(wt[,9]),    # z axis
                      pch=".",
                      main="WT UV1",
                      xlab=expression("log"[2]*"(WT UV1 Repeat 1)"),
                      ylab=expression("log"[2]*"(WT UV1 Repeat 2)"),
                      zlab=expression("log"[2]*"(WT UV1 Repeat 3)")
        )
        
        # WT C4
        scatterplot3d(log2(wt[,4]),   # x axis
                      log2(wt[,5]),     # y axis
                      log2(wt[,6]),    # z axis
                      pch=".",
                      main="WT C4",
                      xlab=expression("log"[2]*"(WT C4 Repeat 1)"),
                      ylab=expression("log"[2]*"(WT C4 Repeat 2)"),
                      zlab=expression("log"[2]*"(WT C4 Repeat 3)")
        )
        # WT UV4
        scatterplot3d(log2(wt[,10]),   # x axis
                      log2(wt[,11]),     # y axis
                      log2(wt[,12]),    # z axis
                      pch=".",
                      main="WT UV4",
                      xlab=expression("log"[2]*"(WT UV4 Repeat 1)"),
                      ylab=expression("log"[2]*"(WT UV4 Repeat 2)"),
                      zlab=expression("log"[2]*"(WT UV4 Repeat 3)")
        )
        
}


annotateExisting <- function(genelist, output.name="Output", sheet.name="Sheet", namespace="molecular_function", col.names=TRUE, row.names=FALSE, append=FALSE){
        # #for testing
# 
#         data <- transcriptome1hr_up
#         output.name = "GO grouping test DELETE ME"
#         col.names=TRUE
#         row.names=FALSE
#         append=FALSE
#         namespace="molecular_function"

        library(xlsx)
        library(biomaRt)
        library(dplyr)
        
        if (output.name=="Output" | sheet.name=="Sheet"){
                cat("WARNING: output.name and/or sheet.name not specified.\nUsing default value(s) instead.\n\n")
        }
        
        xlsx_output <- paste0(Sys.Date(), " ",output.name, " GO_annotation.xlsx")
        ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
        
        summarised_go_names <- getBM(attributes= c('hgnc_symbol','go_id','name_1006', 'description'), filters = c('hgnc_symbol', 'go_parent_name'), values = list(genelist$Gene, namespace), mart = ensembl)
        
        summarised_go_names <- summarised_go_names %>%
                group_by(hgnc_symbol) %>%
                summarize(go_names = do.call(paste, c(as.list(name_1006), sep=",")))
        options(dplyr.width = Inf)
         
        matched <- match(genelist$Gene, summarised_go_names$hgnc_symbol)
        genelist <- cbind(genelist, summarised_go_names[matched,2])
       
        write.xlsx2(genelist, xlsx_output, sheetName=sheet.name, col.names=col.names, row.names=row.names, append=append)
        
        cat("Annotated",nrow(summarised_go_names), "out of", nrow(genelist),"Genes\n",
            "Results successfully written to", xlsx_output)
        
}

getGeneName <- function(genelist, output){
       
        library(biomaRt)
        library(dplyr)
        ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
        
       
        Gene_description <- getBM(attributes= c('hgnc_symbol','description'), filters = 'hgnc_symbol', values = genelist$Gene, mart = ensembl)
        
        matched <- match(genelist$Gene, Gene_description$hgnc_symbol)
        output <<- cbind(as.data.frame(genelist), Gene_description[matched,2])
        rm(matched,Gene_description)        
        
        
        cat("Obtained",nrow(Gene_description),"gene descriptions out of", nrow(genelist),"from input list.\n",
            "Results outputted to the global environment")
        
}
