##########################################################################################
# This script performs differential gene expression analysis with DESeq2 #
##########################################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")

BiocManager::install("DESeq2")
install.packages("pheatmap")
install.packages("RColorBrewer")
BiocManager::install("vsn")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("IHW")
install.packages("tidyverse")
install.packages("ashr")
BiocManager::install("EnhancedVolcano")
install.packages('digest', repos='http://cran.us.r-project.org')
install.packages('rlang', repos='http://cran.us.r-project.org')

### LOAD REQUIRED LIBRARIES
library(digest)
library(rlang)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(vsn)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(genefilter)
library(biomaRt)
library(IHW)
library(tidyverse)
library(ashr)
library(EnhancedVolcano)

### SET WORKING DIRECTORY
# note: this directory should be populated with the raw counts file
setwd("/Users/brianamullins/Desktop/fastqs")

### Import count table and details on experimental design
# NB: Make sure column names in the sample(table) file and counts file are exactly the same and in the same order
CountTable <- read.csv("ALL_features.csv", header=TRUE, row.names=1)
samples <- read.table("samples.txt", header=TRUE)
rowsCountTable = as.data.frame(colnames(CountTable))
input <- CountTable[6:11]
#condition should always be last
feat_data <- DESeqDataSetFromMatrix(countData = input, colData=samples, design=~condition)

### PRELIMINARY ANALYSES ###
# The first steps in your analysis should focus on better understanding the relationship of the datasets being studied. 
# This can be simply achieved by generating a PCA plot showing the relationship of your samples.
# First we transform our raw count data using a variance stabilizing transformation (VST) that roughly mirrors how DeSeq2 models the data.
vsd1 <- varianceStabilizingTransformation(feat_data, blind=FALSE)

# Then we plot a PCA, grouping and coloring our datasets according to batch
plotPCA(vsd1, "condition")
### note that you can attach additional information based on the column headers in your sample table
plotPCA(vsd1, intgroup=c("condition")) + ggtitle("Features PCA")

# we can also attempt to replicate the batch effect correction performed by DeSeq2 using the limma::removeBatchEffect function
vsd2 <- varianceStabilizingTransformation(feat_data, blind=FALSE)
assay(vsd2) <- limma::removeBatchEffect(assay(vsd2), vsd2$batch)
plotPCA(vsd2, "condition") + ggtitle("Features PCA")

### We can also calculate and plot sample distances using either the batch corrected (vsd2) or uncorrected (vsd1) data. 

# uncorrected
sampleDists <- dist( t( assay(vsd1) ) )
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)

# corrected
sampleDistsCorr <- dist( t( assay(vsd2) ) )
sampleDistsCorr
sampleDistCorrMatrix <- as.matrix( sampleDistsCorr )
colnames(sampleDistCorrMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistCorrMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)

# At this stage, you should have a good sense of how your samples cluster and the effect of batch correction (if used)
# In simple RNA-Seq situations (control vs treatment, 3-5 bioreps each), this is all that should be required. 
# For more complex situations, you will need to dive deep into working of DeSeq2
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

### BASIC DGE ANALYSIS USING DESEQ2 ###
#feat_data <- DESeqDataSetFromMatrix(countData = CountTable, colData=samples, design=~condition) # incorporates batch information (if available)

# Run DESEQ and generate a simple plot showing the distribution of regulated and unregulated genes
DatasetProcessed <- DESeq(feat_data) # runs DESEQ
par(mfrow=c(1,1))
plotMA(DatasetProcessed, main="Features DESeq2 MLE", ylim=c(-5,5))

# Next we perform a contrast analysis to produce a list of differentially regulated genes between our two conditions
# First we set CTRL dataset as baseline
feat_data$condition <- relevel(feat_data$condition, "Ctrl")

# Next we create our results object while performing shrinkage of effect size 
# (this reduces the impact of apparent gross changes in low expressed genes)
install.packages("ashr")
library("ashr")
res1 <- lfcShrink(DatasetProcessed, contrast=c("condition","dsDNA","Ctrl"), type = "ashr")
DESeq2::plotMA(res1, main="Features DESeq2 MAP", ylim = c(-6,6))

# Here we modify our output data to include two additional columns that contain the baseMeans (a proxy for counts)
# This is useful for downstream filtering of lowly expressed genes
baseMeanCtrl = rowMeans(counts(DatasetProcessed,normalized=TRUE)[,DatasetProcessed$condition == "Ctrl"])
baseMeanCnot = rowMeans(counts(DatasetProcessed,normalized=TRUE)[,DatasetProcessed$condition == "Cnot"])
res1 = cbind(as.data.frame(res1), baseMeanCtrl, baseMeanCnot)

# Here we add two further columns, the gene symbol (common name) and entrez ID - both of which may be useful downstream
res1$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res1), column="SYMBOL", keytype="ENSEMBL", multiVals="first") # MAPS GENE IDs
res1$entrez <- mapIds(org.Hs.eg.db, keys=row.names(res1), column="ENTREZID", keytype="ENSEMBL", multiVals="first")

# Finally we write the complete results object to an outfile
write.csv(res1, "Features_DEseq_analysis.csv", row.names=TRUE)

############## PLOTTING #######################
#filter DESeq2 results to find which and how many upregulated and downregulated genes
res1_new <- read.csv(file = "Features_DEseq_analysis.csv")
###############################################
nrow(res1)
res1_stat_sig <- res1 %>% filter(pvalue < .00001)
nrow(res1_stat_sig)
res1_upreg <- res1_stat_sig %>% filter(log2FoldChange > 1)
nrow(res1_upreg)
res1_downreg <- res1_stat_sig %>% filter(log2FoldChange < -1)
nrow(res1_downreg)
res1_upreg_downreg <- rbind(res1_upreg,res1_downreg)
###############################################

### Plot heatmap - sorting by padj and log2fc
# PICK ALL GENES WITH pADJ < 0.00001 AND THEN SUBSET FOR THOSE WITH Log2FC > 1 THEN PICK TOP 25 HITS
res1_upreg <- res1_upreg[order(-res1_upreg$log2FoldChange),]
sigGenes <- res1_upreg$...1[1:25]

rows <- match(sigGenes, row.names(vsd1))
mat <- assay(vsd1)[rows,]
mat <- mat-rowMeans(mat) ##plot distances from the mean, makes heatmap clearer
mart <- useMart("ensembl","hsapiens_gene_ensembl") ## assuming human
gns <- getBM(c("hgnc_symbol","ensembl_gene_id"), "ensembl_gene_id", row.names(mat), mart)
row.names(mat)[match(gns[,2], row.names(mat))] <- gns[,1]
df <- as.data.frame(colData(DatasetProcessed))
genes_of_interest <- rownames(mat)
pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df,) 

### order by significance
res1_upreg <- res1_upreg[order(res1_upreg$padj),]
sigGenes <- res1_upreg$...1[1:25]

rows <- match(sigGenes, row.names(vsd1))
mat <- assay(vsd1)[rows,]
mat <- mat-rowMeans(mat) ##plot distances from the mean, makes heatmap clearer
mart <- useMart("ensembl","hsapiens_gene_ensembl") ## assuming human
gns <- getBM(c("hgnc_symbol","ensembl_gene_id"), "ensembl_gene_id", row.names(mat), mart)
row.names(mat)[match(gns[,2], row.names(mat))] <- gns[,1]
df <- as.data.frame(colData(DatasetProcessed))
genes_of_interest <- rownames(mat)
pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df,) 


############## PLOTTING ############## 
#Volcano plot with default parameters (p<.00001, FC > 1, FC < -1)

EnhancedVolcano(res1,
                lab = res1$symbol,
                x = 'log2FoldChange',
                y= 'padj',
                pCutoff = 0.00001,
                FCcutoff = 1,
                xlim = c(-4,4),
                ylim = c(0,300),
                title = "Features")

### VOLCANO PLOTS
feat_data$condition <- relevel(feat_data$condition, "Ctrl")
res1 <- lfcShrink(DatasetProcessed, contrast=c("condition","dsDNA","Ctrl"), type = "ashr")
with(res1, plot(log2FoldChange, -log10(pvalue), pch=16, cex=1.5, main="Volcano plot", xlim=c(-5,5)))
with(subset(res1, padj>0.01), points(log2FoldChange, -log10(pvalue), pch=16, cex=1.5, col="gray"))
with(subset(res1, padj<0.01 & log2FoldChange>0), points(log2FoldChange, -log10(pvalue), pch=16, cex=1.5, col="red"))
with(subset(res1, padj<0.01 & log2FoldChange<0), points(log2FoldChange, -log10(pvalue), pch=16, cex=1.5, col="blue"))
### ADD BELLS AND WHISTLES
pval = 0.01
abline(h = -log10(pval), col = "black", lty = 2, lwd=4)
mtext(paste("pval = 0.01", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)












