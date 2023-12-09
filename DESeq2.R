library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(DataEditR)
library(pheatmap)
library(ggrepel)
library(GGally)

setwd('/Users/juanjovel/jj/data_analysis/walterLopez/Carolina_Tomato_Nematode/new_results_paired-end')

data_file <- 'all_samples_counts_b.tsv'
metadata_file <- 'metadata_b.tsv'

data_df <- read.table(data_file, header = TRUE, sep = '\t', row.names = 1)
metadata <- read.table(metadata_file, header = TRUE, sep = '\t', row.names = 1)

# Filter only genes with 10 or more reads on average
data_df_10up <- subset(data_df, rowMeans(data_df) >= 10)
matrix <- as.matrix(round(data_df_10up, digits = 0))

##### Generate dot plots of pair samples counts
# Import data matrices into a DESeq2 objects
dds <- DESeqDataSetFromMatrix(countData = matrix, 
                                 colData = metadata,
                                 design =~ group)

dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

# Convert to data frame
norm_counts_df <- as.data.frame(norm_counts)

# Adding gene names as row names (if needed)
rownames(norm_counts_df) <- rownames(norm_counts)

# Customize ggpairs plot
ggpairs(norm_counts_df,
        lower = list(continuous = wrap("points", size = 0.1, colour = "blue")),  # Customizing lower triangle
        upper = list(continuous = wrap("cor", size = 3)),  # Customizing upper triangle
        diag = list(continuous = wrap("blankDiag", size = 1)),  # Customizing diagonal plots
        title = "Pairwise Plots of Normalized Counts") +
  theme(text = element_text(size = 8))  # Reducing the font size


# Subset samples groups
is_control_bact <-  rownames(metadata[metadata$group %in% c("01_TC", "03_TB"), ])
is_control_nem <-  rownames(metadata[metadata$group %in% c("01_TC", "02_TN"), ])
is_control_bactNem <-  rownames(metadata[metadata$group %in% c("01_TC", "04_TBN"), ])

meta_control_bact <- metadata[is_control_bact,]
meta_control_nem <- metadata[is_control_nem,]
meta_control_bactNem <- metadata[is_control_bactNem,]

data_control_bact    <- data_df[,is_control_bact]
data_control_bact_10up <- subset(data_control_bact, rowMeans(data_control_bact) >= 10)
matrix1 <- as.matrix(round(data_control_bact_10up, digits = 0))
data_control_nem     <- data_df[,is_control_nem]
data_control_nem_10up <- subset(data_control_nem, rowMeans(data_control_nem) >= 10)
matrix2 <- as.matrix(round(data_control_nem_10up, digits = 0))
data_control_bactNem <- data_df[,is_control_bactNem]
data_control_bactNem_10up <- subset(data_control_bactNem, rowMeans(data_control_bactNem) >= 10)
matrix3 <- as.matrix(round(data_control_bactNem_10up, digits = 0))

# Import data matrices into a DESeq2 objects
dds_cb <- DESeqDataSetFromMatrix(countData = matrix1, 
                              colData = meta_control_bact,
                              design =~ group)

dds_cn <- DESeqDataSetFromMatrix(countData = matrix2, 
                                 colData = meta_control_nem,
                                 design =~ group)

dds_cbn <- DESeqDataSetFromMatrix(countData = matrix3, 
                                 colData = meta_control_bactNem,
                                 design =~ group)

# Calculate size factors 
dds_cb <- estimateSizeFactors(dds_cb)
dds_cn <- estimateSizeFactors(dds_cn)
dds_cbn <- estimateSizeFactors(dds_cbn)

# normalize data
# Regularized logarithmic transformation
rld_cb  <- rlogTransformation(dds_cb)
rld_cn  <- rlogTransformation(dds_cn)
rld_cbn <- rlogTransformation(dds_cbn)

# Hierarchical clustering
makeHCheatmap <- function(rld, metadata){
  dist = dist(t(assay(rld)))
  dist_matrix <- as.matrix(dist)
  row.names(dist_matrix) <- metadata$group
  # Retrieve the name of the dataframe
  df_name <- deparse(substitute(metadata))
  
  # Use gsub to create the new name
  prefix <- gsub("meta_", "", df_name)
  file_name <- paste(prefix, 'HC_plot.png', sep = '_')
  png(file_name)
  pheatmap(dist_matrix, 
           color=colorRampPalette(brewer.pal(n=9,name = "RdYlBu"))(255), 
           clustering_distance_cols = dist, clustering_distance_rows = dist
  )
  dev.off()
}

makeHCheatmap(rld_cb, meta_control_bact)
makeHCheatmap(rld_cn, meta_control_nem)
makeHCheatmap(rld_cbn, meta_control_bactNem)


# PCA plots
makePCA <- function(rld, metadata, ccolor, tcolor) {
  df_name <- deparse(substitute(metadata))
  prefix <- gsub("meta_", "", df_name)
  
  # Use gsub to create the new name
  file_name <- paste(prefix, "PCA_plot.png", sep = '_')
  
  data <- plotPCA(rld, intgroup = c("group"), returnData = TRUE)
  percentVar <- round(100 * attr(data, "percentVar"))
  
  # Define custom color palette for groups
  custom_colors <- c(ccolor, tcolor)
  
  PCA_EucDist <- ggplot(data, aes(x = PC1, y = PC2, color = group)) +
    xlab(paste("PC1 :", percentVar[1], "% variance")) +
    ylab(paste("PC2 :", percentVar[2], "% variance")) +
    ggtitle(paste(prefix, "PCA", sep = ' ')) +
    geom_point(aes(fill = group), shape = 21, size = 12, color = "black") + # Add filled points with black border
    scale_color_manual(values = custom_colors) + # Define point colors
    scale_fill_manual(values = custom_colors) +  # Define fill colors
    theme_bw() +
    theme(legend.position = "right")  # Adjust legend position
  
  # Add encircling ellipses to each group
  PCA_EucDist <- PCA_EucDist +
    stat_ellipse(aes(fill = group), geom = "polygon", level = 0.95, alpha = 0.2, show.legend = FALSE) 
    
  # Save the plot as a PNG file
  png(file_name, width = 600, height = 600)  # Specify width and height as needed
  print(PCA_EucDist)
  dev.off()
}

makePCA(rld_cb, meta_control_bact, "gray", "limegreen")
makePCA(rld_cn, meta_control_nem, "gray", "firebrick1")
makePCA(rld_cbn, meta_control_bactNem, "gray", "dodgerblue2")

# Conduct Diff. Exp. Analyses
diffExp <- function(dds, comp_name){
  dds <- DESeq(dds)
  res <- results(dds)
  resSig <- subset(res, padj < 0.05)
  outfile_name <- paste(comp_name, "q0.05.tsv", sep = '_')
  write.table(resSig, file = outfile_name, sep = "\t", quote = F)
  
  allResfile <- paste0(comp_name, "_allRes.tsv")
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Transcript"
  resdata <- cbind(geneName=rownames(resdata), resdata)
  write.table(resdata, file=allResfile, sep="\t", quote = F, row.names = F)
  
}

diffExp(dds_cb, "control_vs_bacteria")
diffExp(dds_cn, "control_vs_nematode")
diffExp(dds_cbn, "control_vs_bacteria-nematode")


volcanoplot <- function(res, lfcthresh=2, sigthresh=0.1, main='', legendpos="topright", labelsig=TRUE, textcx=1) {
  max_lfc <- max(abs(res$log2FoldChange))
  max_padj <- max(-log10(res$padj))
  
  with(res, plot(log2FoldChange, -log10(padj), pch=20, main=main, cex=0.8, ylim=c(0, 15), xlim=c(-15, 15)))
  with(subset(res, padj<sigthresh), points(log2FoldChange, -log10(padj), pch=20, col="dodgerblue", cex = 1))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(padj), pch=20, col="darkgrey", cex = 1))
  with(subset(res, padj<sigthresh & log2FoldChange > lfcthresh), points(log2FoldChange, -log10(padj), pch=20, col="red", cex = 1))
  with(subset(res, padj<sigthresh & log2FoldChange < -lfcthresh), points(log2FoldChange, -log10(padj), pch=20, col="forestgreen", cex = 1))
  mycol <- c("black", "dodgerblue", "darkgray", "red", "forestgreen")
  legend("topleft", legend=c("Padj >= 0.05", "Padj < 0.05; FC < 2", "Padj > 0.0.5; FC > 2", "Padj < 0.05; FC > 2", "Padj < 0.05; FC < -2"),
         pch=c(19, 19), col=mycol, bg=mycol, pt.cex=1, bty="n", cex=1, text.col="black", horiz=FALSE, inset=c(0.02, 0.02))
  
}

dds_cb <- DESeq(dds_cb)
res_cb <- results(dds_cb)
res_cb_sign <- subset(res_cb, padj < 0.05)

dds_cn <- DESeq(dds_cn)
res_cn <- results(dds_cn)
res_cn_sign <- subset(res_cn, padj < 0.05)

dds_cbn <- DESeq(dds_cbn)
res_cbn <- results(dds_cbn)
res_cbn_sign <- subset(res_cbn, padj < 0.05)

# Make volcano plots
volcanoplot(res_cb, lfcthresh=1, sigthresh=0.05)
volcanoplot(res_cn, lfcthresh=1, sigthresh=0.05)
volcanoplot(res_cbn, lfcthresh=1, sigthresh=0.05)

# Create upset plots
library(UpSetR)

# Function to create a binary matrix for UpSetR
createBinaryMatrix <- function(..., logFC.threshold = 1, upregulated = TRUE) {
  df.list <- list(...)
  names(df.list) <- paste0("Contrast", seq_along(df.list))
  
  # Merge all dataframes into a single binary matrix
  binary.matrix <- Reduce(function(x, y) {
    merge(x, y, by = "gene", all = TRUE)
  }, lapply(names(df.list), function(contrast) {
    df <- as.data.frame(df.list[[contrast]])
    df$gene <- rownames(df)
    if (upregulated) {
      df$present <- ifelse(df$log2FoldChange >= logFC.threshold, 1, 0)
    } else {
      df$present <- ifelse(df$log2FoldChange <= -logFC.threshold, 1, 0)
    }
    df <- df[, c("gene", "present")]
    colnames(df)[2] <- contrast
    return(df)
  }))
  
  # Replace NAs with 0
  binary.matrix[is.na(binary.matrix)] <- 0
  
  return(binary.matrix)
}


# UpSet plots
binary.matrix.up <- createBinaryMatrix(res_cb_sign, res_cn_sign, res_cbn_sign, upregulated = TRUE)
binary.matrix.down <- createBinaryMatrix(res_cb_sign, res_cn_sign, res_cbn_sign, upregulated = FALSE)


pdf("upregulated_features.pdf", width = 11, height = 8.5)
upset(inputTable, sets = samples,
      sets.bar.color = "dodgerblue", matrix.color='black', sets.x.label = "Set Size", point.size = 4,
      line.size = 0.9, mb.ratio = c(0.7, 0.3),  show.numbers = "yes", main.bar.color = barColor,
      number.angles = 20, shade.color = "gray", shade.alpha = 0.5, boxplot.summary = NULL,
      text.scale = 1.2, set_size.angles = 0)
dev.off()


pdf("downregulated_features.pdf", width = 11, height = 8.5)
upset(binary.matrix.down, sets = colnames(binary.matrix)[-1], 
      main.bar.color = "forestgreen", number.angles = 30, point.size = 5, cex =2, line.size = 1)
dev.off()

