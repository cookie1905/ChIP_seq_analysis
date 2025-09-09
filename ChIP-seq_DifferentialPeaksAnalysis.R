library("DESeq2")
library(ggplot2)
library(ComplexHeatmap)
library(pheatmap)

setwd("C:/ChIP_seq_analysis")

# peaks x samples

raw_counts <- readRDS("data/ChIP-seq_RawCounts.Rds")
# head(raw_counts)

colData <- data.frame(
  condition = factor(c(rep("Normal",6), rep("Tumor",6), rep("Encode",6))),
  row.names = colnames(raw_counts)
)

levels(colData$condition) 

# in alphabetical order of conditions: "Encode" "Normal" "Tumor" 
# as a result Encode will be used as reference

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = colData,
                              design = ~ condition)

dds <- DESeq(dds)

resultsNames(dds) # gives: "Intercept" "condition_Normal_vs_Encode" "condition_Tumor_vs_Encode" 

# set Normal as reference
colData$condition <- relevel(colData$condition, ref = "Normal")

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds)
resultsNames(dds)


# Variance-stabilizing transformation
vsdata <- vst(dds, blind = FALSE)

# PCA
plotPCA(vsdata, intgroup="condition")+ theme_grey()

# Tumor vs Normal
res <- results(dds, contrast=c("condition","Tumor","Normal"))

# or
# res <- results(dds, name = "condition_Tumor_vs_Normal")

any(is.na(res)) # TRUE

res <- na.omit(res)
res.df <- as.data.frame(res)


sigs.df <- res.df[
    res.df$padj < 0.05 &                        
    abs(res.df$log2FoldChange) > 1 ,
]



res.df$significant <- "No"
res.df$significant[res.df$padj < 0.05 & res.df$log2FoldChange > 1] <- "Up"
res.df$significant[res.df$padj < 0.05 & res.df$log2FoldChange < -1] <- "Down"


ggplot(res.df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("No"="grey", "Up"="red", "Down"="blue")) +
  theme_minimal() +
  labs(title="Volcano Plot: Tumor vs Normal", x="Log2 Fold Change", y="-log10(FDR)")


mat <- counts(dds, normalized = TRUE)
cor_mat <- cor(mat, method = "pearson")  # or method = "spearman"

pheatmap(cor_mat,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         annotation_col = colData,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 8,
         fontsize_col = 8,
         border_color = NA)

# Extract matrix for significant peaks only

mat <- counts(dds, normalized = TRUE)[rownames(sigs.df), ]
mat.z <- t(apply(mat, 1, scale)) # Compute row variance
colnames(mat.z) <- colnames(mat) # apply() column names often get replaced with numbers (1,2,...).


pheatmap(mat.z,
         color = colorRampPalette(c("blue","white","red"))(100),
         cluster_rows = TRUE,       
         cluster_cols = TRUE,       
         show_rownames = FALSE,    
         annotation_col = colData)

# show top peaks

top_peaks <- rownames(sigs.df[order(abs(sigs.df$log2FoldChange), decreasing = TRUE), ])[1:30]
mat_top <- counts(dds, normalized = TRUE)[top_peaks, ]
mat_top.z <- t(apply(mat_top, 1, scale))
colnames(mat_top.z) <- colnames(mat_top)

pheatmap(mat_top.z,
         color = colorRampPalette(c("blue","white","red"))(100),
         cluster_rows = TRUE,       
         cluster_cols = TRUE,       
         show_rownames = TRUE,    
         annotation_col = colData,
         border_color = NA,
         fontsize_row = 6,
         fontsize_col = 6)   

