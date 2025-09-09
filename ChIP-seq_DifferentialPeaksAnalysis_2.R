library("DESeq2")
library(ggplot2)
library(pheatmap)

setwd("C:/ChIP_seq_analysis")

# peaks x samples

raw_counts <- readRDS("data/ChIP-seq_RawCounts.Rds")
raw_counts <- raw_counts[, -((ncol(raw_counts)-5):ncol(raw_counts))] 

# head(raw_counts)

cell_type = factor(c(rep("squamous",3),rep("adeno", 3),rep("squamous",3),rep("adeno", 3)))

colData <- data.frame(
  condition = factor(c(rep("Normal",6), rep("Tumor",6))),
  cellType = factor(cell_type),
  row.names = colnames(raw_counts)
)

# in alphabetical order
# By default, the first level of each factor is taken as the reference.
# Reference for condition = "Normal"
# Reference for cellType = "adeno"

levels(colData$condition)
levels(colData$cellType)

dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = colData,
  design = ~ condition + cellType
)

dds <- DESeq(dds)

resultsNames(dds) 

# condition_Tumor_vs_Normal: how Tumor differs from Normal on average, regardless of whether samples are squamous or adeno.
# cellType_squamous_vs_adeno: how squamous differs from adeno on average, regardless of whether samples are Normal or Tumor.


# Variance-stabilizing transformation
vsdata <- vst(dds, blind = FALSE)

# PCA
plotPCA(vsdata, intgroup="condition")+ theme_grey()

# Tumor vs Normal
res_condition <- results(dds, contrast = c("condition", "Tumor", "Normal"))

# or use results(dds, name = "condition_Tumor_vs_Normal")

any(is.na(res_condition)) # TRUE

res_condition <- na.omit(res_condition)
res_condition.df <- as.data.frame(res_condition)


sigs_condition.df <- res_condition.df[
  res_condition.df$padj < 0.05 &                        
    abs(res_condition.df$log2FoldChange) > 1 ,
]

# Volcano plot

res_condition.df$significant <- "No"
res_condition.df$significant[res_condition.df$padj < 0.05 & res_condition.df$log2FoldChange > 1] <- "Up"
res_condition.df$significant[res_condition.df$padj < 0.05 & res_condition.df$log2FoldChange < -1] <- "Down"


ggplot(res_condition.df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("No"="grey", "Up"="red", "Down"="blue")) +
  theme_minimal() +
  labs(title="Volcano Plot: Tumor vs Normal", x="Log2 Fold Change", y="-log10(FDR)")

# Squamous vs Adeno:
res_celltype <- results(dds, contrast = c("cellType", "squamous", "adeno"))

any(is.na(res_celltype))

res_celltype <- na.omit(res_celltype)
res_celltype.df <- as.data.frame(res_celltype)

sigs_celltype.df <- res_celltype.df[
  res_celltype.df$padj < 0.05 & abs(res_celltype.df$log2FoldChange) > 1 ,
]

res_celltype.df$significant <- "No"
res_celltype.df$significant[res_celltype.df$padj < 0.05 & res_celltype.df$log2FoldChange > 1] <- "Up"
res_celltype.df$significant[res_celltype.df$padj < 0.05 & res_celltype.df$log2FoldChange < -1] <- "Down"

ggplot(res_celltype.df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("No"="grey", "Up"="red", "Down"="blue")) +
  theme_minimal() +
  labs(title="Volcano Plot: Squamous vs Adeno", x="Log2 Fold Change", y="-log10(FDR)")


# Correlation heatmap

mat <- counts(dds, normalized = TRUE)
cor_mat <- cor(mat, method = "pearson")  

pheatmap(cor_mat,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         annotation_col = colData,
         cluster_rows = T,
         cluster_cols = T,
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

