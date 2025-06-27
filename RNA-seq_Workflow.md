# RNA-seq Differential Expression and Enrichment Analysis Workflow

## Load Libraries
```r
library(DESeq2)
library(edgeR)
library(limma)
library(dplyr)
library(enrichR)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(tidyverse)
library(pheatmap)
library(qs)
library(clusterProfiler)
library(estimate)
library(ggpubr)
```

## Set Working Directory
```r
setwd("/path/to/project/directory")
```

## Load Raw Count Matrix and Metadata
```r
count_data <- read.table("Matrix_count.txt", header = TRUE, row.names = 1)
count_data <- round(count_data, digits = 0)
sample_info <- data.frame(condition = rep(c('WT','SH'), times = c(3,5)),
                          row.names = colnames(count_data))
sample_info$condition <- factor(sample_info$condition, levels = c("WT", "SH"))
```

## Create DESeq2 Dataset and Run Differential Expression
```r
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~ condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "SH", "WT"),
               independentFiltering = TRUE, alpha = 0.05)
```

## Save DEGs Results
```r
qsave(res, "DEGs.result_SH-WT.qs")
write.csv(as.data.frame(res), file = "DEGs_SH-WT.csv")
```

## Volcano Plot for DEGs
```r
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$significance <- with(res_df, ifelse(pvalue < 0.05 & abs(log2FoldChange) > 1,
                                           ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"),
                                           "Not significant"))
res_df$significance <- factor(res_df$significance,
                              levels = c("Upregulated", "Downregulated", "Not significant"))

label_genes <- c("CCL11","CCL4","CCL5","CCL8","CCR9","CXCL10","CXCL11","CXCL12","CXCL5","XCL1")

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significance)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_point(data = subset(res_df, gene %in% label_genes),
             aes(x = log2FoldChange, y = -log10(pvalue)),
             colour = "#E4C755", size = 2) +
  scale_color_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", "Not significant" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano Plot", x = "log2(FoldChange)", y = "-log10(p-value)") +
  theme_classic(base_size = 14) +
  geom_label_repel(data = subset(res_df, gene %in% label_genes),
                   aes(label = gene), size = 4)
```

---

## KEGG Pathway Enrichment Analysis
```r
library(clusterProfiler)
library(org.Mm.eg.db)

# Extract significantly downregulated genes
down_genes <- rownames(subset(res, padj < 0.05 & log2FoldChange < -1))

# Convert gene symbols to Entrez IDs
down_entrez <- bitr(down_genes, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)

# KEGG enrichment
kegg_result <- enrichKEGG(gene = down_entrez$ENTREZID,
                          organism = "mmu",
                          pvalueCutoff = 0.05)

# Make results human-readable
kegg_result <- setReadable(kegg_result, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

# Barplot of top KEGG pathways
barplot(kegg_result, showCategory = 10, title = "KEGG Pathway Enrichment") +
  theme_classic(base_size = 14)

# Save results
write.csv(as.data.frame(kegg_result), file = "kegg_downregulated.csv")
```

---

## GSEA (Gene Set Enrichment Analysis)
```r
library(msigdbr)
library(enrichplot)

# Prepare gene list
gene_list <- res$log2FoldChange
names(gene_list) <- rownames(res)
gene_list <- sort(gene_list, decreasing = TRUE)

# Load gene sets (e.g., Hallmark)
m_df <- msigdbr(species = "Mus musculus", category = "H")
m_df$gs_name <- gsub("^HALLMARK_", "", m_df$gs_name)
TERM2GENE <- m_df[, c("gs_name", "gene_symbol")]

# Run GSEA
gsea_result <- GSEA(geneList = gene_list,
                    TERM2GENE = TERM2GENE,
                    pvalueCutoff = 0.05,
                    minGSSize = 15)

# Plot top GSEA terms
dotplot(gsea_result, showCategory = 20, title = "GSEA - Hallmark Pathways") +
  theme_classic(base_size = 14)

# Save GSEA results
write.csv(as.data.frame(gsea_result), file = "GSEA_Hallmark.csv")
```

---

## Heatmap of Immune-Related Genes
```r
library(pheatmap)

# Extract expression matrix
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
expr_matrix <- assay(vsd)

# Define immune-related genes
genes_of_interest <- c("Cd3e","Cd4","Cd8a","Ifng","Tnf","Tbx21",
                       "Ccl2","Ccl5","Cxcl10","S100a8","S100a9",
                       "Arg1","Nos2","Il10","Stat3","Cd274")

# Filter and plot
filtered_expr <- expr_matrix[rownames(expr_matrix) %in% genes_of_interest, ]
pheatmap(filtered_expr,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = sample_info,
         show_rownames = TRUE,
         main = "Immune-Related Gene Expression")
```

---
