setwd("PATH")

library(DESeq2)
library(stats)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)


##### Differential gene expression paired analysis #####

### Tumor ###

countTable = read.table("Tumor_pb_matrix.csv",header=TRUE, row.names=1, check.names = FALSE,sep = ',')
designMatrix = read.table('Tumor_metadata.csv', header=TRUE, row.names=1,sep = ',')
libType = c("oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor")
condition = factor(designMatrix$Condition)
patient = factor(designMatrix$Patient)
condition <- relevel(condition, "Pre")
experiment_design=data.frame(row.names = colnames(countTable),condition,patient,libType)

colnames(countTable) == rownames(experiment_design)

cds <- DESeqDataSetFromMatrix(countData =countTable, colData=experiment_design, design=~patient + condition)
cds_DESeqED <- DESeq(cds,parallel = TRUE)
res <- results(cds_DESeqED,parallel = TRUE,alpha = 0.05, pAdjustMethod = "BH")
write.table(res,file = "Tumor_Post_vs_Pre.differentialExpression.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
resSig <- subset(res, padj < 0.05)
write.table(resSig,file = "Tumor_Post_vs_Pre.sigDEGs.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
normalizedReadCounts = counts(cds_DESeqED,normalized=TRUE)
write.table(normalizedReadCounts,file = "Tumor_Post_vs_Pre.normalizedCounts.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")


## Volcano plot with labels

Tumor_Post_vs_Pre <- read.table("Tumor_Post_vs_Pre.differentialExpression.txt", header=TRUE, sep = '\t')
Tumor_Post_vs_Pre <- Tumor_Post_vs_Pre[complete.cases(Tumor_Post_vs_Pre), ]

Tumor_Post_vs_Pre$Type <- "NONE"
Tumor_Post_vs_Pre$Type[Tumor_Post_vs_Pre$padj < 0.05 & Tumor_Post_vs_Pre$log2FoldChange < -1] <- "DOWN"
Tumor_Post_vs_Pre$Type[Tumor_Post_vs_Pre$padj < 0.05 & Tumor_Post_vs_Pre$log2FoldChange > 1] <- "UP"
cols <- c("NONE" = "#474657", "DOWN" = "#157ded", "UP" = "#f2800f")

Tumor_Post_vs_Pre$genelabels <- ""

# To make a column with labels according to the threshold
Tumor_Post_vs_Pre$label <- ifelse(
  Tumor_Post_vs_Pre$padj < 0.05,
  Tumor_Post_vs_Pre$X,
  NA
)

png("DESeq2_Tumor_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)

ggplot(Tumor_Post_vs_Pre, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 1, alpha = 0.75, na.rm = TRUE) +
  geom_text_repel(aes(label = label), size = 8, max.overlaps = 30, na.rm = TRUE) +
  theme_bw(base_size = 18) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20)) +
  xlab(expression(log[2]("Fold Change"))) + 
  ylab(expression(-log[10]("padj"))) +
  geom_hline(yintercept = 1.3, colour = "red4", linetype = "dashed") + 
  geom_vline(xintercept = 0.0, colour = "#474657", linetype = "dashed") +
  scale_y_continuous(trans = "log1p") +
  ggtitle("DESeq2 Volcano Plot - Tumor cells Surgery vs Biopsy")

dev.off()





### T cells ###


countTable = read.table("T_cells_pb_matrix.csv",header=TRUE, row.names=1, check.names = FALSE,sep = ',')
designMatrix = read.table('T_cells_metadata.csv', header=TRUE, row.names=1,sep = ',')
libType = c("oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor")
condition = factor(designMatrix$Condition)
patient = factor(designMatrix$Patient)
condition <- relevel(condition, "Pre")
experiment_design=data.frame(row.names = colnames(countTable),condition,patient,libType)

colnames(countTable) == rownames(experiment_design)

cds <- DESeqDataSetFromMatrix(countData =countTable, colData=experiment_design, design=~patient + condition)
cds_DESeqED <- DESeq(cds,parallel = TRUE)
res <- results(cds_DESeqED,parallel = TRUE,alpha = 0.05, pAdjustMethod = "BH")
write.table(res,file = "T_cells_Post_vs_Pre.differentialExpression.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
resSig <- subset(res, padj < 0.05)
write.table(resSig,file = "T_cells_Post_vs_Pre.sigDEGs.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
normalizedReadCounts = counts(cds_DESeqED,normalized=TRUE)
write.table(normalizedReadCounts,file = "T_cells_Post_vs_Pre.normalizedCounts.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")


## volcano plot with labels

T_cells_Post_vs_Pre <- read.table("T_cells_Post_vs_Pre.differentialExpression.txt", header=TRUE, sep = '\t')
T_cells_Post_vs_Pre <- T_cells_Post_vs_Pre[complete.cases(T_cells_Post_vs_Pre), ]

T_cells_Post_vs_Pre$Type <- "NONE"
T_cells_Post_vs_Pre$Type[T_cells_Post_vs_Pre$padj < 0.05 & T_cells_Post_vs_Pre$log2FoldChange < -1] <- "DOWN"
T_cells_Post_vs_Pre$Type[T_cells_Post_vs_Pre$padj < 0.05 & T_cells_Post_vs_Pre$log2FoldChange > 1] <- "UP"
cols <- c("NONE" = "#474657", "DOWN" = "#157ded", "UP" = "#f2800f")

T_cells_Post_vs_Pre$genelabels <- ""

# To make a column with labels according to the threshold
T_cells_Post_vs_Pre$label <- ifelse(
  T_cells_Post_vs_Pre$padj < 0.05,
  T_cells_Post_vs_Pre$X,
  NA
)

png("DESeq2_T_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)

ggplot(T_cells_Post_vs_Pre, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 1, alpha = 0.75, na.rm = TRUE) +
  geom_text_repel(aes(label = label), size = 8, max.overlaps = 30, na.rm = TRUE) +
  theme_bw(base_size = 18) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20)) +
  xlab(expression(log[2]("Fold Change"))) + 
  ylab(expression(-log[10]("padj"))) +
  geom_hline(yintercept = 1.3, colour = "red4", linetype = "dashed") + 
  geom_vline(xintercept = 0.0, colour = "#474657", linetype = "dashed") +
  scale_y_continuous(trans = "log1p") +
  ggtitle("DESeq2 Volcano Plot - T cells Surgery vs Biopsy")

dev.off()





### Fibroblasts ###


countTable = read.table("Fibroblasts_pb_matrix.csv",header=TRUE, row.names=1, check.names = FALSE,sep = ',')
designMatrix = read.table('Fibroblasts_metadata.csv', header=TRUE, row.names=1,sep = ',')
libType = c("oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor")
condition = factor(designMatrix$Condition)
patient = factor(designMatrix$Patient)
condition <- relevel(condition, "Pre")
experiment_design=data.frame(row.names = colnames(countTable),condition,patient,libType)

colnames(countTable) == rownames(experiment_design)

cds <- DESeqDataSetFromMatrix(countData =countTable, colData=experiment_design, design=~patient + condition)
cds_DESeqED <- DESeq(cds,parallel = TRUE)
res <- results(cds_DESeqED,parallel = TRUE,alpha = 0.05, pAdjustMethod = "BH")
write.table(res,file = "Fibroblasts_Post_vs_Pre.differentialExpression.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
resSig <- subset(res, padj < 0.05)
write.table(resSig,file = "Fibroblasts_Post_vs_Pre.sigDEGs.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
normalizedReadCounts = counts(cds_DESeqED,normalized=TRUE)
write.table(normalizedReadCounts,file = "Fibroblasts_Post_vs_Pre.normalizedCounts.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")


## Volcano plot with labels

Fibroblasts_Post_vs_Pre <- read.table("Fibroblasts_Post_vs_Pre.differentialExpression.txt", header=TRUE, sep = '\t')
Fibroblasts_Post_vs_Pre <- Fibroblasts_Post_vs_Pre[complete.cases(Fibroblasts_Post_vs_Pre), ]

Fibroblasts_Post_vs_Pre$Type <- "NONE"
Fibroblasts_Post_vs_Pre$Type[Fibroblasts_Post_vs_Pre$padj < 0.05 & Fibroblasts_Post_vs_Pre$log2FoldChange < -1] <- "DOWN"
Fibroblasts_Post_vs_Pre$Type[Fibroblasts_Post_vs_Pre$padj < 0.05 & Fibroblasts_Post_vs_Pre$log2FoldChange > 1] <- "UP"
cols <- c("NONE" = "#474657", "DOWN" = "#157ded", "UP" = "#f2800f")

Fibroblasts_Post_vs_Pre$genelabels <- ""

# To make a column with labels according to the threshold
Fibroblasts_Post_vs_Pre$label <- ifelse(
  Fibroblasts_Post_vs_Pre$padj < 0.05,
  Fibroblasts_Post_vs_Pre$X,
  NA
)

png("DESeq2_Fibroblasts_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)

ggplot(Fibroblasts_Post_vs_Pre, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 1, alpha = 0.75, na.rm = TRUE) +
  geom_text_repel(aes(label = label), size = 8, max.overlaps = 30, na.rm = TRUE) +
  theme_bw(base_size = 18) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20)) +
  xlab(expression(log[2]("Fold Change"))) + 
  ylab(expression(-log[10]("padj"))) +
  geom_hline(yintercept = 1.3, colour = "red4", linetype = "dashed") + 
  geom_vline(xintercept = 0.0, colour = "#474657", linetype = "dashed") +
  scale_y_continuous(trans = "log1p") +
  ggtitle("DESeq2 Volcano Plot - Fibroblasts cells Surgery vs Biopsy") 

dev.off()



### Myeloid ###


countTable = read.table("Myeloid_pb_matrix.csv",header=TRUE, row.names=1, check.names = FALSE,sep = ',')
designMatrix = read.table('Myeloid_metadata.csv', header=TRUE, row.names=1,sep = ',')
libType = c("oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor")
condition = factor(designMatrix$Condition)
patient = factor(designMatrix$Patient)
condition <- relevel(condition, "Pre")
experiment_design=data.frame(row.names = colnames(countTable),condition,patient,libType)

colnames(countTable) == rownames(experiment_design)

cds <- DESeqDataSetFromMatrix(countData =countTable, colData=experiment_design, design=~patient + condition)
cds_DESeqED <- DESeq(cds,parallel = TRUE)
res <- results(cds_DESeqED,parallel = TRUE,alpha = 0.05, pAdjustMethod = "BH")
write.table(res,file = "Myeloid_Post_vs_Pre.differentialExpression.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
resSig <- subset(res, padj < 0.05)
write.table(resSig,file = "Myeloid_Post_vs_Pre.sigDEGs.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
normalizedReadCounts = counts(cds_DESeqED,normalized=TRUE)
write.table(normalizedReadCounts,file = "Myeloid_Post_vs_Pre.normalizedCounts.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")


## Volcano plot with labels

Myeloid_Post_vs_Pre <- read.table("Myeloid_Post_vs_Pre.differentialExpression.txt", header=TRUE, sep = '\t')
Myeloid_Post_vs_Pre <- Myeloid_Post_vs_Pre[complete.cases(Myeloid_Post_vs_Pre), ]

Myeloid_Post_vs_Pre$Type <- "NONE"
Myeloid_Post_vs_Pre$Type[Myeloid_Post_vs_Pre$padj < 0.05 & Myeloid_Post_vs_Pre$log2FoldChange < -1] <- "DOWN"
Myeloid_Post_vs_Pre$Type[Myeloid_Post_vs_Pre$padj < 0.05 & Myeloid_Post_vs_Pre$log2FoldChange > 1] <- "UP"
cols <- c("NONE" = "#474657", "DOWN" = "#157ded", "UP" = "#f2800f")

Myeloid_Post_vs_Pre$genelabels <- ""

# To make a column with labels according to the threshold
Myeloid_Post_vs_Pre$label <- ifelse(
  Myeloid_Post_vs_Pre$padj < 0.05,
  Myeloid_Post_vs_Pre$X,
  NA
)

png("DESeq2_Myeloid_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)

ggplot(Myeloid_Post_vs_Pre, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 1, alpha = 0.75, na.rm = TRUE) +
  geom_text_repel(aes(label = label), size = 8, max.overlaps = 30, na.rm = TRUE) +
  theme_bw(base_size = 18) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20)) +
  xlab(expression(log[2]("Fold Change"))) + 
  ylab(expression(-log[10]("padj"))) +
  geom_hline(yintercept = 1.3, colour = "red4", linetype = "dashed") + 
  geom_vline(xintercept = 0.0, colour = "#474657", linetype = "dashed") +
  scale_y_continuous(trans = "log1p") +
  ggtitle("DESeq2 Volcano Plot - Myeloid cells Surgery vs Biopsy")

dev.off()



### Endothelial cells ###


countTable = read.table("Endothelial_pb_matrix.csv",header=TRUE, row.names=1, check.names = FALSE,sep = ',')
designMatrix = read.table('Endothelial_metadata.csv', header=TRUE, row.names=1,sep = ',')
libType = c("oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor", "oneFactor")
condition = factor(designMatrix$Condition)
patient = factor(designMatrix$Patient)
condition <- relevel(condition, "Pre")
experiment_design=data.frame(row.names = colnames(countTable),condition,patient,libType)

colnames(countTable) == rownames(experiment_design)

cds <- DESeqDataSetFromMatrix(countData =countTable, colData=experiment_design, design=~patient + condition)
cds_DESeqED <- DESeq(cds,parallel = TRUE)
res <- results(cds_DESeqED,parallel = TRUE,alpha = 0.05, pAdjustMethod = "BH")
write.table(res,file = "Endothelial_Post_vs_Pre.differentialExpression.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
resSig <- subset(res, padj < 0.05)
write.table(resSig,file = "Endothelial_Post_vs_Pre.sigDEGs.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
normalizedReadCounts = counts(cds_DESeqED,normalized=TRUE)
write.table(normalizedReadCounts,file = "Endothelial_Post_vs_Pre.normalizedCounts.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")


## Volcano plot with labels

Endothelial_Post_vs_Pre <- read.table("Endothelial_Post_vs_Pre.differentialExpression.txt", header=TRUE, sep = '\t')
Endothelial_Post_vs_Pre <- Endothelial_Post_vs_Pre[complete.cases(Endothelial_Post_vs_Pre), ]

Endothelial_Post_vs_Pre$Type <- "NONE"
Endothelial_Post_vs_Pre$Type[Endothelial_Post_vs_Pre$padj < 0.05 & Endothelial_Post_vs_Pre$log2FoldChange < -1] <- "DOWN"
Endothelial_Post_vs_Pre$Type[Endothelial_Post_vs_Pre$padj < 0.05 & Endothelial_Post_vs_Pre$log2FoldChange > 1] <- "UP"
cols <- c("NONE" = "#474657", "DOWN" = "#157ded", "UP" = "#f2800f")

Endothelial_Post_vs_Pre$genelabels <- ""

# To make a column with labels according to the threshold
Endothelial_Post_vs_Pre$label <- ifelse(
  Endothelial_Post_vs_Pre$padj < 0.05,
  Endothelial_Post_vs_Pre$X,
  NA
)

png("DESeq2_Endothelial_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)

ggplot(Endothelial_Post_vs_Pre, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 1, alpha = 0.75, na.rm = TRUE) +
  geom_text_repel(aes(label = label), size = 8, max.overlaps = 30, na.rm = TRUE) +
  theme_bw(base_size = 18) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20)) +
  xlab(expression(log[2]("Fold Change"))) + 
  ylab(expression(-log[10]("padj"))) +
  geom_hline(yintercept = 1.3, colour = "red4", linetype = "dashed") + 
  geom_vline(xintercept = 0.0, colour = "#474657", linetype = "dashed") +
  scale_y_continuous(trans = "log1p") +
  ggtitle("DESeq2 Volcano Plot - Endothelial cells Surgery vs Biopsy")

dev.off()


### B cells ###


countTable = read.table("B_cell_pb_matrix.csv",header=TRUE, row.names=1, check.names = FALSE,sep = ',')
designMatrix = read.table('B_cell_metadata.csv', header=TRUE, row.names=1,sep = ',')
libType = c("oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor")
condition = factor(designMatrix$Condition)
patient = factor(designMatrix$Patient)
condition <- relevel(condition, "Pre")
experiment_design=data.frame(row.names = colnames(countTable),condition,patient,libType)

colnames(countTable) == rownames(experiment_design)

cds <- DESeqDataSetFromMatrix(countData =countTable, colData=experiment_design, design=~patient + condition)
cds_DESeqED <- DESeq(cds,parallel = TRUE)
res <- results(cds_DESeqED,parallel = TRUE,alpha = 0.05, pAdjustMethod = "BH")
write.table(res,file = "Bcell_Post_vs_Pre.differentialExpression.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
resSig <- subset(res, padj < 0.05)
write.table(resSig,file = "Bcell_Post_vs_Pre.sigDEGs.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
normalizedReadCounts = counts(cds_DESeqED,normalized=TRUE)
write.table(normalizedReadCounts,file = "Bcell_Post_vs_Pre.normalizedCounts.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")


## Volcano plot with labels

Bcell_Post_vs_Pre <- read.table("/Users/dvalcarcel/Documents/TFM/Analysis_16um/RCTD_deconvolution/RCTC_snRNAseq/v2/PB_final_annotation/Bcell_Post_vs_Pre.differentialExpression.txt", header=TRUE, sep = '\t')
Bcell_Post_vs_Pre <- Bcell_Post_vs_Pre[complete.cases(Bcell_Post_vs_Pre), ]

Bcell_Post_vs_Pre$Type <- "NONE"
Bcell_Post_vs_Pre$Type[Bcell_Post_vs_Pre$padj < 0.05 & Bcell_Post_vs_Pre$log2FoldChange < -1] <- "DOWN"
Bcell_Post_vs_Pre$Type[Bcell_Post_vs_Pre$padj < 0.05 & Bcell_Post_vs_Pre$log2FoldChange > 1] <- "UP"
cols <- c("NONE" = "#474657", "DOWN" = "#157ded", "UP" = "#f2800f")

Bcell_Post_vs_Pre$genelabels <- ""

# To make a column with labels according to the threshold
Bcell_Post_vs_Pre$label <- ifelse(
  Bcell_Post_vs_Pre$padj < 0.05,
  Bcell_Post_vs_Pre$X,
  NA
)

png("DESeq2_B_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)

ggplot(Bcell_Post_vs_Pre, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 1, alpha = 0.75, na.rm = TRUE) +
  geom_text_repel(aes(label = label), size = 8, max.overlaps = 30, na.rm = TRUE) +
  theme_bw(base_size = 18) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20)) +
  xlab(expression(log[2]("Fold Change"))) + 
  ylab(expression(-log[10]("padj"))) +
  geom_hline(yintercept = 1.3, colour = "red4", linetype = "dashed") + 
  geom_vline(xintercept = 0.0, colour = "#474657", linetype = "dashed") +
  scale_y_continuous(trans = "log1p") +
  ggtitle("DESeq2 Volcano Plot - B cells Surgery vs Biopsy") 

dev.off()



### Adipocytes ###


countTable = read.table("Adipocytes_pb_matrix.csv",header=TRUE, row.names=1, check.names = FALSE,sep = ',')
designMatrix = read.table('Adipocytes_metadata.csv', header=TRUE, row.names=1,sep = ',')
libType = c("oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor")
condition = factor(designMatrix$Condition)
patient = factor(designMatrix$Patient)
condition <- relevel(condition, "Pre")
experiment_design=data.frame(row.names = colnames(countTable),condition,patient,libType)

colnames(countTable) == rownames(experiment_design)

cds <- DESeqDataSetFromMatrix(countData =countTable, colData=experiment_design, design=~patient + condition)
cds_DESeqED <- DESeq(cds,parallel = TRUE)
res <- results(cds_DESeqED,parallel = TRUE,alpha = 0.05, pAdjustMethod = "BH")
write.table(res,file = "Adypocites_Post_vs_Pre.differentialExpression.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
resSig <- subset(res, padj < 0.05)
write.table(resSig,file = "Adypocites_Post_vs_Pre.sigDEGs.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
normalizedReadCounts = counts(cds_DESeqED,normalized=TRUE)
write.table(normalizedReadCounts,file = "Adypocites_Post_vs_Pre.normalizedCounts.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")


## Volcano plot with labels

Adypocites_Post_vs_Pre <- read.table("Adypocites_Post_vs_Pre.differentialExpression.txt", header=TRUE, sep = '\t')
Adypocites_Post_vs_Pre <- Adypocites_Post_vs_Pre[complete.cases(Adypocites_Post_vs_Pre), ]

Adypocites_Post_vs_Pre$Type <- "NONE"
Adypocites_Post_vs_Pre$Type[Adypocites_Post_vs_Pre$padj < 0.05 & Adypocites_Post_vs_Pre$log2FoldChange < -1] <- "DOWN"
Adypocites_Post_vs_Pre$Type[Adypocites_Post_vs_Pre$padj < 0.05 & Adypocites_Post_vs_Pre$log2FoldChange > 1] <- "UP"
cols <- c("NONE" = "#474657", "DOWN" = "#157ded", "UP" = "#f2800f")

Adypocites_Post_vs_Pre$genelabels <- ""

# To make a column with labels according to the threshold
Adypocites_Post_vs_Pre$label <- ifelse(
  Adypocites_Post_vs_Pre$padj < 0.05,
  Adypocites_Post_vs_Pre$X,
  NA
)



png("DESeq2_Adypocites_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)

ggplot(Adypocites_Post_vs_Pre, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 1, alpha = 0.75, na.rm = TRUE) +
  geom_text_repel(aes(label = label), size = 8, max.overlaps = 30, na.rm = TRUE) +
  theme_bw(base_size = 18) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20)) +
  xlab(expression(log[2]("Fold Change"))) + 
  ylab(expression(-log[10]("padj"))) +
  geom_hline(yintercept = 1.3, colour = "red4", linetype = "dashed") + 
  geom_vline(xintercept = 0.0, colour = "#474657", linetype = "dashed") +
  scale_y_continuous(trans = "log1p") +
  ggtitle("DESeq2 Volcano Plot - Adypocites Surgery vs Biopsy") 

dev.off()



### Lymphatic ###

countTable = read.table("Lymphatic_pb_matrix.csv",header=TRUE, row.names=1, check.names = FALSE,sep = ',')
designMatrix = read.table('Lymphatic_metadata.csv', header=TRUE, row.names=1,sep = ',')
libType = c("oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor")
condition = factor(designMatrix$Condition)
patient = factor(designMatrix$Patient)
condition <- relevel(condition, "Pre")
experiment_design=data.frame(row.names = colnames(countTable),condition,patient,libType)

colnames(countTable) == rownames(experiment_design)

cds <- DESeqDataSetFromMatrix(countData =countTable, colData=experiment_design, design=~patient + condition)
cds_DESeqED <- DESeq(cds,parallel = TRUE)
res <- results(cds_DESeqED,parallel = TRUE,alpha = 0.05, pAdjustMethod = "BH")
write.table(res,file = "Lymphatic_Post_vs_Pre.differentialExpression.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
resSig <- subset(res, padj < 0.05)
write.table(resSig,file = "Lymphatic_Post_vs_Pre.sigDEGs.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
normalizedReadCounts = counts(cds_DESeqED,normalized=TRUE)
write.table(normalizedReadCounts,file = "Lymphatic_Post_vs_Pre.normalizedCounts.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")


## Volcano plot with labels

Lymphatic_Post_vs_Pre <- read.table("Lymphatic_Post_vs_Pre.differentialExpression.txt", header=TRUE, sep = '\t')
Lymphatic_Post_vs_Pre <- Lymphatic_Post_vs_Pre[complete.cases(Lymphatic_Post_vs_Pre), ]

Lymphatic_Post_vs_Pre$Type <- "NONE"
Lymphatic_Post_vs_Pre$Type[Lymphatic_Post_vs_Pre$padj < 0.05 & Lymphatic_Post_vs_Pre$log2FoldChange < -1] <- "DOWN"
Lymphatic_Post_vs_Pre$Type[Lymphatic_Post_vs_Pre$padj < 0.05 & Lymphatic_Post_vs_Pre$log2FoldChange > 1] <- "UP"
cols <- c("NONE" = "#474657", "DOWN" = "#157ded", "UP" = "#f2800f")

Lymphatic_Post_vs_Pre$genelabels <- ""

# To make a column with labels according to the threshold
Lymphatic_Post_vs_Pre$label <- ifelse(
  Lymphatic_Post_vs_Pre$padj < 0.05,
  Lymphatic_Post_vs_Pre$X,
  NA
)


png("DESeq2_Lymphatic_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)

ggplot(Lymphatic_Post_vs_Pre, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 1, alpha = 0.75, na.rm = TRUE) +
  geom_text_repel(aes(label = label), size = 8, max.overlaps = 30, na.rm = TRUE) +
  theme_bw(base_size = 18) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20)) +
  xlab(expression(log[2]("Fold Change"))) + 
  ylab(expression(-log[10]("padj"))) +
  geom_hline(yintercept = 1.3, colour = "red4", linetype = "dashed") + 
  geom_vline(xintercept = 0.0, colour = "#474657", linetype = "dashed") +
  scale_y_continuous(trans = "log1p") +
  ggtitle("DESeq2 Volcano Plot - Lymphatic Surgery vs Biopsy") 

dev.off()





### Plasmablasts ###

countTable = read.table("Plasmablasts_pb_matrix.csv",header=TRUE, row.names=1, check.names = FALSE,sep = ',')
designMatrix = read.table('Plasmablasts_metadata.csv', header=TRUE, row.names=1,sep = ',')
libType = c("oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor")
condition = factor(designMatrix$Condition)
patient = factor(designMatrix$Patient)
condition <- relevel(condition, "Pre")
experiment_design=data.frame(row.names = colnames(countTable),condition,patient,libType)

colnames(countTable) == rownames(experiment_design)

cds <- DESeqDataSetFromMatrix(countData =countTable, colData=experiment_design, design=~patient + condition)
cds_DESeqED <- DESeq(cds,parallel = TRUE)
res <- results(cds_DESeqED,parallel = TRUE,alpha = 0.05, pAdjustMethod = "BH")
write.table(res,file = "Plasmablasts_Post_vs_Pre.differentialExpression.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
resSig <- subset(res, padj < 0.05)
write.table(resSig,file = "Plasmablasts_Post_vs_Pre.sigDEGs.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
normalizedReadCounts = counts(cds_DESeqED,normalized=TRUE)
write.table(normalizedReadCounts,file = "Plasmablasts_Post_vs_Pre.normalizedCounts.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")


## Volcano plot with labels

Plasmablasts_Post_vs_Pre <- read.table("Plasmablasts_Post_vs_Pre.differentialExpression.txt", header=TRUE, sep = '\t')
Plasmablasts_Post_vs_Pre <- Plasmablasts_Post_vs_Pre[complete.cases(Plasmablasts_Post_vs_Pre), ]

Plasmablasts_Post_vs_Pre$Type <- "NONE"
Plasmablasts_Post_vs_Pre$Type[Plasmablasts_Post_vs_Pre$padj < 0.05 & Plasmablasts_Post_vs_Pre$log2FoldChange < -1] <- "DOWN"
Plasmablasts_Post_vs_Pre$Type[Plasmablasts_Post_vs_Pre$padj < 0.05 & Plasmablasts_Post_vs_Pre$log2FoldChange > 1] <- "UP"
cols <- c("NONE" = "#474657", "DOWN" = "#157ded", "UP" = "#f2800f")

Plasmablasts_Post_vs_Pre$genelabels <- ""

# To make a column with labels according to the threshold
Plasmablasts_Post_vs_Pre$label <- ifelse(
  Plasmablasts_Post_vs_Pre$padj < 0.05,
  Plasmablasts_Post_vs_Pre$X,
  NA
)



png("DESeq2_Plasmablasts_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)

ggplot(Plasmablasts_Post_vs_Pre, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 1, alpha = 0.75, na.rm = TRUE) +
  geom_text_repel(aes(label = label), size = 8, max.overlaps = 30, na.rm = TRUE) +
  theme_bw(base_size = 18) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20)) +
  xlab(expression(log[2]("Fold Change"))) + 
  ylab(expression(-log[10]("padj"))) +
  geom_hline(yintercept = 1.3, colour = "red4", linetype = "dashed") + 
  geom_vline(xintercept = 0.0, colour = "#474657", linetype = "dashed") +
  scale_y_continuous(trans = "log1p") +
  ggtitle("DESeq2 Volcano Plot - Plasmablasts Surgery vs Biopsy") 

dev.off()






### Mast cells ###

countTable = read.table("Mast_cell_pb_matrix.csv",header=TRUE, row.names=1, check.names = FALSE,sep = ',')
designMatrix = read.table('Mast_cell_metadata.csv', header=TRUE, row.names=1,sep = ',')
libType = c("oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor")
condition = factor(designMatrix$Condition)
patient = factor(designMatrix$Patient)
condition <- relevel(condition, "Pre")
experiment_design=data.frame(row.names = colnames(countTable),condition,patient,libType)

colnames(countTable) == rownames(experiment_design)

cds <- DESeqDataSetFromMatrix(countData =countTable, colData=experiment_design, design=~patient + condition)
cds_DESeqED <- DESeq(cds,parallel = TRUE)
res <- results(cds_DESeqED,parallel = TRUE,alpha = 0.05, pAdjustMethod = "BH")
write.table(res,file = "Mast_Post_vs_Pre.differentialExpression.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
resSig <- subset(res, padj < 0.05)
write.table(resSig,file = "Mast_Post_vs_Pre.sigDEGs.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
normalizedReadCounts = counts(cds_DESeqED,normalized=TRUE)
write.table(normalizedReadCounts,file = "Mast_Post_vs_Pre.normalizedCounts.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")



## Volcano plot with labels

Mast_Post_vs_Pre <- read.table("Mast_Post_vs_Pre.differentialExpression.txt", header=TRUE, sep = '\t')
Mast_Post_vs_Pre <- Mast_Post_vs_Pre[complete.cases(Mast_Post_vs_Pre), ]

Mast_Post_vs_Pre$Type <- "NONE"
Mast_Post_vs_Pre$Type[Mast_Post_vs_Pre$padj < 0.05 & Mast_Post_vs_Pre$log2FoldChange < -1] <- "DOWN"
Mast_Post_vs_Pre$Type[Mast_Post_vs_Pre$padj < 0.05 & Mast_Post_vs_Pre$log2FoldChange > 1] <- "UP"
cols <- c("NONE" = "#474657", "DOWN" = "#157ded", "UP" = "#f2800f")

Mast_Post_vs_Pre$genelabels <- ""

# To make a column with labels according to the threshold
Mast_Post_vs_Pre$label <- ifelse(
  Mast_Post_vs_Pre$padj < 0.05,
  Mast_Post_vs_Pre$X,
  NA
)



png("DESeq2_Mast_cell_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)

ggplot(Mast_Post_vs_Pre, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 1, alpha = 0.75, na.rm = TRUE) +
  geom_text_repel(aes(label = label), size = 8, max.overlaps = 30, na.rm = TRUE) +
  theme_bw(base_size = 18) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20)) +
  xlab(expression(log[2]("Fold Change"))) + 
  ylab(expression(-log[10]("padj"))) +
  geom_hline(yintercept = 1.3, colour = "red4", linetype = "dashed") + 
  geom_vline(xintercept = 0.0, colour = "#474657", linetype = "dashed") +
  scale_y_continuous(trans = "log1p") +
  ggtitle("DESeq2 Volcano Plot - Mast cells Surgery vs Biopsy") 

dev.off()




### Normal breast ###

countTable = read.table("Normal_pb_matrix.csv",header=TRUE, row.names=1, check.names = FALSE,sep = ',')
designMatrix = read.table('Normal_metadata.csv', header=TRUE, row.names=1,sep = ',')
libType = c("oneFactor","oneFactor","oneFactor","oneFactor")
condition = factor(designMatrix$Condition)
patient = factor(designMatrix$Patient)
condition <- relevel(condition, "Pre")
experiment_design=data.frame(row.names = colnames(countTable),condition,patient,libType)

colnames(countTable) == rownames(experiment_design)

cds <- DESeqDataSetFromMatrix(countData =countTable, colData=experiment_design, design=~patient + condition)
cds_DESeqED <- DESeq(cds,parallel = TRUE)
res <- results(cds_DESeqED,parallel = TRUE,alpha = 0.05, pAdjustMethod = "BH")
write.table(res,file = "Normal_Post_vs_Pre.differentialExpression.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
resSig <- subset(res, padj < 0.05)
write.table(resSig,file = "Normal_Post_vs_Pre.sigDEGs.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
normalizedReadCounts = counts(cds_DESeqED,normalized=TRUE)
write.table(normalizedReadCounts,file = "Normal_Post_vs_Pre.normalizedCounts.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")


## Volcano plot with labels

Normal_breast_Post_vs_Pre <- read.table("Normal_Post_vs_Pre.differentialExpression.txt", header=TRUE, sep = '\t')
Normal_breast_Post_vs_Pre <- Normal_breast_Post_vs_Pre[complete.cases(Normal_breast_Post_vs_Pre), ]

Normal_breast_Post_vs_Pre$Type <- "NONE"
Normal_breast_Post_vs_Pre$Type[Normal_breast_Post_vs_Pre$padj < 0.05 & Normal_breast_Post_vs_Pre$log2FoldChange < -1] <- "DOWN"
Normal_breast_Post_vs_Pre$Type[Normal_breast_Post_vs_Pre$padj < 0.05 & Normal_breast_Post_vs_Pre$log2FoldChange > 1] <- "UP"
cols <- c("NONE" = "#474657", "DOWN" = "#157ded", "UP" = "#f2800f")

Normal_breast_Post_vs_Pre$genelabels <- ""


# To make a column with labels according to the threshold
Normal_Post_vs_Pre$label <- ifelse(
  Normal_Post_vs_Pre$padj < 0.05,
  Normal_Post_vs_Pre$X,
  NA
)


png("DESeq2_Normal_breast_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)
ggplot(Normal_Post_vs_Pre, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 1, alpha = 0.75, na.rm = TRUE) +
  geom_text_repel(aes(label = label), size = 8, max.overlaps = 30, na.rm = TRUE) +
  theme_bw(base_size = 18) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20)) +
  xlab(expression(log[2]("Fold Change"))) + 
  ylab(expression(-log[10]("padj"))) +
  geom_hline(yintercept = 1.3, colour = "red4", linetype = "dashed") + 
  geom_vline(xintercept = 0.0, colour = "#474657", linetype = "dashed") +
  scale_y_continuous(trans = "log1p") +
  ggtitle("DESeq2 Volcano Plot - Normal breast Surgery vs Biopsy") 

dev.off()







### Perivascular ###

countTable = read.table("Perivascular_pb_matrix.csv",header=TRUE, row.names=1, check.names = FALSE,sep = ',')
designMatrix = read.table('Perivascular_metadata.csv', header=TRUE, row.names=1,sep = ',')
libType = c("oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor")
condition = factor(designMatrix$Condition)
patient = factor(designMatrix$Patient)
condition <- relevel(condition, "Pre")
experiment_design=data.frame(row.names = colnames(countTable),condition,patient,libType)

colnames(countTable) == rownames(experiment_design)

cds <- DESeqDataSetFromMatrix(countData =countTable, colData=experiment_design, design=~patient + condition)
cds_DESeqED <- DESeq(cds,parallel = TRUE)
res <- results(cds_DESeqED,parallel = TRUE,alpha = 0.05, pAdjustMethod = "BH")
write.table(res,file = "Perivascular_Post_vs_Pre.differentialExpression.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
resSig <- subset(res, padj < 0.05)
write.table(resSig,file = "Perivascular_Post_vs_Pre.sigDEGs.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
normalizedReadCounts = counts(cds_DESeqED,normalized=TRUE)
write.table(normalizedReadCounts,file = "Perivascular_Post_vs_Pre.normalizedCounts.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")


## Volcano plot with labels

Perivascular_Post_vs_Pre <- read.table("Perivascular_Post_vs_Pre.differentialExpression.txt", header=TRUE, sep = '\t')
Perivascular_Post_vs_Pre <- Perivascular_Post_vs_Pre[complete.cases(Perivascular_Post_vs_Pre), ]

Perivascular_Post_vs_Pre$Type <- "NONE"
Perivascular_Post_vs_Pre$Type[Perivascular_Post_vs_Pre$padj < 0.05 & Perivascular_Post_vs_Pre$log2FoldChange < -1] <- "DOWN"
Perivascular_Post_vs_Pre$Type[Perivascular_Post_vs_Pre$padj < 0.05 & Perivascular_Post_vs_Pre$log2FoldChange > 1] <- "UP"
cols <- c("NONE" = "#474657", "DOWN" = "#157ded", "UP" = "#f2800f")

Perivascular_Post_vs_Pre$genelabels <- ""

# To make a column with labels according to the threshold
Perivascular_Post_vs_Pre$label <- ifelse(
  Perivascular_Post_vs_Pre$padj < 0.05,
  Perivascular_Post_vs_Pre$X,
  NA
)


png("DESeq2_Perivascular_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)
ggplot(Perivascular_Post_vs_Pre, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 1, alpha = 0.75, na.rm = TRUE) +
  geom_text_repel(aes(label = label), size = 8, max.overlaps = 30, na.rm = TRUE) +
  theme_bw(base_size = 18) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20)) +
  xlab(expression(log[2]("Fold Change"))) + 
  ylab(expression(-log[10]("padj"))) +
  geom_hline(yintercept = 1.3, colour = "red4", linetype = "dashed") + 
  geom_vline(xintercept = 0.0, colour = "#474657", linetype = "dashed") +
  scale_y_continuous(trans = "log1p") +
  ggtitle("DESeq2 Volcano Plot - Perivascular Surgery vs Biopsy") 

dev.off()

