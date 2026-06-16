setwd("PATH")

library(Seurat)
library(DESeq2)
library(stats)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(SeuratExtend)
library(fgsea)
library(tidyr)
library(grid)
library(circlize)


###### Patient_38 ###### 

load("Patient_38_RCTD.rda")

### To get differentially expressed genes for each cell type ### 

Patient_38_RCTD$final_annotation_Denosumab <-paste(Patient_38_RCTD$final_annotation, Patient_38_RCTD$Denosumab, sep = "_")

Idents(Patient_38_RCTD) <- Patient_38_RCTD$final_annotation_Denosumab

Myeloid_Denosumab_marker <- FindMarkers(Patient_38_RCTD, ident.1 = "Myeloid_Post", ident.2 = "Myeloid_Pre")
write.table(Myeloid_Denosumab_marker, file = "Patient_38_Myeloid_Denosumab_DEGs.txt", sep = '\t')

B_Denosumab_marker_38 <- FindMarkers(Patient_38_RCTD, ident.1 = "B_cells_Post", ident.2 = "B_cells_Pre")
write.table(B_Denosumab_marker_38, file = "Patient_38_B_Denosumab_DEGs.txt", sep = '\t')

T_Denosumab_marker_38 <- FindMarkers(Patient_38_RCTD, ident.1 = "T_cells_Post", ident.2 = "T_cells_Pre")
write.table(T_Denosumab_marker_38, file = "Patient_38_T_Denosumab_DEGs.txt", sep = '\t')

Fibroblast_Denosumab_marker_38 <- FindMarkers(Patient_38_RCTD, ident.1 = "Fibroblast_Post", ident.2 = "Fibroblast_Pre")
write.table(Fibroblast_Denosumab_marker_38, file = "Patient_38_Fibroblast_Denosumab_DEGs.txt", sep = '\t')

Lymphatic_Denosumab_marker_38 <- FindMarkers(Patient_38_RCTD, ident.1 = "Lymphatic_Post", ident.2 = "Lymphatic_Pre")
write.table(Lymphatic_Denosumab_marker_38, file = "Patient_38_Lymphatic_Denosumab_DEGs.txt", sep = '\t')

Adipocytes_Denosumab_marker_38 <- FindMarkers(Patient_38_RCTD, ident.1 = "Adipocytes_Post", ident.2 = "Adipocytes_Pre")
write.table(Adipocytes_Denosumab_marker_38, file = "Patient_38_Adipocytes_Denosumab_DEGs.txt", sep = '\t')

Endothelial_Denosumab_marker_38 <- FindMarkers(Patient_38_RCTD, ident.1 = "Endothelial_cells_Post", ident.2 = "Endothelial_cells_Pre")
write.table(Endothelial_Denosumab_marker_38, file = "Patient_38_Endothelial_Denosumab_DEGs.txt", sep = '\t')

plasmablasts_Denosumab_marker_38 <- FindMarkers(Patient_38_RCTD, ident.1 = "Plasmablasts_Post", ident.2 = "Plasmablasts_Pre")
write.table(plasmablasts_Denosumab_marker_38, file = "Patient_38_plasmablasts_Denosumab_DEGs.txt", sep = '\t')

ductal_invasion_Denosumab_marker_38 <- FindMarkers(Patient_38_RCTD, ident.1 = "ductal_invasion_Post", ident.2 = "ductal_invasion_Pre")
write.table(Endothelial_Denosumab_marker_38, file = "Patient_38_ductal_invasion_Denosumab_DEGs.txt", sep = '\t')

Mast_Denosumab_marker_38 <- FindMarkers(Patient_38_RCTD, ident.1 = "Mast_Post", ident.2 = "Mast_Pre")
write.table(Mast_Denosumab_marker_38, file = "Patient_38_Mast_Denosumab_DEGs.txt", sep = '\t')

normal_breast_tissue_Denosumab_marker_38 <- FindMarkers(Patient_38_RCTD, ident.1 = "normal_breast_tissue_Post", ident.2 = "normal_breast_tissue_Pre")
write.table(normal_breast_tissue_Denosumab_marker_38, file = "Patient_38_normal_breast_tissue_Denosumab_DEGs.txt", sep = '\t')

tumor_Denosumab_marker_38 <- FindMarkers(Patient_38_RCTD, ident.1 = "tumor_Post", ident.2 = "tumor_Pre")
write.table(tumor_Denosumab_marker_38, file = "Patient_38_tumor_Denosumab_DEGs.txt", sep = '\t')




### To perform fgsea using preranked DEGs ### 

# We load the gmt with Hallmarks, Reactome and KEGG signatures

sig_file <- "h.all.v2026.1.Hs.symbols_plus_reactome_kegg.txt"

sig <- read.delim(
  sig_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  fill = TRUE,
  quote = "",
  comment.char = ""
)

genesets <- setNames(
  lapply(seq_len(nrow(sig)), function(i) {
    genes <- unlist(sig[i, -1], use.names = FALSE)
    genes <- genes[!is.na(genes) & genes != ""]
    unique(as.character(genes))
  }),
  sig[[1]]
)

# Myeloid
Myeloid_Denosumab_marker <- read.table("Patient_38_Myeloid_Denosumab_DEGs.txt", header = TRUE, sep = "\t")

ranks <- Myeloid_Denosumab_marker$avg_log2FC
names(ranks) <- rownames(Myeloid_Denosumab_marker)


Myeloid_38 <-fgseaMultilevel(
  pathways = genesets,
  ranks,
  sampleSize = 101,
  minSize = 1,
  maxSize = length(ranks)- 1,
  eps = 0.0,
  scoreType = "std",
  nproc = 0,
  gseaParam = 1,
  BPPARAM = NULL,
  nPermSimple = 1000,
  absEps = NULL
)

Myeloid_38$ident <- "Myeloid"

Myeloid_38$leadingEdge <- sapply(Myeloid_38$leadingEdge, function(x) paste(x, collapse = ", "))
write.table(Myeloid_38, file = "Patient_38_Myeloid_Denosumab_GSEA.txt", sep = '\t')


# B_cells
B_cells_Denosumab_marker <- read.table("Patient_38_B_Denosumab_DEGs.txt", header=TRUE, sep = '\t')

ranks_B <- B_cells_Denosumab_marker$avg_log2FC
names(ranks_B) <- rownames(B_cells_Denosumab_marker)


B_cells_38 <-fgseaMultilevel(
  pathways = genesets,
  ranks_B,
  sampleSize = 101,
  minSize = 1,
  maxSize = length(ranks)- 1,
  eps = 0.0,
  scoreType = "std",
  nproc = 0,
  gseaParam = 1,
  BPPARAM = NULL,
  nPermSimple = 1000,
  absEps = NULL
)

B_cells_38$ident <- "B_cells"


B_cells_38$leadingEdge <- sapply(B_cells_38$leadingEdge, function(x) paste(x, collapse = ", "))
write.table(B_cells_38, file = "Patient_38_B_cells_Denosumab_GSEA.txt", sep = '\t')


# T_cells #

T_cells_Denosumab_marker <- read.table("Patient_38_T_Denosumab_DEGs.txt", header=TRUE, sep = '\t')

ranks_T <- T_cells_Denosumab_marker$avg_log2FC
names(ranks_T) <- rownames(T_cells_Denosumab_marker)


T_cells_38 <-fgseaMultilevel(
  pathways = genesets,
  ranks_T,
  sampleSize = 101,
  minSize = 1,
  maxSize = length(ranks)- 1,
  eps = 0.0,
  scoreType = "std",
  nproc = 0,
  gseaParam = 1,
  BPPARAM = NULL,
  nPermSimple = 1000,
  absEps = NULL
)

T_cells_38$ident <- "T_cells"

T_cells_38$leadingEdge <- sapply(T_cells_38$leadingEdge, function(x) paste(x, collapse = ", "))
write.table(T_cells_38, file = "Patient_38_T_cells_Denosumab_GSEA.txt", sep = '\t')


# Fibroblast #

Fibroblast_Denosumab_marker <- read.table("Patient_38_Fibroblast_Denosumab_DEGs.txt", header=TRUE, sep = '\t')

ranks_Fibroblast <- Fibroblast_Denosumab_marker$avg_log2FC
names(ranks_Fibroblast) <- rownames(Fibroblast_Denosumab_marker)


Fibroblast_38 <-fgseaMultilevel(
  pathways = genesets,
  ranks_Fibroblast,
  sampleSize = 101,
  minSize = 1,
  maxSize = length(ranks)- 1,
  eps = 0.0,
  scoreType = "std",
  nproc = 0,
  gseaParam = 1,
  BPPARAM = NULL,
  nPermSimple = 1000,
  absEps = NULL
)

Fibroblast_38$ident <- "Fibroblast"


Fibroblast_38$leadingEdge <- sapply(Fibroblast_38$leadingEdge, function(x) paste(x, collapse = ", "))
write.table(Fibroblast_38, file = "Patient_38_Fibroblast_Denosumab_GSEA.txt", sep = '\t')


# Endothelial_cells #

Endothelial_cells_Denosumab_marker <- read.table("Patient_38_Endothelial_Denosumab_DEGs.txt", header=TRUE, sep = '\t')

ranks_endothelial <- Endothelial_cells_Denosumab_marker$avg_log2FC
names(ranks_endothelial) <- rownames(Endothelial_cells_Denosumab_marker)


endothelial_38 <-fgseaMultilevel(
  pathways = genesets,
  ranks_endothelial,
  sampleSize = 101,
  minSize = 1,
  maxSize = length(ranks)- 1,
  eps = 0.0,
  scoreType = "std",
  nproc = 0,
  gseaParam = 1,
  BPPARAM = NULL,
  nPermSimple = 1000,
  absEps = NULL
)

endothelial_38$ident <- "Endothelial_cells"

endothelial_38$leadingEdge <- sapply(endothelial_38$leadingEdge, function(x) paste(x, collapse = ", "))
write.table(endothelial_38, file = "Patient_38_Endothelial_Denosumab_GSEA.txt", sep = '\t')



# normal_breast_tissue #

normal_breast_tissue_Denosumab_marker <- read.table("Patient_38_normal_breast_tissue_Denosumab_DEGs.txt", header=TRUE, sep = '\t')

ranks_normal <- normal_breast_tissue_Denosumab_marker$avg_log2FC
names(ranks_normal) <- rownames(normal_breast_tissue_Denosumab_marker)


normal_38 <-fgseaMultilevel(
  pathways = genesets,
  ranks_normal,
  sampleSize = 101,
  minSize = 1,
  maxSize = length(ranks)- 1,
  eps = 0.0,
  scoreType = "std",
  nproc = 0,
  gseaParam = 1,
  BPPARAM = NULL,
  nPermSimple = 1000,
  absEps = NULL
)

normal_38$ident <- "normal_breast_tissue"


normal_38$leadingEdge <- sapply(normal_38$leadingEdge, function(x) paste(x, collapse = ", "))
write.table(normal_38, file = "Patient_38_normal_breast_Denosumab_GSEA.txt", sep = '\t')



# tumor #

tumor_Denosumab_marker <- read.table("Patient_38_tumor_Denosumab_DEGs.txt", header=TRUE, sep = '\t')


ranks_tumor <- tumor_Denosumab_marker$avg_log2FC
names(ranks_tumor) <- rownames(tumor_Denosumab_marker)


tumor_38 <-fgseaMultilevel(
  pathways = genesets,
  ranks_tumor,
  sampleSize = 101,
  minSize = 1,
  maxSize = length(ranks)- 1,
  eps = 0.0,
  scoreType = "std",
  nproc = 0,
  gseaParam = 1,
  BPPARAM = NULL,
  nPermSimple = 1000,
  absEps = NULL
)

tumor_38$ident <- "tumor"



tumor_38$leadingEdge <- sapply(tumor_38$leadingEdge, function(x) paste(x, collapse = ", "))
write.table(tumor_38, file = "Patient_38_tumor_Denosumab_GSEA.txt", sep = '\t')


Patient_38_GSEA <- bind_rows(B_cells_38, endothelial_38, Myeloid_38, normal_38, Fibroblast_38, T_cells_38, tumor_38)
Patient_38_GSEA$patient <- "38"

write.table(Patient_38_GSEA, file = "Patient_38_filtered_Denosumab_GSEA.txt", sep = '\t')


