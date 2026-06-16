setwd("PATH")

library(SeuratExtend)
library(dplyr)
library(fgsea)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(grid)
library(ComplexHeatmap)
library(circlize)


load("Patient_38_RCTD.rda")



#### Markers by subclone ####

Patient_38_RCTD <- Patient_38_RCTD[, !is.na(Patient_38_RCTD$tumor_subclones_ptx38)]



# Tumor subclone 1 #


Patient_38_RCTD$tumor_subclones_ptx38[
  Patient_38_RCTD$tumor_subclones_ptx38 %in% c("2", "3")
] <- "2_3"

Idents(Patient_38_RCTD) <- Patient_38_RCTD$tumor_subclones_ptx38

Tumor_subclone_1_marker <- FindMarkers(Patient_38_RCTD, ident.1 = "1", ident.2 = "2_3")
write.table(Tumor_subclone_1_marker, file = "Patient_38_Tumor_subclone_1_marker.txt", sep = '\t')



rm(Patient_38_RCTD)
load("Patient_38_RCTD.rda")
Patient_38_RCTD <- Patient_38_RCTD[, !is.na(Patient_38_RCTD$tumor_subclones_ptx38)]



# Tumor subclone 2 #


Patient_38_RCTD$tumor_subclones_ptx38[
  Patient_38_RCTD$tumor_subclones_ptx38 %in% c("1", "3")
] <- "1_3"

Idents(Patient_38_RCTD) <- Patient_38_RCTD$tumor_subclones_ptx38

Tumor_subclone_2_marker <- FindMarkers(Patient_38_RCTD, ident.1 = "2", ident.2 = "1_3")
write.table(Tumor_subclone_2_marker, file = "Patient_38_Tumor_subclone_2_marker.txt", sep = '\t')


rm(Patient_38_RCTD)
load("Patient_38_RCTD.rda")
Patient_38_RCTD <- Patient_38_RCTD[, !is.na(Patient_38_RCTD$tumor_subclones_ptx38)]



# Tumor subclone 3 #


Patient_38_RCTD$tumor_subclones_ptx38[
  Patient_38_RCTD$tumor_subclones_ptx38 %in% c("1", "2")
] <- "1_2"

Idents(Patient_38_RCTD) <- Patient_38_RCTD$tumor_subclones_ptx38

Tumor_subclone_3_marker <- FindMarkers(Patient_38_RCTD, ident.1 = "3", ident.2 = "1_2")
write.table(Tumor_subclone_3_marker, file = "Patient_38_Tumor_subclone_3_marker.txt", sep = '\t')




#### fgsea for each tumor subclone ####


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


# Tumor subclone 1
Tumor_subclone_1_marker <- read.table("Patient_38_Tumor_subclone_1_marker.txt", header = TRUE, sep = "\t")

ranks <- Tumor_subclone_1_marker$avg_log2FC
names(ranks) <- rownames(Tumor_subclone_1_marker)


Tumor_1_38 <-fgseaMultilevel(
  pathways = genesets,
  ranks,
  sampleSize = 101,
  minSize = 1,
  maxSize = length(ranks)- 1,
  eps = 0.0,
  scoreType = "std",
  nproc = 0,
  fgseaParam = 1,
  BPPARAM = NULL,
  nPermSimple = 1000,
  absEps = NULL
)

Tumor_1_38$ident <- "Tumor_subclone_1"


Tumor_1_38$leadingEdge <- sapply(Tumor_1_38$leadingEdge, function(x) paste(x, collapse = ", "))
write.table(Tumor_1_38, file = "Patient_38_Tumor_subclone_1_fgsea.txt", sep = '\t')



# Tumor subclone 2
Tumor_subclone_2_marker <- read.table("Patient_38_Tumor_subclone_2_marker.txt", header = TRUE, sep = "\t")

ranks <- Tumor_subclone_2_marker$avg_log2FC
names(ranks) <- rownames(Tumor_subclone_2_marker)


Tumor_2_38 <-fgseaMultilevel(
  pathways = genesets,
  ranks,
  sampleSize = 101,
  minSize = 1,
  maxSize = length(ranks)- 1,
  eps = 0.0,
  scoreType = "std",
  nproc = 0,
  fgseaParam = 1,
  BPPARAM = NULL,
  nPermSimple = 1000,
  absEps = NULL
)

Tumor_2_38$ident <- "Tumor_subclone_2"

Tumor_2_38$leadingEdge <- sapply(Tumor_2_38$leadingEdge, function(x) paste(x, collapse = ", "))
write.table(Tumor_2_38, file = "Patient_38_Tumor_subclone_2_fgsea.txt", sep = '\t')




# Tumor subclone 3
Tumor_subclone_3_marker <- read.table("Patient_38_Tumor_subclone_3_marker.txt", header = TRUE, sep = "\t")

ranks <- Tumor_subclone_3_marker$avg_log2FC
names(ranks) <- rownames(Tumor_subclone_3_marker)


Tumor_3_38 <-fgseaMultilevel(
  pathways = genesets,
  ranks,
  sampleSize = 101,
  minSize = 1,
  maxSize = length(ranks)- 1,
  eps = 0.0,
  scoreType = "std",
  nproc = 0,
  fgseaParam = 1,
  BPPARAM = NULL,
  nPermSimple = 1000,
  absEps = NULL
)

Tumor_3_38$ident <- "Tumor_subclone_3"

Tumor_3_38$leadingEdge <- sapply(Tumor_3_38$leadingEdge, function(x) paste(x, collapse = ", "))
write.table(Tumor_3_38, file = "Patient_38_Tumor_subclone_3_fgsea.txt", sep = '\t')


Patient_38_fgsea <- bind_rows(Tumor_1_38, Tumor_2_38, Tumor_3_38)
Patient_38_fgsea$patient <- "38"

write.table(Patient_38_fgsea, file = "Patient_38_tumor_subclones_fgsea.txt", sep = '\t')







