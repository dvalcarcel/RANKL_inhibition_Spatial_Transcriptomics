setwd("PATH")

library(SCEVAN)
library(Seurat)
library(Matrix)


load("Patient_38_RCTD.rda")

#### Inferring clonality from RNA-Seq ####

counts_matrix <- Patient_38_RCTD@assays$Spatial@layers$counts

rownames(counts_matrix) <- rownames(Patient_38_RCTD)   # genes
colnames(counts_matrix) <- colnames(Patient_38_RCTD)   # spots

keep_groups <- c("tumor", "T_cells", "B_cells", "Myeloid", "Endothelial_cells", "Fibroblast")

annotations <- data.frame(
  spot_id = colnames(Patient_38_RCTD),
  cluster = Idents(Patient_38_RCTD)
)

# Filtering the cell identities required for the analysis
spots_to_keep <- annotations$spot_id[annotations$cluster %in% keep_groups]

# Subset the matrix
counts_subset <- counts_matrix[, spots_to_keep]


# To define vector of normal cells
norm_cell <- colnames(Patient_38_RCTD)[
  Patient_38_RCTD$final_annotation %in%
    c( "T_cells", "B_cells", "Myeloid", "Endothelial_cells", "Fibroblast")]



results <- pipelineCNA(
  counts_subset,
  sample = "Patient_38",
  par_cores = 6,
  norm_cell = norm_cell,
  SUBCLONES = TRUE,
  beta_vega = 0.5,
  ClonalCN = TRUE,
  plotTree = TRUE,
  AdditionalGeneSets = NULL,
  SCEVANsignatures = TRUE,
  organism = "human",
  ngenes_chr = 5,
  perc_genes = 10,
  FIXED_NORMAL_CELLS = TRUE
)

save(results, file = "SCEVAN_38.rda")


# To add the metadata into the Seurat object
tumor_subclones_ptx38 <- as.data.frame(results$subclone)
rownames(tumor_subclones_ptx38) <- rownames(results)
colnames(tumor_subclones_ptx38) <- "tumor_subclones_ptx38"
Patient_38_RCTD <- AddMetaData(Patient_38_RCTD, tumor_subclones_ptx38)

save(Patient_38_RCTD, file = "Patient_38_RCTD.rda")




## SpatialDimPlot to spatially visualize the tumour subclones ##


cols <- c(
  "1"   = "#A6CEE3",  
  "2"   = "#1F78B4",  
  "3"   = "#B2DF8A",
  "NA"  = "#7F7F7F"  
)



png("Patient_38_tumor_subclones.png", width = 28, height = 10, units = "in", res = 1000)
SpatialDimPlot(Patient_38_RCTD, group.by ="tumor_subclones_ptx38", cols = cols, image.alpha=0, pt.size.factor=5.5) +
  theme(
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 14)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 6),  
      ncol = 2
    ))
dev.off()



## Normalised Proportion barplot between biopsy and surgery ##


Patient_38_RCTD$final_annotation_Denosumab <-paste(Patient_38_RCTD$final_annotation, Patient_38_RCTD$Denosumab, sep = "_")

df_tumor <- Patient_38_RCTD@meta.data %>%
  filter(grepl("^tumor_P", final_annotation_Denosumab))

df_tumor <- df_tumor %>%
  mutate(
    Denosumab = case_when(
      grepl("_Pre$", final_annotation_Denosumab)  ~ "Pre",
      grepl("_Post$", final_annotation_Denosumab) ~ "Post",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Denosumab))


df_prop <- df_tumor %>%
  count(Denosumab, tumor_subclones_ptx38) %>%
  group_by(Denosumab) %>%
  mutate(
    proportion = n / sum(n)
  )

df_prop$Denosumab <- factor(
  df_prop$Denosumab,
  levels = c("Pre", "Post")
)

png("Patient_38_tumor_subclones_proportion.png", width = 15, height = 10, units = "in", res = 600)
ggplot(df_prop,
       aes(x = tumor_subclones_ptx38, y = proportion, fill = Denosumab)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_fill_manual(values = c("Pre" = "#1f77b4", "Post" = "#ff7f0e")) +
  theme_minimal() +
  ylab("Proportion (normalized by number of tumor spots)") +
  xlab("Tumor subclones") +
  theme(
    legend.title = element_text(size = 22, face = "bold"),
    legend.text  = element_text(size = 22),
    axis.text.x  = element_text(size = 22),  
    axis.text.y  = element_text(size = 22),  
    axis.title.x = element_text(size = 22),  
    axis.title.y = element_text(size = 22)
  )
dev.off()
