library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sctransform)
library(glmGamPoi)
library(clustree)
library(harmony)


setwd("PATH")

load("Patient_4_RCTD_snRNAseq.rda")
load("Patient_06_RCTD_snRNAseq.rda")
load("Patient_38_RCTD_snRNAseq.rda")
load("Patient_41_67_RCTD_snRNAseq.rda")
load("Patient_61_RCTD_snRNAseq.rda")

Patient_4_RCTD <- Patient_4_RCTD[, !is.na(Patient_4_RCTD$Block)]
Patient_06_RCTD <- Patient_06_RCTD[, !is.na(Patient_06_RCTD$Block)]
Patient_38_RCTD <- Patient_38_RCTD[, !is.na(Patient_38_RCTD$Block)]
Patient_41_67_RCTD <- Patient_41_67_RCTD[, !is.na(Patient_41_67_RCTD$Block)]
Patient_61_RCTD <- Patient_61_RCTD[, !is.na(Patient_61_RCTD$Block)]


Patient_4_RCTD <- Patient_4_RCTD[, !is.na(Patient_4_RCTD$final_annotation)]
Patient_06_RCTD <- Patient_06_RCTD[, !is.na(Patient_06_RCTD$final_annotation)]
Patient_38_RCTD <- Patient_38_RCTD[, !is.na(Patient_38_RCTD$final_annotation)]
Patient_41_67_RCTD <- Patient_41_67_RCTD[, !is.na(Patient_41_67_RCTD$final_annotation)]
Patient_61_RCTD <- Patient_61_RCTD[, !is.na(Patient_61_RCTD$final_annotation)]



ST_merged  <- merge(x = Patient_06_RCTD, y = c(Patient_38_RCTD, Patient_4_RCTD, Patient_61_RCTD, Patient_41_67_RCTD), add.cell.ids = c("A1","A2","B1","B2","C1"), merge.data = TRUE)

ST_merged <- JoinLayers(ST_merged,  assay = "Spatial")
ST_merged <- subset(ST_merged, subset = nCount_Spatial > 50)

ST_merged <- NormalizeData(ST_merged)
ST_merged <- FindVariableFeatures(ST_merged)
ST_merged <- ScaleData(ST_merged)
ST_merged <- RunPCA(ST_merged)

pdf('elbow_ST_merged.pdf')
ElbowPlot(ST_merged,ndims=50)
dev.off()

ST_merged <- FindNeighbors(ST_merged, dims = 1:30)
ST_merged <- FindClusters(ST_merged,resolution = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

pdf('clustree_ST_merged.pdf',width = 14, height = 14)
clustree(ST_merged)
dev.off()


## Harmony integration:

ST_integrated <- RunHarmony(ST_merged, "Block",kmeans_init_nstart=20, kmeans_init_iter_max=100)

ST_integrated <- RunUMAP(ST_integrated, reduction="harmony", dims=1:30)

ST_integrated <- FindNeighbors(ST_integrated, dims = 1:30)
ST_integrated <- FindClusters(ST_integrated,resolution = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

pdf('clustree_ST_integrated.pdf',width = 14, height = 14)
clustree(ST_integrated)
dev.off()

Idents(ST_integrated) <- ST_integrated$final_annotation

save(ST_integrated, file = "ST_ER_RCTD_integrated.rda")


# UMAPs

# Defining the same color for each annotation as SpatialDimPlots

cols <- c(
  ductal_invasion        = "#4A0000",  
  tumor                  = "#7F0000",  
  tumor_mixed            = "#F70000",
  "in situ carcinoma"    = "#2A0000",
  normal_breast_tissue   = "#1A9850", 
  normal_breast_mixed    = "#66C2A5",  
  Fibroblast             = "#BCA9FC",  
  artery                 = "#4126AC",
  Endothelial_cells      = "#2166AC",  
  Lymphatic              = "#4393C3", 
  Perivascular           = "#92C5DE",  
  B_cells                = "#F46D43",  
  T_cells                = "#FDAE61",  
  Plasmablasts           = "#1CECFC",  
  Myeloid                = "#ED91E5",  
  Mast                   = "#FCE99D",  
  Adipocytes             = "#FFFF33"   
)


p1 <- DimPlot(object = ST_integrated, reduction = "umap", raster=FALSE, group.by = "Patient") +
  ggtitle("Patient") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 14),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16)) 
p2 <- DimPlot(object = ST_integrated, reduction = "umap", raster=FALSE, group.by = "Denosumab") +
  ggtitle("Denosumab") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 14),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16)) 
p3 <- DimPlot(object = ST_integrated, reduction = "umap", raster=FALSE, cols = cols, group.by = "Pathologist_annotation") +
  ggtitle("Pathologist annotation") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 14),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16))
p4 <- DimPlot(object = ST_integrated, reduction = "umap", raster=FALSE, cols = cols, group.by = "final_annotation") +
  ggtitle("Final annotation") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 14),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16))

png("UMAPs_ST_integrated.png", width = 12000, height = 10000, res = 650)
p1+p2+p3+p4
dev.off()


### FindAllmarkers of all cell idents ###

Idents(ST_integrated) <- ST_integrated$final_annotation

ST_integrated.markers <- FindAllMarkers(ST_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(ST_integrated.markers, "ST_integrated_markers.txt", sep = '\t')

# Top 5 markers dotplot #

ST_integrated.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5_final_dominant_celltype

clusters <- sort(unique(ST_integrated.markers$cluster), decreasing = TRUE)

Idents(ST_integrated) <- factor(Idents(ST_integrated), levels = clusters)

ST_integrated <- ST_integrated[, !(ST_integrated$final_annotation == "artifact")]

png("Pseudobulk_final_annotation_top5markers_dotplot.png", width = 13500, height = 4000, res = 450)

DotPlot(ST_integrated, features = unique(top5_final_dominant_celltype$gene)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
dev.off()


# Counts matrix with the final annotation and treatment timepoint for each patient #

ST_integrated$final_annotation_Denosumab <-paste(ST_integrated$final_annotation, ST_integrated$Denosumab, sep = "_")

agg_expr <- AggregateExpression(
  object = ST_integrated,
  group.by = c("final_annotation_Denosumab", "Patient"),
  assay = "Spatial",
  slot = "counts"  
)

pb_matrix <- as.matrix(agg_expr$Spatial)
write.table(pb_matrix, file = "ST_integrated_pb_matrix_final_annotation.txt", sep = '\t')
