### Visium HD analysis from: https://satijalab.org/seurat/articles/visiumhd_analysis_vignette ###

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clustree)
library(loupeR)


local_dir <- "PATH/square_016um/"
Patient_38 <- Load10X_Spatial(data.dir = local_dir,filename = "filtered_feature_bc_matrix.h5",assay = "Spatial")

# Setting default assay changes between 8um and 16um binning
Assays(Patient_38)
DefaultAssay(Patient_38) <- "Spatial.016um"

pdf("Patient_38_nCount_16um.pdf", height = 7, width = 14)
vln.plot <- VlnPlot(Patient_38, features = "nCount_Spatial", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(Patient_38, features = "nCount_Spatial",pt.size.factor = 3.2) + theme(legend.position = "right")

vln.plot | count.plot
dev.off()

Patient_38 <- subset(Patient_38, subset = nCount_Spatial > 200)

## Normalize data

Patient_38 <- NormalizeData(Patient_38)
Patient_38 <- FindVariableFeatures(Patient_38)
Patient_38 <- ScaleData(Patient_38)
Patient_38 <- RunPCA(Patient_38)
Patient_38 <- RunUMAP(Patient_38, dims = 1:50)


Patient_38 <- FindNeighbors(Patient_38, dims = 1:50)
Patient_38 <- FindClusters(Patient_38, resolution = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

pdf("clustree_Patient_38_16um.pdf",width = 10, height = 14)
clustree(Patient_38)
dev.off()

Idents(Patient_38) <- Patient_38$Spatial_snn_res.0.7


## loupe R

create_loupe_from_seurat(Patient_38, output_name = "Patient_38_16um.cloupe")
Patient_38_metadata <- read.csv("PATH/metadata_blockA2.csv")
Patient_38_16um_bc <- as.data.frame(colnames(Patient_38))
Patient_38_16um_bc$`colnames(Patient_38)` <- gsub("016um", "008um", Patient_38_16um_bc$`colnames(Patient_38)`)
Patient_38_16um_metadata <- Patient_38_metadata[Patient_38_metadata$Barcode %in% Patient_38_16um_bc$`colnames(Patient_38)`, ]

Patient_38_16um_metadata$Barcode <- gsub("008um", "016um", Patient_38_16um_metadata$Barcode)

## Add metadata

Patient_38 <- AddMetaData(object = Patient_38, metadata = Patient_38_16um_metadata$Block, col.name = "Block")
Patient_38 <- AddMetaData(object = Patient_38, metadata = Patient_38_16um_metadata$Patient, col.name = "Patient")
Patient_38 <- AddMetaData(object = Patient_38, metadata = Patient_38_16um_metadata$Denosumab, col.name = "Denosumab")

save(Patient_38, file = "Patient_38_16um.rda")

