library(Seurat)
library(harmony)
library(hdf5r)
library(clustree)
library(dplyr)
library(ggplot2)

setwd("PATH/GSE234817/")


#Loading datasets
snRNAseq_1.data <- Read10X_h5(filename = "GSM7475145_hbca_n01_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_2.data <- Read10X_h5(filename = "GSM7475146_hbca_n02_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_3.data <- Read10X_h5(filename = "GSM7475147_hbca_n03_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_4.data <- Read10X_h5(filename = "GSM7475148_hbca_n04_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_5.data <- Read10X_h5(filename = "GSM7475149_hbca_n05_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_6.data <- Read10X_h5(filename = "GSM7475150_hbca_n06_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_7.data <- Read10X_h5(filename = "GSM7475151_hbca_n07_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_8.data <- Read10X_h5(filename = "GSM7475152_hbca_n08_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_9.data <- Read10X_h5(filename = "GSM7475153_hbca_n09_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_10.data <- Read10X_h5(filename = "GSM7475154_hbca_n10_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_11.data <- Read10X_h5(filename = "GSM7475155_hbca_n11_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_12.data <- Read10X_h5(filename = "GSM7475156_hbca_n12_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_13.data <- Read10X_h5(filename = "GSM7475157_hbca_n13_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_14.data <- Read10X_h5(filename = "GSM7475158_hbca_n14_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_15.data <- Read10X_h5(filename = "GSM7475159_hbca_n15_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_16.data <- Read10X_h5(filename = "GSM7475160_hbca_n16_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_17.data <- Read10X_h5(filename = "GSM7475161_hbca_n17_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_18.data <- Read10X_h5(filename = "GSM7475162_hbca_n18_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_19.data <- Read10X_h5(filename = "GSM7475163_hbca_n19_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_20.data <- Read10X_h5(filename = "GSM7475164_hbca_n20_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_21.data <- Read10X_h5(filename = "GSM7475165_hbca_n21_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_22.data <- Read10X_h5(filename = "GSM7475166_hbca_n22_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_23.data <- Read10X_h5(filename = "GSM7475167_hbca_n23_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)
snRNAseq_24.data <- Read10X_h5(filename = "GSM7475168_hbca_n24_filtered_feature_bc_matrix.h5", 
                         use.names = TRUE, unique.features = TRUE)

#Create Seurat objects

snRNAseq_1.data <- CreateSeuratObject(counts = snRNAseq_1.data, project = "snRNAseq_1")
snRNAseq_2.data <- CreateSeuratObject(counts = snRNAseq_2.data, project = "snRNAseq_2")
snRNAseq_3.data <- CreateSeuratObject(counts = snRNAseq_3.data, project = "snRNAseq_3")
snRNAseq_4.data <- CreateSeuratObject(counts = snRNAseq_4.data, project = "snRNAseq_4")
snRNAseq_5.data <- CreateSeuratObject(counts = snRNAseq_5.data, project = "snRNAseq_5")
snRNAseq_6.data <- CreateSeuratObject(counts = snRNAseq_6.data, project = "snRNAseq_6")
snRNAseq_7.data <- CreateSeuratObject(counts = snRNAseq_7.data, project = "snRNAseq_7")
snRNAseq_8.data <- CreateSeuratObject(counts = snRNAseq_8.data, project = "snRNAseq_8")
snRNAseq_9.data <- CreateSeuratObject(counts = snRNAseq_9.data, project = "snRNAseq_9")
snRNAseq_10.data <- CreateSeuratObject(counts = snRNAseq_10.data, project = "snRNAseq_10")
snRNAseq_11.data <- CreateSeuratObject(counts = snRNAseq_11.data, project = "snRNAseq_11")
snRNAseq_12.data <- CreateSeuratObject(counts = snRNAseq_12.data, project = "snRNAseq_12")
snRNAseq_13.data <- CreateSeuratObject(counts = snRNAseq_13.data, project = "snRNAseq_13")
snRNAseq_14.data <- CreateSeuratObject(counts = snRNAseq_14.data, project = "snRNAseq_14")
snRNAseq_15.data <- CreateSeuratObject(counts = snRNAseq_15.data, project = "snRNAseq_15")
snRNAseq_16.data <- CreateSeuratObject(counts = snRNAseq_16.data, project = "snRNAseq_16")
snRNAseq_17.data <- CreateSeuratObject(counts = snRNAseq_17.data, project = "snRNAseq_17")
snRNAseq_18.data <- CreateSeuratObject(counts = snRNAseq_18.data, project = "snRNAseq_18")
snRNAseq_19.data <- CreateSeuratObject(counts = snRNAseq_19.data, project = "snRNAseq_19")
snRNAseq_20.data <- CreateSeuratObject(counts = snRNAseq_20.data, project = "snRNAseq_20")
snRNAseq_21.data <- CreateSeuratObject(counts = snRNAseq_21.data, project = "snRNAseq_21")
snRNAseq_22.data <- CreateSeuratObject(counts = snRNAseq_22.data, project = "snRNAseq_22")
snRNAseq_23.data <- CreateSeuratObject(counts = snRNAseq_23.data, project = "snRNAseq_23")
snRNAseq_24.data <- CreateSeuratObject(counts = snRNAseq_24.data, project = "snRNAseq_24")

## Cell QC Analysis ##

# Determine % of mitochondrial genes 

snRNAseq_1.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_1.data, pattern = "^MT-")
snRNAseq_2.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_2.data, pattern = "^MT-")
snRNAseq_3.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_3.data, pattern = "^MT-")
snRNAseq_4.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_4.data, pattern = "^MT-")
snRNAseq_5.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_5.data, pattern = "^MT-")
snRNAseq_6.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_6.data, pattern = "^MT-")
snRNAseq_7.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_7.data, pattern = "^MT-")
snRNAseq_8.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_8.data, pattern = "^MT-")
snRNAseq_9.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_9.data, pattern = "^MT-")
snRNAseq_10.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_10.data, pattern = "^MT-")
snRNAseq_11.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_11.data, pattern = "^MT-")
snRNAseq_12.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_12.data, pattern = "^MT-")
snRNAseq_13.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_13.data, pattern = "^MT-")
snRNAseq_14.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_14.data, pattern = "^MT-")
snRNAseq_15.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_15.data, pattern = "^MT-")
snRNAseq_16.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_16.data, pattern = "^MT-")
snRNAseq_17.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_17.data, pattern = "^MT-")
snRNAseq_18.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_18.data, pattern = "^MT-")
snRNAseq_19.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_19.data, pattern = "^MT-")
snRNAseq_20.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_20.data, pattern = "^MT-")
snRNAseq_21.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_21.data, pattern = "^MT-")
snRNAseq_22.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_22.data, pattern = "^MT-")
snRNAseq_23.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_23.data, pattern = "^MT-")
snRNAseq_24.data[["percent.mt"]] <- PercentageFeatureSet(snRNAseq_24.data, pattern = "^MT-")

# Load housekeeping genes
HK_genes <- read.table("HK_genes_human.txt") 
HK_genes <- as.vector(HK_genes$V1)

# Identify the housekeeping genes that are in the dataset
HK_genes_snRNAseq_1.data <- HK_genes[HK_genes %in% rownames(snRNAseq_1.data)]
HK_genes_snRNAseq_2.data <- HK_genes[HK_genes %in% rownames(snRNAseq_2.data)]
HK_genes_snRNAseq_3.data <- HK_genes[HK_genes %in% rownames(snRNAseq_3.data)]
HK_genes_snRNAseq_4.data <- HK_genes[HK_genes %in% rownames(snRNAseq_4.data)]
HK_genes_snRNAseq_5.data <- HK_genes[HK_genes %in% rownames(snRNAseq_5.data)]
HK_genes_snRNAseq_6.data <- HK_genes[HK_genes %in% rownames(snRNAseq_6.data)]
HK_genes_snRNAseq_7.data <- HK_genes[HK_genes %in% rownames(snRNAseq_7.data)]
HK_genes_snRNAseq_8.data <- HK_genes[HK_genes %in% rownames(snRNAseq_8.data)]
HK_genes_snRNAseq_9.data <- HK_genes[HK_genes %in% rownames(snRNAseq_9.data)]
HK_genes_snRNAseq_10.data <- HK_genes[HK_genes %in% rownames(snRNAseq_10.data)]
HK_genes_snRNAseq_11.data <- HK_genes[HK_genes %in% rownames(snRNAseq_11.data)]
HK_genes_snRNAseq_12.data <- HK_genes[HK_genes %in% rownames(snRNAseq_12.data)]
HK_genes_snRNAseq_13.data <- HK_genes[HK_genes %in% rownames(snRNAseq_13.data)]
HK_genes_snRNAseq_14.data <- HK_genes[HK_genes %in% rownames(snRNAseq_14.data)]
HK_genes_snRNAseq_15.data <- HK_genes[HK_genes %in% rownames(snRNAseq_15.data)]
HK_genes_snRNAseq_16.data <- HK_genes[HK_genes %in% rownames(snRNAseq_16.data)]
HK_genes_snRNAseq_17.data <- HK_genes[HK_genes %in% rownames(snRNAseq_17.data)]
HK_genes_snRNAseq_18.data <- HK_genes[HK_genes %in% rownames(snRNAseq_18.data)]
HK_genes_snRNAseq_19.data <- HK_genes[HK_genes %in% rownames(snRNAseq_19.data)]
HK_genes_snRNAseq_20.data <- HK_genes[HK_genes %in% rownames(snRNAseq_20.data)]
HK_genes_snRNAseq_21.data <- HK_genes[HK_genes %in% rownames(snRNAseq_21.data)]
HK_genes_snRNAseq_22.data <- HK_genes[HK_genes %in% rownames(snRNAseq_22.data)]
HK_genes_snRNAseq_23.data <- HK_genes[HK_genes %in% rownames(snRNAseq_23.data)]
HK_genes_snRNAseq_24.data <- HK_genes[HK_genes %in% rownames(snRNAseq_24.data)]

# Count the number of housekeeping genes expressed per cell
snRNAseq_1.data$HK_genes <- colSums(GetAssayData(snRNAseq_1.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_1.data, ] > 0)
snRNAseq_2.data$HK_genes <- colSums(GetAssayData(snRNAseq_2.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_2.data, ] > 0)
snRNAseq_3.data$HK_genes <- colSums(GetAssayData(snRNAseq_3.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_3.data, ] > 0)
snRNAseq_4.data$HK_genes <- colSums(GetAssayData(snRNAseq_4.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_4.data, ] > 0)
snRNAseq_5.data$HK_genes <- colSums(GetAssayData(snRNAseq_5.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_5.data, ] > 0)
snRNAseq_6.data$HK_genes <- colSums(GetAssayData(snRNAseq_6.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_6.data, ] > 0)
snRNAseq_7.data$HK_genes <- colSums(GetAssayData(snRNAseq_7.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_7.data, ] > 0)
snRNAseq_8.data$HK_genes <- colSums(GetAssayData(snRNAseq_8.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_8.data, ] > 0)
snRNAseq_9.data$HK_genes <- colSums(GetAssayData(snRNAseq_9.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_9.data, ] > 0)
snRNAseq_10.data$HK_genes <- colSums(GetAssayData(snRNAseq_10.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_10.data, ] > 0)
snRNAseq_11.data$HK_genes <- colSums(GetAssayData(snRNAseq_11.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_11.data, ] > 0)
snRNAseq_12.data$HK_genes <- colSums(GetAssayData(snRNAseq_12.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_12.data, ] > 0)
snRNAseq_13.data$HK_genes <- colSums(GetAssayData(snRNAseq_13.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_13.data, ] > 0)
snRNAseq_14.data$HK_genes <- colSums(GetAssayData(snRNAseq_14.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_14.data, ] > 0)
snRNAseq_15.data$HK_genes <- colSums(GetAssayData(snRNAseq_15.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_15.data, ] > 0)
snRNAseq_16.data$HK_genes <- colSums(GetAssayData(snRNAseq_16.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_16.data, ] > 0)
snRNAseq_17.data$HK_genes <- colSums(GetAssayData(snRNAseq_17.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_17.data, ] > 0)
snRNAseq_18.data$HK_genes <- colSums(GetAssayData(snRNAseq_18.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_18.data, ] > 0)
snRNAseq_19.data$HK_genes <- colSums(GetAssayData(snRNAseq_19.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_19.data, ] > 0)
snRNAseq_20.data$HK_genes <- colSums(GetAssayData(snRNAseq_20.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_20.data, ] > 0)
snRNAseq_21.data$HK_genes <- colSums(GetAssayData(snRNAseq_21.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_21.data, ] > 0)
snRNAseq_22.data$HK_genes <- colSums(GetAssayData(snRNAseq_22.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_22.data, ] > 0)
snRNAseq_23.data$HK_genes <- colSums(GetAssayData(snRNAseq_23.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_23.data, ] > 0)
snRNAseq_24.data$HK_genes <- colSums(GetAssayData(snRNAseq_24.data, assay = "RNA", slot = "counts")[HK_genes_snRNAseq_24.data, ] > 0)




## Visualization of QC criteria
VlnPlot(snRNAseq_1.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_2.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_3.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_4.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_5.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_6.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_7.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_8.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_9.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_10.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_11.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_12.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_13.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_14.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_15.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_16.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_17.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_18.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_19.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_20.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_21.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_22.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_23.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)
VlnPlot(snRNAseq_24.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "HK_genes"), ncol = 4)


# Merging the samples
snRNAseq_merged <- merge(x = snRNAseq_1.data, y = c(snRNAseq_2.data,
                                                    snRNAseq_3.data,
                                                    snRNAseq_4.data,
                                                    snRNAseq_5.data,
                                                    snRNAseq_6.data,
                                                    snRNAseq_7.data,
                                                    snRNAseq_8.data,
                                                    snRNAseq_9.data,
                                                    snRNAseq_10.data,
                                                    snRNAseq_11.data,
                                                    snRNAseq_12.data,
                                                    snRNAseq_13.data,
                                                    snRNAseq_14.data,
                                                    snRNAseq_15.data,
                                                    snRNAseq_16.data,
                                                    snRNAseq_17.data,
                                                    snRNAseq_18.data,
                                                    snRNAseq_19.data,
                                                    snRNAseq_20.data,
                                                    snRNAseq_21.data,
                                                    snRNAseq_22.data,
                                                    snRNAseq_23.data,
                                                    snRNAseq_24.data),
                         add.cell.ids = c("1","2", "3", "4", "5","6","7","8", "9", "10",
                                          "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                          "20", "21", "22", "23", "24"))

snRNAseq_merged <- JoinLayers(snRNAseq_merged)


#Adding percent.ribo
snRNAseq_merged[["percent.ribo"]]<- PercentageFeatureSet(snRNAseq_merged, pattern = "^RP[SL][[:digit:]]")


# Filtering genes, UMIS, MT 
snRNAseq_merged_filtered <- subset(snRNAseq_merged, subset = nCount_RNA > 500 & nCount_RNA < 20000 & nFeature_RNA > 150 & nFeature_RNA < 5000 & percent.mt < 10 & percent.ribo < 50)

# Repeat preprocessing

snRNAseq_merged_filtered <- NormalizeData(snRNAseq_merged_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
snRNAseq_merged_filtered <- FindVariableFeatures(snRNAseq_merged_filtered, selection.method = "vst", nfeatures = 2000)
snRNAseq_merged_filtered <- ScaleData(snRNAseq_merged_filtered)
snRNAseq_merged_filtered <- RunPCA(snRNAseq_merged_filtered, features = VariableFeatures(object = snRNAseq_merged_filtered))

#Saving merged filtered rda
save(snRNAseq_merged_filtered, file = "snRNAseq_merged_filtered.rda")

load("snRNAseq_merged_filtered.rda")

# Visualize datasets
DimPlot(snRNAseq_merged_filtered, reduction = "pca", group.by = "orig.ident")

#UMAP in snRNAseq_merged
snRNAseq_merged_filtered <- RunUMAP(snRNAseq_merged_filtered, dims = 1:20)

DimPlot(snRNAseq_merged_filtered, reduction = "umap",group.by = "orig.ident")

# Correct batch effect by integration with harmony

snRNAseq_merged_integrated <- RunHarmony(snRNAseq_merged_filtered, "orig.ident")
snRNAseq_merged_integrated <- snRNAseq_merged_integrated %>% RunUMAP(reduction = "harmony",  dims = 1:20)


# Visualize datasets
pdf("UMAP_integrated.pdf", width = 12, height = 7)
DimPlot(object = snRNAseq_merged_integrated, reduction = "umap",group.by = "orig.ident")
dev.off()


#Saving merged filtered rda
save(snRNAseq_merged_integrated, file = "snRNAseq_merged_integrated.rda")



## VISION

library(VISION)

snRNAseq_projection <- snRNAseq_merged_integrated@reductions$umap@cell.embeddings
snRNAseq_vis <- Vision(snRNAseq_merged_integrated@assays$RNA$data,signatures=c("./snRNAseq_top_markers_Kumar_et_al.gmt"), meta = snRNAseq_merged_integrated@meta.data, pool=FALSE)
snRNAseq_vis <- addProjection(snRNAseq_vis, "UMAP", snRNAseq_projection)
snRNAseq_vis <- analyze(snRNAseq_vis)

save("snRNAseq_vis", file = "snRNAseq_vis.rda")

signature_scores <- as.data.frame(snRNAseq_vis@SigScores)

snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$Basal, col.name = "Basal")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$LumHR, col.name = "LumHR")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$LumSec, col.name = "LumSec")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$Fibroblasts, col.name = "Fibroblasts")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$Lymphatic, col.name = "Lymphatic")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$Mast, col.name = "Mast")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$Vascular, col.name = "Vascular")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$Perivascular, col.name = "Perivascular")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$Adipocytes, col.name = "Adipocytes")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$T_cells, col.name = "T_cells")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$Myeloid, col.name = "Myeloid")

basal <- FeaturePlot(snRNAseq_merged_integrated, features = "Basal", reduction = "umap")
pdf("snRNAseq_Basal.pdf")
basal + scale_fill_gradientn(colors = c("blue", "white", "red"))  
dev.off() 
lumHR <- FeaturePlot(snRNAseq_merged_integrated, features = "LumHR", reduction = "umap")
pdf("snRNAseq_LumHR.pdf")
lumHR + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
lumSec <- FeaturePlot(snRNAseq_merged_integrated, features = "LumSec", reduction = "umap")
pdf("snRNAseq_LumSec.pdf")
lumSec + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
fibroblasts <- FeaturePlot(snRNAseq_merged_integrated, features = "Fibroblasts", reduction = "umap")
pdf("snRNAseq_Fibroblasts.pdf")
fibroblasts + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
lymphatic <- FeaturePlot(snRNAseq_merged_integrated, features = "Lymphatic", reduction = "umap")
pdf("snRNAseq_Lymphatic.pdf")
lymphatic + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
mast <- FeaturePlot(snRNAseq_merged_integrated, features = "Mast", reduction = "umap")
pdf("snRNAseq_Mast.pdf")
mast + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
vascular <- FeaturePlot(snRNAseq_merged_integrated, features = "Vascular", reduction = "umap")
pdf("snRNAseq_Vascular.pdf")
vascular + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
perivascular <- FeaturePlot(snRNAseq_merged_integrated, features = "Perivascular", reduction = "umap")
pdf("snRNAseq_Perivascular.pdf")
perivascular + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
adipocytes <- FeaturePlot(snRNAseq_merged_integrated, features = "Adipocytes", reduction = "umap")
pdf("snRNAseq_Adipocytes.pdf")
adipocytes + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
t_cells <- FeaturePlot(snRNAseq_merged_integrated, features = "T_cells", reduction = "umap")
pdf("snRNAseq_T_cells.pdf")
t_cells + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
myeloid <- FeaturePlot(snRNAseq_merged_integrated, features = "Myeloid", reduction = "umap")
pdf("snRNAseq_Myeloid.pdf")
myeloid + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

save("snRNAseq_merged_integrated", file = "snRNAseq_merged_integrated.rda")


# Cluster analysis of integrated dataset

snRNAseq_merged_integrated <- snRNAseq_merged_integrated %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = c(0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) 
clustree(snRNAseq_merged_integrated)
Idents(snRNAseq_merged_integrated) <- snRNAseq_merged_integrated$RNA_snn_res.0.1

pdf("UMAP_integrated_cluster.pdf", width = 12, height = 7)
DimPlot(object = snRNAseq_merged_integrated, reduction = "umap")
dev.off()

save(snRNAseq_merged_integrated, file = "snRNAseq_merged_integrated.rda")

load("snRNAseq_merged_integrated.rda")

#Findallmarkers

# Cluster markers

snRNAseq_merged_integrated <- JoinLayers(snRNAseq_merged_integrated)
save(snRNAseq_merged_integrated, file = "snRNAseq_merged_integrated.rda")


snRNAseq_0.1_markers <- FindAllMarkers(snRNAseq_merged_integrated, only.pos = TRUE)


# Visualize datasets
pdf("UMAP_integrated_labelled.pdf", width = 12, height = 7)
DimPlot(object = snRNAseq_merged_integrated, reduction = "umap",group.by = "RNA_snn_res.0.1", label = TRUE)
dev.off()


# To determine cut off for Mast and Perivascular clusters

pdf("Mast_violin.pdf", width = 12, height = 7)
Mast_violin <- VlnPlot(object = snRNAseq_merged_integrated, features = "Mast")
dev.off()

pdf("Perivascular_violin.pdf")
VlnPlot(snRNAseq_merged_integrated, features = "Perivascular")
dev.off()

snRNAseq_subset <- subset(x = snRNAseq_merged_integrated, idents = "11", invert = TRUE)

pdf("UMAP_integrated_subset.pdf", width = 12, height = 7)
DimPlot(object = snRNAseq_subset, reduction = "umap",group.by = "RNA_snn_res.0.1", label = TRUE)
dev.off()

snRNAseq_subset$Mast_high <- ifelse(snRNAseq_subset$Mast > 4, "high", "low")
snRNAseq_subset$Perivascular_high <- ifelse(snRNAseq_subset$Perivascular > 3, "high", "low")


new.cluster.ids <- c("Fibroblast", "Basal", "LumHR", "LumSec", "Vascular", "Myeloid", "Adipocytes", "T_cells",
                     "Lymphatic", "LumHR", "LumHR","11")
names(new.cluster.ids) <- levels(snRNAseq_subset)
snRNAseq_subset <- RenameIdents(snRNAseq_subset, new.cluster.ids)
snRNAseq_subset$cell_type <- Idents(snRNAseq_subset)


names(new.cluster.ids) <- levels(snRNAseq_merged_integrated)
snRNAseq_merged_integrated <- RenameIdents(snRNAseq_merged_integrated, new.cluster.ids)
snRNAseq_merged_integrated$cell_type <- Idents(snRNAseq_merged_integrated)


# We need to change Myeloid (Pal vs Kumar) name to perform RCTD

current.cluster.ids <- c("Lymphatic", "Adipocytes", "Vascular", "Myeloid", "Fibroblast", "Perivascular", 
                         "Mast", "T_cells", "Basal", "LumSec", "LumHR")


new.cluster.ids <- c("Lymphatic", "Adipocytes", "Vascular", "Myel", "Fibroblast", "Perivascular", 
                         "Mast", "T_cells", "Basal", "LumSec", "LumHR")

Idents(snRNAseq_subset) <- plyr::mapvalues(x = Idents(snRNAseq_subset), from = current.cluster.ids, to = new.cluster.ids)


save(snRNAseq_subset, file = "snRNAseq_subset.rda")

pdf("UMAP_subset.pdf")
DimPlot(snRNAseq_subset, reduction = "umap", group.by = "RNA_snn_res.0.1", label = TRUE)
dev.off()


# Visualize the distribution of populations
plot_data <- as.data.frame(table(Idents(snRNAseq_merged_integrated), snRNAseq_merged_integrated$orig.ident))
colnames(plot_data) <- c("Cluster", "Sample", "Count")

ggplot(plot_data, aes(x = Sample, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(labels = scales::percent) + 
  labs(fill = "Cell population"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),  
    axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.title = element_blank(),  
    axis.text.y = element_blank()  
  )


