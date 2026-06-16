library(Seurat)
library(spacexr)
library(SingleCellExperiment)
library(Matrix)
library(dplyr)



####### Perform RCTD using Pal et al. scRNA-seq ER+ dataset #######

# Step 1: Load scRNA Seurat object

load("PATH/Pal_2021.processed.seurat.RData")  

counts <- processed.seurat@assays$RNA@counts
rownames(counts) <- rownames(processed.seurat)
colnames(counts) <- colnames(processed.seurat)


cell_types <- processed.seurat$celltype
names(cell_types) <- colnames(processed.seurat) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- processed.seurat$nCount_RNA; names(nUMI) <- colnames(processed.seurat) # create nUMI named list

reference <- Reference(counts, cell_types, nUMI)


# Step 2: Load VisiumHD Seurat object

load("PATH/Patient_38_16um.rda")  

# Get tissue coordinates using Seurat v5 function
coords <- GetTissueCoordinates(Patient_38)

# Extract counts matrix
counts <- as.matrix(Patient_38@assays$Spatial$counts)

# Sanity check
stopifnot(all(rownames(coords) == colnames(counts)))

# Create SpatialRNA object for RCTD
spatial <- SpatialRNA(coords = coords[, c("x", "y")], counts = counts)



# Step 3: Run RCTD

rctd <- create.RCTD(spatialRNA = spatial, reference = reference, max_cores = 4)

# Run the model (doublet mode is recommended for complex tissues)
rctd <- run.RCTD(rctd, doublet_mode = 'full')
Patient_38_RCTD <- subset(Patient_38, cells = rownames(rctd@results$weights))

# Step 4: Add results to Seurat object

results <- rctd@results$weights  # cell type proportions per spot
Bhupinder_celltype <- as.data.frame(results)

# Add as metadata to Seurat object
Patient_38_RCTD <- AddMetaData(Patient_38_RCTD, metadata = Bhupinder_celltype)
dominant_celltype <- colnames(results)[apply(results, 1, which.max)]
dominant_df <- data.frame(dominant_celltype)
rownames(dominant_df) <- rownames(results)
Patient_38_RCTD <- AddMetaData(Patient_38_RCTD, metadata = dominant_df)

pdf("Patient_38_Pal_celltype.pdf")
SpatialDimPlot(Patient_38_RCTD, group.by ="dominant_celltype", image.alpha=0, pt.size.factor=3.2)
dev.off()

save(Patient_38_RCTD, file="Patient_38_RCTD.rda")



####### Perform RCTD using Kumar et al. snRNA-seq normal breast dataset #######

# Step 1: Load snRNAseq Seurat object

load("snRNAseq_subset.rda")

counts <- snRNAseq_subset@assays$RNA@layers$counts
rownames(counts) <- rownames(snRNAseq_subset)
colnames(counts) <- colnames(snRNAseq_subset)


cell_types <- snRNAseq_subset$cell_type
names(cell_types) <- colnames(snRNAseq_subset) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- snRNAseq_subset$nCount_RNA; names(nUMI) <- colnames(snRNAseq_subset) # create nUMI named list

reference <- Reference(counts, cell_types, nUMI)


# Step 2: Load VisiumHD Seurat object

load("Patient_38_ER_RCTD.rda")

# Get tissue coordinates using Seurat v5 function
coords <- GetTissueCoordinates(Patient_38_RCTD)

# Extract counts matrix
counts <- as.matrix(Patient_38_RCTD@assays$Spatial$counts)

# Sanity check
stopifnot(all(rownames(coords) == colnames(counts)))

# Create SpatialRNA object for RCTD
spatial <- SpatialRNA(coords = coords[, c("x", "y")], counts = counts)


# Step 3: Run RCTD

rctd <- create.RCTD(spatialRNA = spatial, reference = reference, max_cores = 4)

# Run the model (doublet mode is recommended for complex tissues)
rctd <- run.RCTD(rctd, doublet_mode = 'full')
Patient_38_RCTD <- subset(Patient_38_RCTD, cells = rownames(rctd@results$weights))


# Step 4: Add results to Seurat object

results <- rctd@results$weights  # cell type proportions per spot
Kumar_celltype <- as.data.frame(results)

# Add as metadata to Seurat object
Patient_38_RCTD <- AddMetaData(Patient_38_RCTD, metadata = Kumar_celltype)
Kumar_dominant_celltype <- colnames(results)[apply(results, 1, which.max)]
Kumar_dominant_df <- data.frame(Kumar_dominant_celltype)
rownames(Kumar_dominant_df) <- rownames(results)
Patient_38_RCTD <- AddMetaData(Patient_38_RCTD, metadata = Kumar_dominant_df)

pdf("Patient_38_Kumar_celltype.pdf")
SpatialDimPlot(Patient_38_RCTD, group.by ="Kumar_dominant_celltype", image.alpha=0, pt.size.factor=3.2)
dev.off()

save(Patient_38_RCTD, file="Patient_38_RCTD_snRNAseq.rda")




####### Adding metadata according to Pal and Kumar weights #######

setwd("PATH")

## Patient 38

load("Patient_38_RCTD_snRNAseq.rda")


results <- data.frame(
  B_cells = Patient_38_RCTD$`B cells`,
  Endothelial_cells = Patient_38_RCTD$`Endothelial cells`,
  Epithelial_cells = Patient_38_RCTD$`Epithelial cells`,
  Malignant = Patient_38_RCTD$malignant,
  Myeloid_Pal = Patient_38_RCTD$Myeloid,
  Plasmablasts = Patient_38_RCTD$plasmablasts,
  Stromal = Patient_38_RCTD$Stromal,
  T_cells_Pal = Patient_38_RCTD$`T cells`,
  Lymphatic = Patient_38_RCTD$Lymphatic,
  Adipocytes = Patient_38_RCTD$Adipocytes,
  Vascular = Patient_38_RCTD$Vascular,
  Myeloid_Kumar = Patient_38_RCTD$Myel_Kumar,
  Fibroblast = Patient_38_RCTD$Fibroblast,
  Perivascular = Patient_38_RCTD$Perivascular,
  Mast = Patient_38_RCTD$Mast,
  T_cells_Kumar = Patient_38_RCTD$T_cells,
  Basal = Patient_38_RCTD$Basal,
  LumSec = Patient_38_RCTD$LumSec,
  LumHR = Patient_38_RCTD$LumHR
)

# Add as metadata to Seurat object
Final_dominant_celltype <- colnames(results)[apply(results, 1, which.max)]
Final_dominant_df <- data.frame(Final_dominant_celltype)
rownames(Final_dominant_df) <- rownames(results)
Patient_38_RCTD <- AddMetaData(Patient_38_RCTD, metadata = Final_dominant_df)


pdf("Patient_38_final_dominat_celltype.pdf")
SpatialDimPlot(Patient_38_RCTD, group.by ="Final_dominant_celltype", image.alpha=0, pt.size.factor=3.2)
dev.off()

Idents(Patient_38_RCTD) <- Patient_38_RCTD$Final_dominant_celltype

Final_dominant_celltype_ptx38 <- FindAllMarkers(Patient_38_RCTD, only.pos = TRUE)

write.table(Final_dominant_celltype_ptx38, file = "Final_dominant_celltype_ptx38_markers.txt", sep = "\t")

save(Patient_38_RCTD, file="Patient_38_RCTD_snRNAseq.rda")





####### Final annotation #######

# We decided to maintain the epithelial annotation of the pathologist and annotate the stroma by the RCTD's weights

dominant_cell_type_curated <- data.frame(Patient_38_RCTD$Final_dominant_celltype)
colnames(dominant_cell_type_curated) <- "dominant_cell_type_curated"

Patient_38_RCTD <- AddMetaData(Patient_38_RCTD, metadata = dominant_cell_type_curated)

Idents(Patient_38_RCTD) <- Patient_38_RCTD$dominant_cell_type_curated

# Rename values
Patient_38_RCTD@meta.data$dominant_cell_type_curated[
  Patient_38_RCTD@meta.data$dominant_cell_type_curated == "Myeloid_Kumar"
] <- "Myeloid"

Patient_38_RCTD@meta.data$dominant_cell_type_curated[
  Patient_38_RCTD@meta.data$dominant_cell_type_curated == "Myeloid_Pal"
] <- "Myeloid"

Patient_38_RCTD@meta.data$dominant_cell_type_curated[
  Patient_38_RCTD@meta.data$dominant_cell_type_curated == "T_cells_Kumar"
] <- "T_cells"

Patient_38_RCTD@meta.data$dominant_cell_type_curated[
  Patient_38_RCTD@meta.data$dominant_cell_type_curated == "T_cells_Pal"
] <- "T_cells"

Patient_38_RCTD@meta.data$dominant_cell_type_curated[
  Patient_38_RCTD@meta.data$dominant_cell_type_curated == "Stromal"
] <- "Fibroblast"

Patient_38_RCTD@meta.data$dominant_cell_type_curated[
  Patient_38_RCTD@meta.data$dominant_cell_type_curated == "Vascular"
] <- "Endothelial_cells"


Idents(Patient_38_RCTD) <- Patient_38_RCTD$dominant_cell_type_curated

save(Patient_38_RCTD, file = "Patient_38_RCTD_snRNAseq.rda")


# Create a new metadata to annotate the stroma

Patient_38_RCTD@meta.data$Pathologist_annotation_refined <- as.character(Patient_38_RCTD@meta.data$Pathologist_annotation)

# To identify spots annotated as "stroma" by the pathologist
stromal_idx <- which(Patient_38_RCTD@meta.data$Pathologist_annotation %in% "stroma")

# To replace this annotation by those of RCTD (Pal and Kumar annotations)
Patient_38_RCTD@meta.data$Pathologist_annotation_refined[stromal_idx] <- as.character(Patient_38_RCTD$dominant_cell_type_curated[stromal_idx])

save(Patient_38_RCTD, file = "Patient_38_RCTD_snRNAseq.rda")


####### Epithelial mixed annotation ####### 

# We detect that some stromal spots were annotated as epithelial and were located near to tumor or normal breast spots.

# We decided to refine the annotation by adding two new categories: tumor_mixed and normal_breast_mixed.


# To create a new metadata based on the pathologist’s annotation and the differentiation of the stroma (Pal and Kumar):
Patient_38_RCTD@meta.data$final_annotation <- as.character(Patient_38_RCTD@meta.data$Pathologist_annotation_refined)

# To identify epithelial cells annotated in the stroma:
basal_idx <- which(Patient_38_RCTD@meta.data$final_annotation %in% "Basal")
epithelial_idx <- which(Patient_38_RCTD@meta.data$final_annotation %in% "Epithelial_cells")
LumHR_idx <- which(Patient_38_RCTD@meta.data$final_annotation %in% "LumHR")
LumSec_idx <- which(Patient_38_RCTD@meta.data$final_annotation %in% "LumSec")
malignant_idx <- which(Patient_38_RCTD@meta.data$final_annotation %in% "Malignant")


# To replace those values with ‘tumour_mixed’ (Pal et al.) or ‘normal_breast_mixed’ (Kumar et al.):
Patient_38_RCTD@meta.data$final_annotation[basal_idx] <- "normal_breast_mixed"
Patient_38_RCTD@meta.data$final_annotation[LumSec_idx] <- "normal_breast_mixed"
Patient_38_RCTD@meta.data$final_annotation[LumHR_idx] <- "normal_breast_mixed"
Patient_38_RCTD@meta.data$final_annotation[epithelial_idx] <- "tumor_mixed"
Patient_38_RCTD@meta.data$final_annotation[malignant_idx] <- "tumor_mixed"

Idents(Patient_38_RCTD) <- Patient_38_RCTD$final_annotation

save(Patient_38_RCTD, file = "Patient_38_RCTD_snRNAseq.rda")










