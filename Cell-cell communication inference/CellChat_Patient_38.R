setwd("PATH")

library(Seurat)
library(CellChat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(jsonlite)


load("Patient_38_RCTD.rda")

# Remove spots where Denosumab is NA
Patient_38_RCTD <- Patient_38_RCTD[, !is.na(Patient_38_RCTD$Denosumab)]
Patient_38_RCTD <- Patient_38_RCTD[, !is.na(Patient_38_RCTD$final_annotation)]
Patient_38_RCTD$final_annotation_Denosumab <- paste(Patient_38_RCTD$final_annotation,Patient_38_RCTD$Denosumab,sep="_")
Idents(Patient_38_RCTD) <- Patient_38_RCTD$final_annotation_Denosumab

# Extract expression data and metadata
data.input <- GetAssayData(Patient_38_RCTD, assay = "Spatial", slot = "data")
meta <- data.frame(labels = Idents(Patient_38_RCTD), row.names = names(Idents(Patient_38_RCTD)))
spatial.locs = GetTissueCoordinates(Patient_38_RCTD, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs$cell <- NULL

scalefactors = fromJSON(txt = file.path("/Patient_38", 'scalefactors_json.json'))
spot.size = 16 # the theoretical spot size (um) in 10X Visium
conversion.factor = spot.size/scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

d.spatial <- computeCellDistance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)
min(d.spatial[d.spatial!=0]) 
#[1] 16.0005
#This means that the minimum distance between two adjacent spots is 16 microns

spatial.locs <- as.matrix(spatial.locs)


# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)

# Set the human ligand-receptor interaction database
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# Preprocessing
cellchat <- subsetData(cellchat)  # subset to signaling genes
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Part II: Inference of cell-cell communication network
# Compute communication probabilities with spatial constraint
# when using scale.distance=0.01 got an error with this message: Please increase the value of `scale.distance` and use a value that is slighly smaller than 4.4


cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, distance.use = TRUE, interaction.range = 50, scale.distance = 4.3, 
                              contact.dependent = TRUE, contact.range = 16)

cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)
write.table(df.net, file= "Patient_38_CellChat_Interactions.txt", sep='\t')

# Compute communication at pathway level and aggregate
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


# Visualization

# Heatmap of communication
pdf("Patient_38_CellChat_Heatmap.pdf")
netVisual_heatmap(cellchat)
dev.off()

# Save CellChat object
saveRDS(cellchat, file = "Patient_38_CellChat.rds")


# Circle plots

levels(cellchat@idents) <- sub("_Pre$", "_Biopsy", levels(cellchat@idents))
levels(cellchat@idents) <- sub("_Post$", "_Surgery", levels(cellchat@idents))

idents <- levels(cellchat@idents)

biopsy  <- sort(idents[grepl("_Biopsy$", idents)])
surgery <- sort(idents[grepl("_Surgery$", idents)])

new_idents <- c(biopsy, surgery)

cellchat@idents <- factor(cellchat@idents, levels = new_idents)

idents_to_filter <- c("B_cells_Biopsy", "Endothelial_cells_Biopsy", "Myeloid_Biopsy",
                      "Fibroblast_Biopsy",  "T_cells_Biopsy","tumour_Biopsy", "normal_breast_tissue_Biopsy",
                      "B_cells_Surgery", "Endothelial_cells_Surgery", "Myeloid_Surgery",
                      "Fibroblast_Surgery", "T_cells_Surgery", "tumour_Surgery", "normal_breast_tissue_Surgery")

biopsy_idents  <- sort(idents_to_filter[grepl("_Biopsy$", idents_to_filter)])
surgery_idents <- sort(idents_to_filter[grepl("_Surgery$", idents_to_filter)])

new_idents <- sort(idents)
cellchat@idents <- factor(cellchat@idents, levels = new_idents)

cols <- c(
  ductal_invasion_Biopsy       = "#4A0000",
  ductal_invasion_Surgery      = "#4A0000",
  tumour_Biopsy                 = "#7F0000",
  tumour_Surgery                = "#7F0000",
  tumour_mixed_Biopsy           = "#F70000",
  tumour_mixed_Surgery          = "#F70000",
  "in situ carcinoma_Biopsy"   = "#2A0000",
  "in situ carcinoma_Surgery"  = "#2A0000",
  normal_breast_tissue_Biopsy  = "#1A9850",
  normal_breast_tissue_Surgery = "#1A9850",
  normal_breast_mixed_Biopsy   = "#66C2A5",
  normal_breast_mixed_Surgery  = "#66C2A5",
  Fibroblast_Biopsy            = "#BCA9FC",
  Fibroblast_Surgery           = "#BCA9FC",
  artery_Biopsy                = "#4126AC",
  artery_Surgery               = "#4126AC",
  Endothelial_cells_Biopsy     = "#2166AC",
  Endothelial_cells_Surgery    = "#2166AC",
  Lymphatic_Biopsy             = "#4393C3",
  Lymphatic_Surgery            = "#4393C3",
  Perivascular_Biopsy          = "#92C5DE",
  Perivascular_Surgery         = "#92C5DE",
  B_cells_Biopsy               = "#F46D43",
  B_cells_Surgery              = "#F46D43",
  T_cells_Biopsy               = "#FDAE61",
  T_cells_Surgery              = "#FDAE61",
  Plasmablasts_Biopsy          = "#1CECFC",
  Plasmablasts_Surgery         = "#1CECFC",
  Myeloid_Biopsy               = "#ED91E5",
  Myeloid_Surgery              = "#ED91E5",
  Mast_Biopsy                  = "#FCE99D",
  Mast_Surgery                 = "#FCE99D",
  Adipocytes_Biopsy            = "#FFFF33",
  Adipocytes_Surgery           = "#FFFF33"
)


colnames(cellchat@net$weight) <- sub("_Pre$", "_Biopsy", colnames(cellchat@net$weight))
rownames(cellchat@net$weight) <- sub("_Pre$", "_Biopsy", rownames(cellchat@net$weight))

colnames(cellchat@net$weight) <- sub("_Post$", "_Surgery", colnames(cellchat@net$weight))
rownames(cellchat@net$weight) <- sub("_Post$", "_Surgery", rownames(cellchat@net$weight))


cols_use <- cols[colnames(cellchat@net$weight)]
setdiff(colnames(cellchat@net$weight), names(cols))


# Circle plot - Biopsy

biopsy_nodes <- intersect(biopsy_idents, rownames(cellchat@net$weight))
biopsy_nodes <- sort(biopsy_nodes)

mat_biopsy <- cellchat@net$weight[biopsy_nodes, biopsy_nodes, drop = FALSE]

cols_biopsy <- cols[biopsy_nodes]


png("Patient_38_CellChat_Circle_Biopsy.png", width = 5500, height = 5500, res = 800)
netVisual_circle(mat_biopsy, vertex.weight = rowSums(mat_biopsy), weight.scale = TRUE, label.edge = FALSE,
                 targets.use = biopsy_nodes, sources.use = biopsy_nodes, color.use = cols_biopsy, remove.isolate = TRUE, 
                 vertex.label.cex = 1.5, edge.label.cex = 1.5, arrow.size = 0.5, margin = 0.4)
dev.off()


# Circle plot - Surgery

surgery_nodes <- intersect(surgery_idents, rownames(cellchat@net$weight))

surgery_nodes <- sort(surgery_nodes)

mat_surgery <- cellchat@net$weight[surgery_nodes, surgery_nodes, drop = FALSE]

cols_surgery <- cols[surgery_nodes]


png("Patient_38_CellChat_Circle_Surgery.png", width = 5500, height = 5500, res = 800)
netVisual_circle(mat_surgery, vertex.weight = rowSums(mat_surgery), weight.scale = TRUE, label.edge = FALSE,
                 targets.use = surgery_nodes, sources.use = surgery_nodes, color.use = cols_surgery, remove.isolate = TRUE, 
                 vertex.label.cex = 1.5, edge.label.cex = 1.5, arrow.size = 0.5, margin = 0.4)
dev.off()


