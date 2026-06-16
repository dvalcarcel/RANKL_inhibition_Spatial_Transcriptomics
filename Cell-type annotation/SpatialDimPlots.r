library(Seurat)
library(ggplot2)
library(patchwork)

setwd("PATH")

# Defining the color for each annotation

cols <- c(
  ductal_invasion        = "#4A0000",  
  tumor                  = "#7F0000",  
  tumor_mixed            = "#F70000",
  "in situ carcinoma"    = "#2A0000",
  normal_breast_tissue   = "#1A9850", 
  normal_breast_mixed    = "#66C2A5",  
  Stromal                = "#BCA9FC",  
  stroma                 = "#BCA9FC",
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


########## Patient 4 and 26 ################
load("Patient_4_26_RCTD_snRNAseq.rda")


Patient_4_26_RCTD <- Patient_4_26_RCTD[, !is.na(Patient_4_RCTD$final_annotation)]


png("Patient_4_26_SpatialDimPlot_final_annotation.png", width = 10000, height = 4000, res = 850)
SpatialDimPlot(
  Patient_4_26_RCTD,
  group.by = "final_annotation",
  cols = cols,
  image.alpha = 0,
  pt.size.factor = 3.2
) +
  theme(
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 5),  
      ncol = 2
    )
  )
dev.off()



########## Patient 06 ################
load("Patient_06_RCTD.rda")


Patient_06_RCTD <- Patient_06_RCTD[, !is.na(Patient_06_RCTD$final_annotation)]


png("Patient_06_SpatialDimPlot_final_annotation.png", width = 10000, height = 4000, res = 850)
SpatialDimPlot(
  Patient_06_RCTD,
  group.by = "final_annotation",
  cols = cols,
  image.alpha = 0,
  pt.size.factor = 3.2
) +
  theme(
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 5),  
      ncol = 2
    )
  )

dev.off()




########## Patient 38 ################
load("Patient_38_RCTD.rda")


Patient_38_RCTD <- Patient_38_RCTD[, !is.na(Patient_38_RCTD$final_annotation)]

png("Patient_38_SpatialDimPlot_final_annotation.png", width = 10000, height = 4000, res = 850)

SpatialDimPlot(
  Patient_38_RCTD,
  group.by = "final_annotation",
  cols = cols,
  image.alpha = 0,
  pt.size.factor = 3.2
) +
  theme(
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 5),  
      ncol = 2
    )
  )

dev.off()


########## Patient 41 and 67 ################
load("Patient_41_67_RCTD_snRNAseq.rda")


Patient_41_67_RCTD <- Patient_41_67_RCTD[, !is.na(Patient_41_67_RCTD$final_annotation)]

png("Patient_41_67_SpatialDimPlot_final_annotation.png", width = 10000, height = 4000, res = 850)

SpatialDimPlot(
  Patient_41_67_RCTD,
  group.by = "final_annotation",
  cols = cols,
  image.alpha = 0,
  pt.size.factor = 3.2
) +
  theme(
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 5),  
      ncol = 2
    )
  )

dev.off()



########## Patient 61 ################
load("Patient_61_RCTD.rda")


Patient_61_RCTD <- Patient_61_RCTD[, !is.na(Patient_61_RCTD$final_annotation)]

png("Patient_61_SpatialDimPlot_final_annotation.png", width = 10000, height = 4000, res = 850)

SpatialDimPlot(
  Patient_61_RCTD,
  group.by = "final_annotation",
  cols = cols,
  image.alpha = 0,
  pt.size.factor = 3.2
) +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 5),  
      ncol = 2
    )
  )

dev.off()


