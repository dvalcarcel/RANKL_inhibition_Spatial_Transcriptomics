library(dplyr)
library(ggplot2)
library(Hmisc)
library(scales)
library(stats)
library(RColorBrewer)
library(tidyr)
library(Seurat)
library(STDistance)
library(patchwork)
library(ggpubr)


setwd("PATH")


######### Patient 38 #############


load("Patient_38_RCTD.rda")

# Extract coordinates from the Seurat object
coords_38 <- GetTissueCoordinates(Patient_38_RCTD)

#Create a data frame type “tissue_posi” from the STDistance tutorial
tissue_posi_38 <- data.frame(
  barcode = rownames(coords_38),
  in_tissue = 1,  
  array_row = NA,
  array_col = NA,
  pxl_row_in_fullres = coords_38$y,   
  pxl_col_in_fullres = coords_38$x,   
  Sample = "Patient_38",
  Sampleid = 38,
  Newbarcode = paste0(rownames(coords_38), "_38")
)

# To create metadata from Seurat
meta_38 <- Patient_38_RCTD@meta.data
meta_38$Newbarcode <- paste0(rownames(meta_38), "_38")


# Merging coordinates with metadata
posi_38 <- merge(
  x = tissue_posi_38,
  y = meta_38,
  by = "Newbarcode",
  all.y = TRUE
)

posi_38 <- posi_38[!is.na(posi_38$final_annotation), ]


posi_38_Pre <- posi_38 %>% filter(Denosumab %in% "Pre")
posi_38_Post <- posi_38 %>% filter(Denosumab %in% "Post")




########## Tumor cells as center ##################

########## Myeloid cells as targets ###############

# To calculate nearest distances to Myeloid cells

distance_results_38_TumortoMyeloid_Pre <- calculate_nearest_distances(
  posi_38_Pre,
  reference_type = "tumor",   
  target_types = "Myeloid",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)
distance_results_38_TumortoMyeloid_Pre$Condition  <- "Before"
distance_results_38_TumortoMyeloid_Pre$Patient  <- "38"


distance_results_38_TumortoMyeloid_Post <- calculate_nearest_distances(
  posi_38_Post,
  reference_type = "tumor",   
  target_types = "Myeloid",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)
distance_results_38_TumortoMyeloid_Post$Condition  <- "After"
distance_results_38_TumortoMyeloid_Post$Patient  <- "38"



df_38 <- bind_rows(
  distance_results_38_TumortoMyeloid_Pre,
  distance_results_38_TumortoMyeloid_Post
) 

head(df_38)



########## T cells as targets ###############

# To calculate nearest distances to T cells

distance_results_38_TumortoT_Pre <- calculate_nearest_distances(
  posi_38_Pre,
  reference_type = "tumor",   
  target_types = "T_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)

distance_results_38_TumortoT_Pre$Condition  <- "Before"
distance_results_38_TumortoT_Pre$Patient  <- "38"

distance_results_38_TumortoT_Post <- calculate_nearest_distances(
  posi_38_Post,
  reference_type = "tumor",   
  target_types = "T_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)

distance_results_38_TumortoT_Post$Condition  <- "After"
distance_results_38_TumortoT_Post$Patient  <- "38"

distance_results_38_TumortoT <- bind_rows(
  distance_results_38_TumortoT_Pre,
  distance_results_38_TumortoT_Post
) 

df_final_38 <- df_38 %>%
  left_join(distance_results_38_TumortoT,
            by = c("Newbarcode", "Condition", "Patient"))



########## B cells as targets ###############

# To calculate nearest distances to B cells

distance_results_38_TumortoB_Pre <- calculate_nearest_distances(
  posi_38_Pre,
  reference_type = "tumor",   
  target_types = "B_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)

distance_results_38_TumortoB_Post <- calculate_nearest_distances(
  posi_38_Post,
  reference_type = "tumor",   
  target_types = "B_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)

distance_results_38_TumortoB_Pre$Condition  <- "Before"
distance_results_38_TumortoB_Pre$Patient  <- "38"

distance_results_38_TumortoB_Post$Condition  <- "After"
distance_results_38_TumortoB_Post$Patient  <- "38"

distance_results_38_TumortoB <- bind_rows(
  distance_results_38_TumortoB_Post,
  distance_results_38_TumortoB_Pre
) 

df_final_38 <- df_final_38 %>%
  left_join(distance_results_38_TumortoB,
            by = c("Newbarcode", "Condition", "Patient"))







######### Patient 61 #############


load("Patient_61_RCTD.rda")

# Extract coordinates from the Seurat object
coords_61 <- GetTissueCoordinates(Patient_61_RCTD)

#Create a data frame type “tissue_posi” from the STDistance tutorial
tissue_posi_61 <- data.frame(
  barcode = rownames(coords_61),
  in_tissue = 1,  
  array_row = NA,
  array_col = NA,
  pxl_row_in_fullres = coords_61$y,   
  pxl_col_in_fullres = coords_61$x,   
  Sample = "Patient_61",
  Sampleid = 61,
  Newbarcode = paste0(rownames(coords_61), "_61")
)

# To create metadata from Seurat
meta_61 <- Patient_61_RCTD@meta.data
meta_61$Newbarcode <- paste0(rownames(meta_61), "_61")


# Merging coordinates with metadata
posi_61 <- merge(
  x = tissue_posi_61,
  y = meta_61,
  by = "Newbarcode",
  all.y = TRUE
)

posi_61 <- posi_61[!is.na(posi_61$final_annotation), ]


posi_61_Pre <- posi_61 %>% filter(Denosumab %in% "Pre")
posi_61_Post <- posi_61 %>% filter(Denosumab %in% "Post")




########## Tumor cells as center ##################

########## Myeloid cells as targets ###############

# To calculate nearest distances to Myeloid cells

distance_results_61_TumortoMyeloid_Pre <- calculate_nearest_distances(
  posi_61_Pre,
  reference_type = "tumor",   
  target_types = "Myeloid",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)
distance_results_61_TumortoMyeloid_Pre$Condition  <- "Before"
distance_results_61_TumortoMyeloid_Pre$Patient  <- "61"


distance_results_61_TumortoMyeloid_Post <- calculate_nearest_distances(
  posi_61_Post,
  reference_type = "tumor",   
  target_types = "Myeloid",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)
distance_results_61_TumortoMyeloid_Post$Condition  <- "After"
distance_results_61_TumortoMyeloid_Post$Patient  <- "61"



df_61 <- bind_rows(
  distance_results_61_TumortoMyeloid_Pre,
  distance_results_61_TumortoMyeloid_Post
) 

head(df_61)



########## T cells as targets ###############

# To calculate nearest distances to T cells

distance_results_61_TumortoT_Pre <- calculate_nearest_distances(
  posi_61_Pre,
  reference_type = "tumor",   
  target_types = "T_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)

distance_results_61_TumortoT_Pre$Condition  <- "Before"
distance_results_61_TumortoT_Pre$Patient  <- "61"

distance_results_61_TumortoT_Post <- calculate_nearest_distances(
  posi_61_Post,
  reference_type = "tumor",   
  target_types = "T_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)

distance_results_61_TumortoT_Post$Condition  <- "After"
distance_results_61_TumortoT_Post$Patient  <- "61"

distance_results_61_TumortoT <- bind_rows(
  distance_results_61_TumortoT_Pre,
  distance_results_61_TumortoT_Post
) 

df_final_61 <- df_61 %>%
  left_join(distance_results_61_TumortoT,
            by = c("Newbarcode", "Condition", "Patient"))



########## B cells as targets ###############

# To calculate nearest distances to B cells

distance_results_61_TumortoB_Pre <- calculate_nearest_distances(
  posi_61_Pre,
  reference_type = "tumor",   
  target_types = "B_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)

distance_results_61_TumortoB_Post <- calculate_nearest_distances(
  posi_61_Post,
  reference_type = "tumor",   
  target_types = "B_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)

distance_results_61_TumortoB_Pre$Condition  <- "Before"
distance_results_61_TumortoB_Pre$Patient  <- "61"

distance_results_61_TumortoB_Post$Condition  <- "After"
distance_results_61_TumortoB_Post$Patient  <- "61"

distance_results_61_TumortoB <- bind_rows(
  distance_results_61_TumortoB_Post,
  distance_results_61_TumortoB_Pre
) 

df_final_61 <- df_final_61 %>%
  left_join(distance_results_61_TumortoB,
            by = c("Newbarcode", "Condition", "Patient"))




######### Patient 06 #############


load("Patient_06_RCTD.rda")

# Extract coordinates from the Seurat object
coords_06 <- GetTissueCoordinates(Patient_06_RCTD)

#Create a data frame type “tissue_posi” from the STDistance tutorial
tissue_posi_06 <- data.frame(
  barcode = rownames(coords_06),
  in_tissue = 1,  
  array_row = NA,
  array_col = NA,
  pxl_row_in_fullres = coords_06$y,   
  pxl_col_in_fullres = coords_06$x,   
  Sample = "Patient_06",
  Sampleid = 06,
  Newbarcode = paste0(rownames(coords_06), "_06")
)

# To create metadata from Seurat
meta_06 <- Patient_06_RCTD@meta.data
meta_06$Newbarcode <- paste0(rownames(meta_06), "_06")


# Merging coordinates with metadata
posi_06 <- merge(
  x = tissue_posi_06,
  y = meta_06,
  by = "Newbarcode",
  all.y = TRUE
)

posi_06 <- posi_06[!is.na(posi_06$final_annotation), ]


posi_06_Pre <- posi_06 %>% filter(Denosumab %in% "Pre")
posi_06_Post <- posi_06 %>% filter(Denosumab %in% "Post")




########## Tumor cells as center ##################

########## Myeloid cells as targets ###############

# To calculate nearest distances to Myeloid cells

distance_results_06_TumortoMyeloid_Pre <- calculate_nearest_distances(
  posi_06_Pre,
  reference_type = "tumor",   
  target_types = "Myeloid",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)
distance_results_06_TumortoMyeloid_Pre$Condition  <- "Before"
distance_results_06_TumortoMyeloid_Pre$Patient  <- "06"


distance_results_06_TumortoMyeloid_Post <- calculate_nearest_distances(
  posi_06_Post,
  reference_type = "tumor",   
  target_types = "Myeloid",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)
distance_results_06_TumortoMyeloid_Post$Condition  <- "After"
distance_results_06_TumortoMyeloid_Post$Patient  <- "06"



df_06 <- bind_rows(
  distance_results_06_TumortoMyeloid_Pre,
  distance_results_06_TumortoMyeloid_Post
) 

head(df_06)



########## T cells as targets ###############

# To calculate nearest distances to T cells

distance_results_06_TumortoT_Pre <- calculate_nearest_distances(
  posi_06_Pre,
  reference_type = "tumor",   
  target_types = "T_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)

distance_results_06_TumortoT_Pre$Condition  <- "Before"
distance_results_06_TumortoT_Pre$Patient  <- "06"

distance_results_06_TumortoT_Post <- calculate_nearest_distances(
  posi_06_Post,
  reference_type = "tumor",   
  target_types = "T_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)

distance_results_06_TumortoT_Post$Condition  <- "After"
distance_results_06_TumortoT_Post$Patient  <- "06"

distance_results_06_TumortoT <- bind_rows(
  distance_results_06_TumortoT_Pre,
  distance_results_06_TumortoT_Post
) 

df_final_06 <- df_06 %>%
  left_join(distance_results_06_TumortoT,
            by = c("Newbarcode", "Condition", "Patient"))



########## B cells as targets ###############

# To calculate nearest distances to B cells

distance_results_06_TumortoB_Pre <- calculate_nearest_distances(
  posi_06_Pre,
  reference_type = "tumor",   
  target_types = "B_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)

distance_results_06_TumortoB_Post <- calculate_nearest_distances(
  posi_06_Post,
  reference_type = "tumor",   
  target_types = "B_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)

distance_results_06_TumortoB_Pre$Condition  <- "Before"
distance_results_06_TumortoB_Pre$Patient  <- "06"

distance_results_06_TumortoB_Post$Condition  <- "After"
distance_results_06_TumortoB_Post$Patient  <- "06"

distance_results_06_TumortoB <- bind_rows(
  distance_results_06_TumortoB_Post,
  distance_results_06_TumortoB_Pre
) 

df_final_06 <- df_final_06 %>%
  left_join(distance_results_06_TumortoB,
            by = c("Newbarcode", "Condition", "Patient"))




df_final_all <- bind_rows(df_final_38, df_final_61, df_final_06)

write.table(df_final_all,file = "STDistance_All_Patients.txt",row.names = FALSE,col.names = TRUE,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")




### Raincloud plots ###

df_final_all <- read.table(file = "STDistance_All_Patients.txt", sep = "\t", header = TRUE)

df_final_all <- df_final_all %>%
  mutate(Condition = factor(Condition, levels = c("Before", "After")))


df_final_all <- df_final_all %>%
  mutate(Patient = factor(Patient, levels = c("6", "38", "61")))


pd <- position_dodge(width = 0.8)

# p value by patient + BH adjust
p_df <- df_final_all %>%
  group_by(Patient) %>%
  filter(all(c("Before","After") %in% Condition)) %>%
  summarise(
    p = wilcox.test(Myeloid ~ Condition, data = cur_data())$p.value,
    y.position = max(Myeloid, na.rm = TRUE) * 1.10,
    .groups = "drop"
  ) %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    p.adj.label = paste0("p.adj=", scales::pvalue(p.adj, accuracy = 0.001))
  )

png("STDistance_All_patients_TumorToMyeloid_Raincloud.png", width = 12000, height = 8000, res = 650)

ggplot(df_final_all, aes(x = Patient, y = Myeloid, fill = Condition)) +
  geom_violin(
    position = pd,
    alpha = 0.25,
    width = 1.5,
    trim = TRUE,
    color = NA
  ) +
  geom_point(
    aes(color = Condition),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
    alpha = 0.45,
    size = 1.6
  ) +
  gghalves::geom_half_boxplot(
    side = "r",
    position = pd,
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.55
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 3,
    color = "red",
    position = pd
  ) +
  stat_summary(
    fun.data = function(x) {
      m <- mean(x, na.rm = TRUE)
      se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
      data.frame(y = m, ymin = m - se, ymax = m + se)
    },
    geom = "errorbar",
    width = 0.12,
    color = "darkred",
    position = pd
  ) +
  geom_text(
    data = p_df,
    aes(x = Patient, y = y.position, label = p.adj.label),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold"
  ) +
  
  scale_fill_manual(values = c("Before" = "#1f77b4", "After" = "#ff7f0e")) +
  scale_color_manual(values = c("Before" = "#1f77b4", "After" = "#ff7f0e"), guide = "none") +
  
  labs(
    x = "Patient",
    y = "Nearest Neighbor Distance from Tumor to Myeloid cells (µm)",
    title = "STDistance - Tumor as reference and Myeloid cells as targets",
    subtitle = "Raincloud plot with mean ± SEM, Wilcoxon test by patient (p.adj BH)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 16, face = "bold"),
    axis.text.x  = element_text(size = 16, face = "bold"),
    axis.text.y  = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold")
  )
dev.off()





