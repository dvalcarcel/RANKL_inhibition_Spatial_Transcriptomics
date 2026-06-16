setwd("PATH")

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
library(pheatmap)
library(reshape2)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)



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


# Merging coordinates with metada
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



merged_38_IFN_Pre <- merge(
  posi_38_Pre,
  distance_results_38_TumortoMyeloid_Pre,
  by = "Newbarcode",
  all.x = FALSE 
)

result_correlation_38_TumortoMyeloid_Pre <- calculate_correlations(
  spatial_data = merged_38_IFN_Pre,
  distance_results = distance_results_38_TumortoMyeloid_Pre,
  spatial_feature = "ONLY_EXP_UP_IFN_GAMMA",
  distance_metric = "Myeloid",
  method = "pearson",
  plot = TRUE,
  plot_title = "Before Denosumab"
) 

merged_38_IFN_Post <- merge(
  posi_38_Post,
  distance_results_38_TumortoMyeloid_Post,
  by = "Newbarcode",
  all.x = FALSE 
)

result_correlation_38_TumortoMyeloid_Post <- calculate_correlations(
  spatial_data = merged_38_IFN_Post,
  distance_results = distance_results_38_TumortoMyeloid_Post,
  spatial_feature = "ONLY_EXP_UP_IFN_GAMMA",
  distance_metric = "Myeloid",
  method = "pearson",
  plot = TRUE,
  plot_title = "After Denosumab"
) 


# Create a dataframe with correlation results for pre-treatment
correlation_38_TumortoMyeloid_Pre <- data.frame(
  Patient = "38",
  Treatment = "Pre",
  Target = "Myeloid",
  Correlation = result_correlation_38_TumortoMyeloid_Pre$estimate,
  P_value = result_correlation_38_TumortoMyeloid_Pre$p_value
)

# Create a dataframe with correlation results for post-treatment
correlation_38_TumortoMyeloid_Post <- data.frame(
  Patient = "38",
  Treatment = "Post",
  Target = "Myeloid",
  Correlation = result_correlation_38_TumortoMyeloid_Post$estimate,
  P_value = result_correlation_38_TumortoMyeloid_Post$p_value
)

# Combine the dataframes for "Pre" and "Post" into a single dataframe
correlation_df_myeloid <- rbind(correlation_38_TumortoMyeloid_Pre, correlation_38_TumortoMyeloid_Post)




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


merged_38_IFN_Pre <- merge(
  posi_38_Pre,
  distance_results_38_TumortoT_Pre,
  by = "Newbarcode",
  all.x = FALSE 
)

result_correlation_38_TumortoT_Pre <- calculate_correlations(
  spatial_data = merged_38_IFN_Pre,
  distance_results = distance_results_38_TumortoT_Pre,
  spatial_feature = "ONLY_EXP_UP_IFN_GAMMA",
  distance_metric = "T_cells",
  method = "pearson",
  plot = TRUE,
  plot_title = "Before Denosumab"
  
  
) 

merged_38_IFN_Post <- merge(
  posi_38_Post,
  distance_results_38_TumortoT_Post,
  by = "Newbarcode",
  all.x = FALSE 
)

result_correlation_38_TumortoT_Post <- calculate_correlations(
  spatial_data = merged_38_IFN_Post,
  distance_results = distance_results_38_TumortoT_Post,
  spatial_feature = "ONLY_EXP_UP_IFN_GAMMA",
  distance_metric = "T_cells",
  method = "pearson",
  plot = TRUE,
  plot_title = "After Denosumab"
) 


# Create a dataframe with correlation results for pre-treatment
correlation_38_TumortoT_Pre <- data.frame(
  Patient = "38",
  Treatment = "Pre",
  Target = "T_cells",
  Correlation = result_correlation_38_TumortoT_Pre$estimate,
  P_value = result_correlation_38_TumortoT_Pre$p_value
)

# Create a dataframe with correlation results for post-treatment
correlation_38_TumortoT_Post <- data.frame(
  Patient = "38",
  Treatment = "Post",
  Target = "T_cells",
  Correlation = result_correlation_38_TumortoT_Post$estimate,
  P_value = result_correlation_38_TumortoT_Post$p_value
)

# Combine the dataframes for "Pre" and "Post" into a single dataframe

correlation_df_T <- rbind(correlation_38_TumortoT_Pre, correlation_38_TumortoT_Post)



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



merged_38_IFN_Pre <- merge(
  posi_38_Pre,
  distance_results_38_TumortoB_Pre,
  by = "Newbarcode",
  all.x = FALSE 
)

result_correlation_38_TumortoB_Pre <- calculate_correlations(
  spatial_data = posi_38_Pre,
  distance_results = distance_results_38_TumortoB_Pre,
  spatial_feature = "ONLY_EXP_UP_IFN_GAMMA",
  distance_metric = "B_cells",
  method = "pearson",
  plot = TRUE,
  plot_title = "Before Denosumab"
  
  
) 

merged_38_IFN_Post <- merge(
  posi_38_Post,
  distance_results_38_TumortoB_Post,
  by = "Newbarcode",
  all.x = FALSE 
)

result_correlation_38_TumortoB_Post <- calculate_correlations(
  spatial_data = posi_38_Post,
  distance_results = distance_results_38_TumortoB_Post,
  spatial_feature = "ONLY_EXP_UP_IFN_GAMMA",
  distance_metric = "B_cells",
  method = "pearson",
  plot = TRUE,
  plot_title = "After Denosumab"
) 


# Create a dataframe with correlation results for pre-treatment
correlation_38_TumortoB_Pre <- data.frame(
  Patient = "38",
  Treatment = "Pre",
  Target = "B_cells",
  Correlation = result_correlation_38_TumortoB_Pre$estimate,
  P_value = result_correlation_38_TumortoB_Pre$p_value
)

# Create a dataframe with correlation results for post-treatment
correlation_38_TumortoB_Post <- data.frame(
  Patient = "38",
  Treatment = "Post",
  Target = "B_cells",
  Correlation = result_correlation_38_TumortoB_Post$estimate,
  P_value = result_correlation_38_TumortoB_Post$p_value
)

# Combine the dataframes for "Pre" and "Post" into a single dataframe

correlation_df_B <- rbind(correlation_38_TumortoB_Pre, correlation_38_TumortoB_Post)


correlation_df_38 <- rbind(correlation_df_myeloid, correlation_df_T, correlation_df_B)

write.table(correlation_df_38,file = "STDistance_correlation_38.txt",row.names = FALSE,col.names = TRUE,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")



## Heatmap ##


correlation_df_06 <- read.table("STDistance_correlation_06.txt", 
                                header = TRUE,  
                                sep = "\t",     
                                stringsAsFactors = FALSE,  
                                na.strings = "NA")

correlation_df_38 <- read.table("STDistance_correlation_38.txt", 
                                header = TRUE,  
                                sep = "\t",     
                                stringsAsFactors = FALSE,  
                                na.strings = "NA")

correlation_df_61 <- read.table("STDistance_correlation_61.txt", 
                                header = TRUE,  
                                sep = "\t",     
                                stringsAsFactors = FALSE,  
                                na.strings = "NA")


correlation_df <- rbind(correlation_df_06, correlation_df_38, correlation_df_61)


setDT(correlation_df)

col_fun <- colorRamp2(
  c(-1, 0, 1),
  c("#3B4CC0", "#F2F2F2", "#B40426")
)

cor_wide <- dcast(correlation_df, Patient ~ Treatment + Target, value.var = "Correlation")
p_wide   <- dcast(correlation_df, Patient ~ Treatment + Target, value.var = "P_value")

cor_mat <- as.matrix(cor_wide[, -1, with = FALSE]); rownames(cor_mat) <- cor_wide$Patient
p_mat   <- as.matrix(p_wide[, -1, with = FALSE]);   rownames(p_mat)   <- p_wide$Patient


pre_cols  <- sort(grep("^Pre_",  colnames(cor_mat), value = TRUE))
post_cols <- sort(grep("^Post_", colnames(cor_mat), value = TRUE))
new_order <- c(pre_cols, post_cols)

cor_mat <- cor_mat[, new_order, drop = FALSE]
p_mat   <- p_mat[,   new_order, drop = FALSE]


rownames(cor_mat) <- paste("Patient", rownames(cor_mat))
rownames(p_mat)   <- rownames(cor_mat)  

# p < .001 = ***, p < .01 = **, p < .05 = *, if else "" 
stars <- cut(
  p_mat,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", ""),
  right = FALSE,
  include.lowest = TRUE
)
stars <- matrix(as.character(stars), nrow = nrow(p_mat), dimnames = dimnames(p_mat))
stars[is.na(p_mat)] <- ""

labels_mat <- matrix(
  sprintf("%.2f%s", round(cor_mat, 2), stars),
  nrow = nrow(cor_mat), dimnames = dimnames(cor_mat)
)

get_txt_col <- function(fill_color) {
  rgb <- col2rgb(fill_color) / 255
  lum <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
  ifelse(lum < 0.5, "white", "black")
}

ht <- Heatmap(
  cor_mat,
  name = "Pearson\ncorrelation",
  col = col_fun,
  width = unit(20, "cm"),
  height = unit(20, "cm"),
  row_names_side = "left",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 18),
  column_names_gp = gpar(fontsize = 18),
  column_title_gp = gpar(fontsize = 20),
  heatmap_legend_param = list(
    title = "Pearson correlation",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 14, fontface = "bold"),  
    labels_gp = gpar(fontsize = 12),
    legend_height = unit(12, "cm"),
    grid_width = unit(1, "cm")
  ),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(
      labels_mat[i, j],
      x = x, y = y,
      gp = gpar(col = get_txt_col(fill), fontsize = 16, fontface = "bold"),
      just = "center",
      default.units = "native"
    )
  }
)


right_spacer <- rowAnnotation(
  spacer = anno_empty(width = unit(18, "mm"), border = FALSE)
)

ht_sp <- ht + right_spacer

png("Correlation_STDistance_Tumor_IFN.png", width = 12500, height = 8000, res = 650)
draw(
  ht_sp,
  heatmap_legend_side = "right",
  column_title = "Correlation heatmap - Distance between Tumor to target cell types and HALLMARK IFN GAMMA expression (Bulk RNA-Seq)",
  column_title_gp = gpar(fontsize = 22, fontface = "bold")
)
dev.off()



