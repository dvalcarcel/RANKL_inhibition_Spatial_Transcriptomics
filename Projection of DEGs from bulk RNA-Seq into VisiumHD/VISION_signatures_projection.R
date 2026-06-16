setwd("PATH")

library(Seurat)
library(VISION)
library(ggplot2)

### Projection of denosumab-modulated genes from the bulk RNA-Seq into Visium HD data ###

# To project DEGs from bulk RNA-Seq obtained in the ssGSEA

load("Patient_38_RCTD.rda")

Patient_38_RCTD_projection <- Patient_38_RCTD@reductions$umap@cell.embeddings


Patient_38_RCTD_vis <- Vision(Patient_38_RCTD@assays$Spatial$data,signatures=c("ONLY_EXP_ssGSEA.gmt"), meta = Patient_38_RCTD@meta.data, pool=FALSE)
Patient_38_RCTD_vis <- addProjection(Patient_38_RCTD_vis, "UMAP", Patient_38_RCTD_projection)
Patient_38_RCTD_vis <- analyze(Patient_38_RCTD_vis)

save("Patient_38_RCTD_vis", file = "Patient_38_RCTD_vis.rda")


# To add the signature scores into the Seurat's object

signature_scores <- as.data.frame(Patient_38_RCTD_vis@SigScores)

all(rownames(signature_scores) %in% colnames(Patient_38_RCTD))

Patient_38_RCTD <- AddMetaData(Patient_38_RCTD, metadata = signature_scores)



### Plots for each signature ### 

# Define output directory
output_dir <- "signature_plots_38/"

# Extract signature names from metadata
signature_names <- colnames(signature_scores)

# Generate and save SpatialFeaturePlots for each signature
for (signature in signature_names) {
  p <- SpatialFeaturePlot(Patient_38_RCTD, features = signature,image.alpha=0,pt.size.factor=3.2) + ggtitle(signature)
  
  # Save the plot as a PDF
  pdf_filename <- paste0(output_dir, signature, ".pdf")
  ggsave(pdf_filename, plot = p, width = 10, height = 7)
  
  print(paste("Saved:", pdf_filename))
}



# Generate and save VlnPlots for each signature
for (signature in signature_names) {
  p <- VlnPlot(Patient_38_RCTD, features = signature)
  
  # Save the plot as a PDF
  pdf_filename <- paste0(output_dir, signature, "_VlnPlot.pdf")
  ggsave(pdf_filename, plot = p, width = 10, height = 7)
  
  print(paste("Saved:", pdf_filename))
}


df <- Patient_38_RCTD@meta.data

for (signature in signature_names) {
  df_long <- df %>%
    pivot_longer(
      cols = signature,
      names_to = "Score",
      values_to = "Value"
    )
  
  df_long$Denosumab <- factor(
    df_long$Denosumab,
    levels = c("Pre", "Post")
  )
  p <- ggplot(df_long, aes(x = final_annotation, y = Value, fill = Denosumab)) +
    geom_half_violin(
      data = subset(df_long, Denosumab == "Pre"),
      side = "l",
      alpha = 0.6
    ) +
    geom_half_violin(
      data = subset(df_long, Denosumab == "Post"),
      side = "r",
      alpha = 0.6
    ) +
    facet_wrap(~Score, scales = "free_y") +
    scale_fill_manual(values = c("Pre" = "#1f77b4", "Post" = "#ff7f0e")) +
    theme_minimal() +
    ylab("Score") +
    xlab("Cell idents") +
    ggtitle(signature) + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
      legend.text  = element_text(size = 18),
      axis.text.x  = element_text(size = 18),  
      axis.text.y  = element_text(size = 18),  
      axis.title.x = element_text(size = 18),  
      axis.title.y = element_text(size = 18),
      strip.text = element_text(size = 18, face = "bold")
    )
  
  # Save the plot as a PDF
  pdf_filename <- paste0(output_dir, signature, "_VlnPlot_Denosumab.pdf")
  ggsave(pdf_filename, plot = p, width = 40, height = 10)
  
  print(paste("Saved:", pdf_filename))
}


save(Patient_38_RCTD, file = "Patient_38_RCTD.rda")


