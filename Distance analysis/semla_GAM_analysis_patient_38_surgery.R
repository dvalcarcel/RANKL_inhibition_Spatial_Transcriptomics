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
library(ggrepel)
library(ggtext)



######### Patient 38 #############

load("Patient_38_RCTD.rda")


# To split object by Denosumab_cluster annotation
objs_list <- SplitObject(Patient_38_RCTD, split.by = "Denosumab_cluster")
Pre1_38  <- objs_list[["Pre1_38"]]
Pre2_38 <- objs_list[["Pre2_38"]]
Post_38 <- objs_list[["Post_38"]]



######### Patient 38 post #############################

Post_38 <- Post_38[, !is.na(Post_38$final_annotation)]

Idents(Post_38) <- Post_38$final_annotation

Post_38_semla <- UpdateSeuratForSemla(Post_38)

Post_38_semla <- LoadImages(Post_38_semla)


pdf("Patient_38_Post_MapLabels.pdf", height = 8, width = 15)
MapLabels(Post_38_semla, column_name = "final_annotation", 
          image_use = "raw", pt_alpha = 0.6, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))
dev.off()

Post_38_semla$tumor <- ifelse(Post_38_semla$final_annotation %in% "tumor", "tumor", NA)
pdf("Patient_38_Post_MapLabels_tumor.pdf", height = 8, width = 15)
MapLabels(Post_38_semla, column_name = "tumor", override_plot_dims = TRUE, 
          image_use = "raw", drop_na = TRUE, pt_size = 1.8) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))
dev.off()

# Radial distances from tumor spots

Post_38_semla <- RadialDistance(Post_38_semla, column_name = "final_annotation", selected_groups = "tumor")

pdf("Patient_38_Post_MapFeatures_tumor.pdf", height = 8, width = 15)
MapFeatures(Post_38_semla, features = "r_dist_tumor", center_zero = TRUE, pt_size = 1.8, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)
dev.off()

#Pixel coordinates 

Post_38_semla$r_dist_tumor <- (100/273)*Post_38_semla$r_dist_tumor

Post_38_semla <- RadialDistance(Post_38_semla, column_name = "final_annotation", 
                                selected_groups = "tumor", convert_to_microns = TRUE)


Post_38_semla$r_dist_tumor_sqrt <- sign(Post_38_semla$r_dist_tumor)*sqrt(abs(Post_38_semla$r_dist_tumor))
pdf("Patient_38_Post_MapFeatures_tumor_2.pdf", height = 8, width = 15)
MapFeatures(Post_38_semla, features = "r_dist_tumor", center_zero = TRUE, pt_size = 1.8, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)
dev.off()


sample_use   <- "Patient 38 - Post"

outdir <- "Patient 38 - GAM"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


p_dist <- MapFeatures(
  Post_38_semla,
  features = "r_dist_tumor",
  center_zero = TRUE,
  pt_size = 1.8,
  colors = brewer.pal(n = 11, name = "RdBu") |> rev(),
  override_plot_dims = TRUE
) +
  ggtitle(paste0(sample_use, " | r_dist_tumor (microns)")) +
  theme(plot.title = element_text(hjust = 0.5))

pdf(file.path(outdir, paste0(sample_use, "_RadialDistance_r_dist_tumor.pdf")), width = 9, height = 9)
print(p_dist)
dev.off()

p_dist_sqrt <- MapFeatures(
  Post_38_semla,
  features = "r_dist_tumor_sqrt",
  center_zero = TRUE,
  pt_size = 1.8,
  colors = brewer.pal(n = 11, name = "RdBu") |> rev(),
  override_plot_dims = TRUE
) +
  ggtitle(paste0(sample_use, " | r_dist_tumor_sqrt")) +
  theme(plot.title = element_text(hjust = 0.5))

pdf(file.path(outdir, paste0(sample_use, "_RadialDistance_r_dist_tumor_sqrt.pdf")), width = 9, height = 9)
print(p_dist_sqrt)
dev.off()


saveRDS(Post_38_semla, file.path(outdir, paste0(sample_use, "_semla_with_radial_distance.rds")))

message("Radial-distance outputs saved in: ", normalizePath(outdir))

#### GAM analysis #### 
## Model: y ~ s(r_dist_1, bs = "cs")
## Window: r_dist_1 in [0, 1000] microns
## Filter: gene expressed (>0) in >= 5% of spots


gam_outdir  <- file.path(outdir, "GAM_Patient 38 - GAM")
dir.create(gam_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_outdir, recursive = TRUE, showWarnings = FALSE)


## Parameters
dist_var <- "r_dist_tumor"
dist_min <- 0
dist_max <- 1000
min_detect_frac <- 0.05
slot_use <- "data"       # log-normalized data
top_n_plots <- 50
n_grid_deriv <- 200


## Helper: mean first derivative

compute_mean_derivative <- function(fit, dist_var, dist_min, dist_max, n_grid = 200) {
  grid <- data.frame(r_dist_tumor = seq(dist_min, dist_max, length.out = n_grid))
  
  pred <- predict(
    fit,
    newdata = grid,
    type = "terms",
    terms = paste0("s(", dist_var, ")"),
    se.fit = FALSE
  )
  
  yhat <- as.numeric(pred[, 1])
  x <- grid[[dist_var]]
  dydx <- diff(yhat) / diff(x)
  
  mean(dydx, na.rm = TRUE)
}

## Build metadata table

meta_df <- Post_38_semla@meta.data %>%
  as_tibble(rownames = "barcode") %>%
  filter(!is.na(.data[[dist_var]])) %>%
  filter(.data[[dist_var]] >= dist_min, .data[[dist_var]] <= dist_max)

message("Spots used for GAM after distance filtering = ", nrow(meta_df))


## Expression matrix aligned to meta_df

expr_mat <- GetAssayData(Post_38_semla, slot = slot_use)
expr_mat <- expr_mat[, meta_df$barcode, drop = FALSE]

message("Genes in expression matrix = ", nrow(expr_mat))


## Gene filtering

detect_frac <- Matrix::rowMeans(expr_mat > 0)
keep_genes <- names(detect_frac)[detect_frac >= min_detect_frac]

message("Genes passing detection filter (>= ", min_detect_frac * 100, "%) = ", length(keep_genes))

expr_mat_filt <- expr_mat[keep_genes, , drop = FALSE]


## Fit GAM per gene

# Define the GAM formula.
# For each gene, we model expression (y) as a smooth function of radial distance.
# bs = "cs" means cubic regression spline with shrinkage,
# which helps prevent overfitting if the smooth term is weak.
form <- as.formula("y ~ s(r_dist_tumor, bs = 'cs')")

# Set GAM fitting controls.
# epsilon = convergence tolerance
# maxit   = maximum number of fitting iterations
ctrl <- mgcv::gam.control(epsilon = 1e-7, maxit = 200)

# Create an empty list to store one result table per gene.
# The length equals the number of genes that passed the expression filter.
results <- vector("list", length = nrow(expr_mat_filt))

# Name each list element by gene name for easier tracking.
names(results) <- rownames(expr_mat_filt)

# Create a progress bar so you can monitor progress during the loop.
pb <- txtProgressBar(min = 0, max = nrow(expr_mat_filt), style = 3)
i <- 0

# Loop over genes one by one
for (gene in rownames(expr_mat_filt)) {
  i <- i + 1
  setTxtProgressBar(pb, i)
  
  # Extract expression values for the current gene across all selected spots.
  # This gives one numeric value per spot.
  y <- as.numeric(expr_mat_filt[gene, ])
  
  # Build a small data frame for GAM fitting.
  # Each row = one spot
  # y        = expression of current gene
  # r_dist_1 = radial distance of that spot from the ROI border
  df_g <- meta_df %>%
    transmute(
      y = y,
      r_dist_tumor = .data[[dist_var]]
    )
  
  # Fit the GAM:
  # y ~ smooth function of radial distance
  # method = "REML" is generally stable and recommended for mgcv
  # tryCatch prevents the whole script from crashing if one gene fails
  fit <- tryCatch(
    mgcv::gam(form, data = df_g, method = "REML", control = ctrl),
    error = function(e) NULL
  )
  
  # If the model fails, save NA values for this gene and continue
  if (is.null(fit)) {
    results[[gene]] <- tibble(
      gene = gene,
      n_spots = nrow(df_g),                       # number of spots used
      detect_frac = as.numeric(detect_frac[gene]),# fraction of spots expressing gene
      edf = NA_real_,                             # effective degrees of freedom.
      # edf ≈ 1 → the relationship is close to linear
      # edf > 1 → the relationship is nonlinear / curved
      # edf < 1 → simple trend
      # higher edf → more wiggliness in the fitted smooth, i.e. more complex pattern across distance
      p_smooth = NA_real_,                        # p-value of smooth term
      mean_derivative = NA_real_,                 # mean slope across distance window
      note = "fit_failed"
    )
    next
  }
  
  # Extract model summary
  sm <- summary(fit)
  
  # s.table contains statistics for smooth terms
  # including edf and significance
  s_tab <- sm$s.table
  
  # Expected name of the distance smooth term
  smooth_name <- "s(r_dist_tumor)"
  
  # Sometimes mgcv labels the smooth term slightly differently,
  # so this checks the row names more flexibly.
  if (!(smooth_name %in% rownames(s_tab))) {
    idx <- grep("^s\\(r_dist_tumor\\)", rownames(s_tab))
    if (length(idx) == 1) smooth_name <- rownames(s_tab)[idx]
  }
  
  # If we still cannot find the distance smooth term,
  # store NA for this gene and continue
  if (!(smooth_name %in% rownames(s_tab))) {
    results[[gene]] <- tibble(
      gene = gene,
      n_spots = nrow(df_g),
      detect_frac = as.numeric(detect_frac[gene]),
      edf = NA_real_,
      p_smooth = NA_real_,
      mean_derivative = NA_real_,
      note = "smooth_term_not_found"
    )
    next
  }
  
  # Compute the mean first derivative of the fitted smooth curve
  # across the selected distance window (here 0 to 1000 µm).
  #
  # Interpretation:
  # positive mean_derivative = expression tends to increase with distance
  # negative mean_derivative = expression tends to decrease with distance
  # near zero                = little overall directional change
  mean_deriv <- tryCatch(
    compute_mean_derivative(
      fit = fit,
      dist_var = dist_var,
      dist_min = dist_min,
      dist_max = dist_max,
      n_grid = n_grid_deriv
    ),
    error = function(e) NA_real_
  )
  
  # Store one summary row for this gene
  results[[gene]] <- tibble(
    gene = gene,
    n_spots = nrow(df_g),
    detect_frac = as.numeric(detect_frac[gene]),
    edf = as.numeric(s_tab[smooth_name, "edf"]),
    # Usually the last column of s.table is the p-value for the smooth term
    p_smooth = as.numeric(s_tab[smooth_name, ncol(s_tab)]),
    mean_derivative = mean_deriv,
    note = "ok"
  )
}

# Close progress bar after loop finishes
close(pb)

# Combine all per-gene results into one data frame
res_df <- bind_rows(results) %>%
  mutate(
    # Adjust p-values across all genes using Benjamini-Hochberg FDR correction
    p_adj = p.adjust(p_smooth, method = "BH"),
    
    # Add a simple direction label based on mean derivative
    direction = case_when(
      is.na(mean_derivative) ~ NA_character_,
      mean_derivative > 0 ~ "up_with_distance",
      mean_derivative < 0 ~ "down_with_distance",
      TRUE ~ "flat"
    )
  ) %>%
  arrange(p_adj, p_smooth)



## Save GAM results

write.table(
  res_df,
  file.path(gam_outdir, "GAM_distance_dependent_genes_full_results_Post_38.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

saveRDS(
  res_df,
  file.path(gam_outdir, "GAM_distance_dependent_genes_full_results_Post_38.rds")
)

message("Saved GAM results to: ", normalizePath(gam_outdir))





## Plot ##

res_df <- read.table("Patient 38 - GAM/GAM_Patient 38 - GAM/GAM_distance_dependent_genes_full_results_Post_38.tsv",
                     header=TRUE, row.names=1, sep = '\t', check.names = FALSE)

res_df$gene <- rownames(res_df)


padj_cut <- 1e-3
label_n_each_side <- 8
caption_pathway_IFN_GAMMA <- "HALLMARK IFN GAMMA" 
caption_pathway_EMT <- "HALLMARK EMT" 


highlight_gene <- "ESR1"

HALLMARK_IFN_GAMMA <- c(
  "ADAR","APOL6","ARID5B","ARL4A","AUTS2","B2M","BANK1","BATF2","BPGM","BST2","BTG1","C1R","C1S",
  "CASP1","CASP3","CASP4","CASP7","CASP8","CCL2","CCL5","CCL7","CD274","CD38","CD40","CD69","CD74",
  "CD86","CDKN1A","CFB","CFH","CIITA","CMKLR1","CMPK2","CMTR1","CSF2RB","CXCL10","CXCL11","CXCL9",
  "DDX60","DHX58","EIF2AK2","EIF4E3","EPSTI1","FAS","FCGR1A","FGL2","FPR1","GBP4","GBP6","GCH1","GPR18",
  "GZMA","HELZ2","HERC6","HIF1A","HLA-A","HLA-B","HLA-DMA","HLA-DQA1","HLA-DRB1","HLA-G","ICAM1","IDO1",
  "IFI27","IFI30","IFI35","IFI44","IFI44L","IFIH1","IFIT1","IFIT2","IFIT3","IFITM2","IFITM3","IFNAR2",
  "IL10RA","IL15","IL15RA","IL18BP","IL2RB","IL4R","IL6","IL7","IRF1","IRF2","IRF4","IRF5","IRF7","IRF8",
  "IRF9","ISG15","ISG20","ISOC1","ITGB7","JAK2","KLRK1","LAP3","LATS2","LCP2","LGALS3BP","LY6E","LYSMD2",
  "MARCHF1","MT2A","MTHFD2","MVP","MX1","MX2","MYD88","NAMPT","NCOA3","NFKB1","NFKBIA","NLRC5","NMI","NOD1",
  "NUP93","OAS2","OAS3","OASL","OGFR","P2RY14","PARP12","PARP14","PDE4B","PELI1","PFKP","PIM1","PLA2G4A",
  "PLSCR1","PML","PNP","PNPT1","PSMA2","PSMA3","PSMB10","PSMB2","PSMB8","PSMB9","PSME1","PSME2","PTGS2",
  "PTPN1","PTPN2","PTPN6","RAPGEF6","RBCK1","RIGI","RIPK1","RIPK2","RNF213","RNF31","RSAD2","RTP4","SAMD9L",
  "SAMHD1","SECTM1","SELP","SERPING1","SLAMF7","SLC25A28","SOCS1","SOCS3","SOD2","SP110","SPPL2A","SRI","SSPN",
  "ST3GAL5","ST8SIA4","STAT1","STAT2","STAT3","STAT4","TAP1","TAPBP","TDRD7","TMT1B","TNFAIP2","TNFAIP3",
  "TNFAIP6","TNFSF10","TOR1B","TRAFD1","TRIM14","TRIM21","TRIM25","TRIM26","TXNIP","UBE2L6","UPP1","USP18",
  "VAMP5","VAMP8","VCAM1","WARS1","XAF1","XCL1","ZBP1","ZNFX1"
)

HALLMARK_EMT <- c(
  "ABI3BP","ACTA2","ADAM12","ANPEP","APLP1","AREG","BASP1","BDNF","BGN","BMP1",
  "CADM1","CALD1","CALU","CAP2","CAPG","CCN1","CCN2","CD44","CD59","CDH11",
  "CDH2","CDH6","COL11A1","COL12A1","COL16A1","COL1A1","COL1A2","COL3A1",
  "COL4A1","COL4A2","COL5A1","COL5A2","COL5A3","COL6A2","COL6A3","COL7A1",
  "COL8A2","COLGALT1","COMP","COPA","CRLF1","CTHRC1","CXCL1","CXCL12",
  "CXCL6","CXCL8","DAB2","DCN","DKK1","DPYSL3","DST","ECM1","ECM2","EDIL3",
  "EFEMP2","ELN","EMP3","ENO2","FAP","FAS","FBLN1","FBLN2","FBLN5","FBN1",
  "FBN2","FERMT2","FGF2","FLNA","FMOD","FN1","FOXC2","FSTL1","FSTL3",
  "FUCA1","FZD8","GADD45A","GADD45B","GAS1","GEM","GJA1","GLIPR1","GPC1",
  "GPX7","GREM1","HTRA1","ID2","IGFBP2","IGFBP3","IGFBP4","IL15","IL32",
  "IL6","INHBA","ITGA2","ITGA5","ITGAV","ITGB1","ITGB3","ITGB5","JUN",
  "LAMA1","LAMA2","LAMA3","LAMC1","LAMC2","LGALS1","LOX","LOXL1","LOXL2",
  "LRP1","LRRC15","LUM","MAGEE1","MATN2","MATN3","MCM7","MEST","MFAP5",
  "MGP","MMP1","MMP14","MMP2","MMP3","MSX1","MXRA5","MYL9","MYLK","NID2",
  "NNMT","NOTCH2","NT5E","NTM","OXTR","P3H1","PCOLCE","PCOLCE2","PDGFRB",
  "PDLIM4","PFN2","PLAUR","PLOD1","PLOD2","PLOD3","PMEPA1","PMP22","POSTN",
  "PPIB","PRRX1","PRSS2","PTHLH","PTX3","PVR","QSOX1","RGS4","RHOB","SAT1",
  "SCG2","SDC1","SDC4","SERPINE1","SERPINE2","SERPINH1","SFRP1","SFRP4",
  "SGCB","SGCD","SGCG","SLC6A8","SLIT2","SLIT3","SNAI2","SNTB1","SPARC",
  "SPOCK1","SPP1","TAGLN","TFPI2","TGFB1","TGFBI","TGFBR3","TGM2","THBS1",
  "THBS2","THY1","TIMP1","TIMP3","TNC","TNFAIP3","TNFRSF11B","TNFRSF12A",
  "TPM1","TPM2","TPM4","VCAM1","VCAN","VEGFA","VEGFC","VIM","WIPF1","WNT5A"
)



# Filtering and ranking

res_sig <- res_df %>%
  filter(!is.na(mean_derivative), !is.na(p_adj), note == "ok", p_adj < padj_cut) %>%
  arrange(mean_derivative) %>%
  mutate(
    rank = row_number(),
    in_label = gene %in% genes_label,
    IFN_related    = gene %in% HALLMARK_IFN_GAMMA,
    EMT_related = gene %in% HALLMARK_EMT,
    sig_group = case_when(
      IFN_related    ~ "IFNγ",
      EMT_related ~ "EMT",
      TRUE           ~ "Other"
    ),
    sig_group = factor(sig_group, levels = c("Other","IFNγ","EMT"))
  )

# To label all the genes
lab_df <- res_sig %>%
  filter(in_label) %>%
  distinct(gene, .keep_all = TRUE)


# Counts for the legend
n_pos <- sum(res_sig$mean_derivative > 0, na.rm = TRUE)
n_neg <- sum(res_sig$mean_derivative < 0, na.rm = TRUE)


idx_IFN <- res_sig$gene %in% HALLMARK_IFN_GAMMA
n_pos_IFN <- sum(res_sig$mean_derivative[idx_IFN] > 0, na.rm = TRUE)

idx_EMT <- res_sig$gene %in% HALLMARK_EMT
n_pos_EMT <- sum(res_sig$mean_derivative[idx_EMT] > 0, na.rm = TRUE)


sig_cols <- c(
  "IFNγ"         = "#63B4A4",
  "EMT" = "#EB4CF7"
)



p <- ggplot(res_sig, aes(rank, mean_derivative)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.6) +
  
  geom_point(
    data = subset(res_sig, sig_group == "Other"),
    color = "grey90",
    alpha = 0.08,
    size  = 0.8
  ) +
  
  geom_point(
    data = subset(res_sig, sig_group != "Other"),
    aes(fill = sig_group),
    alpha  = 1,
    size   = 3.5,
    shape  = 21,
    color  = "white",
    stroke = 0.35
  ) +
  scale_fill_manual(values = sig_cols) +
  guides(fill = "none") +
  
  scale_y_continuous(labels = label_scientific(digits = 1)) +
  labs(
    title = "Patient 38 - Surgery \nGAM derivatives across genes",
    subtitle = sprintf("(adj. P < %.3g)", padj_cut),
    x = "Genes ranked by mean derivative (low -> high)",
    y = "Mean Derivative"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, face = "italic"),
    axis.title = element_text(face = "bold")
  ) 


rich_label <- sprintf(
  "<span style='color:#63B4A4'>● %s &gt; 0, %d genes</span><br>
   <span style='color:#EB4CF7'>● %s &gt; 0, %d genes</span><br>
   Mean derivative &gt; 0, %d genes<br>
   Mean derivative &lt; 0, %d genes",
  caption_pathway_IFN_GAMMA,
  n_pos_IFN,
  caption_pathway_EMT,
  n_pos_EMT,
  n_pos, n_neg
)

p <- p +
  ggtext::geom_richtext(
    aes(x = Inf, y = -Inf, label = rich_label),
    inherit.aes = FALSE,
    hjust = 1.02, vjust = -0.1,
    size = 4.5,
    fill = NA, label.color = NA
  )



if (!is.na(highlight_gene) && highlight_gene %in% res_sig$gene) {
  
  dx <- -0.03 * max(res_sig$rank, na.rm = TRUE)
  dy <- 0.06 * diff(range(res_sig$mean_derivative, na.rm = TRUE))
  
  
  p <- p +
    geom_point(
      data = res_sig %>% dplyr::filter(gene == highlight_gene),
      shape = 17, color = "black", size = 4
    ) +
    geom_text_repel(
      data = res_sig %>% dplyr::filter(gene == highlight_gene),
      aes(label = gene),
      color = "#7F0000",
      size = 5,
      box.padding = 0.35,
      point.padding = 0.25,
      max.overlaps = Inf,
      nudge_x = dx,
      nudge_y = dy,
      force = 2,
      min.segment.length = 0
    )
}

png("Patient 38 - GAM/Patient_38_GAM_Surgery_tumor_signatures.png", width = 4500, height = 4500, res = 450)
p
dev.off()
