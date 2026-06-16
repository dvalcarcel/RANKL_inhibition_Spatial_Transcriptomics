setwd("PATH")


library(dplyr)
library(purrr)
library(tibble)

### To analyse signature scores and perform downsampling analysis ###

median_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(list(n=0, q1=NA_real_, med=NA_real_, q3=NA_real_, iqr=NA_real_))
  q <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  list(n=length(x), q1=q[1], med=q[2], q3=q[3], iqr=q[3]-q[1])
}

rank_biserial <- function(W, n1, n2) {
  U <- as.numeric(W) - n1*(n1+1)/2
  (2*U)/(n1*n2) - 1
}

downsample_mw <- function(x_pre, x_post, B = 500, seed = 123) {
  x_pre  <- x_pre[is.finite(x_pre)]
  x_post <- x_post[is.finite(x_post)]
  
  n_pre  <- length(x_pre)
  n_post <- length(x_post)
  n_ds   <- min(n_pre, n_post)
  
  if (n_ds < 3) {
    return(list(
      pre_n = n_pre, post_n = n_post, n_ds = n_ds,
      p_median = NA_real_, p_q25 = NA_real_, p_q75 = NA_real_, prop_005 = NA_real_,
      delta_median_ds = NA_real_,
      effect_rbc_median = NA_real_
    ))
  }
  
  set.seed(seed)
  
  pvals  <- numeric(B)
  deltas <- numeric(B)
  rbcs   <- numeric(B)
  
  for (b in seq_len(B)) {
    s_pre  <- sample(x_pre,  n_ds, replace = FALSE)
    s_post <- sample(x_post, n_ds, replace = FALSE)
    
    wt <- suppressWarnings(stats::wilcox.test(s_pre, s_post, exact = FALSE))
    pvals[b] <- wt$p.value
    deltas[b] <- stats::median(s_post) - stats::median(s_pre)
    rbcs[b] <- rank_biserial(wt$statistic, n_ds, n_ds)
  }
  
  list(
    pre_n = n_pre, post_n = n_post, n_ds = n_ds,
    p_median = stats::median(pvals, na.rm = TRUE),
    p_q25 = stats::quantile(pvals, 0.25, na.rm = TRUE),
    p_q75 = stats::quantile(pvals, 0.75, na.rm = TRUE),
    prop_005 = mean(pvals < 0.05, na.rm = TRUE),
    delta_median_ds = stats::median(deltas, na.rm = TRUE),
    effect_rbc_median = stats::median(rbcs, na.rm = TRUE)
  )
}


# Main function

signature_analysis <- function(
    seu,
    signatures,
    celltype_col = "final_annotation",
    time_col = "Denosumab",
    pre_label = "Pre",
    post_label = "Post",
    B = 500,
    seed = 123,
    min_cells = 3
) {
  
  if (!inherits(seu, "Seurat")) stop("The object doesn't appear to be a Seurat object.")
  md <- seu@meta.data
  
  needed_cols <- c(celltype_col, time_col)
  missing_needed <- setdiff(needed_cols, colnames(md))
  if (length(missing_needed) > 0) {
    stop("These columns are missing: ", paste(missing_needed, collapse = ", "))
  }
  
  sig_vec <- unique(signatures)
  sig_present <- sig_vec[sig_vec %in% colnames(md)]
  sig_missing <- setdiff(sig_vec, sig_present)
  
  if (length(sig_present) == 0) stop("There is no signature in the object's metadata.")
  if (length(sig_missing) > 0) message("Missing signatures (ignored): ", paste(sig_missing, collapse = ", "))
  
  df0 <- md %>%
    mutate(
      cell_type = .data[[celltype_col]],
      time = .data[[time_col]]
    ) %>%
    filter(time %in% c(pre_label, post_label)) %>%
    select(cell_type, time, all_of(sig_present))
  
  cell_types <- sort(unique(df0$cell_type))

  
  df_ct <- map_dfr(cell_types, function(ct) {
    dct <- df0 %>% filter(cell_type == ct)
    
    map_dfr(sig_present, function(sig) {
      x_pre  <- dct %>% filter(time == pre_label)  %>% pull(!!sym(sig))
      x_post <- dct %>% filter(time == post_label) %>% pull(!!sym(sig))
      
      s_pre  <- median_iqr(x_pre)
      s_post <- median_iqr(x_post)
      
      if (s_pre$n < min_cells || s_post$n < min_cells) {
        return(tibble(
          cell_type = ct, signature = sig,
          pre_n = s_pre$n,
          pre_summary = sprintf("%.4f [%.4f–%.4f] (IQR=%.4f)", s_pre$med, s_pre$q1, s_pre$q3, s_pre$iqr),
          post_n = s_post$n,
          post_summary = sprintf("%.4f [%.4f–%.4f] (IQR=%.4f)", s_post$med, s_post$q1, s_post$q3, s_post$iqr),
          delta_median = s_post$med - s_pre$med,
          effect_rbc = NA_real_,
          p_value = NA_real_
        ))
      }
      
      wt <- suppressWarnings(stats::wilcox.test(x_pre, x_post, exact = FALSE))
      rbc <- rank_biserial(wt$statistic, s_pre$n, s_post$n)
      
      tibble(
        cell_type = ct, signature = sig,
        pre_n = s_pre$n,
        pre_summary = sprintf("%.4f [%.4f–%.4f] (IQR=%.4f)", s_pre$med, s_pre$q1, s_pre$q3, s_pre$iqr),
        post_n = s_post$n,
        post_summary = sprintf("%.4f [%.4f–%.4f] (IQR=%.4f)", s_post$med, s_post$q1, s_post$q3, s_post$iqr),
        delta_median = s_post$med - s_pre$med,
        effect_rbc = rbc,
        p_value = wt$p.value
      )
    })
  }) %>%
    group_by(cell_type) %>%
    mutate(p_adj_fdr_within_ct = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(p_adj_fdr_global = p.adjust(p_value, method = "BH"))
  
  
  df_ct_downsample <- map_dfr(cell_types, function(ct) {
    dct <- df0 %>% filter(cell_type == ct)
    
    map_dfr(sig_present, function(sig) {
      x_pre  <- dct %>% filter(time == pre_label)  %>% pull(!!sym(sig))
      x_post <- dct %>% filter(time == post_label) %>% pull(!!sym(sig))
      
      seed_i <- seed + as.integer(factor(paste(ct, sig)))
      
      ds <- downsample_mw(x_pre, x_post, B = B, seed = seed_i)
      
      tibble(
        cell_type = ct,
        signature = sig,
        pre_n = ds$pre_n,
        post_n = ds$post_n,
        n_ds = ds$n_ds,
        p_ds_median = ds$p_median,
        p_ds_q25 = ds$p_q25,
        p_ds_q75 = ds$p_q75,
        p_ds_prop_0.05 = ds$prop_005,
        delta_median_ds = ds$delta_median_ds,
        effect_rbc_ds_median = ds$effect_rbc_median
      )
    })
  }) %>%
    group_by(cell_type) %>%
    mutate(p_ds_adj_fdr_within_ct = p.adjust(p_ds_median, method = "BH")) %>%
    ungroup() %>%
    mutate(p_ds_adj_fdr_global = p.adjust(p_ds_median, method = "BH")) %>%
    mutate(
      ds_summary = ifelse(
        is.na(p_ds_median),
        NA_character_,
        sprintf("p_med=%s [%s–%s]; prop(p<0.05)=%.2f; n_ds=%d; B=%d",
                format.pval(p_ds_median, digits = 3, eps = 1e-300),
                format.pval(p_ds_q25,    digits = 3, eps = 1e-300),
                format.pval(p_ds_q75,    digits = 3, eps = 1e-300),
                p_ds_prop_0.05, n_ds, B)
      )
    )
  

  df_combined <- df_ct %>%
    left_join(
      df_ct_downsample %>%
        select(cell_type, signature,
               n_ds,
               delta_median_ds, effect_rbc_ds_median,
               p_ds_median, p_ds_adj_fdr_within_ct, p_ds_adj_fdr_global),
      by = c("cell_type", "signature")
    ) %>%
    rename(
      `Cell type` = cell_type,
      `Signature` = signature,
      `Biopsy (median [Q1–Q3] IQR)` = pre_summary,
      `Surgery (median [Q1–Q3] IQR)` = post_summary,
      `Δ (Surgery-Biopsy)` = delta_median,
      `Effect size (rbc)` = effect_rbc,
      `p value` = p_value,
      `FDR (within cell type)` = p_adj_fdr_within_ct,
      `FDR (global)` = p_adj_fdr_global,
      `N downsampling` = n_ds,
      `Δ DS (Surgery-Biopsy)` = delta_median_ds,
      `Effect size DS (rbc)` = effect_rbc_ds_median,
      `p DS (median)` = p_ds_median,
      `FDR DS (within cell type)` = p_ds_adj_fdr_within_ct,
      `FDR DS (global)` = p_ds_adj_fdr_global
    ) %>%
    arrange(`Cell type`, `FDR DS (within cell type)`, `FDR (within cell type)`)
  
  list(
    df_ct = df_ct_pretty,
    df_ct_downsample = df_ct_downsample,
    df_combined = df_combined
  )
}


load("/Users/dvalcarcel/Documents/TFM/final_objects/Patient_38_RCTD.rda")


signatures_ssGSEA <- c(
  "ONLY_EXP_REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING",
  "ONLY_EXP_HALLMARK_HYPOXIA",
  "ONLY_EXP_KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "ONLY_EXP_HALLMARK_IL2_STAT5_SIGNALING",
  "ONLY_EXP_HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "ONLY_EXP_HALLMARK_INFLAMMATORY_RESPONSE",
  "ONLY_EXP_REACTOME_CD28_CO_STIMULATION",
  "ONLY_EXP_REACTOME_TNF_SIGNALING",
  "ONLY_EXP_HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "ONLY_EXP_REACTOME_INTERLEUKIN_10_SIGNALING",
  "ONLY_EXP_REACTOME_INTERLEUKIN_17_SIGNALING",
  "ONLY_EXP_KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "ONLY_EXP_REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_CYCLE_GENES",
  "ONLY_EXP_REACTOME_CONSTITUTIVE_SIGNALING_BY_AKT1_E17K_IN_CANCER",
  "ONLY_EXP_HALLMARK_E2F_TARGETS",
  "ONLY_EXP_HALLMARK_TGF_BETA_SIGNALING",
  "ONLY_EXP_REACTOME_RESP_ELECTRON_TRANSPORT_UCP",
  "ONLY_EXP_REACTOME_NICOTINATE_METABOLISM",
  "ONLY_EXP_KEGG_NICOTINATE_AND_NICOTINAMIDE_METABOLISM",
  "ONLY_EXP_REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  "ONLY_EXP_REACTOME_MITOCHONDRIAL_BIOGENESIS",
  "ONLY_EXP_REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING",
  "ONLY_EXP_KEGG_ADIPOCYTOKINE_SIGNALING_PATHWAY",
  "ONLY_EXP_HALLMARK_FATTY_ACID_METABOLISM",
  "ONLY_EXP_HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "ONLY_EXP_HALLMARK_MTORC1_SIGNALING"
)


res_pt38 <- signature_analysis(
  seu = Patient_38_RCTD,
  signatures = signatures_ssGSEA,
  celltype_col = "final_annotation",
  time_col = "Denosumab",
  B = 500,
  seed = 123
)

df_combined_pt38 <- res_pt38$df_combined
write_xlsx(df_combined_pt38, "Patient_38_df_signatures.xlsx")

