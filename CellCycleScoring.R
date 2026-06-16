library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)
library(tibble)

setwd("PATH")


## Patient 38 ##

load("Patient_38_RCTD.rda")

Patient_38_RCTD <- UpdateSeuratObject(Patient_38_RCTD)

Idents(Patient_38_RCTD) <- Patient_38_RCTD$final_annotation
Patient_38 <- subset(Patient_38_RCTD, idents = c("tumor", "T_cells", "B_cells"))
Patient_38 <- Patient_38[, !is.na(Patient_38$final_annotation)]
Patient_38 <- Patient_38[, !is.na(Patient_38$Denosumab)]

## Cell cycle scoring

s.genes  <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

Patient_38 <- CellCycleScoring(
  Patient_38,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE
)

# To merge S + G2M phases and get the proportions

meta <- Patient_38@meta.data %>%
  mutate(
    timepoint = factor(Denosumab, levels = c("Pre","Post")),
    Phase = factor(Phase, levels = c("G1","S","G2M")),
    Phase2 = ifelse(Phase == "G1", "G1", "S+G2M"),
    Phase2 = factor(Phase2, levels = c("G1", "S+G2M"))
  )

counts_long <- meta %>%
  count(final_annotation, timepoint, Phase, name = "n") %>%
  group_by(final_annotation, timepoint) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


# Fisher test

global_test_one_annotation <- function(df_cells_ann) {
  tab <- table(df_cells_ann$Phase2, df_cells_ann$timepoint)
  
  if (ncol(tab) < 2) {
    return(tibble(
      method = NA_character_,
      p_value = NA_real_,
      OR = NA_real_,
      pre_prop_prolif = NA_real_,
      post_prop_prolif = NA_real_,
      delta_prolif = NA_real_,
      note = "Biopsy or surgery missing"
    ))
  }
  
  ft <- fisher.test(tab)
  
  pre_prop_prolif  <- mean(df_cells_ann$Phase2[df_cells_ann$timepoint == "Pre"]  == "S+G2M")
  post_prop_prolif <- mean(df_cells_ann$Phase2[df_cells_ann$timepoint == "Post"] == "S+G2M")
  
  tibble(
    method = "Fisher's Exact Test",
    p_value = ft$p.value,
    OR = unname(ft$estimate),
    pre_prop_prolif = pre_prop_prolif,
    post_prop_prolif = post_prop_prolif,
    delta_prolif = post_prop_prolif - pre_prop_prolif,
    note = NA_character_
  )
}

global_results <- meta %>%
  group_by(final_annotation) %>%
  group_modify(~ global_test_one_annotation(.x)) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))



# Plot labels

counts_long_lbl <- counts_long %>%
  mutate(lbl = ifelse(!is.na(prop) & prop >= 0.07, as.character(n), ""))

ann_order <- counts_long_lbl %>%
  distinct(final_annotation) %>%
  pull(final_annotation) %>%
  as.character()

tp_levels <- c("Pre", "Post")

base <- seq(1, by = 3, length.out = length(ann_order))

xmap <- tidyr::expand_grid(final_annotation = ann_order, timepoint = tp_levels) %>%
  mutate(
    x = rep(base, each = 2) + ifelse(timepoint == "Post", 1, 0)
  )

ann_pos <- tibble(
  final_annotation = ann_order,
  x_mid = base + 0.5
)

counts_long_lbl2 <- counts_long_lbl %>%
  mutate(
    final_annotation = as.character(final_annotation),
    timepoint = as.character(timepoint)
  ) %>%
  left_join(xmap, by = c("final_annotation", "timepoint"))

# Plot p values

use_adj <- FALSE   

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    "NA",
    ifelse(
      p < 1e-3,
      "<0.001",
      formatC(p, format = "f", digits = 3)
    )
  )
}

pcol_global <- if (use_adj) "p_adj" else "p_value"


p_labels <- global_results %>%
  mutate(final_annotation = as.character(final_annotation)) %>%
  left_join(ann_pos, by = "final_annotation") %>%
  mutate(
    label = paste0(
      "S+G2M fraction p",
      if (use_adj) "(adj)=" else "=",
      fmt_p(.data[[pcol_global]])
    ),
    y = 1.08
  )

ann_labels_df <- ann_pos %>%
  mutate(
    y = -0.06,
    label = final_annotation
  )

# Plot

phase_cols <- c(
  "G1" = "#4C78A8",
  "S" = "#F58518",
  "G2M" = "#54A24B"
)

p1 <- ggplot(counts_long_lbl2, aes(x = x, y = prop, fill = Phase)) +
  geom_col(width = 0.9, color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = lbl),
    position = position_stack(vjust = 0.5),
    size = 8, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(values = phase_cols, breaks = c("G1","S","G2M"), drop = FALSE) +
  scale_x_continuous(
    breaks = xmap$x,
    labels = xmap$timepoint,
    expand = expansion(add = 0.6)
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  coord_cartesian(ylim = c(-0.12, 1.22), clip = "off") +
  labs(x = NULL, y = "Normalized proportion", fill = "Phase") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 14),
    plot.margin = margin(10, 10, 25, 10)
  ) +
  geom_text(
    data = ann_labels_df,
    aes(x = x_mid, y = y, label = label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 7
  ) +
  geom_text(
    data = p_labels,
    aes(x = x_mid, y = y, label = label),
    inherit.aes = FALSE,
    size = 7.0,
    lineheight = 0.95
  ) +
  ggtitle("Patient 38 – Cell Cycle Scoring") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  )

png("Patient_38_CellCycleScoring.png", width = 7000, height = 5000, res = 450)
p1
dev.off()



