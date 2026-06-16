setwd("PATH")

library(ggplot2)
library(tidyverse)
library(pheatmap)
library(dplyr)
library(stringr)
library(readxl)
library(ggsci)


# Previously, we have merged all patients interactions in a single file
df <- read_excel("All_Patients_CellChat_Interactions_50μm.xls")
df <- as.data.frame(df)

# To add denosumab time-point
df_mod <- df |> 
  mutate(Condition = ifelse(grepl("Pre", source), "Pre", "Post"))

# To delete pre and post endings from cell type annotations
df_mod <- df_mod %>%
  mutate(source = str_remove(source, "(_Pre|_Post)$"))

df_mod <- df_mod %>%
  mutate(target = str_remove(target, "(_Pre|_Post)$"))


# To add TotalProb column by summing all the interaction probabilities by pathway
df_long_adj <- df_mod %>%
  mutate(prob = as.numeric(prob),
         CellType = target) %>%
  group_by(Patient, source, CellType, Condition, interaction_name, pathway_name) %>%
  summarise(TotalProb = sum(prob, na.rm = TRUE), .groups = "drop")

write.table(df_long_adj, file="All_Patients_CellChat_TotalProb_by_Pathway_50μm.txt", sep = '\t')




#### Defining the color for each pathway #####
pathway_colors <- c(
  APP="#e41a1c", COLLAGEN="#ff7f0e", LAMININ="#4daf4a", PECAM1="#17becf",
  SPP1="#9467bd", FN1="#8c564b", GAP="#f781bf", THBS="#999999",
  Netrin="#66c2a5", CDH="#98df8a", MIF="#8da0cb", SEMA3="#e78ac3",
  COMPLEMENT="#a6d854", CXCL="#ffd92f", ADIPONECTIN="#e5c494",
  CCL="#b3b3b3", IL16="#1b9e77", SELL="#d95f02", GRN="#a9961d",
  PTPR="#e7298a", PECAM2="#66a61e", TENASCIN="#e6ab02", ADGRG="#a6761d",
  CD96="#666666", CD99="#8dd3c7", IGFBP="#ffffb3", CDH5="#bebada",
  ESAM="#fb8072", NOTCH="#80b1d3", ANGPTL="#fdb462", CD46="#b3de69",
  EPHB="#fccde5", GAS="#d9d9d9", JAM="#bc80bd", MK="#ccebc5",
  PDGF="#ffed6f", VEGF="#1f78b4", PERIOSTIN="#88a02c", PTPRM="#fb9a99",
  ADGRA="#a6cee3", SEMA4="#b2df8a", SEMA6="#fdbf6f", TGFb="#cab2d6",
  ANNEXIN="#ffff99", PLAU="#6a3d9a", THY1="#ff7f7f", VCAM="#b15928",
  MMP="#e31a1c", VWF="#6a51a3", HSPG="black", AGRN="#e7298a",
  EPHA="#3182bd", ACTIVIN="#09d304", ADGRB="#756bb1", RELN="#636363",
  CSF="#43a2ca", PARs="#fb6a4a", RA="#238b45", TWEAK="#a1d99b",
  Prostaglandin="#bcbd22", ADGRE="#9467bd", PTN="#17becf",
  DESMOSOME="#9edae5", CysLTs="#8c564b", MPZ="#c49c94", CypA="#e377c2",
  SEMA5="#7f7f7f", CDH1="#aec7e8", Desmosterol="#ffbb78", KIT="#98df8a",
  ICAM="#ff9896", NECTIN="#c5b0d5", CEACAM="#f0f0f0", OCLN="#dbdb8d",
  WNT="#9edae5", PCDH="#17becf", GALECTIN="#1f77b4", BAFF="#ff7f0e",
  CD6="#2ca02c", CD39="#fb8072", CALCR="#9467bd", APRIL="#8c564b",
  CADM="#e377c2", BMP="#7f7f7f", CLDN="#bcbd22", Cholesterol="#79becf",
  DHT="#aec7e8", EGF="#ffbb78", VISFATIN="#fc8d62", ncWNT="#ff9896",
  CSPG4="#c5b0d5", "Other pathways"="#e5c494"
)



###### Tumor-tumor interactions ######


df_tumor <- df_long_adj[df_long_adj$source %in% "tumor",  ]


df_tumor <- df_tumor[df_tumor$CellType %in% "tumor",  ]


df_tumor <- df_tumor %>%
  mutate(pathway_name = if_else(TotalProb < 4e-4, "Other pathways", pathway_name))


png("All_patients_tumour-tumour_50µm_interactions.png", width = 9000, height = 4500, res = 800)

df_tumor %>%
  mutate(
    Condition = factor(Condition, levels = c("Biopsy", "Surgery")),
    Patient = factor(Patient, levels = sort(unique(Patient))),
    pathway_name = factor(pathway_name, levels = names(pathway_colors))
  ) %>%
  ggplot(aes(x = Condition,
             y = TotalProb,
             fill = pathway_name)) +
  geom_col(position = "stack") +
  facet_grid(~ Patient, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pathway_colors) +
  labs(
    x = "Condition",
    y = "Sum of interaction probabilities",
    title = "All patients - Tumour-Tumour interactions by pathway (50µm)",
    fill = "Pathway"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    strip.text = element_text(face = "bold", size = 14)
  )
dev.off()




