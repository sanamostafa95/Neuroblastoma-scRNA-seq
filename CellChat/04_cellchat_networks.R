# ==============================================================
# 04_cellchat_networks.R
# Task 4: Cell-Cell Communication & Network Modeling
#
# Goal:
#   1. Run CellChat separately for Pre- and Post-Chemotherapy
#   2. Construct cell-cell communication networks
#   3. Compare communication between Pre and Post
#   4. Prioritize:
#        - NECTIN2 -> TIGIT
#        - CD274 (PD-L1) -> PDCD1 (PD-1)
#   5. Generate:
#        - network circle plots
#        - LR probability heatmaps
#        - target LR result tables
# ==============================================================

# --------------------------------------------------------------
# 1. Load packages
# --------------------------------------------------------------
library(Seurat)
library(CellChat)
library(dplyr)
library(ggplot2)
library(patchwork)
set.seed(42)

# --------------------------------------------------------------
# 2. Create output directories
# --------------------------------------------------------------

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("objects", showWarnings = FALSE)

# --------------------------------------------------------------
# 3. Load final object from Task 3
# --------------------------------------------------------------

obj <- readRDS("objects/immune_seurat_dysfunction_scored.rds")

DefaultAssay(obj) <- "RNA"


# --------------------------------------------------------------
# 4. Check required metadata
# --------------------------------------------------------------

required_metadata <- c(
  "cell_type",
  "treatment_status",
  "patient_id"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(obj@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    "Missing required metadata: ",
    paste(missing_metadata, collapse = ", ")
  )
}
# --------------------------------------------------------------
# 5. Clean CellChat metadata
# --------------------------------------------------------------

obj$cellchat_group <- as.character(obj$cell_type)

obj$cellchat_group[
  is.na(obj$cellchat_group) |
    obj$cellchat_group == ""
] <- "Unknown"

obj$treatment_status <- factor(
  obj$treatment_status,
  levels = c("Pre-Chemo", "Post-Chemo")
)

message("\nCell type distribution:")
print(table(obj$cellchat_group))

message("\nTreatment distribution:")
print(table(obj$treatment_status, useNA = "ifany"))

# --------------------------------------------------------------
# 6. Split object into Pre- and Post-Chemotherapy
# --------------------------------------------------------------

obj_pre <- subset(
  obj,
  subset = treatment_status == "Pre-Chemo"
)

obj_post <- subset(
  obj,
  subset = treatment_status == "Post-Chemo"
)

message("\nPre-Chemo cells: ", ncol(obj_pre))
message("Post-Chemo cells: ", ncol(obj_post))

# --------------------------------------------------------------
# 7. Create CellChat objects
# --------------------------------------------------------------
data_pre <- GetAssayData(
  obj_pre,
  assay = "RNA",
  layer = "data"
)

data_post <- GetAssayData(
  obj_post,
  assay = "RNA",
  layer = "data"
)

cellchat_pre <- createCellChat(
  object = data_pre,
  meta = obj_pre@meta.data,
  group.by = "cellchat_group"
)

cellchat_post <- createCellChat(
  object = data_post,
  meta = obj_post@meta.data,
  group.by = "cellchat_group"
)
# --------------------------------------------------------------
# 8. Set CellChat ligand-receptor database
# --------------------------------------------------------------

CellChatDB <- CellChatDB.human

cellchat_pre@DB <- CellChatDB
cellchat_post@DB <- CellChatDB


# --------------------------------------------------------------
# 9. Preprocess CellChat — Pre-Chemo
# --------------------------------------------------------------

cellchat_pre <- subsetData(cellchat_pre)

cellchat_pre <- identifyOverExpressedGenes(cellchat_pre)

cellchat_pre <- identifyOverExpressedInteractions(cellchat_pre)

cellchat_pre <- computeCommunProb(
  cellchat_pre,
  type = "truncatedMean",
  trim = 0.1
)

cellchat_pre <- filterCommunication(
  cellchat_pre,
  min.cells = 10
)

cellchat_pre <- computeCommunProbPathway(cellchat_pre)

cellchat_pre <- aggregateNet(cellchat_pre)

# --------------------------------------------------------------
# 10. Preprocess CellChat — Post-Chemo
# --------------------------------------------------------------

cellchat_post <- subsetData(cellchat_post)

cellchat_post <- identifyOverExpressedGenes(cellchat_post)

cellchat_post <- identifyOverExpressedInteractions(cellchat_post)

cellchat_post <- computeCommunProb(
  cellchat_post,
  type = "truncatedMean",
  trim = 0.1
)

cellchat_post <- filterCommunication(
  cellchat_post,
  min.cells = 10
)

cellchat_post <- computeCommunProbPathway(cellchat_post)

cellchat_post <- aggregateNet(cellchat_post)




# --------------------------------------------------------------
# 11. Save CellChat objects
# --------------------------------------------------------------

saveRDS(
  cellchat_pre,
  "objects/cellchat_pre.rds"
)

saveRDS(
  cellchat_post,
  "objects/cellchat_post.rds"
)


# --------------------------------------------------------------
# 12. Compare overall communication
# --------------------------------------------------------------

weight_pre <- sum(cellchat_pre@net$weight, na.rm = TRUE)
weight_post <- sum(cellchat_post@net$weight, na.rm = TRUE)

number_pre <- sum(cellchat_pre@net$count, na.rm = TRUE)
number_post <- sum(cellchat_post@net$count, na.rm = TRUE)

overall_summary <- data.frame(
  Condition = c("Pre-Chemo", "Post-Chemo"),
  Total_Communication_Weight = c(
    weight_pre,
    weight_post
  ),
  Number_of_Interactions = c(
    number_pre,
    number_post
  )
)

write.csv(
  overall_summary,
  "results/cellchat_overall_communication.csv",
  row.names = FALSE
)

# --------------------------------------------------------------
# 13. Extract all ligand-receptor interactions
# --------------------------------------------------------------

lr_pre <- subsetCommunication(cellchat_pre)

lr_post <- subsetCommunication(cellchat_post)

write.csv(
  lr_pre,
  "results/cellchat_LR_Pre_Chemo.csv",
  row.names = FALSE
)

write.csv(
  lr_post,
  "results/cellchat_LR_Post_Chemo.csv",
  row.names = FALSE
)



# --------------------------------------------------------------
# 14. Prioritize NECTIN2 -> TIGIT
# --------------------------------------------------------------

nectin2_pre <- lr_pre %>%
  filter(
    ligand == "NECTIN2" &
      receptor == "TIGIT"
  )

nectin2_post <- lr_post %>%
  filter(
    ligand == "NECTIN2" &
      receptor == "TIGIT"
  )

write.csv(
  nectin2_pre,
  "results/NECTIN2_TIGIT_Pre_Chemo.csv",
  row.names = FALSE
)

write.csv(
  nectin2_post,
  "results/NECTIN2_TIGIT_Post_Chemo.csv",
  row.names = FALSE
)


# --------------------------------------------------------------
# 15. Prioritize CD274 -> PDCD1
# --------------------------------------------------------------

pdl1_pre <- lr_pre %>%
  filter(
    ligand == "CD274" &
      receptor == "PDCD1"
  )

pdl1_post <- lr_post %>%
  filter(
    ligand == "CD274" &
      receptor == "PDCD1"
  )

write.csv(
  pdl1_pre,
  "results/CD274_PDCD1_Pre_Chemo.csv",
  row.names = FALSE
)

write.csv(
  pdl1_post,
  "results/CD274_PDCD1_Post_Chemo.csv",
  row.names = FALSE
)


# --------------------------------------------------------------
# 16. Save target interaction summary
# --------------------------------------------------------------

target_LR_summary <- bind_rows(
  nectin2_pre %>%
    mutate(Condition = "Pre-Chemo"),
  
  nectin2_post %>%
    mutate(Condition = "Post-Chemo"),
  
  pdl1_pre %>%
    mutate(Condition = "Pre-Chemo"),
  
  pdl1_post %>%
    mutate(Condition = "Post-Chemo")
)

write.csv(
  target_LR_summary,
  "results/target_LR_summary.csv",
  row.names = FALSE
)

# --------------------------------------------------------------
# 17. Network circle plots
# --------------------------------------------------------------

# Pdf version


pdf(
  "figures/CellChat_network_circle_Pre_Post.pdf",
  width = 12,
  height = 6
)

par(mfrow = c(1, 2))

netVisual_circle(
  cellchat_pre@net$weight,
  vertex.weight = as.numeric(table(cellchat_pre@idents)),
  weight.scale = TRUE,
  title.name = "Pre-Chemo"
)

netVisual_circle(
  cellchat_post@net$weight,
  vertex.weight = as.numeric(table(cellchat_post@idents)),
  weight.scale = TRUE,
  title.name = "Post-Chemo"
)

dev.off()

# PNG - Pre-Chemo

png(
  "figures/CellChat_network_circle_Pre.png",
  width = 7,
  height = 6,
  units = "in",
  res = 300
)

netVisual_circle(
  cellchat_pre@net$weight,
  vertex.weight = as.numeric(table(cellchat_pre@idents)),
  weight.scale = TRUE,
  title.name = "Pre-Chemo"
)

dev.off()

# PNG - Post-Chemo

png(
  "figures/CellChat_network_circle_Post.png",
  width = 7,
  height = 6,
  units = "in",
  res = 300
)

netVisual_circle(
  cellchat_post@net$weight,
  vertex.weight = as.numeric(table(cellchat_post@idents)),
  weight.scale = TRUE,
  title.name = "Post-Chemo"
)

dev.off()

# --------------------------------------------------------------
# 18. LR probability heatmaps
# --------------------------------------------------------------


# Pdf version


pdf(
  "figures/CellChat_LR_probability_heatmaps_Pre_Post.pdf",
  width = 14,
  height = 7
)

par(mfrow = c(1, 2))

netVisual_heatmap(
  cellchat_pre,
  measure = "weight",
  color.heatmap = "Reds",
  title.name = "Pre-Chemo LR Probability"
)

netVisual_heatmap(
  cellchat_post,
  measure = "weight",
  color.heatmap = "Reds",
  title.name = "Post-Chemo LR Probability"
)


dev.off()


# PNG - Pre-Chemo

png(
  "figures/CellChat_LR_probability_heatmap_Pre.png",
  width = 7,
  height = 7,
  units = "in",
  res = 300
)

netVisual_heatmap(
  cellchat_pre,
  measure = "weight",
  color.heatmap = "Reds",
  title.name = "Pre-Chemo LR Probability"
)

dev.off()


# PNG - Post-Chemo

png(
  "figures/CellChat_LR_probability_heatmap_Post.png",
  width = 7,
  height = 7,
  units = "in",
  res = 300
)

netVisual_heatmap(
  cellchat_post,
  measure = "weight",
  color.heatmap = "Reds",
  title.name = "Post-Chemo LR Probability"
)


dev.off()



# --------------------------------------------------------------
# 19. Save final CellChat analysis objects
# --------------------------------------------------------------



saveRDS(
  list(
    pre = cellchat_pre,
    post = cellchat_post,
    overall_summary = overall_summary,
    target_LR_summary = target_LR_summary
  ),
  "objects/Task4_CellChat_final.rds"
)


