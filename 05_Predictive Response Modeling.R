# ==============================================================================
# Script 05: Predictive Response Modeling & Dysfunction Correlation Analysis
# ==============================================================================
library(Seurat)
library(CellChat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)

set.seed(42)

# Create output directories
dir.create("figures/predictive_modeling", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Load Seurat Object & CellChat Target Interactions
# ------------------------------------------------------------------------------
seurat_path <- "objects/immune_seurat_dysfunction_scored.rds"
cellchat_csv <- "results/target_checkpoint_interactions_summary.csv"

if (!file.exists(seurat_path)) {
  stop("Missing 'objects/immune_seurat_dysfunction_scored.rds'. Complete Task 3 first.")
}
if (!file.exists(cellchat_csv)) {
  stop("Missing 'results/target_checkpoint_interactions_summary.csv'. Complete Task 4 first.")
}

message("Loading Seurat object and CellChat target checkpoint data...")
seurat_obj <- readRDS(seurat_path)
target_lr   <- read.csv(cellchat_csv)

# ------------------------------------------------------------------------------
# 2. Extract Sample-Level T/NK Dysfunction & Cytotoxicity Scores
# ------------------------------------------------------------------------------
# Extract metadata containing dysfunction/cytotoxicity module scores
meta_df <- seurat_obj@meta.data

# Identify dysfunction/cytotoxicity score columns
score_cols <- grep("dysfunction|cytotoxicity|exhaustion|score", colnames(meta_df), ignore.case = TRUE, value = TRUE)
message("Found scoring metadata columns: ", paste(score_cols, collapse = ", "))

# Calculate mean dysfunction scores per patient/sample & treatment condition
sample_scores <- meta_df %>%
  group_by(patient_id, treatment_status) %>%
  summarise(
    across(any_of(score_cols), \(x) mean(x, na.rm = TRUE)),
    cell_count = n(),
    .groups = "drop"
  )

# Normalize condition names for joining (e.g., Pre-Chemo vs Pretreatment)
sample_scores <- sample_scores %>%
  mutate(Condition_Clean = ifelse(grepl("pre", treatment_status, ignore.case = TRUE), "Pre-Chemo", "Post-Chemo"))

# ------------------------------------------------------------------------------
# 3. Aggregate Checkpoint LR Interaction Probabilities by Condition & Sender Subcluster
# ------------------------------------------------------------------------------
# Process ligand-receptor probabilities across macrophage subclusters
lr_summary <- target_lr %>%
  group_by(Condition, source, target, ligand, receptor) %>%
  summarise(
    mean_prob = mean(prob, na.rm = TRUE),
    max_prob  = max(prob, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    Condition_Clean = ifelse(grepl("pre", Condition, ignore.case = TRUE), "Pre-Chemo", "Post-Chemo"),
    Axis = paste(ligand, receptor, sep = " -> ")
  )

# Write aggregated LR summary table
write.csv(lr_summary, "results/aggregated_checkpoint_probabilities.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Statistical Correlation: LR Probabilities vs. Cellular Dysfunction
# ------------------------------------------------------------------------------
DefaultAssay(seurat_obj) <- "RNA"

# FIX FOR SEURAT v5: Join layers and normalize to populate 'data' layer
message("Joining layers and normalizing expression data for per-cell extraction...")
seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
seurat_obj <- NormalizeData(seurat_obj, assay = "RNA", verbose = FALSE)

# Extract expression matrix and clean gene names (e.g. 'ENSG...--NECTIN2' -> 'NECTIN2')
expr_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
clean_genes <- sub(".*--", "", rownames(expr_mat))
rownames(expr_mat) <- clean_genes

# Safely extract gene expression vectors
nectin2_idx <- match("NECTIN2", clean_genes)
cd274_idx   <- match("CD274", clean_genes)

seurat_obj$Axis_NECTIN2_TIGIT <- if (!is.na(nectin2_idx)) expr_mat[nectin2_idx, ] else 0
seurat_obj$Axis_CD274_PDCD1   <- if (!is.na(cd274_idx))   expr_mat[cd274_idx, ] else 0

# Fetch cell-level metadata and gene expression
cell_level_data <- FetchData(
  seurat_obj, 
  vars = c("patient_id", "treatment_status", "cell_type_detailed", 
           "Axis_NECTIN2_TIGIT", "Axis_CD274_PDCD1", score_cols)
)

# Subset to Myeloid senders and T/NK receivers
myeloid_tnk_data <- cell_level_data %>%
  filter(grepl("Macrophage|Monocyte|T_Cells|NK_Cells", cell_type_detailed))

# Calculate sample-level ligand expression vs T-cell dysfunction
sample_correlation_df <- myeloid_tnk_data %>%
  group_by(patient_id, treatment_status) %>%
  summarise(
    mean_NECTIN2 = mean(Axis_NECTIN2_TIGIT, na.rm = TRUE),
    mean_CD274   = mean(Axis_CD274_PDCD1, na.rm = TRUE),
    across(any_of(score_cols), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )

write.csv(sample_correlation_df, "results/sample_level_checkpoint_correlations.csv",
          row.names = FALSE)

# ------------------------------------------------------------------------------
# 2. Extract Biological Dysfunction / Cytotoxicity Scores
# ------------------------------------------------------------------------------
meta_df <- seurat_obj@meta.data

# Target biological module score keywords directly (excluding generic QC/doublet metrics)
all_cols <- colnames(meta_df)

# Filter for Task 3 module scores (e.g., dysfunction, cytotoxicity, exhaustion, UCell)
score_cols <- all_cols[grepl("dysfunction|cytotoxicity|exhaustion|UCell|Module", all_cols, ignore.case = TRUE)]

# Exclude doublet, read count, feature count, and mitochondrial metadata
qc_blacklist <- "scDblFinder|doublet|nCount|nFeature|percent|mito|ribo"
score_cols   <- score_cols[!grepl(qc_blacklist, score_cols, ignore.case = TRUE)]

# Fallback check: If no specific name matched, list non-QC numeric columns
if (length(score_cols) == 0) {
  message("Available metadata columns in your object:")
  print(all_cols)
  stop("Could not automatically identify Task 3 dysfunction score column. Please check the column names printed above.")
}

primary_score_col <- score_cols[1]
message("Successfully selected target score column: '", primary_score_col, "'")

# Aggregate sample-level dysfunction scores
sample_scores <- meta_df %>%
  group_by(patient_id, treatment_status) %>%
  summarise(
    across(all_of(score_cols), \(x) mean(x, na.rm = TRUE)),
    cell_count = n(),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 5. Visualizations: Scatter Plots & Statistical Correlation
# ------------------------------------------------------------------------------
if (!exists("primary_score_col") || is.null(primary_score_col)) {
  stop("`primary_score_col` is missing. Ensure Section 2 ran successfully.")
}

message("Generating scatter plots using Y-axis metric: ", primary_score_col)

# Scatter Plot A: NECTIN2 vs Biological Score
p1 <- ggplot(sample_correlation_df, aes(x = mean_NECTIN2, y = .data[[primary_score_col]], color = treatment_status)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw(base_size = 12) +
  labs(
    title = "NECTIN2 Expression vs T/NK Dysfunction",
    x = "Mean NECTIN2 Expression (Myeloid)",
    y = paste("Mean Score:", primary_score_col),
    color = "Treatment"
  )

# Scatter Plot B: CD274 (PD-L1) vs Biological Score
p2 <- ggplot(sample_correlation_df, aes(x = mean_CD274, y = .data[[primary_score_col]], color = treatment_status)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw(base_size = 12) +
  labs(
    title = "CD274 (PD-L1) Expression vs T/NK Dysfunction",
    x = "Mean CD274 Expression (Myeloid)",
    y = paste("Mean Score:", primary_score_col),
    color = "Treatment"
  )

# Save Scatter Plots
pdf("figures/predictive_modeling/05_checkpoint_vs_dysfunction_scatter.pdf", width = 12, height = 5)
print(p1 + p2)
dev.off()

png("figures/predictive_modeling/05_checkpoint_vs_dysfunction_scatter.png", width = 12, height = 5, units = "in", res = 300)
print(p1 + p2)
dev.off()

# ------------------------------------------------------------------------------
# 6. Granular Subcluster Signaling Contribution Barplot
# ------------------------------------------------------------------------------
p3 <- ggplot(lr_summary, aes(x = source, y = mean_prob, fill = Condition_Clean)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Axis, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Granular Myeloid Subcluster Checkpoint Signaling Strength",
    x = "Macrophage / Monocyte Subcluster",
    y = "Mean Interaction Probability",
    fill = "Condition"
  )

pdf("figures/predictive_modeling/05_granular_signaling_strength.pdf", width = 11, height = 6)
print(p3)
dev.off()

png("figures/predictive_modeling/05_granular_signaling_strength.png", width = 11, height = 6, units = "in", res = 300)
print(p3)
dev.off()

# ------------------------------------------------------------------------------
# 7. Select Target Exhaustion TFs & Co-Inhibitory Markers
# ------------------------------------------------------------------------------
target_markers <- c("TOX", "EOMES", "NFATC1", "LAG3", "HAVCR2", "TIGIT", "PDCD1")
available_markers <- target_markers[target_markers %in% clean_genes]

message("Found T-cell exhaustion markers in matrix: ", paste(available_markers, collapse = ", "))

# ------------------------------------------------------------------------------
# 8. Extract Myeloid Ligand Predictors & T-Cell Outcome Genes
# ------------------------------------------------------------------------------
meta_df <- FetchData(seurat_obj, vars = c("patient_id", "treatment_status", "cell_type_detailed"))

# Myeloid cells (Monocytes & Macrophages)
myeloid_ids <- rownames(meta_df[grepl("Monocyte|Macrophage", meta_df$cell_type_detailed, ignore.case = TRUE), ])

myeloid_ligands_df <- data.frame(
  cell_id = myeloid_ids,
  patient_id = seurat_obj$patient_id[myeloid_ids],
  treatment_status = seurat_obj$treatment_status[myeloid_ids],
  NECTIN2 = expr_mat["NECTIN2", myeloid_ids],
  CD274   = expr_mat["CD274", myeloid_ids]
) %>%
  group_by(patient_id, treatment_status) %>%
  summarise(
    myeloid_NECTIN2 = mean(NECTIN2, na.rm = TRUE),
    myeloid_CD274   = mean(CD274, na.rm = TRUE),
    .groups = "drop"
  )

# Receiving T cells
t_cell_ids <- rownames(meta_df[meta_df$cell_type_detailed == "T_Cells", ])

t_cell_marker_mat <- t(expr_mat[available_markers, t_cell_ids, drop = FALSE])
t_cell_markers_df <- data.frame(
  cell_id = t_cell_ids,
  patient_id = seurat_obj$patient_id[t_cell_ids],
  t_cell_marker_mat
) %>%
  group_by(patient_id) %>%
  summarise(
    across(all_of(available_markers), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )

# Merge Myeloid Predictors with T-Cell Marker Outcomes
tf_correlation_df <- inner_join(myeloid_ligands_df, t_cell_markers_df, by = "patient_id")
write.csv(tf_correlation_df, "results/downstream_tf_correlations.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 9. Correlation Matrix Heatmap: Myeloid Ligands vs. T-Cell TF Program
# ------------------------------------------------------------------------------
corr_matrix <- cor(
  tf_correlation_df %>% select(myeloid_NECTIN2, myeloid_CD274),
  tf_correlation_df %>% select(all_of(available_markers)),
  method = "spearman"
)

corr_melted <- melt(corr_matrix)
colnames(corr_melted) <- c("Myeloid_Ligand", "T_Cell_Marker", "Spearman_R")

p_heatmap <- ggplot(corr_melted, aes(x = T_Cell_Marker, y = Myeloid_Ligand, fill = Spearman_R)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Spearman_R)), color = "black", size = 4) +
  scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", midpoint = 0, limit = c(-1, 1)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Myeloid Checkpoint Ligands vs. T-Cell Exhaustion Program",
    x = "T-Cell Transcription Factors & Inhibitory Checkpoints",
    y = "Myeloid Expression",
    fill = "Spearman R"
  )

# ------------------------------------------------------------------------------
# 10. Scatter Plots: NECTIN2 vs Master Regulator TOX and Checkpoint LAG3
# ------------------------------------------------------------------------------
p_tox <- ggplot(tf_correlation_df, aes(x = myeloid_NECTIN2, y = TOX, color = treatment_status)) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw(base_size = 12) +
  labs(
    title = "Myeloid NECTIN2 vs T-Cell TOX Expression",
    x = "Mean NECTIN2 Expression (Myeloid)",
    y = "Mean TOX Expression (T Cells)",
    color = "Treatment"
  )

p_lag3 <- ggplot(tf_correlation_df, aes(x = myeloid_NECTIN2, y = LAG3, color = treatment_status)) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw(base_size = 12) +
  labs(
    title = "Myeloid NECTIN2 vs T-Cell LAG3 Expression",
    x = "Mean NECTIN2 Expression (Myeloid)",
    y = "Mean LAG3 Expression (T Cells)",
    color = "Treatment"
  )

# ------------------------------------------------------------------------------
# 11. Save Visualizations
# ------------------------------------------------------------------------------
ggsave("figures/predictive_modeling/06_tf_correlation_heatmap.png", p_heatmap, width = 8, height = 4, dpi = 300)
ggsave("figures/predictive_modeling/06_nectin2_vs_tox.png", p_tox, width = 6.5, height = 4.5, dpi = 300)
ggsave("figures/predictive_modeling/06_nectin2_vs_lag3.png", p_lag3, width = 6.5, height = 4.5, dpi = 300)

message("\n=== Script 06 execution complete! Results exported to 'figures/predictive_modeling/' and 'results/' ===")