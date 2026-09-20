
list.files(pattern = "\\.rds$")

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

immune_seurat <- readRDS("immune_seurat_step2_clustered.rds")

#######################################################################
# 02_clustering_annotation.R
# Member 2: Immune Sub-clustering & Paper Reproduction Lead
#
# Goal:
#   1. Isolate CD45+ (PTPRC+) immune population from the full dataset
#   2. Re-cluster at high resolution
#   3. Annotate detailed T-cell / NK / myeloid subsets
#   4. Reproduce baseline immune UMAP + marker dotplots for validation
#
# Inputs:
#   - combined_seurat.rds  (output of 01_load_merge_qc.R)
#
# Outputs:
#   - immune_seurat.rds
#   - plots/immune_umap*.pdf/png
#   - plots/immune_marker_dotplot*.pdf/png
#   - tables/immune_cluster_markers*.csv
#
# This dataset's specific quirks that this script handles automatically:
#   - Features are named "ENSEMBL_ID--SYMBOL" (e.g. "ENSG00000081237--PTPRC"),
#     not plain gene symbols.
#   - combined_seurat comes out of merge() with UNJOINED per-sample layers
#     (Seurat v5). Joining the FULL object can exhaust memory (std::bad_alloc),
#     so we extract genes layer-by-layer before subsetting, and only join
#     the much smaller immune subset.
#   - Per-sample min.cells filtering in script 01 means some genes (incl.
#     PTPRC) may be entirely absent from some samples' layers -> treated
#     as zero expression rather than causing an error.
#######################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(42)

dir.create("plots", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

## ------------------------------------------------------------------
## Helper functions
## ------------------------------------------------------------------

# Resolve a gene symbol against feature names that may be either plain
# symbols OR "ENSEMBL_ID--SYMBOL" style names. Case-insensitive.
resolve_gene <- function(symbol, obj, warn_if_missing = TRUE) {
  available <- rownames(obj)
  hit <- available[match(tolower(symbol), tolower(available))]
  if (is.na(hit)) {
    suffix_hits <- available[grepl(paste0("(^|--)", symbol, "$"), available, ignore.case = TRUE)]
    if (length(suffix_hits) > 0) hit <- suffix_hits[1]
  }
  if (is.na(hit) && warn_if_missing) {
    warning(sprintf("Gene '%s' not found (checked exact name and '--SYMBOL' suffix).", symbol))
    return(NA_character_)
  }
  hit
}

resolve_genes <- function(symbols, obj) {
  resolved <- vapply(symbols, resolve_gene, obj = obj, FUN.VALUE = character(1))
  missing <- symbols[is.na(resolved)]
  if (length(missing) > 0) {
    message("  Not found (skipped): ", paste(missing, collapse = ", "))
  }
  unname(resolved[!is.na(resolved)])
}

# Safely pull one gene's expression across a (possibly unjoined, v5) Seurat
# object without ever materializing the full joined matrix. Missing gene
# in a given layer -> treated as zero expression for that layer's cells.
get_gene_vector_safe <- function(obj, gene, assay = "RNA", layer_search = "data") {
  layers <- Layers(obj, assay = assay, search = layer_search)
  if (length(layers) == 0) stop("No layers found matching '", layer_search, "' in assay '", assay, "'")
  pieces <- lapply(layers, function(ly) {
    ld <- LayerData(obj, assay = assay, layer = ly)
    if (gene %in% rownames(ld)) {
      ld[gene, ]
    } else {
      setNames(rep(0, ncol(ld)), colnames(ld))
    }
  })
  vec <- unlist(pieces)
  vec[is.na(vec)] <- 0
  vec
}

# JoinLayers with a retry after garbage collection if memory allocation fails.
safe_join_layers <- function(obj, assay = "RNA") {
  tryCatch({
    JoinLayers(obj, assay = assay)
  }, error = function(e) {
    if (grepl("bad_alloc", conditionMessage(e), ignore.case = TRUE)) {
      message("JoinLayers ran out of memory, retrying after gc()...")
      gc()
      JoinLayers(obj, assay = assay)
    } else {
      stop(e)
    }
  })
}

## ------------------------------------------------------------------
## 0. Load data
## ------------------------------------------------------------------

rds_path <- "combined_seurat.rds"

if (file.exists(rds_path)) {
  combined_seurat <- readRDS(rds_path)
} else if (exists("combined_seurat")) {
  message("No .rds file found at '", rds_path, "' -- using `combined_seurat` already in memory.")
} else {
  stop("Could not find '", rds_path, "' and no `combined_seurat` object exists in memory.\n",
       "  - Check working directory with getwd()\n",
       "  - Or re-run script 01 and ensure it ends with saveRDS(combined_seurat, 'combined_seurat.rds')")
}

DefaultAssay(combined_seurat) <- "RNA"

message("Assays present: ", paste(Assays(combined_seurat), collapse = ", "))
message("Default assay: ", DefaultAssay(combined_seurat))
message("First 10 feature names: ", paste(head(rownames(combined_seurat), 10), collapse = ", "))
message("Layers (RNA, unjoined is normal here): ", paste(Layers(combined_seurat, assay = "RNA"), collapse = ", "))

## ------------------------------------------------------------------
## 1. Isolate CD45+ (PTPRC+) immune population
## ------------------------------------------------------------------
# NOTE: We deliberately do NOT call JoinLayers() on the full combined_seurat
# here -- on a large multi-sample object this can exhaust memory
# (std::bad_alloc). Instead we resolve the gene name once, then pull its
# expression layer-by-layer via get_gene_vector_safe(), which needs almost
# no extra memory since it only touches one gene's row at a time.

ptprc_symbol <- resolve_gene("PTPRC", combined_seurat)
if (is.na(ptprc_symbol)) {
  stop("No gene resembling 'PTPRC' was found in the object at all.\n",
       "  Check: head(rownames(combined_seurat)) and Assays(combined_seurat)")
}
message("Resolved PTPRC as: ", ptprc_symbol)

ptprc_vec <- get_gene_vector_safe(combined_seurat, ptprc_symbol, assay = "RNA", layer_search = "data")

immune_cells <- names(ptprc_vec)[ptprc_vec > 0]

message(sprintf("Identified %d / %d cells as CD45+ (%s+)",
                length(immune_cells), ncol(combined_seurat), ptprc_symbol))

if (length(immune_cells) == 0) {
  stop("Zero cells passed the PTPRC+ threshold. Check summary(ptprc_vec).")
}

immune_seurat <- subset(combined_seurat, cells = immune_cells)

# Now that the object is much smaller, joining layers is safe and needed
# for FindAllMarkers, GetAssayData, DotPlot, etc. downstream.

combined_seurat <- NormalizeData(combined_seurat)
ptprc_vec <- get_gene_vector_safe(combined_seurat, ptprc_symbol, assay = "RNA", layer_search = "data")
immune_cells <- names(ptprc_vec)[ptprc_vec > 0]

# 1. Subset to just the immune (PTPRC+) cells
immune_seurat <- subset(combined_seurat, cells = immune_cells)

# 2. Join layers on this smaller object (safe now, memory-wise)
immune_seurat <- safe_join_layers(immune_seurat, assay = "RNA")

# 3. Save checkpoint
saveRDS(immune_seurat, "immune_seurat_step1_subset.rds")


immune_seurat <- NormalizeData(immune_seurat)
immune_seurat <- FindVariableFeatures(immune_seurat, selection.method = "vst", nfeatures = 2000)

regress_vars <- intersect(c("nCount_RNA", "percent.mt"), colnames(immune_seurat@meta.data))
var_features <- VariableFeatures(immune_seurat)
immune_seurat <- ScaleData(immune_seurat, features = var_features, vars.to.regress = regress_vars)
immune_seurat <- RunPCA(immune_seurat, npcs = 30, features = var_features)

immune_seurat <- FindNeighbors(immune_seurat, dims = 1:20)
immune_seurat <- FindClusters(immune_seurat, resolution = 1.2)
immune_seurat <- RunUMAP(immune_seurat, dims = 1:20)

saveRDS(immune_seurat, "immune_seurat_step2_clustered.rds")



## ------------------------------------------------------------------
## 2. Re-cluster the immune subset at high resolution
## ------------------------------------------------------------------

stopifnot(
  "immune_seurat does not exist -- Step 1 must run before Step 2" = exists("immune_seurat")
)

immune_seurat <- NormalizeData(immune_seurat)
immune_seurat <- FindVariableFeatures(immune_seurat, selection.method = "vst", nfeatures = 2000)

regress_vars <- intersect(c("nCount_RNA", "percent.mt"), colnames(immune_seurat@meta.data))
if (length(regress_vars) < 2) {
  message("Note: not all standard QC columns found for regression. Using: ",
          paste(regress_vars, collapse = ", "))
}

var_features <- VariableFeatures(immune_seurat)
immune_seurat <- ScaleData(immune_seurat, features = var_features, vars.to.regress = regress_vars)
immune_seurat <- RunPCA(immune_seurat, npcs = 30, features = var_features)

ElbowPlot(immune_seurat, ndims = 30)
ggsave("plots/immune_elbow_plot.pdf", width = 6, height = 4)

n_pcs <- 20  # adjust based on elbow plot above

immune_seurat <- FindNeighbors(immune_seurat, dims = 1:n_pcs)
immune_seurat <- FindClusters(immune_seurat, resolution = 1.2)
immune_seurat <- RunUMAP(immune_seurat, dims = 1:n_pcs)

# CHECKPOINT
saveRDS(immune_seurat, "immune_seurat_step2_clustered.rds")
message("Checkpoint saved: immune_seurat_step2_clustered.rds")

## ------------------------------------------------------------------
## 3. Baseline UMAP (cluster-level, for pipeline validation)
## ------------------------------------------------------------------

stopifnot(
  "immune_seurat does not exist -- Step 2 must run before Step 3" = exists("immune_seurat")
)

umap_clusters <- DimPlot(immune_seurat, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("Immune sub-clusters (unannotated)")
print(umap_clusters)
ggsave("plots/immune_umap_clusters.pdf", umap_clusters, width = 7, height = 6)
ggsave("plots/immune_umap_clusters.png", umap_clusters, width = 7, height = 6, dpi = 300)

## ------------------------------------------------------------------
## 4. Marker panel for annotation (resolved against ENSEMBL--SYMBOL names)
## ------------------------------------------------------------------

marker_panel_raw <- list(
  "Pan-immune"        = c("PTPRC"),
  "Pan-T"              = c("CD3D", "CD3E", "CD3G"),
  "CD8+ T"             = c("CD8A", "CD8B"),
  "DUSP4hi CD4+ T"     = c("CD4", "DUSP4"),
  "TPT1hi CD4+ T"      = c("CD4", "TPT1"),
  "Treg"               = c("FOXP3", "IL2RA", "CTLA4"),
  "gamma-delta T"      = c("TRDC", "TRGC1", "TRGC2"),
  "NK"                 = c("NKG7", "GNLY", "KLRD1", "NCAM1"),
  "Pan-myeloid"        = c("LYZ", "CD68", "CD14"),
  "S100A8/A9+ Mono"    = c("S100A8", "S100A9"),
  "IL10hi Macrophage"  = c("IL10"),
  "CCL2hi Macrophage"  = c("CCL2"),
  "APOhi Macrophage"   = c("APOE", "APOC1"),
  "MAFhi Macrophage"   = c("MAF")
)

if (DefaultAssay(immune_seurat) != "RNA" && "RNA" %in% Assays(immune_seurat)) {
  message("DefaultAssay was '", DefaultAssay(immune_seurat), "' -- switching to 'RNA'.")
  DefaultAssay(immune_seurat) <- "RNA"
}

message("\nResolving marker panel against dataset feature names...")
marker_panel <- lapply(names(marker_panel_raw), function(nm) {
  message("  ", nm, ":")
  resolve_genes(marker_panel_raw[[nm]], immune_seurat)
})
names(marker_panel) <- names(marker_panel_raw)
marker_panel <- marker_panel[lengths(marker_panel) > 0]

all_markers <- unique(unlist(marker_panel))
message("\nTotal resolved markers for dotplot: ", length(all_markers))

if (length(all_markers) == 0) {
  stop("None of the curated marker genes were found in `immune_seurat`.\n",
       "  DefaultAssay: ", DefaultAssay(immune_seurat), "\n",
       "  First 20 features: ", paste(head(rownames(immune_seurat), 20), collapse = ", "))
}

## ------------------------------------------------------------------
## 5. Cluster marker discovery (data-driven check on top of curated panel)
## ------------------------------------------------------------------

stopifnot(
  "immune_seurat does not exist -- earlier steps must run first" = exists("immune_seurat")
)

# Layers are already joined (Step 1), so FindAllMarkers works normally here.
Idents(immune_seurat) <- "seurat_clusters"

use_presto <- requireNamespace("presto", quietly = TRUE)
if (use_presto) {
  presto_res <- presto::wilcoxauc(immune_seurat, group_by = "seurat_clusters")
  cluster_markers <- presto_res %>%
    filter(pct_in >= 25, logFC > 0.25) %>%
    transmute(
      cluster = group,
      gene    = feature,
      avg_log2FC = logFC,
      p_val      = pval,
      p_val_adj  = padj,
      pct.1      = pct_in / 100,
      pct.2      = pct_out / 100
    )
} else {
  message("Tip: install `presto` (devtools::install_github('immunogenomics/presto')) for much faster marker finding.")
  cluster_markers <- FindAllMarkers(immune_seurat,
                                    only.pos = TRUE,
                                    min.pct = 0.25,
                                    logfc.threshold = 0.25)
}

if (nrow(cluster_markers) == 0 || !"cluster" %in% colnames(cluster_markers)) {
  stop("No DE genes were identified across clusters. This usually means layers ",
       "were not joined before this step, or clustering produced degenerate clusters.\n",
       "  Layers(immune_seurat): ", paste(Layers(immune_seurat, assay = "RNA"), collapse = ", "))
}

top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(cluster_markers, "tables/immune_cluster_markers_full.csv", row.names = FALSE)
write.csv(top_markers, "tables/immune_cluster_markers_top10.csv", row.names = FALSE)

## ------------------------------------------------------------------
## 6. Marker dotplot (curated panel, for paper reproduction)
## ------------------------------------------------------------------

dotplot <- DotPlot(immune_seurat, features = all_markers, cluster.idents = TRUE) +
  RotatedAxis() +
  ggtitle("Immune lineage marker expression by cluster")
print(dotplot)
ggsave("plots/immune_marker_dotplot.pdf", dotplot, width = 12, height = 7)
ggsave("plots/immune_marker_dotplot.png", dotplot, width = 12, height = 7, dpi = 300)

## ------------------------------------------------------------------
## 7. Manual annotation
## ------------------------------------------------------------------
# IMPORTANT: Inspect `dotplot`, `top_markers`, and `cluster_markers` above
# to decide which numeric cluster corresponds to which cell type, then
# fill in the mapping below. This step is inherently manual/expert-driven.
# The placeholder below assumes 11 clusters (0-10) -- update to match your
# actual cluster count and marker-driven identities before finalizing.

cluster_annotations <- c(
  "0"  = "CD8+ T",
  "1"  = "DUSP4hi CD4+ T",
  "2"  = "TPT1hi CD4+ T",
  "3"  = "Treg",
  "4"  = "gamma-delta T",
  "5"  = "NK",
  "6"  = "S100A8/A9+ Monocyte",
  "7"  = "IL10hi Macrophage",
  "8"  = "CCL2hi Macrophage",
  "9"  = "APOhi Macrophage",
  "10" = "MAFhi Macrophage"
  # add/remove entries to match the actual number of clusters produced
)

present_clusters <- levels(immune_seurat$seurat_clusters)
unmapped <- setdiff(present_clusters, names(cluster_annotations))
if (length(unmapped) > 0) {
  warning("Unannotated clusters detected: ", paste(unmapped, collapse = ", "),
          " -- update `cluster_annotations` before finalizing.")
}

cluster_char <- as.character(immune_seurat$seurat_clusters)
recoded <- cluster_annotations[cluster_char]
recoded[is.na(recoded)] <- cluster_char[is.na(recoded)]  # keep original label if unmapped
immune_seurat$cell_type <- unname(recoded)

Idents(immune_seurat) <- "cell_type"

## ------------------------------------------------------------------
## 8. Annotated UMAP (final figure for paper reproduction)
## ------------------------------------------------------------------

umap_annotated <- DimPlot(immune_seurat, reduction = "umap", group.by = "cell_type",
                          label = TRUE, repel = TRUE) +
  ggtitle("Immune lineage annotation") +
  theme(legend.position = "right")
print(umap_annotated)
ggsave("plots/immune_umap_annotated.pdf", umap_annotated, width = 9, height = 6)
ggsave("plots/immune_umap_annotated.png", umap_annotated, width = 9, height = 6, dpi = 300)

dotplot_annotated <- DotPlot(immune_seurat, features = all_markers, group.by = "cell_type") +
  RotatedAxis() +
  ggtitle("Immune lineage marker expression by annotated cell type")
print(dotplot_annotated)
ggsave("plots/immune_marker_dotplot_annotated.pdf", dotplot_annotated, width = 12, height = 7)
ggsave("plots/immune_marker_dotplot_annotated.png", dotplot_annotated, width = 12, height = 7, dpi = 300)

## ------------------------------------------------------------------
## 9. Save final object
## ------------------------------------------------------------------

saveRDS(immune_seurat, "immune_seurat.rds")

message("Done. Outputs written to plots/, tables/, and immune_seurat.rds")
install.packages("Seurat")
install.packages("patchwork")
install.packages("dyplyr")
install.packages("DESeq2")
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(DESeq2)

set.seed(42)

# Ensure output directories exist
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Load Data & Resolve Feature Names
# ------------------------------------------------------------------------------
seurat_path <- if (file.exists("immune_seurat.rds")) {
  "immune_seurat.rds"
} else if (file.exists("objects/nb_singlet_clustered.rds")) {
  "objects/nb_singlet_clustered.rds"
} else {
  stop("No processed Seurat object found (checked 'immune_seurat.rds' and 'objects/nb_singlet_clustered.rds').")
}

message("Loading object from: ", seurat_path)
obj <- readRDS(seurat_path)
DefaultAssay(obj) <- "RNA"

# Helper function to match raw gene symbols or ENSEMBL--SYMBOL strings
resolve_gene <- function(symbol, seurat_obj) {
  available <- rownames(seurat_obj)
  hit <- available[match(tolower(symbol), tolower(available))]
  if (is.na(hit)) {
    suffix_hits <- available[grepl(paste0("(^|--)", symbol, "$"), available, ignore.case = TRUE)]
    if (length(suffix_hits) > 0) hit <- suffix_hits[1]
  }
  return(hit)
}

resolve_genes <- function(symbols, seurat_obj) {
  resolved <- vapply(symbols, resolve_gene, seurat_obj = seurat_obj, FUN.VALUE = character(1))
  unname(resolved[!is.na(resolved)])
}

# ------------------------------------------------------------------------------
# 2. Extract Longitudinal Metadata (Pre vs. Post Chemotherapy)
# ------------------------------------------------------------------------------
# Extract timepoint or treatment state from sample_id or patient_id if not explicitly present
if (!"treatment_status" %in% colnames(obj@meta.data)) {
  if ("sample_id" %in% colnames(obj@meta.data)) {
    obj$treatment_status <- ifelse(grepl("T2|post|Post|relapse|treated", obj$sample_id, ignore.case = TRUE), 
                                   "Post-Chemo", "Pre-Chemo")
  } else {
    warning("Could not determine treatment status automatically. Defaulting dummy grouping.")
    obj$treatment_status <- "Pre-Chemo"
  }
}

obj$treatment_status <- factor(obj$treatment_status, levels = c("Pre-Chemo", "Post-Chemo"))
message("Treatment distribution:")
print(table(obj$treatment_status))

# Set active identity to cell types
if ("cell_type" %in% colnames(obj@meta.data)) {
  Idents(obj) <- "cell_type"
} else if ("predicted_cell_type" %in% colnames(obj@meta.data)) {
  Idents(obj) <- "predicted_cell_type"
} else {
  Idents(obj) <- "seurat_clusters"
}

# ------------------------------------------------------------------------------
# 3. Pseudobulk Differential Gene Expression (Pre vs. Post Chemotherapy)
# ------------------------------------------------------------------------------
run_pseudobulk_de <- function(seurat_obj, cell_subset_label) {
  message(sprintf("\n--- Running Pseudobulk DE for: %s ---", cell_subset_label))
  
  sub_obj <- subset(seurat_obj, idents = cell_subset_label)
  
  if (length(unique(sub_obj$treatment_status)) < 2) {
    warning(sprintf("Skipping %s: insufficient treatment categories.", cell_subset_label))
    return(NULL)
  }
  
  # Pseudobulk aggregation by patient_id & treatment_status
  ps_obj <- AggregateExpression(
    sub_obj,
    group.by = c("patient_id", "treatment_status"),
    return.seurat = TRUE,
    slot = "counts"
  )
  
  counts_matrix <- GetAssayData(ps_obj, slot = "counts")
  
  # Construct metadata table for DESeq2
  col_data <- data.frame(row.names = colnames(counts_matrix))
  col_data$patient_id <- gsub("_.*", "", rownames(col_data))
  col_data$treatment_status <- factor(
    ifelse(grepl("Post-Chemo", rownames(col_data)), "Post-Chemo", "Pre-Chemo"),
    levels = c("Pre-Chemo", "Post-Chemo")
  )
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_matrix),
    colData = col_data,
    design = ~ treatment_status
  )
  
  # Filter low-count genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("treatment_status", "Post-Chemo", "Pre-Chemo"))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Save results table
  out_csv <- sprintf("results/pseudobulk_DE_%s_post_vs_pre.csv", gsub("[^A-Za-z0-9]", "_", cell_subset_label))
  write.csv(res_df, out_csv, row.names = FALSE)
  
  # Volcano Plot
  res_df$significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.58, "Significant", "Not Sig")
  
  volcano_p <- ggplot(res_df[!is.na(res_df$padj), ], aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("gray60", "red3")) +
    theme_classic() +
    labs(
      title = paste("DE Genes (Post vs Pre):", cell_subset_label),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-Value"
    )
  
  out_png <- sprintf("figures/volcano_DE_%s.png", gsub("[^A-Za-z0-9]", "_", cell_subset_label))
  ggsave(out_png, plot = volcano_p, width = 7, height = 5, dpi = 300)
  
  return(res_df)
}

# Run DE across T and NK populations
target_cell_types <- c("CD8+ T", "NK", "T_cells", "NK_cells", "Treg")
available_targets <- intersect(target_cell_types, unique(Idents(obj)))

de_results_list <- list()
for (ct in available_targets) {
  de_results_list[[ct]] <- run_pseudobulk_de(obj, ct)
}

# ------------------------------------------------------------------------------
# 4. T-Cell Dysfunction & Exhaustion Scoring
# ------------------------------------------------------------------------------
message("\n--- Computing T-Cell Exhaustion Scores ---")

exhaustion_genes_raw <- c("PDCD1", "LAG3", "CTLA4", "TIGIT", "HAVCR2", "TOX", "ENTPD1")
exhaustion_genes <- resolve_genes(exhaustion_genes_raw, obj)

if (length(exhaustion_genes) > 0) {
  obj <- AddModuleScore(
    object = obj,
    features = list(exhaustion_genes),
    name = "T_Cell_Exhaustion_Score"
  )
  
  # Plot Exhaustion Scores
  p_exh <- VlnPlot(
    obj, 
    features = "T_Cell_Exhaustion_Score1", 
    split.by = "treatment_status",
    group.by = if ("cell_type" %in% colnames(obj@meta.data)) "cell_type" else "seurat_clusters",
    pt.size = 0.1
  ) + ggtitle("T-Cell Exhaustion / Dysfunction Score (Pre vs. Post)")
  
  ggsave("figures/tcell_exhaustion_score_vln.png", plot = p_exh, width = 10, height = 6, dpi = 300)
}

# ------------------------------------------------------------------------------
# 5. NK Cell Cytotoxicity & TGF-beta Pathway Scoring
# ------------------------------------------------------------------------------
message("\n--- Computing NK Cytotoxicity & TGF-beta Scores ---")

cytotoxicity_genes_raw <- c("GZMB", "PRF1", "GNLY", "NKG7", "KLRD1")
cytotoxicity_genes <- resolve_genes(cytotoxicity_genes_raw, obj)

if (length(cytotoxicity_genes) > 0) {
  obj <- AddModuleScore(
    object = obj,
    features = list(cytotoxicity_genes),
    name = "NK_Cytotoxicity_Score"
  )
}

tgfbeta_genes_raw <- c("TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4")
tgfbeta_genes <- resolve_genes(tgfbeta_genes_raw, obj)

if (length(tgfbeta_genes) > 0) {
  obj <- AddModuleScore(
    object = obj,
    features = list(tgfbeta_genes),
    name = "TGFbeta_Pathway_Score"
  )
}

# Combine dysfunction & cytotoxicity feature plots
p_scores <- FeaturePlot(
  obj, 
  features = intersect(c("T_Cell_Exhaustion_Score1", "NK_Cytotoxicity_Score1", "TGFbeta_Pathway_Score1"), colnames(obj@meta.data)),
  reduction = "umap",
  ncol = 2
)
ggsave("figures/functional_module_scores_umap.png", plot = p_scores, width = 10, height = 8, dpi = 300)

# Save updated dataset with dysfunction metrics
saveRDS(obj, "objects/immune_seurat_dysfunction_scored.rds")
message("\nCompleted Task 3 execution. Outputs generated in 'results/' and 'figures/'.")
getwd()
list.files()
list.files("figures")
list.files("results")
list.files("objects")
list.files("C:/Users/clinical_1/Documents", 
           pattern = "\\.(png|rds|csv|tsv)$", 
           recursive = TRUE, 
           full.names = TRUE)
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("objects", showWarnings = FALSE)
saveRDS(obj, "objects/immune_seurat_dysfunction_scored.rds")
ls()
exists("obj")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
a
# Objectives:
#   1. Pseudobulk DE analysis (Pre- vs. Post-chemotherapy) for T and NK subsets.
#   2. Compute module dysfunction scores for T cells (exhaustion markers).
#   3. Compute cytotoxicity scores for NK cells (GZMB, PRF1, GNLY) & TGF-beta activity.
#   4. Save DE tables to results/ and plots to figures/
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(DESeq2)

set.seed(42)

# Ensure output directories exist
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Load Data & Resolve Feature Names
# ------------------------------------------------------------------------------
seurat_path <- if (file.exists("immune_seurat.rds")) {
  "immune_seurat.rds"
} else if (file.exists("objects/nb_singlet_clustered.rds")) {
  "objects/nb_singlet_clustered.rds"
} else {
  stop("No processed Seurat object found (checked 'immune_seurat.rds' and 'objects/nb_singlet_clustered.rds').")
}

message("Loading object from: ", seurat_path)
obj <- readRDS(seurat_path)
DefaultAssay(obj) <- "RNA"

# Helper function to match raw gene symbols or ENSEMBL--SYMBOL strings
resolve_gene <- function(symbol, seurat_obj) {
  available <- rownames(seurat_obj)
  hit <- available[match(tolower(symbol), tolower(available))]
  if (is.na(hit)) {
    suffix_hits <- available[grepl(paste0("(^|--)", symbol, "$"), available, ignore.case = TRUE)]
    if (length(suffix_hits) > 0) hit <- suffix_hits[1]
  }
  return(hit)
}

resolve_genes <- function(symbols, seurat_obj) {
  resolved <- vapply(symbols, resolve_gene, seurat_obj = seurat_obj, FUN.VALUE = character(1))
  unname(resolved[!is.na(resolved)])
}

# ------------------------------------------------------------------------------
# 2. Extract Longitudinal Metadata (Pre vs. Post Chemotherapy)
# ------------------------------------------------------------------------------
# Extract timepoint or treatment state from sample_id or patient_id if not explicitly present
if (!"treatment_status" %in% colnames(obj@meta.data)) {
  if ("sample_id" %in% colnames(obj@meta.data)) {
    obj$treatment_status <- ifelse(grepl("T2|post|Post|relapse|treated", obj$sample_id, ignore.case = TRUE), 
                                   "Post-Chemo", "Pre-Chemo")
  } else {
    warning("Could not determine treatment status automatically. Defaulting dummy grouping.")
    obj$treatment_status <- "Pre-Chemo"
  }
}

obj$treatment_status <- factor(obj$treatment_status, levels = c("Pre-Chemo", "Post-Chemo"))
message("Treatment distribution:")
print(table(obj$treatment_status))

# Set active identity to cell types
if ("cell_type" %in% colnames(obj@meta.data)) {
  Idents(obj) <- "cell_type"
} else if ("predicted_cell_type" %in% colnames(obj@meta.data)) {
  Idents(obj) <- "predicted_cell_type"
} else {
  Idents(obj) <- "seurat_clusters"
}

# ------------------------------------------------------------------------------
# 3. Pseudobulk Differential Gene Expression (Pre vs. Post Chemotherapy)
# ------------------------------------------------------------------------------
run_pseudobulk_de <- function(seurat_obj, cell_subset_label) {
  message(sprintf("\n--- Running Pseudobulk DE for: %s ---", cell_subset_label))
  
  sub_obj <- subset(seurat_obj, idents = cell_subset_label)
  
  if (length(unique(sub_obj$treatment_status)) < 2) {
    warning(sprintf("Skipping %s: insufficient treatment categories.", cell_subset_label))
    return(NULL)
  }
  
  # Pseudobulk aggregation by patient_id & treatment_status
  ps_obj <- AggregateExpression(
    sub_obj,
    group.by = c("patient_id", "treatment_status"),
    return.seurat = TRUE,
    slot = "counts"
  )
  
  counts_matrix <- GetAssayData(ps_obj, slot = "counts")
  
  # Construct metadata table for DESeq2
  col_data <- data.frame(row.names = colnames(counts_matrix))
  col_data$patient_id <- gsub("_.*", "", rownames(col_data))
  col_data$treatment_status <- factor(
    ifelse(grepl("Post-Chemo", rownames(col_data)), "Post-Chemo", "Pre-Chemo"),
    levels = c("Pre-Chemo", "Post-Chemo")
  )
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_matrix),
    colData = col_data,
    design = ~ treatment_status
  )
  
  # Filter low-count genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("treatment_status", "Post-Chemo", "Pre-Chemo"))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Save results table
  out_csv <- sprintf("results/pseudobulk_DE_%s_post_vs_pre.csv", gsub("[^A-Za-z0-9]", "_", cell_subset_label))
  write.csv(res_df, out_csv, row.names = FALSE)
  
  # Volcano Plot
  res_df$significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.58, "Significant", "Not Sig")
  
  volcano_p <- ggplot(res_df[!is.na(res_df$padj), ], aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("gray60", "red3")) +
    theme_classic() +
    labs(
      title = paste("DE Genes (Post vs Pre):", cell_subset_label),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-Value"
    )
  
  out_png <- sprintf("figures/volcano_DE_%s.png", gsub("[^A-Za-z0-9]", "_", cell_subset_label))
  ggsave(out_png, plot = volcano_p, width = 7, height = 5, dpi = 300)
  
  return(res_df)
}

# Run DE across T and NK populations
target_cell_types <- c("CD8+ T", "NK", "T_cells", "NK_cells", "Treg")
available_targets <- intersect(target_cell_types, unique(Idents(obj)))

de_results_list <- list()
for (ct in available_targets) {
  de_results_list[[ct]] <- run_pseudobulk_de(obj, ct)
}

# ------------------------------------------------------------------------------
# 4. T-Cell Dysfunction & Exhaustion Scoring
# ------------------------------------------------------------------------------
message("\n--- Computing T-Cell Exhaustion Scores ---")

exhaustion_genes_raw <- c("PDCD1", "LAG3", "CTLA4", "TIGIT", "HAVCR2", "TOX", "ENTPD1")
exhaustion_genes <- resolve_genes(exhaustion_genes_raw, obj)

if (length(exhaustion_genes) > 0) {
  obj <- AddModuleScore(
    object = obj,
    features = list(exhaustion_genes),
    name = "T_Cell_Exhaustion_Score"
  )
  
  # Plot Exhaustion Scores
  p_exh <- VlnPlot(
    obj, 
    features = "T_Cell_Exhaustion_Score1", 
    split.by = "treatment_status",
    group.by = if ("cell_type" %in% colnames(obj@meta.data)) "cell_type" else "seurat_clusters",
    pt.size = 0.1
  ) + ggtitle("T-Cell Exhaustion / Dysfunction Score (Pre vs. Post)")
  
  ggsave("figures/tcell_exhaustion_score_vln.png", plot = p_exh, width = 10, height = 6, dpi = 300)
}

# ------------------------------------------------------------------------------
# 5. NK Cell Cytotoxicity & TGF-beta Pathway Scoring
# ------------------------------------------------------------------------------
message("\n--- Computing NK Cytotoxicity & TGF-beta Scores ---")

cytotoxicity_genes_raw <- c("GZMB", "PRF1", "GNLY", "NKG7", "KLRD1")
cytotoxicity_genes <- resolve_genes(cytotoxicity_genes_raw, obj)

if (length(cytotoxicity_genes) > 0) {
  obj <- AddModuleScore(
    object = obj,
    features = list(cytotoxicity_genes),
    name = "NK_Cytotoxicity_Score"
  )
}

tgfbeta_genes_raw <- c("TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4")
tgfbeta_genes <- resolve_genes(tgfbeta_genes_raw, obj)

if (length(tgfbeta_genes) > 0) {
  obj <- AddModuleScore(
    object = obj,
    features = list(tgfbeta_genes),
    name = "TGFbeta_Pathway_Score"
  )
}

# Combine dysfunction & cytotoxicity feature plots
p_scores <- FeaturePlot(
  obj, 
  features = intersect(c("T_Cell_Exhaustion_Score1", "NK_Cytotoxicity_Score1", "TGFbeta_Pathway_Score1"), colnames(obj@meta.data)),
  reduction = "umap",
  ncol = 2
)
ggsave("figures/functional_module_scores_umap.png", plot = p_scores, width = 10, height = 8, dpi = 300)

# Save updated dataset with dysfunction metrics
saveRDS(obj, "objects/immune_seurat_dysfunction_scored.rds")
message("\nCompleted Task 3 execution. Outputs generated in 'results/' and 'figures/'.")
