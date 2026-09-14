# ============================================================
# Neuroblastoma scRNA-seq - Task 1
# Dataset: GSE218003
# QC, Preprocessing, Integration and Clustering
# ============================================================
install.packages("scDblFinder")
install.packages("harmony")
install.packages("dplyr")
install.packages("SingleCellExperiment")
install.packages("Matrix")

# Load Packages
library(Seurat)
library(Matrix)
library(SingleCellExperiment)
library(scDblFinder)
library(ggplot2)
library(harmony)
library(dplyr)

# Set Random Seed
set.seed(1234)

# Create Output Folders
dir.create("figures",showWarnings=FALSE)
dir.create("results",showWarnings=FALSE)
dir.create("objects",showWarnings=FALSE)

# Set Data Directory
data_dir<-"C:/Users/Ahmed Amer/Downloads/GSE218003_RAW"

# Load Count Files
count_files<-sort(list.files(data_dir,pattern="_counts\\.tsv\\.gz$",full.names=TRUE))
count_ids<-sub("_counts\\.tsv\\.gz$","",basename(count_files))

# Check Count Files
cat("Count files:",length(count_files),"\n")

# Create Seurat Objects
seurat_list<-lapply(seq_along(count_files),function(i){
  counts<-read.delim(gzfile(count_files[i]),header=TRUE,row.names=1,check.names=FALSE)
  stopifnot(!anyNA(counts))
  rownames(counts)<-make.unique(gsub("_","-",sub("^.*__","",rownames(counts))))
  counts<-Matrix::Matrix(as.matrix(counts),sparse=TRUE)
  sample_id<-sub("_scR_CEL-Seq2.*$","",sub("^GSM[0-9]+_","",count_ids[i]))
  patient_id<-sub("_T[0-9]*$","",sample_id)
  obj<-CreateSeuratObject(counts=counts,project="Neuroblastoma_scRNAseq")
  obj$plate_id<-count_ids[i]
  obj$sample_id<-sample_id
  obj$patient_id<-patient_id
  obj
})

# Name Seurat Objects
names(seurat_list)<-count_ids

# Check Loaded Data
cat("Objects:",length(seurat_list),"\n")
cat("NA:",sum(is.na(seurat_list[[1]][["RNA"]]$counts)),"\n")
head(seurat_list[[1]]@meta.data)

# Calculate Mitochondrial Percentage
seurat_list<-lapply(seurat_list,function(x){x[["percent.mt"]]<-PercentageFeatureSet(x,pattern="^MT-");x})

# Merge Seurat Objects
nb<-merge(seurat_list[[1]],y=seurat_list[-1],add.cell.ids=names(seurat_list))

# Check Sample Metadata
cat("Technical units:",length(unique(nb$plate_id)),"\n")
cat("Samples:",length(unique(nb$sample_id)),"\n")
cat("Patients:",length(unique(nb$patient_id)),"\n")

# Save Sample List
sample_manifest<-unique(nb@meta.data[,c("plate_id","sample_id","patient_id")])
write.csv(sample_manifest,"results/sample_manifest.csv",row.names=FALSE)

# QC Summary Before Filtering
summary(nb@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt")])

# Cell Count Before Filtering
cells_before_qc<-ncol(nb)
cat("Cells before QC:",cells_before_qc,"\n")

# QC Violin Plot Before Filtering
qc_violin_before<-VlnPlot(nb,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
qc_violin_before

# QC Scatter Plots Before Filtering
plot1<-FeatureScatter(nb,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(nb,feature1="nCount_RNA",feature2="nFeature_RNA")
qc_scatter_before<-plot1+plot2
qc_scatter_before

# Save QC Plots Before Filtering
ggsave("figures/01_qc_violin_before.png",plot=qc_violin_before,width=12,height=5,dpi=300)
ggsave("figures/02_qc_scatter_before.png",plot=qc_scatter_before,width=10,height=5,dpi=300)

# Cell Filtering
nb_filtered<-subset(nb,subset=nFeature_RNA>=500&nCount_RNA>=1000&nCount_RNA<=150000&percent.mt<10)

# Cell Count After QC
cells_after_qc<-ncol(nb_filtered)
cat("Cells after QC:",cells_after_qc,"\n")

# QC Summary After Filtering
summary(nb_filtered@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt")])

# QC Violin Plot After Filtering
qc_violin_after<-VlnPlot(nb_filtered,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
qc_violin_after

# QC Scatter Plots After Filtering
plot1<-FeatureScatter(nb_filtered,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(nb_filtered,feature1="nCount_RNA",feature2="nFeature_RNA")
qc_scatter_after<-plot1+plot2
qc_scatter_after

# Save QC Plots After Filtering
ggsave("figures/03_qc_violin_after.png",plot=qc_violin_after,width=12,height=5,dpi=300)
ggsave("figures/04_qc_scatter_after.png",plot=qc_scatter_after,width=10,height=5,dpi=300)

# Check Cells Per Plate (Technical QC - identify failed plates)
plate_cell_counts<-table(nb_filtered$plate_id)
sort(plate_cell_counts)

# Document Plates With Likely Technical Failure
failed_plates<-plate_cell_counts[plate_cell_counts<5]
cat("Plates with likely technical failure (<5 cells surviving QC):",length(failed_plates),"\n")
print(failed_plates)
write.csv(as.data.frame(plate_cell_counts),"results/cells_per_plate.csv",row.names=FALSE)

# ============================================================
# Doublet Detection Block
# ============================================================
nb_filtered<-JoinLayers(nb_filtered)

# Check Cells Per Sample
sample_cell_counts<-table(nb_filtered$sample_id)
sort(sample_cell_counts)

# Remove Samples With Too Few Cells for scDblFinder
min_cells_sample<-30
valid_samples<-names(sample_cell_counts[sample_cell_counts>=min_cells_sample])
cat("Samples removed:",length(sample_cell_counts)-length(valid_samples),"\n")
cat("Samples kept:",length(valid_samples),"\n")
nb_filtered<-subset(nb_filtered,subset=sample_id %in% valid_samples)
cat("Cells after removing small samples:",ncol(nb_filtered),"\n")

sce<-as.SingleCellExperiment(nb_filtered)
cat("Doublet groups (samples):",length(unique(nb_filtered$sample_id)),"\n")
sce<-scDblFinder(sce,samples="sample_id")
nb_filtered$doublet_score<-sce$scDblFinder.score
nb_filtered$doublet_class<-sce$scDblFinder.class
table(nb_filtered$doublet_class)
hist(nb_filtered$doublet_score,breaks=50,main="Doublet Score",xlab="Doublet Score")
VlnPlot(nb_filtered,features="doublet_score",group.by="doublet_class")
doublet_scatter<-ggplot(nb_filtered@meta.data,aes(x=nCount_RNA,y=doublet_score,color=doublet_class))+geom_point(size=1,alpha=0.6)+labs(x="nCount_RNA",y="Doublet Score")+theme_classic()
doublet_scatter
nb_singlet<-subset(nb_filtered,subset=doublet_class=="singlet")
cells_after_doublets<-ncol(nb_singlet)
cat("Cells after doublet removal:",cells_after_doublets,"\n")
























# Save Cell Count Summary
cell_count_summary<-data.frame(Stage=c("Before_QC","After_QC","After_Doublet_Removal"),Cells=c(cells_before_qc,cells_after_qc,cells_after_doublets))
cell_count_summary
write.csv(cell_count_summary,"results/cell_count_summary.csv",row.names=FALSE)

# Save QC Object
saveRDS(nb_singlet,"objects/nb_singlet_qc.rds")

# Normalization
nb_singlet<-NormalizeData(nb_singlet,normalization.method="LogNormalize",scale.factor=10000)

# Highly Variable Features
nb_singlet<-FindVariableFeatures(nb_singlet,selection.method="vst",nfeatures=2000)

# Variable Feature Plot
plot1<-VariableFeaturePlot(nb_singlet)
top10<-head(VariableFeatures(nb_singlet),10)
plot2<-LabelPoints(plot=plot1,points=top10,repel=TRUE)
plot1+plot2

# Scale Data
nb_singlet<-ScaleData(nb_singlet)

# PCA
nb_singlet<-RunPCA(nb_singlet,features=VariableFeatures(nb_singlet))

# PCA Plot
DimPlot(nb_singlet,reduction="pca")

# Elbow Plot
ElbowPlot(nb_singlet,ndims=50)

# 20 PCs Selected From Elbow Plot

# UMAP Before Integration
nb_singlet<-RunUMAP(nb_singlet,reduction="pca",dims=1:20,reduction.name="umap_unintegrated",reduction.key="UMAPUNINT_")

# UMAP by Sample Before Integration
umap_before<-DimPlot(nb_singlet,reduction="umap_unintegrated",group.by="sample_id")
umap_before

# Save UMAP Before Integration
ggsave("figures/05_umap_before_integration.png",plot=umap_before,width=10,height=8,dpi=300)

# Harmony Integration by Technical Plate
nb_singlet<-RunHarmony(nb_singlet,group.by.vars="plate_id",reduction.use="pca",dims.use=1:20,verbose=FALSE)

# Check Harmony
Reductions(nb_singlet)

# Find Neighbors
nb_singlet<-FindNeighbors(nb_singlet,reduction="harmony",dims=1:20)

# Clustering
nb_singlet<-FindClusters(nb_singlet,resolution=0.5)

# Check Cluster Number
cat("Clusters:",length(unique(nb_singlet$seurat_clusters)),"\n")
table(nb_singlet$seurat_clusters)

# UMAP After Integration
nb_singlet<-RunUMAP(nb_singlet,reduction="harmony",dims=1:20)

# Cluster UMAP
cluster_umap<-DimPlot(nb_singlet,reduction="umap",label=TRUE)
cluster_umap

# Save Cluster UMAP
ggsave("figures/06_cluster_umap.png",plot=cluster_umap,width=9,height=7,dpi=300)

# UMAP by Sample After Integration (verify batch correction worked)
umap_after_sample<-DimPlot(nb_singlet,reduction="umap",group.by="sample_id")
ggsave("figures/06b_umap_after_integration_by_sample.png",plot=umap_after_sample,width=10,height=8,dpi=300)

# UMAP by Plate After Integration (verify technical effect removed)
umap_after_plate<-DimPlot(nb_singlet,reduction="umap",group.by="plate_id")
ggsave("figures/06c_umap_after_integration_by_plate.png",plot=umap_after_plate,width=10,height=8,dpi=300)

# Set Cluster Identities
Idents(nb_singlet)<-"seurat_clusters"

# Find Cluster Markers
markers<-FindAllMarkers(nb_singlet,only.pos=TRUE,min.pct=0.25,logfc.threshold=0.25)

# Top 10 Markers Per Cluster
top_markers<-markers%>%dplyr::group_by(cluster)%>%dplyr::slice_max(order_by=avg_log2FC,n=10,with_ties=FALSE)

# View Top Markers
print(top_markers, n=Inf)

# Save Marker Tables
write.csv(markers,"results/all_cluster_markers.csv",row.names=FALSE)
write.csv(top_markers,"results/top10_cluster_markers.csv",row.names=FALSE)

# Broad Marker Validation
broad_marker_plot<-DotPlot(nb_singlet,features=c("PHOX2B","HAND2","DBH","TH","CHGB","CD3D","CD3E","IL7R","NKG7","GNLY","PRF1","MS4A1","CD79A","CD19","LST1","TYROBP","C1QA","FOLR2","PECAM1","VWF","KDR","CLDN5","CDH5","COL1A1","COL1A2","DCN","LUM","SFRP2","SOX10","S100B","PLP1","ERBB3","NGFR","HBB","HBA1","GYPA"))+RotatedAxis()
broad_marker_plot

# Save Broad Marker Plot
ggsave("figures/07_broad_marker_dotplot.png",plot=broad_marker_plot,width=15,height=8,dpi=300)




# Save Both Plots
ggsave("figures/08_cluster_umap_numbered.png", plot=cluster_umap_numbers, width=9, height=7, dpi=300)
ggsave("figures/09_cluster_umap_labeled.png", plot=cluster_umap_labeled, width=9, height=7, dpi=300)

# Save Cluster Annotation Table
write.csv(cluster_scores, "results/cluster_predicted_types.csv", row.names=FALSE)

exists("nb_singlet")
exists("cluster_scores")

# ============================================================
# Cluster Annotation Based on Marker Scores (Preliminary)
# ============================================================

marker_sets <- list(
  Malignant_Neuroblast = c("PHOX2B","HAND2","DBH","TH","CHGB"),
  T_cells               = c("CD3D","CD3E","IL7R"),
  NK_cells              = c("NKG7","GNLY","PRF1"),
  B_cells               = c("MS4A1","CD79A","CD19"),
  Myeloid               = c("LST1","TYROBP","C1QA","FOLR2"),
  Endothelial           = c("PECAM1","VWF","KDR","CLDN5","CDH5"),
  Fibroblast            = c("COL1A1","COL1A2","DCN","LUM","SFRP2"),
  Schwann_cells         = c("SOX10","S100B","PLP1","ERBB3","NGFR"),
  Erythroid             = c("HBB","HBA1","GYPA")
)

all_markers <- unlist(marker_sets)
missing_markers <- all_markers[!all_markers %in% rownames(nb_singlet)]
cat("Missing markers:", length(missing_markers), "\n")
print(missing_markers)

marker_sets <- lapply(marker_sets, function(x) x[x %in% rownames(nb_singlet)])

for (ct in names(marker_sets)) {
  nb_singlet <- AddModuleScore(
    nb_singlet,
    features = list(marker_sets[[ct]]),
    name = ct
  )
}

score_cols <- paste0(names(marker_sets), "1")

library(dplyr)
cluster_scores <- nb_singlet@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean))

print(cluster_scores)

cluster_scores$predicted_type <- names(marker_sets)[
  apply(cluster_scores[,score_cols], 1, which.max)
]

print(cluster_scores[,c("seurat_clusters","predicted_type")])

cluster_scores$seurat_clusters <- as.character(cluster_scores$seurat_clusters)

nb_singlet$predicted_cell_type <- cluster_scores$predicted_type[
  match(as.character(nb_singlet$seurat_clusters), cluster_scores$seurat_clusters)
]

table(nb_singlet$predicted_cell_type, useNA="always")

cluster_umap_numbered <- DimPlot(nb_singlet, reduction="umap", label=TRUE, group.by="seurat_clusters") +
  ggtitle("Clusters (Numbered)")

cluster_umap_labeled <- DimPlot(nb_singlet, reduction="umap", label=TRUE, group.by="predicted_cell_type") +
  ggtitle("Clusters (Predicted Cell Type)")

cluster_umap_numbered
cluster_umap_labeled

ggsave("figures/08_cluster_umap_numbered.png", plot=cluster_umap_numbered, width=9, height=7, dpi=300)
ggsave("figures/09_cluster_umap_labeled.png", plot=cluster_umap_labeled, width=9, height=7, dpi=300)

write.csv(cluster_scores, "results/cluster_predicted_types.csv", row.names=FALSE)







# Save Clustered Object
saveRDS(nb_singlet,"objects/nb_singlet_clustered.rds")

# Check Saved Object
file.exists("objects/nb_singlet_clustered.rds")

# Save Session Information
writeLines(capture.output(sessionInfo()),"results/sessionInfo.txt")

# Session Information
sessionInfo()


# How many of the 13 original clusters were classified as Malignant?
malignant_clusters <- cluster_scores$seurat_clusters[cluster_scores$predicted_type == "Malignant_Neuroblast"]
cat("Number of malignant clusters:", length(malignant_clusters), "\n")
print(malignant_clusters)

# How many cells total are malignant?
cat("Total malignant cells:", sum(nb_singlet$predicted_cell_type == "Malignant_Neuroblast", na.rm=TRUE), "\n")

# Malignant cells broken down by patient
table(nb_singlet$patient_id[nb_singlet$predicted_cell_type == "Malignant_Neuroblast"])


unique(nb_singlet$predicted_cell_type)
