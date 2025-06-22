## Leona Task 2 ##

# Downloading raw counts files via CL
# Bash: wget -O raw_count.zip "https://www.ebi.ac.uk/gxa/sc/experiment/E-MTAB-10596/download/zip?fileType=quantification-raw&accessKey="
# Bash: unzip raw_count.zip
# This unzips three files (mtx,cols,rows)

# Installing/loading packs
install.packages('BiocManager')
install.packages('patchwork')
install.packages('Seurat')
BiocManager::install("org.Hs.eg.db")

library(Seurat)
library(patchwork)
library(dplyr)
library(BiocManager)
library(org.Hs.eg.db)
library(ggplot2)

# Non-obligatory packages
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))

# Install the remotes package
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages('Signac')
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)
install.packages("installr")

# Assign file paths to variables (files were downloaded via bash CL from Matrix Market archive and copied to windows)
matrix_file <- "C:\\Users\\leona\\Desktop\\MSc_Bioinformatik\\1_Semester\\07_TrGe_Übung\\Aufgabe_2\\Raw_files\\E-MTAB-10596.aggregated_filtered_counts.mtx"
barcodes_cols_file <- "C:\\Users\\leona\\Desktop\\MSc_Bioinformatik\\1_Semester\\07_TrGe_Übung\\Aufgabe_2\\Raw_files\\E-MTAB-10596.aggregated_filtered_counts.mtx_cols"
features_rows_file <- "C:\\Users\\leona\\Desktop\\MSc_Bioinformatik\\1_Semester\\07_TrGe_Übung\\Aufgabe_2\\Raw_files\\E-MTAB-10596.aggregated_filtered_counts.mtx_rows"

# Checking if the file exists
file.exists(features_rows_file)

# Make Seurat object 
seurat_data <- ReadMtx(mtx = matrix_file, cells = barcodes_cols_file, features = features_rows_file )
seurat_obj <- CreateSeuratObject(counts = seurat_data, project = "TrGeUE_Leona")
seurat_obj

# View first 5 rows of metadata with QC metrics
head(seurat_obj@meta.data, 10)
VlnPlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA", ncol = 2))

# Summary statistics for QC
summary(seurat_obj$nCount_RNA)
summary(seurat_obj$nFeature_RNA)

# Filtering based on summary stats
seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > 500 & 
                       nFeature_RNA < 9000 & 
                       nCount_RNA > 1000 & 
                       nCount_RNA < 120000)
# Scatter plot for QC
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1

# Normalizing data
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Finding highly variable features (high cell-cell variation)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat_obj), 10)

# Plot variable features
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size = 2, xnudge = 0, ynudge = 0)
plot1
plot2

# Scaling all genes (centering, var normalization, linear transformation)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

# PCA on variable features
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)

# PCA visualizations
DimPlot(seurat_obj, reduction = "pca")  + NoLegend()

# Heat maps for PCA1, 1:5 and 1:30
DimHeatmap(seurat_obj, dims = 1, cells = 300, balanced = TRUE)
DimHeatmap(seurat_obj, dims = 1:5, cells = 200, balanced = TRUE)
DimHeatmap(seurat_obj, dims = 1:20, cells = 200, balanced = TRUE)

VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")

ElbowPlot(seurat_obj)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:12)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
head(Idents(seurat_obj), 5)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:12)
DimPlot(seurat_obj, reduction = "umap")

# Finding markers for every cluster
# compared to all other cells (only positive)

seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)

# Filter markers

filter_markers <- seurat_obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
filter_markers

# Filter top 10 markers per cluster with avg_log2FC > 1

top10 <- seurat_obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() 

# Heatmap for top10
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()

# Naming clusters

# Define new cluster identities as numbers 0-6
new.cluster.ids <- as.character(0:6)
names(new.cluster.ids) <- levels(seurat_obj)

# Rename cluster identities
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)

# Plot UMAP with numeric labels
plot3 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 1, label.size = 4) + 
  xlab("UMAP 1") + 
  ylab("UMAP 2") +
  theme(axis.title = element_text(size = 10), legend.text = element_text(size = 10)) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

plot3

# Save the plot
ggsave(filename = "C:/Users/leona/Desktop/MSc_Bioinformatik/1_Semester/07_TrGe_Übung/Aufgabe_2/Umap_Vracar_1.jpeg",
       height = 7, width = 12, plot = plot3, quality = 50)

# Save the Seurat object as an RDS file
saveRDS(seurat_obj, file = "C:/Users/leona/Desktop/MSc_Bioinformatik/1_Semester/07_TrGe_Übung/Aufgabe_2/Seurat_Vracar_1.rds")

