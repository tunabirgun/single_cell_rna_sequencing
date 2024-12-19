# Load R packages for analysis and visualizations
library(Seurat)
library(dplyr)
library(patchwork)

# Set directory
setwd('path')


# Create Seurat Objects for Disease (for two samples)
  ## Disease group first sample

  disease_1_counts <- read.csv('path', header = TRUE, row.names = 1)
  disease_1_seurat <- CreateSeuratObject(counts = disease_1_counts)

    ### Assign metadata to Seurat Object

    disease_1_seurat$metadata_1 <- 'metadata1'
    disease_1_seurat$metadata_2 <- 'metadata2'

  ## Disease group second sample
  
  disease_2_counts <- read.csv('path', header = TRUE, row.names = 1)
  disease_2_seurat <- CreateSeuratObject(counts = disease_2_counts)
  
    ### Assign metadata to Seurat Object
  
    disease_2_seurat$metadata_1 <- 'metadata1'
    disease_2_seurat$metadata_2 <- 'metadata2'

    
# Create Seurat Objects for Control (for two samples)
  ## Control group first sample
    
  control_1_counts <- read.csv('path', header = TRUE, row.names = 1)
  control_1_seurat <- CreateSeuratObject(counts = control_1_counts)
    
    ### Assign metadata to Seurat Object
    
    control_1_seurat$metadata_1 <- 'metadata1'
    control_1_seurat$metadata_2 <- 'metadata2'
    
  ## Disease group second sample
    
  control_2_counts <- read.csv('path', header = TRUE, row.names = 1)
  control_2_seurat <- CreateSeuratObject(counts = control_2_counts)
    
    ### Assign metadata to Seurat Object
    
    control_2_seurat$metadata_1 <- 'metadata1'
    control_2_seurat$metadata_2 <- 'metadata2'

    
# Merge all Seurat Objects
merged_seurat <- merge(disease_1_seurat, y = c(disease_2_seurat, control_1_seurat, control_2_seurat), add.cell.ids = c('disease_1', 'disease_2', 'control_1', 'control_2'), project = 'projectname')

# Save the merged Seurat Object as *.rds file to the current directory
SaveSeuratRds(merged_seurat, file = 'merged_seurat.rds', compress = TRUE)

# Filter the genes by Mitochondrial gene percentage, number of counts, and number of features
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')

filtered_seurat <- subset(merged_seurat, subset = nCount_RNA > 800 &
                                  nFeature_RNA > 200 & nFeature_RNA < 6000 &
                                  mitoPercent < 5)

# Normalize the filtered Seurat Object
normalized_seurat <- NormalizeData(filtered_seurat)

    # Find top 10 variable features of Seurat Object and plot
      vf_seurat <- FindVariableFeatures(normalized_seurat, selection.method = "vst", nfeatures = 1000)
      vf_top10 <- head(VariableFeatures(vf_seurat), 10)
      
        plot1 <- VariableFeaturePlot(VF_MergedSeurat)
        plot2 <- LabelPoints(
          plot = plot1, 
          points = VF_Top10, 
          labels = VF_Top10, 
          repel = FALSE)
      ## plot1 + plot2

  
# Scale the Seurat Object using top thousand variable features
all_genes <- rownames(normalized_seurat)
scaled_seurat <- ScaleData(normalized_seurat, features = all_genes)
vf_scaled_seurat <- FindVariableFeatures(scaled_seurat, selection.method = 'vst', nfeatures = 1000)
vf_scaled_seurat_1000 <- head(VariableFeatures(vf_scaled_seurat), 1000)

# Run PCA using the top thousand features
pca <- RunPCA(scaled_seurat, features = vf_scaled_seurat_1000)

# Save the PCA as Seurat Object as *.rds file to the current directory
SaveSeuratRds(pca, file = 'pca.rds')

# Visualize the elbow plot for PCA to determine the dimensions
ElbowPlot(pca)

pca_neighboors <- FindNeighbors(pca, dims = 1:15)
pca_neighboors_clusters <- FindClusters(pca_neighboors, resolution = 1)

# Run and save UMAP
umap <- RunUMAP(pca_neighboors_clusters, dims = 1:15)
SaveSeuratRds(umap, file = 'umap.rds')

# Visualize the UMAP
DimPlot(umap, reduction = 'umap', label = TRUE)
DimPlot(umap, reduction = 'umap', group.by = 'metadata_1', label = TRUE)
DimPlot(umap, reduction = 'umap', split.by = 'metadata_1', label = TRUE)

# Join the layers of the UMAP
Joined_UMAP <- JoinLayers(UMAP)