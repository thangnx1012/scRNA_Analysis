setwd("/media/hkh/8TB/XUANTHANG/MammaryGland/scRNA_ER/Results")
# Load libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)
library(scater)
library(stringr)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(EnsDb.Mmusculus.v79)
library(tidyverse)
library(RCurl)
library(cowplot)
##load only filtered dataset
load(file="R/MammaryGland_seurat_filtered.RData")
#Clustering workflow
# 1.Normalization, variance stabilization, and regression of unwanted variation (e.g. mitochondrial transcript abundance, cell cycle phase, etc.) for each sample
# 2.Integration of the samples using shared highly variable genes (optional, but recommended to align cells from different samples/conditions if cell types are separating by sample/condition)
# 3.Clustering cells based on top PCs (metagenes)
# 4.Exploration of quality control metrics: determine whether clusters are unbalanced wrt UMIs, genes, cell cycle, mitochondrial content, samples, etc.
# 5.Searching for expected cell types using known cell type-specific gene markers


## 4.Clustering workflow---- 
## 4.1 Cell cycle scoring----
## It is recommended to check the cell cycle phase before performing the sctransform method
# Normalize as defaults parameters
MammaryGland_filtered_phase <- NormalizeData(MammaryGland_filtered)
# Load cell cycle markers: g2m_genes_Mouse and s_genes_Mouse
load("/media/hkh/8TB/XUANTHANG/References/Cell_Cycle_marker/Mouse_Gene_CellCycle.rda")
# Score cells for cell cycle
#If set.ident =TRUE, the Identity will provide the information with cell cycle phase, other will visualize the experiment sample (ctrl/E2)
MammaryGland_filtered_phase <- CellCycleScoring(MammaryGland_filtered_phase, 
                                                g2m.features = g2m_genes_Mouse, 
                                                s.features = s_genes_Mouse,
                                                set.ident = TRUE)
# View cell cycle scores and phases assigned to cells                                 
head(MammaryGland_filtered_phase[[]]) #or View(MammaryGland_filtered_phase@meta.data)   
# Visualize the distribution of cell cycle markers across
RidgePlot(MammaryGland_filtered_phase, 
          features = c("Pcna", "Krt18", "Ccl2", "Ccl7"), 
          ncol = 2)

# Identify the most variable genes
MammaryGland_filtered_phase <- FindVariableFeatures(MammaryGland_filtered_phase, 
                                                    selection.method = "vst",
                                                    verbose = FALSE)
# Scale the counts
MammaryGland_filtered_phase <- ScaleData(MammaryGland_filtered_phase)
# Perform PCA
MammaryGland_filtered_phase <- RunPCA(MammaryGland_filtered_phase) 
                                      # features = c(s.genes_Mouse, g2m.genes_Mouse))
# Plot the PCA colored by cell cycle phase
DimPlot(MammaryGland_filtered_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
# We do not see large differences due to cell cycle phase. Based on this plot, we would not regress out the variation due to cell cycle.
# Plot UMAP
MammaryGland_filtered_phase <- RunUMAP(MammaryGland_filtered_phase, 
                                       dims = 1:40, 
                                       reduction = "pca")
DimPlot(MammaryGland_filtered_phase) 
DimPlot(MammaryGland_filtered_phase, group.by = "sample")
DimPlot(MammaryGland_filtered_phase, split.by = "sample") #split.by = "Phase"
# => Data should be intergrated to perform analysis

# 4.2 SCTransform and Intergrate----
# sctransform method as a more accurate method of normalizing, estimating the variance of the raw filtered data, 
# and identifying the most variable genes. By default, sctransform accounts for cellular sequencing depth, or nUMIs.

# adjust the limit for allowable object sizes within R (Default is 500 * 1024 ^ 2 = 500 Mb)
options(future.globals.maxSize = 4000 * 1024^2)
# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
MammaryGland_split <- SplitObject(MammaryGland_filtered, split.by = "sample")
MammaryGland_split <- MammaryGland_split[c("Ctrl", "E2")]
# SCTransform will rank the genes by residual variance and output the 3000 most variant genes by default
for (i in 1:length(MammaryGland_split)) {
  MammaryGland_split[[i]] <- NormalizeData(MammaryGland_split[[i]], verbose = TRUE)
  MammaryGland_split[[i]] <- CellCycleScoring(MammaryGland_split[[i]], g2m.features=g2m_genes_Mouse, s.features=s_genes_Mouse)
  MammaryGland_split[[i]] <- SCTransform(MammaryGland_split[[i]], vars.to.regress = c("mitoRatio"))
}
# Check which assays are stored in objects
MammaryGland_split$Ctrl@assays



## Integrate samples using shared highly variable genes
# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = MammaryGland_split, 
                                            nfeatures = 3000) 
# Prepare the SCTransform  list object for integration
MammaryGland_split <- PrepSCTIntegration(object.list = MammaryGland_split, 
                                   anchor.features = integ_features)
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = MammaryGland_split, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
# Integrate across conditions
MammaryGland_integrated <- IntegrateData(anchorset = integ_anchors, 
                                         normalization.method = "SCT")
# Save integrated seurat object
saveRDS(MammaryGland_integrated, "R/MammaryGland_integrated_seurat.rds")
save(MammaryGland_integrated, file="R/MammaryGland_integrated_seurat.RData")

# UMAP visualization
# Run PCA and Plot PCA
MammaryGland_integrated <- RunPCA(object = MammaryGland_integrated)
PCAPlot(MammaryGland_integrated, split.by = "sample")  
# Run UMAP and Plot UMAP   
MammaryGland_integrated <- RunUMAP(MammaryGland_integrated, dims = 1:40, reduction = "pca")
DimPlot(MammaryGland_integrated)                             
DimPlot(MammaryGland_integrated, split.by = "sample") 

## 4.3 Clustering cells based on top PCs (metagenes)----
# Explore heatmap of PCs
DimHeatmap(MammaryGland_integrated, 
           dims = 1:10, 
           cells = 500, 
           balanced = TRUE)
# Printing out the most variable genes driving PCs
print(x = MammaryGland_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

# Plot the elbow plot regarding as helpful way to determine how many PCs to use for clustering 
ElbowPlot(object = MammaryGland_integrated, ndims = 40) #change the number of dimension to fit the graph
# Choose the number of dimension is 30

# Determine the K-nearest neighbor graph
MammaryGland_integrated <- FindNeighbors(object = MammaryGland_integrated, 
                                   dims = 1:30)
# Determine the clusters for various resolutions                                
MammaryGland_integrated <- FindClusters(object = MammaryGland_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
# Explore resolutions
MammaryGland_integrated@meta.data %>% View()
# Assign identity of clusters
Idents(object = MammaryGland_integrated) <- "integrated_snn_res.0.8" #often pick something in the middle of the range like 0.6 or 0.8.
# Plot the UMAP
DimPlot(MammaryGland_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

# Visualization both UMAP resolusion
p1 <- DimPlot(MammaryGland_integrated, reduction = "umap", group.by = "sample")
p2 <- DimPlot(MammaryGland_integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(MammaryGland_integrated, reduction = "umap", split.by = "sample") # visualize the two conditions side-by-side

##4.4 Exploration of quality control metrics----
# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(MammaryGland_integrated, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)
View(n_cells) # View table
#Segregation of clusters by sample
# UMAP of cells in each cluster by sample
DimPlot(MammaryGland_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()
# Explore whether clusters segregate by cell cycle phase
DimPlot(MammaryGland_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()
## We do not see much clustering by cell cycle score, so we can proceed with the QC.

#Segregation of clusters by various sources of uninteresting variation
# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")
FeaturePlot(MammaryGland_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
# Sample clustering has been exhibition at group of cluster 0,1,2,4,5,7,9 and 8,11,15


### Exploration of the PCs driving the different clusters
# Defining the information in the seurat object of interest with 15 PC
columns <- c(paste0("PC_", 1:15),
             "ident",
             "UMAP_1", "UMAP_2")


# Extracting this data from the seurat object
pc_data <- FetchData(MammaryGland_integrated, vars = columns)
# Let's take a quick look at the top 16 PCs:
# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(MammaryGland_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:15), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)
# Then, Examine PCA results by gene list driving this PC
print(MammaryGland_integrated[["pca"]], dims = 1:5, nfeatures = 5)


## 4.5 Exploring known cell type markers


# Select the RNA counts slot to be the default assay
DefaultAssay(MammaryGland_integrated) <- "RNA"
# Normalize RNA data for visualization purposes
MammaryGland_integrated <- NormalizeData(MammaryGland_integrated, verbose = FALSE)

#The FeaturePlot() function from seurat makes it easy to visualize a handful of 
#genes using the gene IDs stored in the Seurat object. For example if we were 
#interested in exploring known immune cell markers, such as:
#Macrophage marker
FeaturePlot(MammaryGland_integrated, 
            reduction = "umap", 
            features = c("Marco", "Itgam", "Adgre1"),
            sort.cell = TRUE,
            min.cutoff = 'q10', cols = c("grey", "red"),
            label = TRUE)
#B-cell marker
FeaturePlot(MammaryGland_integrated, 
            reduction = "umap", 
            features = c("Cd79a", "Ms4a1"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', cols = c("grey", "red"),
            label = TRUE)
#T-cell marker #8
FeaturePlot(MammaryGland_integrated, 
            reduction = "umap", 
            features = c("Cd3d"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', cols = c("grey", "red"),
            label = TRUE)














