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
library(AnnotationHub)
library(org.Mm.eg.db)
library(EnsDb.Mmusculus.v79)
library(tidyverse)
library(RCurl)
library(multtest)
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
          features = c("Col1a1", "Krt18", "Ptprc", "Epcam"), 
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
DimPlot(MammaryGland_filtered_phase, group.by = "sample")
DimPlot(MammaryGland_filtered_phase, group.by = "batch")
DimPlot(MammaryGland_filtered_phase, split.by = "sample") #split.by = "Phase" "batch"
# => Data should be intergrated to perform analysis

# 4.2 SCTransform and Intergrate----
# sctransform method as a more accurate method of normalizing, estimating the variance of the raw filtered data, 
# and identifying the most variable genes. By default, sctransform accounts for cellular sequencing depth, or nUMIs.

# adjust the limit for allowable object sizes within R (Default is 500 * 1024 ^ 2 = 500 Mb)
options(future.globals.maxSize = 4000 * 1024^2)
# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
MammaryGland_split <- SplitObject(MammaryGland_filtered, split.by = "batch")
MammaryGland_split <- MammaryGland_split[c("Rep1", "Rep2")]
# SCTransform will rank the genes by residual variance and output the 3000 most variant genes by default
for (i in 1:length(MammaryGland_split)) {
  MammaryGland_split[[i]] <- NormalizeData(MammaryGland_split[[i]], verbose = TRUE)
  MammaryGland_split[[i]] <- CellCycleScoring(MammaryGland_split[[i]], g2m.features=g2m_genes_Mouse, s.features=s_genes_Mouse)
  MammaryGland_split[[i]] <- SCTransform(MammaryGland_split[[i]], vars.to.regress = c("mitoRatio"))
}
# Check which assays are stored in objects
MammaryGland_split$Rep1@assays

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
saveRDS(MammaryGland_integrated, "R/MammaryGland_integrated_seurat_batchCorrection.rds")
save(MammaryGland_integrated, file="R/MammaryGland_integrated_seurat_batchCorrection.RData")

# UMAP visualization
# Run PCA and Plot PCA
MammaryGland_integrated <- RunPCA(object = MammaryGland_integrated)
PCAPlot(MammaryGland_integrated, split.by = "sample")  
PCAPlot(MammaryGland_integrated, split.by = "batch")  
# Run UMAP and Plot UMAP   
MammaryGland_integrated <- RunUMAP(MammaryGland_integrated, dims = 1:40, reduction = "pca")
DimPlot(MammaryGland_integrated)                             
DimPlot(MammaryGland_integrated, split.by = "sample") 
DimPlot(MammaryGland_integrated, split.by = "batch") 

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
# Elbow plot: quantitative approach
ElbowPlot(object = MammaryGland_integrated, ndims = 40) #change the number of dimension to fit the graph
# Choose follow on graph or calculate by formular below
####################################################################################################
# Determine percent of variation associated with each PC
pct <- MammaryGland_integrated[["pca"]]@stdev / sum(MammaryGland_integrated[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# The point where the percent change in variation between the consecutive PCs is less than 0.1%.
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# Minimum of the two calculation
pcs <- min(co1, co2)
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
#Choose dimension 22 for other calculator
####################################################################################################

# Determine the K-nearest neighbor graph
MammaryGland_integrated <- FindNeighbors(object = MammaryGland_integrated, 
                                   dims = 1:22)
# Determine the clusters for various resolutions                                
MammaryGland_integrated <- FindClusters(object = MammaryGland_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
# Assign identity of clusters: #often pick something in the middle of the range like 0.6 or 0.8. whereas seurat chose 0.5
Idents(object = MammaryGland_integrated) <- "integrated_snn_res.0.4" 
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
DimPlot(MammaryGland_integrated, reduction = "umap", split.by = "batch")
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
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
# Sample clustering has been exhibition at group of cluster 0,1,4,5,6 and 7,8,12


### Exploration of the PCs driving the different clusters
# Defining the information in the seurat object of interest with 15 PC
columns <- c(paste0("PC_", 1:12),
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
map(paste0("PC_", 1:12), function(pc){
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


## 4.5 Exploring known cell type markers----
# Select the RNA counts slot to be the default assay
DefaultAssay(MammaryGland_integrated) <- "RNA"
# Normalize RNA data for visualization purposes
MammaryGland_integrated <- NormalizeData(MammaryGland_integrated, verbose = FALSE)

#The FeaturePlot() function from seurat makes it easy to visualize a handful of 
#genes using the gene IDs stored in the Seurat object. For example if we were 
#interested in exploring known immune cell markers, such as:

#Immune_marker #7,8,12
FeaturePlot(MammaryGland_integrated, 
            reduction = "umap", 
            features = c("Ptprc"), 
            order = TRUE,
            min.cutoff = 'q10', cols = c("grey", "red"),
            label = TRUE)
#Luminal marker #2,3,11
FeaturePlot(MammaryGland_integrated, 
            reduction = "umap", 
            features = c("Krt18", "Krt19",  "Epcam"), 
            order = TRUE,
            min.cutoff = 'q10', cols = c("grey", "red"),
            label = TRUE)
#Basal marker Krt14 #9 and Endothelial cell marker Pecam1 #10
FeaturePlot(MammaryGland_integrated, 
            reduction = "umap", 
            features = c("Krt14", "Pecam1"), 
            order = TRUE,
            min.cutoff = 'q10', cols = c("grey", "red"),
            label = TRUE)
#ECM/fibroblasts 0,1,5,4,6,
FeaturePlot(MammaryGland_integrated, 
            reduction = "umap", 
            features = c("Col1a1", "Col1a2", "Col6a1", "Mmp2"), 
            order = TRUE,
            min.cutoff = 'q10', cols = c("grey", "red"),
            label = TRUE)
#Check Esr1 expression: Luminal endothelial cell and ECM/Fibroblast
FeaturePlot(MammaryGland_integrated, 
            reduction = "umap", 
            features = c("Col4a2"), 
            order = TRUE,
            min.cutoff = 'q10', cols = c("grey", "red"),
            label = TRUE)

# 4.6 Single-cell RNA-seq marker identification----

# 4.6.1:  Evaluating marker genes using converse
# Identification of conserved markers in all conditions
DefaultAssay(MammaryGland_integrated) <- "RNA"
# Call for annotations file (that created from annotationHub)
annotations <- read.csv("/media/hkh/8TB/XUANTHANG/References/Annotation_mitogene/annotation_Mouse.csv")
# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(MammaryGland_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}
# Running on multiple samples (0:12) or function across desired clusters c(2,3,7)
conserved_markers <- map_dfr(c(0:10), get_conserved) #cluster11 not enough cell number
# some cases you will have clusters that do not have enough cells for a particular group - and your function will fail. 
# For these clusters you will need to use FindAllMarkers().
# Then check list on genecards.org to indentify cell clustering

## Evaluating marker genes
# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (Ctrl_avg_log2FC + E2_avg_log2FC + PDBE_avg_log2FC) /3) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)
# Visualize top 10 markers per cluster
View(top10)
# Plot interesting marker gene expression for cluster 20
FeaturePlot(object = MammaryGland_integrated, 
            features = c("Ackr3", "Ifi27l2a", "Sfrp4", "Pi16", "C3"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
# Vln plot - also explore the range in expression of specific markers by using violin plots:
VlnPlot(object = MammaryGland_integrated, 
        features = c("Ackr3", "Ifi27l2a", "Sfrp4", "Pi16", "C3"))

FindAllMarkers(MammaryGland_integrated)

# Rename all identities
MammaryGland_integrated <- RenameIdents(object = MammaryGland_integrated, 
                                  "0" = "ECM/fibroblasts cells",
                                  "1" = "ECM/fibroblasts cells",
                                  "2" = "Luminal cells",
                                  "3" = "Luminal cells",
                                  "4" = "ECM/fibroblasts cells",
                                  "5" = "ECM/fibroblasts cells",
                                  "6" = "ECM/fibroblasts cells",
                                  "7" = "Immune cells",
                                  "8" = "Immune cells",
                                  "9" = "Basal cell",
                                  "10" = "Endothelial cells",
                                  "11" = "Luminal cells",
                                  "12" = "Immune cells")


# Plot the UMAP
DimPlot(object = MammaryGland_integrated, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)
# The DotPlot function with the split.by parameter can be useful for viewing conserved cell type markers across conditions
markers.to.plot <- c("Ackr3", "Ifi27l2a", "Sfrp4", "Pi16", "C3")
DotPlot(MammaryGland_integrated, features = rev(markers.to.plot), cols = c("blue", "red", "yellow"), dot.scale = 8, 
        split.by = "sample") + RotatedAxis()



# If we wanted to remove the potentially stressed cells, we could use the subset() function:
# Remove the Immune cells
seurat_subset_labeled <- subset(MammaryGland_integrated,
                                idents = c("ECM/fibroblasts cells", "Immune cells"), invert = TRUE)

# Re-visualize the clusters
DimPlot(object = seurat_subset_labeled, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)

# Save final R object
write_rds(MammaryGland_integrated,
          file = "R/MammaryGland_integrated_labelled.rds")       
read_rds(file = "R/MammaryGland_integrated_labelled.rds")
save(MammaryGland_integrated, file="R/MammaryGland_integrated_lableled.RData")

