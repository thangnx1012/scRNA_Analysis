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
## 1 Load the single cell dataset----
## Create each individual Seurat object for each sample by loop
for (file in c("SRR8452054Gene/filtered", 
               "SRR8452055Gene/filtered",
               "SRR8452056Gene/filtered",
               "SRR8452057Gene/filtered",
               "SRR8452058Gene/filtered",
               "SRR8452059Gene/filtered")){
  seurat_data <- Read10X(data.dir = paste0("STARsolo/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, min.cells = 3, min.features = 100, project = file)
  assign(file, seurat_obj)}

# Check the metadata in the new Seurat objects
head(`SRR8452054Gene/filtered`@meta.data)

## Merging All Seurat Objects (Only E2 and control)
MammaryGland <- merge(x = `SRR8452054Gene/filtered`, 
                      y = c(`SRR8452057Gene/filtered`, `SRR8452055Gene/filtered`, `SRR8452058Gene/filtered`),
                      add.cell.ids = c("Ctrl_1", "Ctrl_2", "E2_1", "E2_2"), 
                      project = "MammaryGland")
MammaryGland
# Check that the merged object has the appropriate sample-specific prefixes
head(colnames(MammaryGland))
tail(colnames(MammaryGland))
# unique(sapply(X = strsplit(colnames(MammaryGland), split = "_"), FUN = "[", 1))
table(MammaryGland$orig.ident)


## 2 Generating Quality control metrics----
## Add number of genes per UMI for each cell to metadata
MammaryGland$log10GenesPerUMI <- log10(MammaryGland$nFeature_RNA) / log10(MammaryGland$nCount_RNA)
# Compute percent mito ratio
MammaryGland$mitoRatio <- PercentageFeatureSet(object = MammaryGland, pattern = "^mt-")
MammaryGland$mitoRatio <- MammaryGland@meta.data$mitoRatio / 100
head(MammaryGland@meta.data, 5)

## Create the metadata dataframe by extracting the meta.data slot from the Seurat object:
metadata <- MammaryGland@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>% dplyr::rename(seq_folder = orig.ident, 
                                       nUMI = nCount_RNA, 
                                       nGene = nFeature_RNA)
# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^Ctrl_"))] <- "Ctrl"
metadata$sample[which(str_detect(metadata$cells, "^E2_"))] <- "E2"
# Add metadata back to Seurat object
MammaryGland@meta.data <- metadata
# Create .RData that was update Seurat object to load at any time
save(MammaryGland, metadata,  file="R/MammaryGland_update_seurat.RData")


## 3. Assessing the quality metrics
load(file="R/MammaryGland_update_seurat.RData")
## Cell counts
# Visualize the number of cell counts per sample
# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells")

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells vs Number of Genes")

##UMIs vs. genes detected
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# Visualize the distribution of mitochondrial gene expression detected per cell
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI: Expected 0.80
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.80)

# Visualize QC metrics as a violin plot
VlnPlot(MammaryGland, features = c("nUMI", "nGene", "mitoRatio"), ncol = 3)

# Filtering
# Filter out low quality reads using selected thresholds - these will change with experiment
MammaryGland_filtered <- subset(x = MammaryGland, 
                                subset= (nUMI >= 500) & (nGene >= 500) & 
                                  (log10GenesPerUMI > 0.80) & (mitoRatio < 0.20))
# Gene-level filtering
# Extract counts
counts <- GetAssayData(object = MammaryGland_filtered, slot = "counts")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
# Only keeping those genes expressed in more than 10 cells in filtered Seurat object
MammaryGland_filtered <- CreateSeuratObject(counts[keep_genes, ], meta.data = MammaryGland_filtered@meta.data)
# Visualize QC metrics as a violin plot
VlnPlot(MammaryGland_filtered, features = c("nUMI", "nGene", "mitoRatio"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships
plot1 <- FeatureScatter(MammaryGland_filtered, feature1 = "nUMI", feature2 = "mitoRatio")
plot2 <- FeatureScatter(MammaryGland_filtered, feature1 = "nUMI", feature2 = "nGene")
plot1 + plot2
# Create filtered .RData object to load at any time
save(MammaryGland_filtered, file="R/MammaryGland_seurat_filtered.RData")
