##################### QC and filtering #########################
library(Seurat)
library(dplyr)
library(ggplot2)
library(scDblFinder)
library(BiocParallel)
library(hdf5r)
options(Seurat.object.assay.version = 'v5')

######## Doublet Finder ##########
# Using scDblFinder to remove doublets
# script is run for each sample individually
# use h5 files after running cellbender

files<-list.files("/mnt/parscratch/users/mdq19mm/nucseq/combined_runs/raw_data/cellbender/",full.names = T)
working_directory <- "/mnt/parscratch/users/mdq19mm/nucseq/combined_runs/doublet_finder"
setwd(working_directory)

# Load data
sample <- "n2"
filepath<-files[grep(toupper(sample),files)] #find filename for specific sample and make seurat object
data<- Read10X_h5(filepath)
n2<- CreateSeuratObject(counts = data, project = "neuroma", min.cells = 3, min.features = 200)
n2 <- n2[!grepl("^DEPRECATED-", rownames(n2)), ] #removing deprecated probes
counts <- GetAssayData(n2, slot = "counts")

# Running scDbltFinder
set.seed(123)
sce <- scDblFinder(counts)

# Adding scDblFinder results to sample metadata and save object
doublet_score <- sce$scDblFinder.score
doublet_class <- sce$scDblFinder.class
n2$doublet_score <- doublet_score
n2$doublet_class <- doublet_class

saveRDS(n2, "n2_scDblFinder.rds")

########### Filtering ###########
# Load Seurat object with scDblFinder results
doublet_file <- "/mnt/parscratch/users/mdq19mm/nucseq/combined_runs/doublet_finder/n2_scDblFinder.rds"
working_directory <- "/mnt/parscratch/users/mdq19mm/nucseq/combined_runs/integration/"
setwd(working_directory)
images <- "n2"

n2 <- readRDS(doublet_file)
features <- length(rownames(n2))

# remove deprecated genes if wasn't done before
n2 <- n2[!grepl("^DEPRECATED-", rownames(n2)), ]
deprecated.genes <- features - length(rownames(n2))
print(paste0("Number of deprecated genes removed:",deprecated.genes))
print(n2)

# Adding QC metrics to metadata
n2 <- PercentageFeatureSet(n2, pattern = "^MT-", col.name = "percent.mt") # mito genes
n2$log10GenesPerUMI <- log10(n2$nFeature_RNA) / log10(n2$nCount_RNA) # novelty score
metadata <- n2@meta.data # get metadata data frame for QC visualisation
metadata$cells <- rownames(metadata)

# Visualize the number UMIs/transcripts per cell
plot3<- metadata %>% 
  ggplot(aes(color=doublet_class, x=nCount_RNA, fill= doublet_class)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell
plot4<- metadata %>% 
  ggplot(aes(color=doublet_class, x=nFeature_RNA, fill= doublet_class)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 400)

# Visualize the distribution of mitochondrial gene expression detected per cell
plot5 <- metadata %>% 
  ggplot(aes(color=doublet_class, x=percent.mt, fill=doublet_class)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 3)

QC3<- plot3+plot4+plot5
ggsave("doublets_vs_singlets_QC_n2.pdf",plot=QC3,path=images, width=12, height=5)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
QC4<-metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, colour=doublet_score)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 400) +
  facet_wrap(~orig.ident)
ggsave("genes_per_UMI_n2.pdf",plot=QC4,path=images)

# Subset singlet only
n2 <- subset(n2,doublet_class=="singlet")
print(n2)

# Visualize QC metrics as a violin plot
QC<-VlnPlot(n2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("QC_singlets_n2.pdf",plot=QC,path=images)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(n2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(n2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
QC2<-plot1 + plot2
ggsave("QC_featureScatter_singlets_n2.pdf",plot=QC2,path=images, width=10, height=5)

# filter according to preset pareameters
n2 <- subset(x = n2, subset= (nCount_RNA >= 500) & (nFeature_RNA >= 250) & (log10GenesPerUMI > 0.80) & (percent.mt < 5))
print(n2)

QC3<-VlnPlot(n2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("QC_filtered_n2.pdf",plot=QC3,path=images)

####### Normalisation #########
## run sctransform
n2 <- SCTransform(n2, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

## Dimensionality reduction
n2 <- RunPCA(n2, verbose = FALSE)
n2 <- RunUMAP(n2, dims = 1:30, verbose = FALSE)
n2 <- FindNeighbors(n2, dims = 1:30, verbose = FALSE)
n2 <- FindClusters(n2, verbose = FALSE, resolution=0.5)
p1<-DimPlot(n2, label = TRUE) + NoLegend()
ggsave("n2_clusters_UMAP_res05.pdf",plot=p1,path=images)

saveRDS(n2,"n2.rds")


