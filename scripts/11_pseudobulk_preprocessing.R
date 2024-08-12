######################## Pseudobulk DE analysis with DESeq2 ########################
############################## All samples Visium 2023 #############################
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleCellExperiment)
library(Matrix)
options(future.globals.maxSize = 1e9)

# this script is to prepare the pseudo-bulk RNA seq dataset from visium samples
# all visium sections are merged in a seurat object 
# barcodes overlaying nerve fascicles (identifying on Loupe) are subset
# aggregated count matrix is obtained

working_directory <- "/mnt/parscratch/users/mdq19mm/visium_analysis/seurat"
setwd(working_directory)

filedir <- "/mnt/parscratch/users/mdq19mm/visium_data"

##### Create list of samples ####
# Loading samples
sample.ids<- list.dirs(filedir,full.names = FALSE, recursive = FALSE)
sample.paths<- paste(filedir,sample.ids,"outs",sep="/")
area <- c("0","A1","B1","C1","D1")
spatial <- list()

############ LN1 ##############
sample<-"LN1"
pain.status<-"P"
sex<-"F"

for (i in 1:5) {
  spatial[[i]] <- Load10X_Spatial(sample.paths[i],filename="filtered_feature_bc_matrix.h5", assay="Spatial",slice=sample.ids[i],filter.matrix=TRUE)
  spatial[[i]]@images[[1]]@coordinates <- spatial[[i]]@images[[1]]@coordinates %>% mutate_all(function(x) as.numeric(as.character(x)))
  spatial[[i]]$orig.ident <- sample.ids[[i]]
  spatial[[i]]$pain.status <- pain.status
  spatial[[i]]$sex <- sex
  spatial[[i]]$sample<- sample
}
for (i in 1:5) {spatial[[i]]$capture.area<- area[i]}

############ LN2 ############
sample<-"LN2"
pain.status<-"NP"
sex<-"M"

for (i in 14:18) {
  spatial[[i]] <- Load10X_Spatial(sample.paths[i],filename="filtered_feature_bc_matrix.h5", assay="Spatial",slice=sample.ids[i],filter.matrix=TRUE)
  spatial[[i]]@images[[1]]@coordinates <- spatial[[i]]@images[[1]]@coordinates %>% mutate_all(function(x) as.numeric(as.character(x)))
  spatial[[i]]$orig.ident <- sample.ids[[i]]
  spatial[[i]]$pain.status <- pain.status
  spatial[[i]]$sex <- sex
  spatial[[i]]$sample<- sample
}
for (i in 1:5) {spatial[[i+13]]$capture.area<- area[i]}

############ LN7 ############
sample<-"LN7"
pain.status<-"NP"
sex<-"F"

for (i in 19:23) {
  spatial[[i]] <- Load10X_Spatial(sample.paths[i],filename="filtered_feature_bc_matrix.h5", assay="Spatial",slice=sample.ids[i],filter.matrix=TRUE)
  spatial[[i]]@images[[1]]@coordinates <- spatial[[i]]@images[[1]]@coordinates %>% mutate_all(function(x) as.numeric(as.character(x)))
  spatial[[i]]$orig.ident <- sample.ids[[i]]
  spatial[[i]]$pain.status <- pain.status
  spatial[[i]]$sex <- sex
  spatial[[i]]$sample<- sample
}
for (i in 1:5) {spatial[[i+18]]$capture.area<- area[i]}

############ LN8 ############
sample<-"LN8"
pain.status<-"NP"
sex<-"M"

for (i in 24:28) {
  spatial[[i]] <- Load10X_Spatial(sample.paths[i],filename="filtered_feature_bc_matrix.h5", assay="Spatial",slice=sample.ids[i],filter.matrix=TRUE)
  spatial[[i]]@images[[1]]@coordinates <- spatial[[i]]@images[[1]]@coordinates %>% mutate_all(function(x) as.numeric(as.character(x)))
  spatial[[i]]$orig.ident <- sample.ids[[i]]
  spatial[[i]]$pain.status <- pain.status
  spatial[[i]]$sex <- sex
  spatial[[i]]$sample<- sample
}
for (i in 1:5) {spatial[[i+23]]$capture.area<- area[i]}

############ LN12 ############
sample<-"LN12"
pain.status<-"P"
sex<-"F"
area <- c("A1","B1")

for (i in 6:7) {
  spatial[[i]] <- Load10X_Spatial(sample.paths[i],filename="filtered_feature_bc_matrix.h5", assay="Spatial",slice=sample.ids[i],filter.matrix=TRUE)
  spatial[[i]]@images[[1]]@coordinates <- spatial[[i]]@images[[1]]@coordinates %>% mutate_all(function(x) as.numeric(as.character(x)))
  spatial[[i]]$orig.ident <- sample.ids[[i]]
  spatial[[i]]$pain.status <- pain.status
  spatial[[i]]$sex <- sex
  spatial[[i]]$sample<- sample
}
for (i in 1:2) {spatial[[i+5]]$capture.area<- area[i]}

############ LN13 ############
sample<-"LN13"
pain.status<-"P"
sex<-"F"
area <- c("C1","D1")

for (i in 8:9) {
  spatial[[i]] <- Load10X_Spatial(sample.paths[i],filename="filtered_feature_bc_matrix.h5", assay="Spatial",slice=sample.ids[i],filter.matrix=TRUE)
  spatial[[i]]@images[[1]]@coordinates <- spatial[[i]]@images[[1]]@coordinates %>% mutate_all(function(x) as.numeric(as.character(x)))
  spatial[[i]]$orig.ident <- sample.ids[[i]]
  spatial[[i]]$pain.status <- pain.status
  spatial[[i]]$sex <- sex
  spatial[[i]]$sample<- sample
}
for (i in 1:2) {spatial[[i+7]]$capture.area<- area[i]}

############ LN15 ############
sample<-"LN15"
pain.status<-"P"
sex<-"F"
area <- c("A1","B1","C1","D1")

for (i in 10:13) {
  spatial[[i]] <- Load10X_Spatial(sample.paths[i],filename="filtered_feature_bc_matrix.h5", assay="Spatial",slice=sample.ids[i],filter.matrix=TRUE)
  spatial[[i]]@images[[1]]@coordinates <- spatial[[i]]@images[[1]]@coordinates %>% mutate_all(function(x) as.numeric(as.character(x)))
  spatial[[i]]$orig.ident <- sample.ids[[i]]
  spatial[[i]]$pain.status <- pain.status
  spatial[[i]]$sex <- sex
  spatial[[i]]$sample<- sample
}
for (i in 1:4) {spatial[[i+9]]$capture.area<- area[i]}

######### SCTransform for normalisation ##########
spatial <- lapply(spatial,SCTransform, assay = "Spatial", verbose = FALSE)

############# QC plots ############
images<- "QC/"
for (i in 1:length(sample.ids)) {
  plot1 <- VlnPlot(spatial[[i]], features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  plot2 <- SpatialFeaturePlot(spatial[[i]], features = "nCount_Spatial") + theme(legend.position = "right")
  plotQC<-plot1+plot2
  plot3 <- SpatialFeaturePlot(spatial[[i]], features = c("SCN9A","MPZ","NRXN1","SLC2A1","CD68","PTPRC","PI16","EGFL7","TNNT1","KRT7","PDGFRB","ACTA2"),ncol = 4)
  ggsave(paste0(sample.ids[i],"_QC.pdf"),p=plotQC,path = images, width = 5, height = 10, limitsize = FALSE)
  ggsave(paste0(sample.ids[i],"_marker_genes.pdf"),p=plot3,path = images, width = 20, height = 15, limitsize = FALSE)
}

############# save spatial combined object ############
saveRDS(spatial, "/mnt/parscratch/users/mdq19mm/visium_analysis/spatial_combined.rds")

############# Adding metadata with nerve classification from Loupe ############
loupe.csv <- "/mnt/parscratch/users/mdq19mm/visium_analysis/DEseq2/loupe_files/"
for (i in 1:length(sample.ids)) {
  nerve.fascicles <- read.csv(paste0(loupe.csv,"Nerve_fascicles_",sample.ids[i],".csv"),row.names = 1) %>% filter(grepl('nerves',Nerve.fascicles))
  spatial[[i]] <- subset(spatial[[i]], cells=row.names(nerve.fascicles))
  spatial[[i]]$loupe.annotation <- nerve.fascicles
}

saveRDS(spatial,"/mnt/parscratch/users/mdq19mm/visium_analysis/spatial_nerve_fascicles.rds")

############# QC plots for nerve fascicles ############
images<- "QC_nerve_fascicles/"

genes <- c("SCN9A","MPZ","NRXN1","SLC2A1","CD68","PTPRC","PI16","EGFL7","TNNT1","KRT7","PDGFRB","ACTA2","MB","MYH7","CKM","HLA.A")
for (i in 1:length(sample.ids)) {
  DefaultAssay(spatial[[i]])<- "Spatial"
  plot1 <- VlnPlot(spatial[[i]], features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  plot2 <- SpatialFeaturePlot(spatial[[i]], features = "nCount_Spatial") + theme(legend.position = "right")
  plotQC<-plot1+plot2
  plot3 <- SpatialFeaturePlot(spatial[[i]], features = genes,ncol = 4)
  ggsave(paste0(sample.ids[i],"_QC.pdf"),p=plotQC,path = images, width = 5, height = 10, limitsize = FALSE)
  ggsave(paste0(sample.ids[i],"_marker_genes.pdf"),p=plot3,path = images, width = 20, height = 15, limitsize = FALSE)
}

############# Merging all samples on the cluster ############
spatial<- readRDS("/mnt/parscratch/users/mdq19mm/visium_analysis/spatial_nerve_fascicles.rds")
spatial<-merge(spatial[[1]], list(spatial[[2]],
                                  spatial[[3]],
                                  spatial[[4]],
                                  spatial[[5]],
                                  spatial[[6]],
                                  spatial[[7]],
                                  spatial[[8]],
                                  spatial[[9]],
                                  spatial[[10]],
                                  spatial[[11]],
                                  spatial[[12]],
                                  spatial[[13]],
                                  spatial[[14]],
                                  spatial[[15]],
                                  spatial[[16]],
                                  spatial[[17]],
                                  spatial[[18]],
                                  spatial[[19]],
                                  spatial[[20]],
                                  spatial[[21]],
                                  spatial[[22]],
                                  spatial[[23]],
                                  spatial[[24]],
                                  spatial[[25]],
                                  spatial[[26]],
                                  spatial[[27]],
                                  spatial[[28]]))


############# Make aggregated count matrix ###########
# Making count matrix
spatial@meta.data <- spatial@meta.data %>% 
  mutate(loupe.annotation = str_replace_all(loupe.annotation,"nerves3","3")) %>% 
  mutate(loupe.annotation = str_replace_all(loupe.annotation,"nerves2","2")) %>% 
  mutate(loupe.annotation = str_replace_all(loupe.annotation,"nerves1","1")) %>%
  mutate(loupe.annotation = str_replace_all(loupe.annotation,"nerves","1")) %>%
  mutate(section=paste0(orig.ident,"_",loupe.annotation))

Idents(spatial)<- "section"

# Extracting counts and metadata
DefaultAssay(spatial) <- "Spatial"
spatial <- JoinLayers(spatial)
counts<-spatial[["Spatial"]]$counts
metadata<-spatial@meta.data

# Aggregating counts
sce<- SingleCellExperiment(assays=list(counts=counts), colData=metadata)
aggr_counts <- aggregate(t(counts(sce)), by = list(spatial$section), FUN = "sum") %>%
  column_to_rownames(var="Group.1") %>% data.matrix()

# Explore output matrix
class(aggr_counts)
dim(aggr_counts)
aggr_counts[1:4, 1:4]

#create metadata
metadata <- colData(sce) %>% as.data.frame() %>% dplyr::select(orig.ident,pain.status,sex,sample,capture.area,section)
metadata <- metadata[!duplicated(metadata), ] # Exclude duplicated rows
rownames(metadata) <- metadata$section # Rename rows

# Number of cells per sample 
t <- table(colData(sce)$section)
metadata$barcodes_count<-as.numeric(t)

# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts <- t(aggr_counts)

# save matrix
write.table(aggr_counts, 'DEseq2/spatial_aggr.txt')
write.table(metadata,"DEseq2/spatial_metadata.txt")
