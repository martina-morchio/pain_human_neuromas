###################### Making Seurat object for all Visium samples ##############
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleCellExperiment)
library(Matrix)
library(here)
library(qs)
library(data.table)
options(future.globals.maxSize = 1e9)

filedir <- "D:/visium"

##### Create list of samples ####
# Loading samples
folders<- list.dirs(filedir,full.names = FALSE, recursive = FALSE)
sample.ids <- grep("^LN", folders, value = TRUE)
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
  spatial[[i]]<-RenameCells(spatial[[i]], new.names = paste0(sample.ids[[i]], '-', colnames(spatial[[i]])))
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
  spatial[[i]]<-RenameCells(spatial[[i]], new.names = paste0(sample.ids[[i]], '-', colnames(spatial[[i]])))
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
  spatial[[i]]<-RenameCells(spatial[[i]], new.names = paste0(sample.ids[[i]], '-', colnames(spatial[[i]])))
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
  spatial[[i]]<-RenameCells(spatial[[i]], new.names = paste0(sample.ids[[i]], '-', colnames(spatial[[i]])))
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
  spatial[[i]]<-RenameCells(spatial[[i]], new.names = paste0(sample.ids[[i]], '-', colnames(spatial[[i]])))
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
  spatial[[i]]<-RenameCells(spatial[[i]], new.names = paste0(sample.ids[[i]], '-', colnames(spatial[[i]])))
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
  spatial[[i]]<-RenameCells(spatial[[i]], new.names = paste0(sample.ids[[i]], '-', colnames(spatial[[i]])))
}
for (i in 1:4) {spatial[[i+9]]$capture.area<- area[i]}

######### SCTransform for normalisation ##########
spatial <- lapply(spatial,SCTransform, assay = "Spatial", verbose = FALSE)

############# QC plots ############
# images<- here('outs','spatial', 'QC')
# if (!dir.exists(images)) {dir.create(images, recursive=T)}
# 
# for (i in 1:length(sample.ids)) {
#   plot1 <- VlnPlot(spatial[[i]], features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
#   plot2 <- SpatialFeaturePlot(spatial[[i]], features = "nCount_Spatial") + theme(legend.position = "right")
#   plotQC<-plot1+plot2
#   plot3 <- SpatialFeaturePlot(spatial[[i]], features = c("SCN9A","MPZ","NRXN1","SLC2A1","CD68","PTPRC","PI16","EGFL7","TNNT1","KRT7","PDGFRB","ACTA2"),ncol = 4)
#   ggsave(paste0(sample.ids[i],"_QC.pdf"),p=plotQC,path = images, width = 5, height = 10, limitsize = FALSE)
#   ggsave(paste0(sample.ids[i],"_marker_genes.pdf"),p=plot3,path = images, width = 20, height = 15, limitsize = FALSE)
# }

# merge
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

############# add annotation from giotto ############
meta<- fread(here('outs','reanalysis', 'giotto', 'metadata_giotto_harmony_integrated_ann2024.csv'))

cluster_ann<- list('1'="Fibro",
                '2'="Endo",
                '3'="SC1",
                '4'="Myo1",
                '5'="Peri",
                '6'="SC2",
                '7'="Myo2",
                '8'="SC3",
                '9'="Myo3",
                '10'="SC4",
                '11'="Myo4",
                '12'="SC5",
                '13'="Bcells",
                '14'="SC6",
                '15'="Macro",
                '16'='NA',
                "17"='NA')

meta<- meta %>% mutate(cell_type = recode(leiden_harmony, !!!cluster_ann))

# first subset with used barcodes only
sp<- subset(spatial, cells = meta$cell_ID)
old.meta <- sp@meta.data %>% rownames_to_column('cell_ID')
new.meta<- inner_join(old.meta, meta, by= 'cell_ID')
rownames(new.meta) <- new.meta$cell_ID
new.meta$cell_ID <- NULL

# remove unnecessary info
new.meta<- new.meta %>% select(-c('nCount_Spatial', 'nFeature_Spatial', 'nCount_SCT', 'nFeature_SCT', 'in_tissue', 'array_row', 'array_col', 'orig.ident', 'leiden_harmony'))

# update metadata
sp@meta.data <- new.meta

# remove cluster 16 and 17
Idents(sp) <- 'cell_type'
sp <- subset(sp, idents = 'NA', invert = TRUE)

# get dimreduc from giotto object
seu<-qread(here('outs','spatial', 'giotto_as_seurat.qs'))
sp@reductions[['pca']]<- seu@reductions[['pca']]
sp@reductions[['umap_harmony']]<- seu@reductions[['umap_harmony']]
sp@reductions[['spatial']]<- seu@reductions[['spatial']]

# save
qsave(sp, here('outs','spatial', 'spatial_annotated.qs'))
#sp<-qread(here('outs','spatial', 'spatial_annotated.qs'))

