######## Cellchat for spatial transcriptomics computation ############
library(CellChat)
library(patchwork)
library(Seurat)
library(tidyverse)
options(stringsAsFactors = FALSE)

working_directory <- "/mnt/parscratch/users/mdq19mm/visium_analysis/cellchat/"
setwd(working_directory)
images<- "images/"

###### Prepare seurat objects ######
# read spatial samples list generated in 11_pseudobulk_preprocessing.R
file <- "/mnt/parscratch/users/mdq19mm/visium_analysis/spatial_combined.rds"
spatial <- readRDS(file)

# read clusters from Giotto generated in 08_sce_object_from_giotto.R
ann<-read.csv("metadata_giotto_harmony_integrated.csv")

# converting leiden clusters into cell types
new_idents <- c("1"="Fibro",
                "2"="Endo",
                "3"="SC1",
                "4"="Myo",
                "5"="Peri",
                "6"="SC",
                "7"="Myo",
                "8"="SC",
                "9"="Myo",
                "10"="SC",
                "11"="Myo",
                "12"="SC",
                "13"="Bcells",
                "14"="SC",
                "15"="Macro")
ann<- ann %>% mutate(cellchat_id=new_idents[as.character(leiden_harmony)])

############## Painful samples - adding metadata to all replicates ##############
sample<- "painful"
sample_names<- c("LN1_0-","LN1_A1-","LN1_B1-","LN1_D1-","LN12_A1-","LN12_B1-","LN13_C1-","LN13_D1-","LN15_A1-","LN15_B1-","LN15_C1-","LN15_D1-")
samples<- list(spatial[[1]],spatial[[2]],spatial[[3]],spatial[[5]],spatial[[6]],spatial[[7]],spatial[[8]],spatial[[9]],spatial[[10]],spatial[[11]],spatial[[12]],spatial[[13]])
for (i in 1:length(sample_names)) {
  ann_sample <- ann %>% filter(str_detect(cell_ID,sample_names[i])) %>% mutate(cell_ID = str_replace(cell_ID,sample_names[i],"")) %>% 
    column_to_rownames("cell_ID") %>% select(list_ID, cellchat_id)
  samples[[i]] <- AddMetaData(samples[[i]], ann_sample)
}


# tidying metadata and collapsing clusters
for (i in 1:length(sample_names)) {
  Idents(samples[[i]]) <- "cellchat_id"
  samples[[i]]$list_ID <- NULL
  if ("Macro" %in% levels(samples[[i]])) {
    samples[[i]]<- subset(samples[[i]], idents=c("SC","Fibro","Endo","Peri","Bcells","Macro"))
  } else {
    samples[[i]]<- subset(samples[[i]], idents=c("SC","Fibro","Endo","Peri","Bcells"))
  }
}

# record number of spots per cluster
clus_distrib<- data.frame()
for (i in 1:length(sample_names)) {
  table<-table(samples[[i]]$cellchat_id)
  clus_distrib<- rbind(clus_distrib,table)
  colnames(clus_distrib)<- names(table)
}
rownames(clus_distrib)<- sample_names
write.csv(clus_distrib, file=paste0(sample,"_cluster_distribution_cellchat.csv"), quote=F)

# image for QC
color.use <- scPalette(nlevels(samples[[1]])); names(color.use) <- levels(samples[[1]])
for (i in 1:length(sample_names)) {
  g<-SpatialDimPlot(samples[[i]], label = T, label.size = 3, cols=color.use)
  ggsave(paste0(sample_names[i],"clusters.pdf"), path=images, plot=g, width=7, height=7)
}

### Prepare input data for CellChat analysis ###
data.input<-list()
for (i in 1:length(sample_names)) {
  data.input[[i]] <- Seurat::GetAssayData(samples[[i]], slot = "data", assay = "SCT") # normalized data matrix
  colnames(data.input[[i]]) <- paste0(sample_names[i], colnames(data.input[[i]]))
}
rownames<- lapply(data.input, rownames)
genes.common <- Reduce(intersect, rownames) # getting all common genes
for (i in 1:length(sample_names)) {
  data.input[[i]]<-data.input[[i]][genes.common,]
}
data.input.all <- do.call(cbind, data.input) #matrix with all counts

# define the meta data:
meta <- list()
sample_names<- gsub("-", "", sample_names)
for (i in 1:length(sample_names)) {
  meta[[i]] <- data.frame(labels = Seurat::Idents(samples[[i]]), samples = sample_names[i], row.names = names(Seurat::Idents(samples[[i]]))) # manually create a dataframe consisting of the cell labels
  meta[[i]]$samples <- factor(meta[[i]]$samples)
}
meta.all <- do.call(rbind, meta)
rownames(meta.all) <- colnames(data.input.all)

# Load spatial transcriptomics information
spatial.locs <- list()
for (i in 1:length(sample_names)) {
  spatial.locs[[i]] <- Seurat::GetTissueCoordinates(samples[[i]], scale = NULL, cols = c("imagerow", "imagecol")) 
}
spatial.locs.all <- do.call(rbind, spatial.locs)
rownames(spatial.locs.all) <- colnames(data.input.all)

json_path <- paste0("/mnt/parscratch/users/mdq19mm/visium_data/",sample_names,"/outs/spatial/")
scalefactors<- list()
for (i in 1:length(sample_names)) {
  scalefactors[[i]] = jsonlite::fromJSON(txt = file.path(json_path[i], 'scalefactors_json.json'))
}
spot.size = 65 # the theoretical spot size (um) in 10X Visium
spatial.factors<- list()
for ( i in 1:length(sample_names)) {
  conversion.factor = spot.size/scalefactors[[i]]$spot_diameter_fullres
  spatial.factors[[i]] = data.frame(ratio = conversion.factor, tol = spot.size/2)
  d.spatial <- computeCellDistance(coordinates = spatial.locs[[i]], ratio = spatial.factors[[i]]$ratio, tol = spatial.factors[[i]]$tol)
  distance<-min(d.spatial[d.spatial!=0])
  print(paste(sample_names[[i]], distance, sep=": ")) 
}
spatial.factors.all <- do.call(rbind, spatial.factors)
rownames(spatial.factors.all) <- sample_names

## Create cellchat object ##
cellchat <- createCellChat(object = data.input.all, meta = meta.all, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs.all, spatial.factors = spatial.factors.all)
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

## Inference of cell-cell communication ##
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# interactions and weight/strengths
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf(paste0(images,"aggregated_cellcell_comm_network.pdf"), width=7, height=7)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

## Save cellchat object ##
saveRDS(cellchat, file = "cellchat_Painful.rds")



############## Non Painful samples - adding metadata to all replicates ##############
sample<- "non_painful"
sample_names<- c("LN2_0-","LN2_A1-","LN2_B1-","LN2_C1-","LN2_D1-","LN7_0-","LN7_A1-","LN7_B1-","LN7_C1-","LN7_D1-","LN8_0-","LN8_A1-","LN8_B1-","LN8_C1-","LN8_D1-")
samples<- list(spatial[[14]],spatial[[15]],spatial[[16]],spatial[[17]],spatial[[18]],spatial[[19]],spatial[[20]],spatial[[21]],spatial[[22]],spatial[[23]],spatial[[24]],spatial[[25]],spatial[[26]],spatial[[27]],spatial[[28]])
for (i in 1:length(sample_names)) {
  ann_sample <- ann %>% filter(str_detect(cell_ID,sample_names[i])) %>% mutate(cell_ID = str_replace(cell_ID,sample_names[i],"")) %>% 
    column_to_rownames("cell_ID") %>% select(list_ID, cellchat_id)
  samples[[i]] <- AddMetaData(samples[[i]], ann_sample)
}


# tidying metadata and collapsing clusters
for (i in 1:length(sample_names)) {
  Idents(samples[[i]]) <- "cellchat_id"
  samples[[i]]$list_ID <- NULL
  samples[[i]]<- subset(samples[[i]], idents=c("SC","Fibro","Endo","Peri","Bcells","Macro"))
}

# record number of spots per cluster
clus_distrib<- data.frame()
for (i in 1:length(sample_names)) {
  table<-table(samples[[i]]$cellchat_id)
  clus_distrib<- rbind(clus_distrib,table)
  colnames(clus_distrib)<- names(table)
}
rownames(clus_distrib)<- sample_names
write.csv(clus_distrib, file=paste0(sample,"_cluster_distribution_cellchat.csv"), quote=F)

# image for QC
color.use <- scPalette(nlevels(samples[[1]])); names(color.use) <- levels(samples[[1]])
for (i in 1:length(sample_names)) {
  g<-SpatialDimPlot(samples[[i]], label = T, label.size = 3, cols=color.use)
  ggsave(paste0(sample_names[i],"clusters.pdf"), path=images, plot=g, width=7, height=7)
}

### Prepare input data for CellChat analysis ###
data.input<-list()
for (i in 1:length(sample_names)) {
  data.input[[i]] <- Seurat::GetAssayData(samples[[i]], slot = "data", assay = "SCT") # normalized data matrix
  colnames(data.input[[i]]) <- paste0(sample_names[i], colnames(data.input[[i]]))
}
rownames<- lapply(data.input, rownames)
genes.common <- Reduce(intersect, rownames) # getting all common genes
for (i in 1:length(sample_names)) {
  data.input[[i]]<-data.input[[i]][genes.common,]
}
data.input.all <- do.call(cbind, data.input) #matrix with all counts

# define the meta data:
meta <- list()
sample_names<- gsub("-", "", sample_names)
for (i in 1:length(sample_names)) {
  meta[[i]] <- data.frame(labels = Seurat::Idents(samples[[i]]), samples = sample_names[i], row.names = names(Seurat::Idents(samples[[i]]))) # manually create a dataframe consisting of the cell labels
  meta[[i]]$samples <- factor(meta[[i]]$samples)
}
meta.all <- do.call(rbind, meta)
rownames(meta.all) <- colnames(data.input.all)

# Load spatial transcriptomics information
spatial.locs <- list()
for (i in 1:length(sample_names)) {
  spatial.locs[[i]] <- Seurat::GetTissueCoordinates(samples[[i]], scale = NULL, cols = c("imagerow", "imagecol")) 
}
spatial.locs.all <- do.call(rbind, spatial.locs)
rownames(spatial.locs.all) <- colnames(data.input.all)

json_path <- paste0("/mnt/parscratch/users/mdq19mm/visium_data/",sample_names,"/outs/spatial/")
scalefactors<- list()
for (i in 1:length(sample_names)) {
  scalefactors[[i]] = jsonlite::fromJSON(txt = file.path(json_path[i], 'scalefactors_json.json'))
}
spot.size = 65 # the theoretical spot size (um) in 10X Visium
spatial.factors<- list()
for ( i in 1:length(sample_names)) {
  conversion.factor = spot.size/scalefactors[[i]]$spot_diameter_fullres
  spatial.factors[[i]] = data.frame(ratio = conversion.factor, tol = spot.size/2)
  d.spatial <- computeCellDistance(coordinates = spatial.locs[[i]], ratio = spatial.factors[[i]]$ratio, tol = spatial.factors[[i]]$tol)
  distance<-min(d.spatial[d.spatial!=0])
  print(paste(sample_names[[i]], distance, sep=": ")) 
}
spatial.factors.all <- do.call(rbind, spatial.factors)
rownames(spatial.factors.all) <- sample_names

## Create cellchat object ##
cellchat <- createCellChat(object = data.input.all, meta = meta.all, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs.all, spatial.factors = spatial.factors.all)
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

## Inference of cell-cell communication ##
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# interactions and weight/strengths
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf(paste0(images,"aggregated_cellcell_comm_network.pdf"), width=7, height=7)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

## Save cellchat object ##
saveRDS(cellchat, file = "cellchat_NonPainful.rds")

