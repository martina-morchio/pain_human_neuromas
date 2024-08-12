############ Cellchat spatial analysis images #############
library(CellChat)
library(Seurat)
library(tidyverse)
library(patchwork)
library(cowplot)
library(ComplexHeatmap)
options(future.globals.maxSize = 1e10)
options(stringsAsFactors = FALSE)

# This script is used to generate images to visualise cell-cell communication with cellchat v2
# the cellchat objects of the painful and non-painful samples are generated in 13_cellchat_object_computation.R

working_directory<- "/Users/martina/Documents/PhD/Main_projects/Visium/analysis/cellchat"
setwd(working_directory)
images <- "painful_vs_nonpainful/"
state <- c("Painful","NonPainful")
file.path <- "/Volumes/shared/boissonade/User/mdq19mm/visium_analysis/cellchat/"
files <- paste0(file.path, "cellchat_", state, ".rds")

# Read files
cellchat.P <- readRDS(files[1])
cellchat.NP <- readRDS(files[2])

# order levels
cell.levels<- c("Fibro","Endo","SC","Peri","Macro","Bcells")
cellchat.P <- updateClusterLabels(cellchat.P, new.order = cell.levels)
cellchat.NP <- updateClusterLabels(cellchat.NP, new.order = cell.levels)

# Merge them together
object.list <- list(NP = cellchat.NP, P = cellchat.P)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Identify altered interactions
ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
g<-gg1 + gg2
ggsave("altered_interactions.pdf", plot=g, path=images, width=3,height=1.5)

# compare differential number of interactions
pdf(paste0(images,"interactions_cells_difference.pdf"), width=5,height=5)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

# comparison of heatmaps
pdf(paste0(images,"difference_heatmaps.pdf"), width=6,height=3)
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2
dev.off()

# comparing number of interactions among different datasets
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf(paste0(images,"interactions_P_vs_NP.pdf"), width=7,height=3.5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

# Identify cell populations with significant changes in sending or receiving signals
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
pdf(paste0(images,"changes_cell_populations_sending_receiving.pdf"), width=8,height=4)
patchwork::wrap_plots(plots = gg)
dev.off()

# identify signalling changes of specific cell populations
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Bcells")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Macro")
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "SC")
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Endo")
gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Fibro")
gg6 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Peri")
pdf(paste0(images,"individual_cell_types_changes.pdf"), width=10,height=10)
patchwork::wrap_plots(plots = list(gg1,gg2,gg3,gg4,gg5,gg6), ncol = 2)
dev.off()

# Signalling groups based on functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
pdf(paste0(images,"groups_functional_similarity.pdf"), width=5,height=5)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 1)
dev.off()

# Pathway distance in the learned joint manifold
pdf(paste0(images,"pathway_distance.pdf"), width=3,height=12)
rankSimilarity(cellchat, type = "functional")
dev.off()

# Compare overall information flow
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", cutoff.pvalue=0.001, sources.use = c("SC","Peri","Macro","Endo"), stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", cutoff.pvalue=0.001, thresh=0.001, targets.use = NULL, stacked = F, do.stat = TRUE)
pdf(paste0(images,"change_in_information_flow.pdf"), width=6,height=14)
gg1 + gg2
dev.off()

# Endo information flow
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", cutoff.pvalue=0.001, thresh=0.001, sources.use = c("Endo"), targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", cutoff.pvalue=0.001, thresh=0.001, sources.use = c("Endo"), targets.use = NULL, stacked = F, do.stat = TRUE)
pdf(paste0(images,"Endo_change_in_information_flow.pdf"), width=7,height=13)
gg1 + gg2
dev.off()

# Bcells information flow
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", cutoff.pvalue=0.001, sources.use = c("Bcells"), targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", cutoff.pvalue=0.001, sources.use = c("Bcells"), targets.use = NULL, stacked = F, do.stat = TRUE)
pdf(paste0(images,"Bcells_change_in_information_flow.pdf"), width=7,height=13)
gg1 + gg2
dev.off()

# SC information flow
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", cutoff.pvalue=0.001, sources.use = c("SC"), targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", cutoff.pvalue=0.001, sources.use = c("SC"), targets.use = NULL, stacked = F, do.stat = TRUE)
pdf(paste0(images,"SC_change_in_information_flow.pdf"), width=7,height=13)
gg1 + gg2
dev.off()

# Fibro information flow
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", cutoff.pvalue=0.001, sources.use = c("Fibro"), targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", cutoff.pvalue=0.001, sources.use = c("Fibro"), targets.use = NULL, stacked = F, do.stat = TRUE)
pdf(paste0(images,"Fibro_change_in_information_flow.pdf"), width=7,height=13)
gg1 + gg2
dev.off()

# Peri information flow
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", cutoff.pvalue=0.001, sources.use = c("Peri"), targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", cutoff.pvalue=0.001, sources.use = c("Peri"), targets.use = NULL, stacked = F, do.stat = TRUE)
pdf(paste0(images,"Peri_change_in_information_flow.pdf"), width=7,height=13)
gg1 + gg2
dev.off()

# Macro information flow
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", cutoff.pvalue=0.001, sources.use = c("Macro"), targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", cutoff.pvalue=0.001, sources.use = c("Macro"), targets.use = NULL, stacked = F, do.stat = TRUE)
pdf(paste0(images,"Macro_change_in_information_flow.pdf"), width=7,height=13)
gg1 + gg2
dev.off()

## Identify upregulated and down-regulated signalling LR pairs
#### DE analysis to identify differential signalling ####
images<-"painful_vs_nonpainful/DE_analysis/"
pos.dataset = "P" # define positive dataset
features.name = paste0(pos.dataset, ".merged")

# calculate de genes
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, 
                                       features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, 
                                       thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)

write.csv(net,paste0(images,"DE_signalling.csv"),quote=F,row.names=F)

net.up <- subsetCommunication(cellchat, net = net, datasets = "P",ligand.logFC = 0.5, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "NP",ligand.logFC = -0.5, receptor.logFC = NULL)

# extract genes
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

# visualise with bubble plot
celltypes<-levels(cellchat.P@idents)
pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

for (i in 1:length(celltypes)) {
  gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
                          sources.use = celltypes[i], comparison = c(1, 2),  
                          angle.x = 90, remove.isolate = T,
                          title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
  ggsave(paste0(celltypes[i],"_up_bubble.pdf"), plot=gg1, path=images, width=5, height=20)
  
  gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, 
                          sources.use = celltypes[i], comparison = c(1, 2),  
                          angle.x = 90, remove.isolate = T,
                          title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
  ggsave(paste0(celltypes[i],"_down_bubble.pdf"), plot=gg2, path=images, width=5, height=20)
}

## images for use
# endothelial
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
                        sources.use = "Endo", targets.use = c("Endo","SC","Peri"), comparison = c(1, 2),  
                        angle.x = 90, remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pdf(paste0(images,"final_up_endo_bubble.pdf"), width=5, height=7)
gg1
dev.off()

# visualise with chord diagram
par(mfrow = c(1,2), xpd=TRUE)
pdf(paste0(images,"chord_endo_genes.pdf"), width=10, height=10)
netVisual_chord_gene(object.list[[2]], sources.use = "Endo", 
                     slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
dev.off()

##### Visualise individual pathways #####
images<- "painful_vs_nonpainful/pathways_oi/"
pathways.P<- cellchat.P@netP$pathways
pathways.NP<- cellchat.NP@netP$pathways
pathways.oi <- c("ICAM", "HGF", "VISTA", "PTPR", "SLITRK", "IL16", "NECTIN", "CDH", "CNTN", 
                 "TRAIL", "Desmosterol", "CysLTs", "NGF", "CD45", "OCLN", "CADM", "SELPLG", 
                 "COMPLEMENT", "CHEMERIN", "SPP1", "NEGR", "SN", "CD96", "SELL", "CD39", 
                 "VEGI", "Adenosine", "TNF", "CSF3", "NGL", "LIGHT", "ANNEXIN", "NRG", "ApoE")

pathways.oi <- c("TNF")

for (i in 1:length(pathways.oi)) {
  pdf(paste0("netVisual_aggregate_chord_",pathways.oi[i],".pdf"),width=7, height=7)
  for (j in 1:length(object.list)) {
    netVisual_aggregate(object.list[[j]], signaling = pathways.oi[i], layout = "chord", signaling.name = paste(pathways.oi[i], names(object.list)[j]))
  }
  dev.off()
  p1<-netAnalysis_contribution(cellchat.P, signaling = pathways.oi[i], title = paste("Contribution of each L-R pair",pathways.oi[i]," in painful"))
  p2<-netAnalysis_contribution(cellchat.NP, signaling = pathways.oi[i], title = paste("Contribution of each L-R pair",pathways.oi[i], " in non-painful"))
  ggsave(paste0("LR_contribution_",pathways.oi[i],".pdf"),plot=p1+p2,path=images,width=7, height=5)

  pdf(paste0(images,"netVisual_heatmap_",pathways.oi[i],".pdf"),width=3, height=2.5)
  p1<-netVisual_heatmap(cellchat.P, signaling = pathways.oi[i], color.heatmap = "Reds", title.name = paste(pathways.oi[i], " in painful"))
  draw(p1)
  p2<-netVisual_heatmap(cellchat.NP, signaling = pathways.oi[i], color.heatmap = "Reds", title.name = paste(pathways.oi[i], " in non-painful"))
  draw(p2)
  dev.off()

  pdf(paste0(images,"chord_celltypes_",pathways.oi[i],".pdf"),width=9, height=4.5)
  par(mfrow = c(1,2), xpd=TRUE)
  for (j in 1:length(object.list)) {
    netVisual_aggregate(object.list[[j]], signaling = pathways.oi[i], layout = "chord",
                        signaling.name = paste(pathways.oi[i], names(object.list)[j]))
  }
  dev.off()

  p1<-netVisual_bubble(cellchat.P, signaling = pathways.oi[i], remove.isolate = FALSE) +NoLegend() +ggtitle(paste0("Significant interactions in painful"))
  p2<-netVisual_bubble(cellchat.NP, signaling = pathways.oi[i], remove.isolate = FALSE) +ggtitle(paste0("Significant interactions in non-painful"))
  ggsave(paste0("netVisual_bubble_",pathways.oi[i],".pdf"),plot=p1+p2,path=images,width=10, height=4)

  pdf(paste0(images,"netVisual_chord_gene_",pathways.oi[i],".pdf"),width=5, height=4)
  netVisual_chord_gene(cellchat.P, signaling = pathways.oi[i], lab.cex = 0.5,legend.pos.y = 50,legend.pos.x = 10)
  netVisual_chord_gene(cellchat.NP, signaling = pathways.oi[i], lab.cex = 0.5,legend.pos.y = 50,legend.pos.x = 10)
  dev.off()
}

####### explore with cellchat explorer ######
runCellChatApp(cellchat.NP)
runCellChatApp(cellchat.P)

# Save them for later
save(object.list, file = paste0(file.path,"cellchat_object.list_neuromas.RData"))
save(cellchat, file = paste0(file.path,"cellchat_merged_neuromas.RData"))

