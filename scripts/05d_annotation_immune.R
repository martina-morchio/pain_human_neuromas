####### snRNAseq immune annotation ###########
library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
library(enrichR)
library(writexl)
library(ggrepel)
library(qs)
library(here)
options(future.globals.maxSize = 1e10)
options(Seurat.object.assay.version = "v5")

here::i_am('pain_human_neuromas.Rproj')

img.path<- here('outs', 'reanalysis', 'immune', 'images')
tbl.path<- here('outs', 'reanalysis', 'immune', 'tables')

paths<- c(img.path,tbl.path)
for (path in paths) {
  if (!dir.exists(path)) {dir.create(path, recursive = T)}
}

# making images for visualisation of snRNAseq data
# load annotated R object
nerves.comb <- readRDS(here('data','nerves_annotated.rds'))

# let's look at marker expression suggested by reviewer
# endo
immune<- subset(nerves.comb, idents = c("Lympho", "Macro"))
DefaultAssay(immune)<- "RNA"
immune<- immune %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()

# try reclustering properly
# split the dataset into a list of two seurat objects (nerve and neuroma)
immune.list <- SplitObject(immune, split.by = "type")

# SCT normalization of each sample type
immune.list.SCT <- lapply(X = immune.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = immune.list.SCT, nfeatures = 3000)
immune.list.SCT <- PrepSCTIntegration(object.list = immune.list.SCT, anchor.features = features)
immune.list.SCT <- lapply(X = immune.list.SCT, FUN = RunPCA, features = features)

# reciprocal PCA integration with SCT normalization
immune.anchors <- FindIntegrationAnchors(object.list = immune.list.SCT, normalization.method = "SCT",
                                        anchor.features = features, dims = 1:30, reduction = "rpca")
immune.comb.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30)

# Run the standard workflow for visualization and clustering
immune.comb.sct <- immune.comb.sct %>%
  ScaleData() %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)

res<- 0.3
immune.comb.sct <- FindClusters(immune.comb.sct, resolution=res, cluster.name = paste0('reclustering_',res))

# check violin plots
Idents(immune.comb.sct) <- paste0('reclustering_',res)
DefaultAssay(immune.comb.sct) <- 'RNA'
genes <- c("CD163", "CD68", "CXCR4", "PTPRC", "CLEC9A", "XCR1", "BATF3", "CLEC4C", 
           "IRF8", "FCER1A", "CD1C", "CLEC10A", "CD19", "CD79B", "JCHAIN", "IGHG1",  "LYZ", 
           "CD14", "MRC1", "MS4A7","CX3CR1", "TREM2", "GPNMB", "SPP1", "APOE")

violin<-VlnPlot(immune.comb.sct, features = genes, pt.size = 0, ncol = 1) & 
  theme(axis.title = element_blank(), text = element_text(size = 24)) 


for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave(paste0("immune_genes_violin_reclustering_",res,".png"), p=violin, path = img.path, width = 5, height = 30, limitsize = FALSE)

#umap
umap<-FeaturePlot(immune.comb.sct, reduction = "umap", features=genes, label = TRUE, repel = TRUE, label.size = 3) + NoAxes() + NoLegend()
ggsave(paste0("immune_genes_umap_reclustering_",res,".png"), p=umap, path=img.path, width = 20, height = 15, limitsize = FALSE)

# look at proportion
table(immune.comb.sct$reclustering_0.3, immune.comb.sct$type)

# annotation based on subtypes suggested by reviewer
ann<- list('0'='Macro',
           '1'='Macro',
           '2'='Lympho',
           '3'= 'Macro',
           '4'='Macro',
           '5'='Macro',
           '6'="Macro",
           '7'='Macro')

meta <- immune.comb.sct@meta.data
new.meta <- meta %>% mutate(new_clustering=recode(reclustering_0.3, !!!ann))

immune.comb.sct@meta.data <- new.meta

# check violin plots
Idents(immune.comb.sct) <- 'new_clustering'
genes <- c("CD163", "CD68", "CXCR4", "PTPRC", "CLEC9A", "XCR1", "BATF3", "CLEC4C", 
           "IRF8", "FCER1A", "CD1C", "CLEC10A", "CD19", "CD79B", "JCHAIN", "IGHG1",  "LYZ", 
           "CD14", "MRC1", "MS4A7","CX3CR1", "TREM2", "GPNMB", "SPP1", "APOE")

violin<-VlnPlot(immune.comb.sct, features = genes, pt.size = 0, ncol = 1) & 
  theme(axis.title = element_blank(), text = element_text(size = 24)) 


for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave(paste0("immune_genes_violin_new.png"), p=violin, path = img.path, width = 2, height = 45, limitsize = FALSE)

#umap
umap<-FeaturePlot(immune.comb.sct, reduction = "umap", features=genes, label = TRUE, repel = TRUE, label.size = 3) + NoAxes() + NoLegend()
ggsave(paste0("immune_genes_umap_new.png"), p=umap, path=img.path, width = 20, height = 20, limitsize = FALSE)

# look at proportion
table(immune.comb.sct$new_clustering, immune.comb.sct$type)


## further analysis
# heatmap for immune cells
immune_clusters <- c('Macro', 'Lympho')
immune.comb.sct <- PrepSCTFindMarkers(immune.comb.sct)
immune.comb.sct.markers<- FindAllMarkers(immune.comb.sct, only.pos = TRUE,
                                        logfc.threshold = 0.25,  
                                        min.pct = 0.25, assay= 'SCT')

signif.markers <- immune.comb.sct.markers %>%
  filter(p_val_adj < 0.05, 
         avg_log2FC > 0.5,    
         pct.1 > 0.50) 

immune.comb.sct.markers.clusters<- list()

for (i in 1:length(levels(immune.comb.sct))) {
  immune.comb.sct.markers.clusters[[i]] <- signif.markers %>% 
    filter(cluster==levels(immune.comb.sct)[[i]]) %>% arrange(desc(avg_log2FC))
}

write_xlsx(immune.comb.sct.markers.clusters, path=here(tbl.path,"immune_markers_clusters_seuratV4.xlsx"))

# images
top5_genes_immune <- c()
for ( i in 1:2 ) {
  top5_genes_cluster <- immune.comb.sct.markers.clusters[[i]]$gene[1:5]
  top5_genes_immune <- append(top5_genes_immune,top5_genes_cluster)
}

immune.comb.sct <- ScaleData(immune.comb.sct, features=top5_genes_immune)
heatmap <- DoHeatmap(immune.comb.sct, features = top5_genes_immune, size=3, slot = "scale.data",disp.min = -2,disp.max = 2) + 
  NoLegend() + theme(plot.margin= unit(c(1, 2, 1, 1), "cm"))
ggsave("immune_markers_heatmap.png",p=heatmap,path = img.path, width = 5, height = 3.5, limitsize = FALSE)

# Distribution across sample type
id_count_table <- as.data.frame(table(immune.comb.sct@meta.data$new_clustering, immune.comb.sct@meta.data$type)) %>% filter(Var1 %in% immune_clusters)
colnames(id_count_table) <- c("cluster","type","frequency")

id_count_table_nerve <- id_count_table %>% filter(type=="nerve") %>% mutate(perc = frequency/sum(frequency))
id_count_table_neuroma <- id_count_table %>% filter(type=="neuroma") %>% mutate(perc = frequency/sum(frequency))
id_count_table <- bind_rows(id_count_table_nerve,id_count_table_neuroma)

histogram<- ggplot(id_count_table, aes(x=type, y=perc, fill=cluster))+
  geom_bar(stat="identity", color="black") +theme_classic() + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  labs(y="Frequency") 

ggsave("immune_clusters_histogram.png", plot=histogram, width= 3, height=7, path=img.path)


# umap
umap<-DimPlot(immune.comb.sct, reduction = "umap", label = TRUE, repel = TRUE, label.size = 3) + NoAxes() + NoLegend()
ggsave("immune_UMAP.png", p=umap, path=img.path, width = 4, height = 4, limitsize = FALSE)


# save object
qsave(immune.comb.sct, file =here('outs', 'reanalysis', 'immune', 'immune.qs'))
immune.comb.sct<- qread(here('outs', 'reanalysis', 'immune', 'immune.qs'))

