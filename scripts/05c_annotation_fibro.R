####### snRNAseq fibro annotation ###########
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

img.path<- here('outs', 'reanalysis', 'fibro', 'images')
tbl.path<- here('outs', 'reanalysis', 'fibro', 'tables')

paths<- c(img.path,tbl.path)
for (path in paths) {
  if (!dir.exists(path)) {dir.create(path, recursive = T)}
}

# making images for visualisation of snRNAseq data
# load annotated R object
nerves.comb <- readRDS(here('data','nerves_annotated.rds'))

# let's look at marker expression suggested by reviewer
# endo
fibro<- subset(nerves.comb, idents = c("EndoF", "PeriF_1", "PeriF_2", "PeriF_3", "PeriF_4", "PFF"))
DefaultAssay(fibro)<- "RNA"
fibro<- fibro %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()

# try reclustering properly
# split the dataset into a list of two seurat objects (nerve and neuroma)
fibro.list <- SplitObject(fibro, split.by = "type")

# SCT normalization of each sample type
fibro.list.SCT <- lapply(X = fibro.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = fibro.list.SCT, nfeatures = 3000)
fibro.list.SCT <- PrepSCTIntegration(object.list = fibro.list.SCT, anchor.features = features)
fibro.list.SCT <- lapply(X = fibro.list.SCT, FUN = RunPCA, features = features)

# reciprocal PCA integration with SCT normalization
fibro.anchors <- FindIntegrationAnchors(object.list = fibro.list.SCT, normalization.method = "SCT",
                                          anchor.features = features, dims = 1:30, reduction = "rpca")
fibro.comb.sct <- IntegrateData(anchorset = fibro.anchors, normalization.method = "SCT", dims = 1:30)

# Run the standard workflow for visualization and clustering
fibro.comb.sct <- fibro.comb.sct %>%
  ScaleData() %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)

res<- 0.3
fibro.comb.sct <- FindClusters(fibro.comb.sct, resolution=res, cluster.name = paste0('reclustering_',res))

# check violin plots
Idents(fibro.comb.sct) <- paste0('reclustering_',res)
DefaultAssay(fibro.comb.sct) <- 'RNA'
genes <- c("VIM","COL1A1","PRRX1","PI16","ABCA10","NGFR","SLC2A1","PTGDS",
                 "ABCA10","ABCA6","SOX9","CSPG4","IGFBP6","THBS1","NGFR",
                 "SLC2A1","FOSB","PTGDS","COL1A1","FOSB","PRRX1","FBLN1", 
                 "CCBE1", "COMP", "FXYD5", "ALPL", "CRABP2", "OGN")

violin<-VlnPlot(fibro.comb.sct, features = genes, pt.size = 0, ncol = 1) & 
  theme(axis.title = element_blank(), text = element_text(size = 24)) 


for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave(paste0("fibro_genes_violin_reclustering_",res,".png"), p=violin, path = img.path, width = 5, height = 30, limitsize = FALSE)

#umap
umap<-FeaturePlot(fibro.comb.sct, reduction = "umap", features=genes, label = TRUE, repel = TRUE, label.size = 3) + NoAxes() + NoLegend()
ggsave(paste0("fibro_genes_umap_reclustering_",res,".png"), p=umap, path=img.path, width = 20, height = 15, limitsize = FALSE)

# look at proportion
table(fibro.comb.sct$reclustering_0.3, fibro.comb.sct$type)

# annotation based on subtypes suggested by reviewer
ann<- list('0'='PFF',
           '1'='PeriF_1',
           '2'='EndoF',
           '3'= 'PeriF_2',
           '4'='PeriF_1',
           '5'='PeriF_1',
           '6'="PeriF_1",
           '7'='EndoF',
           '8'='PFF')

meta <- fibro.comb.sct@meta.data
new.meta <- meta %>% mutate(new_clustering=recode(reclustering_0.3, !!!ann))

fibro.comb.sct@meta.data <- new.meta

# check violin plots
Idents(fibro.comb.sct) <- 'new_clustering'
genes <- c("VIM","COL1A1","PRRX1","PI16","ABCA10","NGFR","SLC2A1","PTGDS",
           "ABCA10","ABCA6","PI16","SOX9","CSPG4","IGFBP6","THBS1","NGFR",
           "SLC2A1","FOSB","PTGDS","COL1A1","FOSB","PRRX1","FBLN1", 
           "CCBE1", "COMP", "FXYD5", "ALPL", "CRABP2", "OGN")

violin<-VlnPlot(fibro.comb.sct, features = genes, pt.size = 0, ncol = 1) & 
  theme(axis.title = element_blank(), text = element_text(size = 24)) 


for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave(paste0("fibro_genes_violin_new.png"), p=violin, path = img.path, width = 3, height = 45, limitsize = FALSE)

#umap
umap<-FeaturePlot(fibro.comb.sct, reduction = "umap", features=genes, label = TRUE, repel = TRUE, label.size = 3) + NoAxes() + NoLegend()
ggsave(paste0("fibro_genes_umap_new.png"), p=umap, path=img.path, width = 20, height = 20, limitsize = FALSE)

# look at proportion
table(fibro.comb.sct$new_clustering, fibro.comb.sct$type)


## further analysis
# heatmap for fibro cells
fibro.comb.sct <- PrepSCTFindMarkers(fibro.comb.sct)
fibro.comb.sct.markers<- FindAllMarkers(fibro.comb.sct, only.pos = TRUE,
                                          logfc.threshold = 0.25,  
                                          min.pct = 0.25, assay= 'SCT')

signif.markers <- fibro.comb.sct.markers %>%
  filter(p_val_adj < 0.05, 
         avg_log2FC > 0.5,    
         pct.1 > 0.50) 

fibro.comb.sct.markers.clusters<- list()

for (i in 1:length(levels(fibro.comb.sct))) {
  fibro.comb.sct.markers.clusters[[i]] <- signif.markers %>% 
    filter(cluster==levels(fibro.comb.sct)[[i]]) %>% arrange(desc(avg_log2FC))
}

write_xlsx(fibro.comb.sct.markers.clusters, path=here(tbl.path,"fibro_markers_clusters_seuratV4.xlsx"))

# images
top5_genes_fibro <- c()
for ( i in 1:4 ) {
  top5_genes_cluster <- fibro.comb.sct.markers.clusters[[i]]$gene[1:5]
  top5_genes_fibro <- append(top5_genes_fibro,top5_genes_cluster)
}

fibro.comb.sct <- ScaleData(fibro.comb.sct, features=top5_genes_fibro)
heatmap <- DoHeatmap(fibro.comb.sct, features = top5_genes_fibro, size=3, slot = "scale.data",disp.min = -2,disp.max = 2) + 
  NoLegend() + theme(plot.margin= unit(c(1, 2, 1, 1), "cm"))
ggsave("fibro_markers_heatmap.png",p=heatmap,path = img.path, width = 5, height = 5, limitsize = FALSE)

# Distribution across sample type
id_count_table <- as.data.frame(table(fibro.comb.sct@meta.data$new_clustering, fibro.comb.sct@meta.data$type)) %>% filter(Var1 %in% fibro_clusters)
colnames(id_count_table) <- c("cluster","type","frequency")

id_count_table_nerve <- id_count_table %>% filter(type=="nerve") %>% mutate(perc = frequency/sum(frequency))
id_count_table_neuroma <- id_count_table %>% filter(type=="neuroma") %>% mutate(perc = frequency/sum(frequency))
id_count_table <- bind_rows(id_count_table_nerve,id_count_table_neuroma)

histogram<- ggplot(id_count_table, aes(x=type, y=perc, fill=cluster))+
  geom_bar(stat="identity", color="black") +theme_classic() + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  labs(y="Frequency") 

ggsave("fibro_clusters_histogram.png", plot=histogram, width= 3, height=7, path=img.path)


# umap
umap<-DimPlot(fibro.comb.sct, reduction = "umap", label = TRUE, repel = TRUE, label.size = 3) + NoAxes() + NoLegend()
ggsave("fibro_UMAP.png", p=umap, path=img.path, width = 4, height = 4, limitsize = FALSE)


# save object
qsave(fibro.comb.sct, file =here('outs', 'reanalysis', 'fibro', 'fibro.qs'))
fibro.comb.sct<- qread(here('outs', 'reanalysis', 'fibro', 'fibro.qs'))

