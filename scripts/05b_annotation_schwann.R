####### snRNAseq schwann annotation ###########
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

img.path<- here('outs', 'reanalysis', 'schwann', 'images')
tbl.path<- here('outs', 'reanalysis', 'schwann', 'tables')

paths<- c(img.path,tbl.path)
for (path in paths) {
  if (!dir.exists(path)) {dir.create(path, recursive = T)}
}

# making images for visualisation of snRNAseq data
# load annotated R object
nerves.comb <- readRDS(here('data','nerves_annotated.rds'))

# let's look at marker expression suggested by reviewer
# endo
schwann<- subset(nerves.comb, idents = c('MSC_1', 'MSC_2', 'MSC_3', 'NMSC'))
DefaultAssay(schwann)<- "RNA"
schwann<- schwann %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()

# try reclustering properly
# split the dataset into a list of two seurat objects (nerve and neuroma)
schwann.list <- SplitObject(schwann, split.by = "type")

# SCT normalization of each sample type
schwann.list.SCT <- lapply(X = schwann.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = schwann.list.SCT, nfeatures = 3000)
schwann.list.SCT <- PrepSCTIntegration(object.list = schwann.list.SCT, anchor.features = features)
schwann.list.SCT <- lapply(X = schwann.list.SCT, FUN = RunPCA, features = features)

# reciprocal PCA integration with SCT normalization
schwann.anchors <- FindIntegrationAnchors(object.list = schwann.list.SCT, normalization.method = "SCT",
                                           anchor.features = features, dims = 1:30, reduction = "rpca")
schwann.comb.sct <- IntegrateData(anchorset = schwann.anchors, normalization.method = "SCT", dims = 1:30)

# Run the standard workflow for visualization and clustering
schwann.comb.sct <- schwann.comb.sct %>%
  ScaleData() %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)

res<- 0.3
schwann.comb.sct <- FindClusters(schwann.comb.sct, resolution=res, cluster.name = paste0('reclustering_',res))

# check violin plots
Idents(schwann.comb.sct) <- paste0('reclustering_',res)
DefaultAssay(schwann.comb.sct) <- 'RNA'
genes<- c("S100B","PRX","SOX10","DHH","L1CAM","NGFR",
          "ERBB3","PMP2","MPZ","MBP","NRXN1","SCN7A","PRIMA1", 
          "GAP43", "NCAM1","L1CAM", "BDNF", "GDNF","RUNX2" ,"ATF3","JUN",
          "EGR1", "FOS")

violin<-VlnPlot(schwann.comb.sct, features = genes, pt.size = 0, ncol = 1) & 
  theme(axis.title = element_blank(), text = element_text(size = 24)) 


for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave(paste0("schwann_genes_violin_reclustering_",res,".pdf"), p=violin, path = img.path, width = 5, height = 30, limitsize = FALSE)

#umap
umap<-FeaturePlot(schwann.comb.sct, reduction = "umap", features=genes, label = TRUE, repel = TRUE, label.size = 3) + NoAxes() + NoLegend()
ggsave(paste0("schwann_genes_umap_reclustering_",res,".pdf"), p=umap, path=img.path, width = 20, height = 15, limitsize = FALSE)

# look at proportion
table(schwann.comb.sct$reclustering_0.3, schwann.comb.sct$type)

# annotation based on subtypes suggested by reviewer
ann<- list('0'='mSC_1',
           '1'='mSC_2',
           '2'='nmSC',
           '3'= 'mSC_1',
           '4'='nmSC',
           '5'='damage_mSC',
           '6'="mSC_2",
           '7'='damage_nmSC')

meta <- schwann.comb.sct@meta.data
new.meta <- meta %>% mutate(new_clustering=recode(reclustering_0.3, !!!ann))

schwann.comb.sct@meta.data <- new.meta

# check violin plots
Idents(schwann.comb.sct) <- 'new_clustering'
genes<- c("S100B","PRX","SOX10","DHH","L1CAM","NGFR",
          "ERBB3","PMP2","MPZ","MBP","NRXN1","SCN7A","PRIMA1", 
          "NCAM1","L1CAM","ATF3","JUN","EGR1", "FOS")

violin<-VlnPlot(schwann.comb.sct, features = genes, pt.size = 0, ncol = 1) & 
  theme(axis.title = element_blank(), text = element_text(size = 24)) 


for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave(paste0("schwann_genes_violin_new.pdf"), p=violin, path = img.path, width = 3, height = 45, limitsize = FALSE)

#umap
umap<-FeaturePlot(schwann.comb.sct, reduction = "umap", features=genes, label = TRUE, repel = TRUE, label.size = 3) + NoAxes() + NoLegend()
ggsave(paste0("schwann_genes_umap_new.pdf"), p=umap, path=img.path, width = 20, height = 20, limitsize = FALSE)

# look at proportion
table(schwann.comb.sct$new_clustering, schwann.comb.sct$type)

## further analysis
# heatmap for schwann cells
schwann_clusters <- c("mSC_1","mSC_2", "nmSC", "damage_mSC", "damage_nmSC")

schwann.comb.sct <- PrepSCTFindMarkers(schwann.comb.sct)
schwann.comb.sct.markers<- FindAllMarkers(schwann.comb.sct, only.pos = TRUE,
                                          logfc.threshold = 0.25,  
                                          min.pct = 0.25, assay= 'SCT')

signif.markers <- schwann.comb.sct.markers %>%
  filter(p_val_adj < 0.05, 
         avg_log2FC > 0.5,    
         pct.1 > 0.50) 

schwann.comb.sct.markers.clusters<- list()

for (i in 1:length(levels(schwann.comb.sct))) {
  schwann.comb.sct.markers.clusters[[i]] <- signif.markers %>% 
    filter(cluster==levels(schwann.comb.sct)[[i]]) %>% arrange(desc(avg_log2FC))
}

write_xlsx(schwann.comb.sct.markers.clusters, path=here(tbl.path,"schwann_markers_clusters_seuratV4.xlsx"))

# images
top5_genes_schwann <- c()
for ( i in 1:4 ) {
  top5_genes_cluster <- schwann.comb.sct.markers.clusters[[i]]$gene[1:5]
  top5_genes_schwann <- append(top5_genes_schwann,top5_genes_cluster)
}

schwann.comb.sct <- ScaleData(schwann.comb.sct, features=top5_genes_schwann)
heatmap <- DoHeatmap(schwann.comb.sct, features = top5_genes_schwann, size=3, slot = "scale.data",disp.min = -2,disp.max = 2) + 
  NoLegend() + theme(plot.margin= unit(c(1, 2, 1, 1), "cm"))
ggsave("schwann_markers_heatmap.png",p=heatmap,path = img.path, width = 5, height = 5, limitsize = FALSE)

# Distribution across sample type
id_count_table <- as.data.frame(table(schwann.comb.sct@meta.data$new_clustering, schwann.comb.sct@meta.data$type)) %>% filter(Var1 %in% schwann_clusters)
colnames(id_count_table) <- c("cluster","type","frequency")

id_count_table_nerve <- id_count_table %>% filter(type=="nerve") %>% mutate(perc = frequency/sum(frequency))
id_count_table_neuroma <- id_count_table %>% filter(type=="neuroma") %>% mutate(perc = frequency/sum(frequency))
id_count_table <- bind_rows(id_count_table_nerve,id_count_table_neuroma)

histogram<- ggplot(id_count_table, aes(x=type, y=perc, fill=cluster))+
  geom_bar(stat="identity", color="black") +theme_classic() + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  labs(y="Frequency") 

ggsave("schwann_clusters_histogram.pdf", plot=histogram, width= 3, height=5, path=img.path)


# umap
umap<-DimPlot(schwann.comb.sct, reduction = "umap", label = TRUE, repel = TRUE, label.size = 3) + NoAxes() + NoLegend()
ggsave("schwann_UMAP.pdf", p=umap, path=img.path, width = 4, height = 4, limitsize = FALSE)

# save object
qsave(schwann.comb.sct, file =here('outs', 'reanalysis', 'schwann', 'schwann.qs'))
schwann.comb.sct<- qread(here('outs', 'reanalysis', 'schwann', 'schwann.qs'))

