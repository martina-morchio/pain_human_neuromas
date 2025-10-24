####### snRNAseq vascular annotation ###########
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

img.path<- here('outs', 'reanalysis', 'vascular', 'images')
tbl.path<- here('outs', 'reanalysis', 'vascular', 'tables')

paths<- c(img.path,tbl.path)
for (path in paths) {
  if (!dir.exists(path)) {dir.create(path)}
}

# making images for visualisation of snRNAseq data
# load annotated R object
nerves.comb <- readRDS(here('data','nerves_annotated.rds'))

# endo
vascular<- subset(nerves.comb, idents = c('Endo_1', 'Endo_2', 'Endo_3', 'Endo_4', 'Mural'))
DefaultAssay(vascular)<- "RNA"
vascular<- vascular %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()

# check violin plots
genes<- c('PECAM1', 'TIE1', 'PROX1', 'SELE', 'ACTA2', 'MYH11')

violin<-VlnPlot(vascular, features = genes, pt.size = 0, ncol = 1) & 
  theme(axis.title = element_blank(), text = element_text(size = 24)) #& 
geom_boxplot(width=0.1, fill="white",position=position_dodge(1))


for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave("vascular_genes_violin.png", p=violin, path = img.path, width = 5, height = 20, limitsize = FALSE)

#umap
umap<-FeaturePlot(vascular, reduction = "umap", features=genes, label = TRUE, repel = TRUE, label.size = 3) + NoAxes() + NoLegend()
ggsave(paste0("vascular_genes_umap.png"), p=umap, path=img.path, width = 8, height = 10, limitsize = FALSE)

# try reclustering properly
# split the dataset into a list of two seurat objects (nerve and neuroma)
vascular.list <- SplitObject(vascular, split.by = "type")

# SCT normalization of each sample type
vascular.list.SCT <- lapply(X = vascular.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = vascular.list.SCT, nfeatures = 3000)
vascular.list.SCT <- PrepSCTIntegration(object.list = vascular.list.SCT, anchor.features = features)
vascular.list.SCT <- lapply(X = vascular.list.SCT, FUN = RunPCA, features = features)

# reciprocal PCA integration with SCT normalization
vascular.anchors <- FindIntegrationAnchors(object.list = vascular.list.SCT, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "rpca")
vascular.comb.sct <- IntegrateData(anchorset = vascular.anchors, normalization.method = "SCT", dims = 1:30)

# Run the standard workflow for visualization and clustering
vascular.comb.sct <- vascular.comb.sct %>%
  ScaleData() %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)

res<- 0.3
vascular.comb.sct <- FindClusters(vascular.comb.sct, resolution=res, cluster.name = paste0('reclustering_',res))


# check violin plots
Idents(vascular.comb.sct) <- paste0('reclustering_',res)
DefaultAssay(vascular.comb.sct) <- 'RNA'
genes<- c('PECAM1', 'TIE1', 'PROX1', 'SELE', 'ACTA2', 'MYH11')

violin<-VlnPlot(vascular.comb.sct, features = genes, pt.size = 0, ncol = 1) & 
  theme(axis.title = element_blank(), text = element_text(size = 24)) #& 
geom_boxplot(width=0.1, fill="white",position=position_dodge(1))


for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave(paste0("vascular_genes_violin_reclustering_",res,".png"), p=violin, path = img.path, width = 10, height = 20, limitsize = FALSE)

#umap
umap<-FeaturePlot(vascular.comb.sct, reduction = "umap", features=genes, label = TRUE, repel = TRUE, label.size = 3) + NoAxes() + NoLegend()
ggsave(paste0("vascular_genes_umap_reclustering_",res,".png"), p=umap, path=img.path, width = 8, height = 10, limitsize = FALSE)

# look at proportion
table(vascular.comb.sct$reclustering_0.3, vascular.comb.sct$type)

# annotation based on subtypes suggested by reviewer
ann<- list('0'='Infl_Endo',
           '1'='Endo',
           '2'='Infl_Endo',
           '3'= 'SMC',
           '4'='Endo',
           '5'='Pericytes',
           '6'='SMC',
           '7'='Pericytes',
           '8'= 'Lymph_Endo')

meta <- vascular.comb.sct@meta.data
new.meta <- meta %>% mutate(new_clustering=recode(reclustering_0.3, !!!ann))

vascular.comb.sct@meta.data <- new.meta
DefaultAssay(vascular.comb.sct) <- 'RNA'

# check violin plots
Idents(vascular.comb.sct) <- 'new_clustering'
genes<- c('PECAM1', 'TIE1', 'PROX1', 'SELE', 'ACTA2', 'MYH11')

violin<-VlnPlot(vascular.comb.sct, features = genes, pt.size = 0, ncol = 1) & 
  theme(axis.title = element_blank(), text = element_text(size = 24)) #& 
geom_boxplot(width=0.1, fill="white",position=position_dodge(1))


for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave(paste0("vascular_genes_violin_new.png"), p=violin, path = img.path, width = 6, height = 20, limitsize = FALSE)

#umap
umap<-FeaturePlot(vascular.comb.sct, reduction = "umap", features=genes, label = TRUE, repel = TRUE, label.size = 3) + NoAxes() + NoLegend()
ggsave(paste0("vascular_genes_umap_new.png"), p=umap, path=img.path, width = 8, height = 10, limitsize = FALSE)

# look at proportion
table(vascular.comb.sct$new_clustering, vascular.comb.sct$type)

## further analysis
# heatmap for vascular cells
vascular_clusters <- c("Infl_Endo","Endo","Lymph_Endo","SMC","Pericytes")

vascular.comb.sct <- PrepSCTFindMarkers(vascular.comb.sct)
vascular.comb.sct.markers<- FindAllMarkers(vascular.comb.sct, only.pos = TRUE,
                                           logfc.threshold = 0.25,  
                                           min.pct = 0.25, assay= 'SCT')

signif.markers <- vascular.comb.sct.markers %>%
  filter(p_val_adj < 0.05, 
         avg_log2FC > 0.5,    
         pct.1 > 0.50) 
vascular.comb.sct.markers.clusters<- list()

for (i in 1:length(levels(vascular.comb.sct))) {
  vascular.comb.sct.markers.clusters[[i]] <- signif.markers %>% 
    filter(cluster==levels(vascular.comb.sct)[[i]]) %>% arrange(desc(avg_log2FC))
}

write_xlsx(vascular.comb.sct.markers.clusters, path=here(tbl.path,"vascular_markers_clusters_seuratV4.xlsx"))

# images
top5_genes_vascular <- c()
for ( i in 1:5 ) {
  top5_genes_cluster <- vascular.comb.sct.markers.clusters[[i]]$gene[1:5]
  top5_genes_vascular <- append(top5_genes_vascular,top5_genes_cluster)
}

vascular.comb.sct <- ScaleData(vascular.comb.sct, features=top5_genes_vascular)
heatmap <- DoHeatmap(vascular.comb.sct, assay = 'RNA', features = top5_genes_vascular, size=3, ) + 
  theme(plot.margin= unit(c(1, 2, 1, 1), "cm")) + NoLegend()
ggsave("vascular_markers_heatmap.png",p=heatmap,path = img.path, width = 5, height = 5, limitsize = FALSE)

# Distribution across sample type
id_count_table <- as.data.frame(table(vascular.comb.sct@meta.data$new_clustering, vascular.comb.sct@meta.data$type)) %>% filter(Var1 %in% vascular_clusters)
colnames(id_count_table) <- c("cluster","type","frequency")

id_count_table_nerve <- id_count_table %>% filter(type=="nerve") %>% mutate(perc = frequency/sum(frequency))
id_count_table_neuroma <- id_count_table %>% filter(type=="neuroma") %>% mutate(perc = frequency/sum(frequency))
id_count_table <- bind_rows(id_count_table_nerve,id_count_table_neuroma)

histogram<- ggplot(id_count_table, aes(x=type, y=perc, fill=cluster))+
  geom_bar(stat="identity", color="black") +theme_classic() + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  labs(y="Frequency") 

ggsave("vascular_clusters_histogram.pdf", plot=histogram, width= 3, height=5, path=img.path)

# umap
umap<-DimPlot(vascular.comb.sct, reduction = "umap", label = TRUE, repel = TRUE, label.size = 3) + NoAxes() + NoLegend()
ggsave("vascular_UMAP.png", p=umap, path=img.path, width = 4, height = 4, limitsize = FALSE)


# save object
qsave(vascular.comb.sct, file =here('outs', 'reanalysis', 'vascular', 'vascular.qs'))
vascular.comb.sct<- qread(here('outs', 'reanalysis', 'vascular', 'vascular.qs'))
