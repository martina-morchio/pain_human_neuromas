####### Images snRNAseq ###########

library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
library(enrichR)
library(writexl)
library(ggrepel)
options(future.globals.maxSize = 1e10)
options(Seurat.object.assay.version = "v5")

# making images for visualisation of snRNAseq data
# load annotated R object
nerves.comb <- readRDS("/Volumes/shared/boissonade/User/mdq19mm/snRNAseq/combined_runs/nerves_annotated.rds")
working_directory <- "/Users/martina/Documents/PhD/PAPER/images/nucseq/R"
setwd(working_directory)

# UMAPs
DefaultAssay(nerves.comb) <- "integrated"
Idents(nerves.comb) <- "annotation"
p1 <- DimPlot(nerves.comb, reduction = "umap", group.by = "type", order = c("neuroma","nerve")) + NoAxes() + ggtitle("")
p2 <- DimPlot(nerves.comb, reduction = "umap", label = TRUE, repel = TRUE, label.box = TRUE, label.size = 3) + NoAxes() + NoLegend()
umap<-cowplot::plot_grid(plotlist=list(p1,p2),nrow = 1)
ggsave("all_UMAP.pdf", p=umap, path=working_directory, width = 12, height = 6, limitsize = FALSE)
p3 <- DimPlot(nerves.comb, reduction = "umap", split.by = "type") + NoAxes()
ggsave("splitbytype_UMAP.pdf", p=p3, path=working_directory, width = 12, height = 6, limitsize = FALSE)

# QC by cluster
QC<-VlnPlot(nerves.comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "annotation")
ggsave("QC.pdf",plot=QC,width = 20, height = 5, path=working_directory)

# Distribution across samples
id_count_table <- as.data.frame(table(nerves.comb@meta.data$annotation, nerves.comb@meta.data$type))
colnames(id_count_table) <- c("cluster","type","frequency")

id_count_table_nerve <- id_count_table %>% filter(type=="nerve") %>% mutate(perc = frequency/sum(frequency))
id_count_table_neuroma <- id_count_table %>% filter(type=="neuroma") %>% mutate(perc = frequency/sum(frequency))
id_count_table <- bind_rows(id_count_table_nerve,id_count_table_neuroma)

histogram<- ggplot(id_count_table, aes(x=type, y=perc, fill=cluster))+
  geom_bar(stat="identity", color="black") +theme_classic() + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  labs(y="Frequency") 

ggsave("nerves_clusters_histogram.pdf", plot=histogram, width= 5, height=7, path=working_directory)

# dotplot
markers.to.plot <- c("SOX10","PMP22","S100B","MBP","NCAM1","NRXN1","L1CAM","GFAP","PLP1","COL1A1","DCN","PI16","ABCA9","COL15A1",
                     "SLC2A1","CLDN1","SFRP2","CRABP2","PRRX1","SELE","ICAM1","PECAM1","EGFL7","ACTA2","PDGFRB","CSF1R","CD163","PTPRC","IL7R","TRAC",
                     "MUC5B","KRT7","KRT14","TNNT1")
DefaultAssay(nerves.comb) <- "RNA"

p1<-DotPlot(nerves.comb, features = markers.to.plot, assay="RNA", cols = c("blue", "deeppink"), dot.scale = 8, split.by = "type") +
  RotatedAxis() + theme(axis.text.y = element_text(size=30),axis.text.x = element_text(size=15))
p2<-DotPlot(nerves.comb, features = markers.to.plot, assay="RNA", dot.scale = 8) + RotatedAxis() + 
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_blank())
p3<-DotPlot(nerves.comb, features = markers.to.plot, assay="RNA", cols = c("blue", "darkturquoise","blueviolet","deeppink"), dot.scale = 8, split.by = "id") +
  RotatedAxis()+ theme(axis.text.y = element_text(size=30),axis.text.x = element_text(size=15))

ggsave("markers_split_by_type_dotplot.pdf",p=p1,path = working_directory, width = 15, height = 16, limitsize = FALSE)
ggsave("markers_dotplot.pdf",p=p2,path = working_directory, width = 17, height = 10, limitsize = FALSE)
ggsave("markers_split_by_id_dotplot.pdf",p=p3,path = working_directory, width = 15, height = 32, limitsize = FALSE)


#### Subtype characterization ####
# load top DE genes for each cluster calculated with 04_finding_markers.R to generate heatmaps
library(readxl)
excel<- "/Users/martina/Documents/PhD/Main_projects/snRNAseq/combined/integration/annotated_nerves_markers_clusters_seuratV4.xlsx"
nerves.comb.markers.clusters <- list()
for (i in 1:24) {
  nerves.comb.markers.clusters[[i]]<- read_excel(excel, sheet=i)
}

for (i in 1:length(nerves.comb.markers.clusters)) {
  names(nerves.comb.markers.clusters)[i] <- nerves.comb.markers.clusters[[i]]$cluster[1]
}

#### Schwann cells #####
schwann <- readRDS("/Volumes/shared/boissonade/User/mdq19mm/snRNAseq/combined_runs/schwann.rds")
working_directory <- "schwann/"
DefaultAssay(schwann)<- "RNA"
schwann<- NormalizeData(schwann)
schwann<- ScaleData(schwann)

#Dim plot UMAP
umap<-DimPlot(schwann, reduction = "umap", label = TRUE, repel = TRUE, label.size = 3) + NoAxes() + NoLegend()
ggsave("schwann_UMAP.pdf", p=umap, path=working_directory, width = 4, height = 4, limitsize = FALSE)

# heatmap for schwann cells
top5_genes_schwann <- c()
schwann_clusters <- c("mSC_1","mSC_2","mSC_3","nmSC")
for ( i in schwann_clusters ) {
  top5_genes_cluster <- nerves.comb.markers.clusters[[i]]$gene[1:5]
  top5_genes_schwann <- append(top5_genes_schwann,top5_genes_cluster)
}
schwann <- ScaleData(schwann, features=top5_genes_schwann)
heatmap <- DoHeatmap(schwann, features = top5_genes_schwann,size=3) + NoLegend()
ggsave("schwann_markers_heatmap.pdf",p=heatmap,path = working_directory, width = 5, height = 5, limitsize = FALSE)

# Distribution across sample type
id_count_table <- as.data.frame(table(schwann@meta.data$annotation, schwann@meta.data$type)) %>% filter(Var1 %in% c("mSC_1","mSC_2","mSC_3","nmSC"))
colnames(id_count_table) <- c("cluster","type","frequency")

id_count_table_nerve <- id_count_table %>% filter(type=="nerve") %>% mutate(perc = frequency/sum(frequency))
id_count_table_neuroma <- id_count_table %>% filter(type=="neuroma") %>% mutate(perc = frequency/sum(frequency))
id_count_table <- bind_rows(id_count_table_nerve,id_count_table_neuroma)

histogram<- ggplot(id_count_table, aes(x=type, y=perc, fill=cluster))+
  geom_bar(stat="identity", color="black") +theme_classic() + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  labs(y="Frequency") 

ggsave("schwann_clusters_histogram.pdf", plot=histogram, width= 3, height=7, path=working_directory)

# Looking at particular genes
DefaultAssay(schwann) <- "RNA"

schwann_genes <- c("S100B","PRX","SOX10","DHH","L1CAM","NGFR","ERBB3","PMP2","MPZ","MBP","NRXN1","SCN7A","PRIMA1")
violin<-VlnPlot(schwann, features = schwann_genes, pt.size = 0, ncol = 1) & 
  theme(axis.title = element_blank(), text = element_text(size = 24)) #& 
  geom_boxplot(width=0.1, fill="white",position=position_dodge(1))


for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave("schwann_genes_violin.pdf", p=violin, path = working_directory, width = 5, height = 30, limitsize = FALSE)

#### Fibroblasts ####
fibro <- readRDS("/Volumes/shared/boissonade/User/mdq19mm/snRNAseq/combined_runs/fibro.rds")
working_directory <- "fibro/"
DefaultAssay(fibro)<- "RNA"
fibro<- NormalizeData(fibro)
fibro<- ScaleData(fibro)

#Dim plot UMAP
umap<-DimPlot(fibro, reduction = "umap", label = TRUE, repel = TRUE, label.size = 3) + NoAxes() + NoLegend()
ggsave("fibro_UMAP.pdf", p=umap, path=working_directory, width = 4, height = 4, limitsize = FALSE)

# heatmap for fibro cells
top5_genes_fibro <- c()
fibro_clusters <- c("EndoF","PeriF_1","PeriF_2","PeriF_3","PeriF_4","PFF")
for ( i in fibro_clusters ) {
  top5_genes_cluster <- nerves.comb.markers.clusters[[i]]$gene[1:5]
  top5_genes_fibro <- append(top5_genes_fibro,top5_genes_cluster)
}
fibro <- ScaleData(fibro, features=top5_genes_fibro)
heatmap <- DoHeatmap(fibro, features = top5_genes_fibro, size=3) + NoLegend() 
ggsave("fibro_markers_heatmap.pdf",p=heatmap,path = working_directory, width = 5, height = 5, limitsize = FALSE)

# Distribution across sample type
id_count_table <- as.data.frame(table(fibro@meta.data$annotation, fibro@meta.data$type)) %>% filter(Var1 %in% fibro_clusters)
colnames(id_count_table) <- c("cluster","type","frequency")

id_count_table_nerve <- id_count_table %>% filter(type=="nerve") %>% mutate(perc = frequency/sum(frequency))
id_count_table_neuroma <- id_count_table %>% filter(type=="neuroma") %>% mutate(perc = frequency/sum(frequency))
id_count_table <- bind_rows(id_count_table_nerve,id_count_table_neuroma)

histogram<- ggplot(id_count_table, aes(x=type, y=perc, fill=cluster))+
  geom_bar(stat="identity", color="black") +theme_classic() + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  labs(y="Frequency") 

ggsave("fibro_clusters_histogram.pdf", plot=histogram, width= 3, height=7, path=working_directory)

# Looking at particular genes
DefaultAssay(fibro) <- "RNA"
fibro_genes <- c("VIM","COL1A1","PRRX1","PI16","ABCA10","NGFR","SLC2A1","PTGDS")
violin<-VlnPlot(fibro, features = fibro_genes, pt.size = 0, ncol = 1) 

for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave("fibro_genes_violin.pdf", p=violin, path = working_directory, width = 5, height = 20, limitsize = FALSE)

# Endoneurial
DefaultAssay(fibro) <- "RNA"
fibro_genes <- c("ABCA10","ABCA6","PI16","CSPG4")
violin<-VlnPlot(fibro, features = fibro_genes, pt.size = 0, ncol = 1) 

for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave("endoF_genes_violin.pdf", p=violin, path = working_directory, width = 5, height = 10, limitsize = FALSE)

# Perineurial
DefaultAssay(fibro) <- "RNA"
fibro_genes <- c("IGFBP6","THBS1","NGFR","SLC2A1","FOSB","PTGDS")
violin<-VlnPlot(fibro, features = fibro_genes, pt.size = 0, ncol = 3) 

for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave("periF_genes_violin.pdf", p=violin, path = working_directory, width = 15, height = 6, limitsize = FALSE)

# Pro-fibrotic
DefaultAssay(fibro) <- "RNA"
fibro_genes <- c("COL1A1","FOSB","PRRX1","FBLN1")
violin<-VlnPlot(fibro, features = fibro_genes, pt.size = 0, ncol = 2) 

for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave("prrx1_fibro_genes_violin.pdf", p=violin, path = working_directory, width = 10, height = 6, limitsize = FALSE)

#### Immune ####
immune <- readRDS("/Volumes/shared/boissonade/User/mdq19mm/snRNAseq/combined_runs/immune.rds")
working_directory <- "immune/"
DefaultAssay(immune)<- "RNA"
immune<- NormalizeData(immune)
immune<- ScaleData(immune)

# heatmap for immune cells
top5_genes_immune <- c()
immune_clusters <- c("Macro","Lympho")
for ( i in immune_clusters ) {
  top5_genes_cluster <- nerves.comb.markers.clusters[[i]]$gene[1:5]
  top5_genes_immune <- append(top5_genes_immune,top5_genes_cluster)
}
immune <- ScaleData(immune, features=top5_genes_immune)

heatmap <- DoHeatmap(immune, features = top5_genes_immune, size=3) + NoLegend() 
ggsave("immune_markers_heatmap.pdf",p=heatmap,path = working_directory, width = 5, height = 2.5, limitsize = FALSE)

# Distribution across sample type
Idents(immune) <- "annotation"
id_count_table <- as.data.frame(table(immune@meta.data$annotation, immune@meta.data$type)) %>% filter(Var1 %in% immune_clusters)
colnames(id_count_table) <- c("cluster","type","frequency")

id_count_table_nerve <- id_count_table %>% filter(type=="nerve") %>% mutate(perc = frequency/sum(frequency))
id_count_table_neuroma <- id_count_table %>% filter(type=="neuroma") %>% mutate(perc = frequency/sum(frequency))
id_count_table <- bind_rows(id_count_table_nerve,id_count_table_neuroma)

histogram<- ggplot(id_count_table, aes(x=type, y=perc, fill=cluster))+
  geom_bar(stat="identity", color="black") +theme_classic() + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  labs(y="Frequency") 

ggsave("immune_clusters_histogram.pdf", plot=histogram, width= 3, height=7, path=working_directory)

## Particular genes
# Lymphocytes
DefaultAssay(immune) <- "RNA"
immune_genes <- c("CXCR4","CD69","PTPRC","TRBC2","CSF1R","CD163","CLEC7A","CD68")
violin<-VlnPlot(immune, features = immune_genes, pt.size = 0, ncol = 2) 

#+ scale_fill_manual("cluster", values = c("Macro"='blueviolet',"Lympho"='darkolivegreen')) to change colour
for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave("lymphocytes_genes_violin.pdf", p=violin, path = working_directory, width = 6, height = 10, limitsize = FALSE)

# Macrophages
DefaultAssay(immune) <- "RNA"
immune_genes <- c("CSF1R","CD163","CLEC7A","CD68")
#immune_genes <- c("CD68","OSM")
violin<-VlnPlot(immune, features = immune_genes, pt.size = 0, ncol = 1) 

for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + scale_fill_manual("cluster", values = c("Macro"='blueviolet',"Lympho"='darkolivegreen')) + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave("macrophages_genes_violin.pdf", p=violin, path = working_directory, width = 3, height = 10, limitsize = FALSE)

#### Vascular ####
vascular <- readRDS("/Volumes/shared/boissonade/User/mdq19mm/snRNAseq/combined_runs/vascular.rds")
working_directory <- "vascular/"
DefaultAssay(schwann)<- "RNA"
vascular<- NormalizeData(vascular)
vascular<- ScaleData(vascular)

# heatmap for vascular cells
top5_genes_vascular <- c()
vascular_clusters <- c("Endo_1","Endo_2","Endo_3","Endo_4","Mural")
for ( i in vascular_clusters ) {
  top5_genes_cluster <- nerves.comb.markers.clusters[[i]]$gene[1:5]
  top5_genes_vascular <- append(top5_genes_vascular,top5_genes_cluster)
}
vascular <- ScaleData(vascular, features=top5_genes_vascular)
heatmap <- DoHeatmap(vascular, features = top5_genes_vascular, size=3) + NoLegend() 
ggsave("vascular_markers_heatmap.pdf",p=heatmap,path = working_directory, width = 5, height = 5, limitsize = FALSE)

# Distribution across sample type
id_count_table <- as.data.frame(table(vascular@meta.data$annotation, vascular@meta.data$type)) %>% filter(Var1 %in% vascular_clusters)
colnames(id_count_table) <- c("cluster","type","frequency")

id_count_table_nerve <- id_count_table %>% filter(type=="nerve") %>% mutate(perc = frequency/sum(frequency))
id_count_table_neuroma <- id_count_table %>% filter(type=="neuroma") %>% mutate(perc = frequency/sum(frequency))
id_count_table <- bind_rows(id_count_table_nerve,id_count_table_neuroma)

histogram<- ggplot(id_count_table, aes(x=type, y=perc, fill=cluster))+
  geom_bar(stat="identity", color="black") +theme_classic() + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  labs(y="Frequency") 

ggsave("vascular_clusters_histogram.pdf", plot=histogram, width= 3, height=7, path=working_directory)

# specific genes
DefaultAssay(vascular) <- "RNA"

vascular_genes<- c("CLDN5","EGFL7","IL6","SELE","ICAM1","ACTA2","PDGFRB")
violin<-VlnPlot(vascular, features = vascular_genes, pt.size = 0, ncol = 1) 

for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave("vascular_genes_violin.pdf", p=violin, path = working_directory, width = 5, height = 17, limitsize = FALSE)

#### supplementary ####
neuroma.spec <- subset(nerves.comb, idents = c("SGC_1","SGC_2","SGC_3","SGC_4","Myo"))
nerves.spec <- subset(nerves.comb, idents = c("Astro", "Oligo", "MenF_1", "MenF_2"))
saveRDS(neuroma.spec, "/Volumes/shared/boissonade/User/mdq19mm/snRNAseq/combined_runs/neuroma.spec.rds")
saveRDS(nerves.spec, "/Volumes/shared/boissonade/User/mdq19mm/snRNAseq/combined_runs/nerves.spec.rds")

#### neuromas specific clusters ####
neuroma.spec <- readRDS("/Volumes/shared/boissonade/User/mdq19mm/snRNAseq/combined_runs/neuroma.spec.rds")
working_directory <- "other/"
DefaultAssay(neuroma.spec)<- "RNA"
neuroma.spec<- NormalizeData(neuroma.spec)
neuroma.spec<- ScaleData(neuroma.spec)
neuroma.spec_clusters <- c("SGC_1","SGC_2","SGC_3","SGC_4","Myo")

# Distribution across sample type
id_count_table <- as.data.frame(table(neuroma.spec@meta.data$annotation, neuroma.spec@meta.data$id)) %>% filter(Var1 %in% neuroma.spec_clusters)
colnames(id_count_table) <- c("cluster","id","frequency")

id_count_table_nerve <- id_count_table %>% filter(id=="n1") %>% mutate(perc = frequency/sum(frequency))
id_count_table_neuroma <- id_count_table %>% filter(id=="n2") %>% mutate(perc = frequency/sum(frequency))
id_count_table <- bind_rows(id_count_table_nerve,id_count_table_neuroma)

histogram<- ggplot(id_count_table, aes(x=id, y=perc, fill=cluster))+
  geom_bar(stat="identity", color="black") +theme_classic() + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  labs(y="Frequency") 

ggsave("neuroma.spec_clusters_histogram.pdf", plot=histogram, width= 3, height=7, path=working_directory)

# specific genes
DefaultAssay(neuroma.spec) <- "RNA"

neuroma.spec_genes<- c("MUC5B", "AQP5", "KRT19", "KRT17", "KRT14","MYL1", "TNNT1", "TNNT3", "TNNI1","ACTA2")
violin<-VlnPlot(neuroma.spec, features = neuroma.spec_genes, pt.size = 0, ncol = 1) 

for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave("neuroma.spec_genes_violin.pdf", p=violin, path = working_directory, width = 5, height = length(neuroma.spec_genes)*3, limitsize = FALSE)

#### nerves specific clusters ####
nerves.spec <- readRDS("/Volumes/shared/boissonade/User/mdq19mm/snRNAseq/combined_runs/nerves.spec.rds")
working_directory <- "other/"
DefaultAssay(nerves.spec)<- "RNA"
nerves.spec<- NormalizeData(nerves.spec)
nerves.spec<- ScaleData(nerves.spec)
nerves.spec_clusters <- c("Astro", "Oligo", "MenF_1", "MenF_2")

# Distribution across sample type
id_count_table <- as.data.frame(table(nerves.spec@meta.data$annotation, nerves.spec@meta.data$id)) %>% filter(Var1 %in% nerves.spec_clusters)
colnames(id_count_table) <- c("cluster","id","frequency")

id_count_table_nerve <- id_count_table %>% filter(id=="tg1") %>% mutate(perc = frequency/sum(frequency))
id_count_table_neuroma <- id_count_table %>% filter(id=="tg2") %>% mutate(perc = frequency/sum(frequency))
id_count_table <- bind_rows(id_count_table_nerve,id_count_table_neuroma)

histogram<- ggplot(id_count_table, aes(x=id, y=perc, fill=cluster))+
  geom_bar(stat="identity", color="black") +theme_classic() + 
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  labs(y="Frequency") 

ggsave("nerves.spec_clusters_histogram.pdf", plot=histogram, width= 3, height=7, path=working_directory)

# specific genes
DefaultAssay(nerves.spec) <- "RNA"

nerves.spec_genes<- c("FXYD5", "ALPL", "CRABP2", "OGN", "PTGDS", "GFAP", "OLIG1", "OLIG2", "MOG", "CNP", "PLP1")
violin<-VlnPlot(nerves.spec, features = nerves.spec_genes, pt.size = 0, ncol = 1) 

for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + theme(axis.title = element_blank(), text = element_text(size = 24)) + geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave("nerves.spec_genes_violin.pdf", p=violin, path = working_directory, width = 5, height = length(nerves.spec_genes)*3, limitsize = FALSE)

######### DE Genes #############
# generating dotplot of genes differentially expressed in painful and non-painful neuromas identified with spatial transcriptomics
de.genes <- read.csv("/Users/martina/Documents/Visium/analysis/DEseq2/De_tables/P_vs_NP_all_genes.csv")

Idents(nerves.comb) <- "type"
neuromas <- subset(nerves.comb, idents = "neuroma")

sig.genes <- de.genes %>% filter(padj<0.05)
filt.genes <- de.genes %>% filter(padj<0.05, abs(log2FoldChange) > 2)
Idents(neuromas) <- "annotation"
Idents(nerves.comb) <- "annotation"
genes <- c("CXCL2","CXCL8","SLC52A3","MMP19","ADGRG1","JMJD1C","GSTM3","TNFAIP","NLRC5")
dotplot<-DotPlot(neuromas, features=genes, dot.scale = 8) +
  RotatedAxis() + theme(axis.text.y = element_text(size=30),axis.text.x = element_text(size=30))
ggsave("dotplot_DE_neuromas_selected_genes.pdf",plot=dotplot, path = working_directory, height=12, width=13)


