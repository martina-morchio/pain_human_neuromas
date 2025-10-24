####### snRNAseq annotation updated ###########
library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
library(enrichR)
library(writexl)
library(readxl)
library(ggrepel)
library(qs)
library(here)
options(future.globals.maxSize = 1e10)
options(Seurat.object.assay.version = "v5")

here::i_am('pain_human_neuromas.Rproj')

img.path<- here('outs', 'reanalysis', 'paper_figures')

paths<- c(img.path,tbl.path)
for (path in paths) {
  if (!dir.exists(path)) {dir.create(path, recursive = T)}
}

# load objects
nerves.comb <- qread(here('outs', 'reanalysis', 'all', 'nerves.comb.qs'))
immune<- qread(here('outs', 'reanalysis', 'immune', 'immune.qs'))
fibro<- qread(here('outs', 'reanalysis', 'fibro', 'fibro.qs'))
schwann<- qread(here('outs', 'reanalysis', 'schwann', 'schwann.qs'))
vascular<- qread(here('outs', 'reanalysis', 'vascular', 'vascular.qs'))

# fig 1B
DefaultAssay(nerves.comb) <- "integrated"
p2 <- DimPlot(nerves.comb, reduction = "umap", label = TRUE, repel = TRUE, label.box = TRUE, label.size = 3) + 
  NoAxes() + 
  NoLegend() &
  theme(panel.background = element_blank(),
        plot.background = element_blank())

ggsave("all_UMAP.png", p=p2, path=img.path, width = 6, height = 6, limitsize = FALSE, bg = "transparent")

# fig 1C dotplot
markers.to.plot <- c("SOX10","PMP22","S100B","MBP","SCN7A","NRXN1","L1CAM", "FOSB", "ATF3","GFAP","PLP1","COL1A1","DCN","PI16","ABCA9","COL15A1",
                     "SLC2A1","CLDN1","SFRP2","CRABP2","PRRX1","PECAM1","EGFL7","TIE1", "SELE","PROX1","ACTA2","MYH11","CSF1R","CD163","PTPRC","IL7R","TRAC",
                     "MUC5B","KRT7","KRT14","TNNT1")
DefaultAssay(nerves.comb) <- "RNA"

p2<-DotPlot(nerves.comb, features = markers.to.plot, assay="RNA", dot.scale = 8) + RotatedAxis() + 
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_blank()) + 
  theme(axis.text.y = element_text(size=20), 
        axis.text.x = element_text(size=20, face = "italic"), 
        axis.title = element_blank())

ggsave("markers_dotplot.png",p=p2,path = img.path, width = 17, height = 10, limitsize = FALSE)

# fig 1D Schwann
# violin plots
genes<- c("SOX10","MBP", "PRX", "PMP2",  "SCN7A", "PRIMA1", "FOSB", "ATF3")

violin<-VlnPlot(schwann, features = genes, pt.size = 0, ncol = 4) & 
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + 
    theme(axis.title = element_blank(),
          plot.title = element_text(size = 24, face = 'bold.italic'), 
          text = element_text(size = 24, face='plain'),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent")) + 
    geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave(paste0("schwann_genes_violin_new.png"), p=violin, path = img.path, width = 15, height = 5, limitsize = FALSE, bg = "transparent")

# proportion
id_count_table <- as.data.frame(table(schwann@meta.data$new_clustering, schwann@meta.data$type)) %>% filter(Var1 %in% schwann_clusters)
colnames(id_count_table) <- c("cluster","type","frequency")

id_count_table_nerve <- id_count_table %>% filter(type=="nerve") %>% mutate(perc = frequency/sum(frequency))
id_count_table_neuroma <- id_count_table %>% filter(type=="neuroma") %>% mutate(perc = frequency/sum(frequency))
id_count_table <- bind_rows(id_count_table_nerve,id_count_table_neuroma)

histogram <- ggplot(id_count_table, aes(x = type, y = perc, fill = cluster)) +
  geom_bar(stat = "identity", color = "black") + 
  theme_classic() + 
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, size=14, hjust = 1),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent")) +
  scale_x_discrete(labels = c('Nerve', 'Neuroma')) + 
  labs(fill='Cluster')

ggsave("schwann_clusters_histogram.png", plot = histogram, width = 3, height = 5.5, path = img.path, bg = "transparent")

# heatmap
sheets <- excel_sheets(here('outs', 'reanalysis', 'schwann', 'tables',"schwann_markers_clusters_seuratV4.xlsx"))
top5_genes <- data.frame()
for (i in 1:length(sheets)) {
  markers.clusters <- read_xlsx(here('outs', 'reanalysis', 'schwann', 'tables',"schwann_markers_clusters_seuratV4.xlsx"), sheet = sheets[[i]])
  top5_genes_cluster <- data.frame(genes = markers.clusters$gene[1:5],
                                   name = markers.clusters$cluster[1:5])
  top5_genes <- rbind(top5_genes,top5_genes_cluster)
}

schwann <- ScaleData(schwann, features=top5_genes$genes)
heatmap <- DoHeatmap(schwann, features = top5_genes$genes, size=3, slot = "scale.data",disp.min = -2,disp.max = 2) + 
  NoLegend() + theme(plot.margin= unit(c(1, 2, 1, 1), "cm"), axis.text.y = element_text(face = "italic"),
                     panel.background = element_rect(fill = "transparent"),
                     plot.background = element_rect(fill = "transparent"))
ggsave("schwann_markers_heatmap.png",p=heatmap,path = img.path, width = 6, height = 5, limitsize = FALSE, bg='transparent')


# fig 1E Fibro
# violin plots
fibro_clusters<- c("PFF","PeriF_1","EndoF","PeriF_2")
genes<- c("VIM", "PI16", "PRRX1", "SLC2A1", "COL1A1", "NGFR",  "THBS1","PTGDS")

violin<-VlnPlot(fibro, features = genes, pt.size = 0, ncol = 4) & 
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + 
    theme(axis.title = element_blank(),
          plot.title = element_text(size = 24, face = 'bold.italic'), 
          text = element_text(size = 24, face='plain'),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent")) + 
    geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave(paste0("fibro_genes_violin_new.png"), p=violin, path = img.path, width = 15, height = 5, limitsize = FALSE, bg = "transparent")

# proportion
id_count_table <- as.data.frame(table(fibro@meta.data$new_clustering, fibro@meta.data$type)) %>% filter(Var1 %in% fibro_clusters)
colnames(id_count_table) <- c("cluster","type","frequency")
id_count_table_nerve <- id_count_table %>% filter(type=="nerve") %>% mutate(perc = frequency/sum(frequency))
id_count_table_neuroma <- id_count_table %>% filter(type=="neuroma") %>% mutate(perc = frequency/sum(frequency))
id_count_table <- bind_rows(id_count_table_nerve,id_count_table_neuroma)

histogram <- ggplot(id_count_table, aes(x = type, y = perc, fill = cluster)) +
  geom_bar(stat = "identity", color = "black") + 
  theme_classic() + 
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, size=14, hjust = 1),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent")) +
  scale_x_discrete(labels = c('nerve'='Nerve', 'neuroma'='Neuroma')) + 
  labs(fill='Cluster')

ggsave("fibro_clusters_histogram.png", plot = histogram, width = 2.5, height = 5.5, path = img.path, bg = "transparent")

# heatmap
sheets <- excel_sheets(here('outs', 'reanalysis', 'fibro', 'tables',"fibro_markers_clusters_seuratV4.xlsx"))
top5_genes <- data.frame()
for (i in 1:length(sheets)) {
  markers.clusters <- read_xlsx(here('outs', 'reanalysis', 'fibro', 'tables',"fibro_markers_clusters_seuratV4.xlsx"), sheet = sheets[[i]])
  top5_genes_cluster <- data.frame(genes = markers.clusters$gene[1:5],
                                   name = markers.clusters$cluster[1:5])
  top5_genes <- rbind(top5_genes,top5_genes_cluster)
}

fibro <- ScaleData(fibro, features=top5_genes$genes)
heatmap <- DoHeatmap(fibro, features = top5_genes$genes, size=3, slot = "scale.data",disp.min = -2,disp.max = 2) + 
  NoLegend() + theme(plot.margin= unit(c(1, 2, 1, 1), "cm"), axis.text.y = element_text(face = "italic"),
                     panel.background = element_rect(fill = "transparent"),
                     plot.background = element_rect(fill = "transparent"))
ggsave("fibro_markers_heatmap.png",p=heatmap,path = img.path, width = 6, height = 5, limitsize = FALSE, bg='transparent')


# fig 1F immune
# violin plots
immune_clusters<-  c("Lympho", "Macro")
genes<- c("CXCR4", "CD69", "CSF1R", "CD163", "PTPRC", "TRBC2", "CLEC7A", "CD68")

violin<-VlnPlot(immune, features = genes, pt.size = 0, ncol = 4) & 
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + 
    theme(axis.title = element_blank(),
          plot.title = element_text(size = 24, face = 'bold.italic'), 
          text = element_text(size = 24, face='plain'),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent")) + 
    geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave(paste0("immune_genes_violin_new.png"), p=violin, path = img.path, width = 10, height = 5, limitsize = FALSE, bg = "transparent")

# proportion
id_count_table <- as.data.frame(table(immune@meta.data$new_clustering, immune@meta.data$type)) %>% filter(Var1 %in% immune_clusters)
colnames(id_count_table) <- c("cluster","type","frequency")
id_count_table_nerve <- id_count_table %>% filter(type=="nerve") %>% mutate(perc = frequency/sum(frequency))
id_count_table_neuroma <- id_count_table %>% filter(type=="neuroma") %>% mutate(perc = frequency/sum(frequency))
id_count_table <- bind_rows(id_count_table_nerve,id_count_table_neuroma)

histogram <- ggplot(id_count_table, aes(x = type, y = perc, fill = cluster)) +
  geom_bar(stat = "identity", color = "black") + 
  theme_classic() + 
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, size=14, hjust = 1),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent")) +
  scale_x_discrete(labels = c('nerve'='Nerve', 'neuroma'='Neuroma')) + 
  labs(fill='Cluster')

ggsave("immune_clusters_histogram.png", plot = histogram, width = 2.5, height = 5.5, path = img.path, bg = "transparent")

# heatmap
sheets <- excel_sheets(here('outs', 'reanalysis', 'immune', 'tables',"immune_markers_clusters_seuratV4.xlsx"))
top5_genes <- data.frame()
for (i in 1:length(sheets)) {
  markers.clusters <- read_xlsx(here('outs', 'reanalysis', 'immune', 'tables',"immune_markers_clusters_seuratV4.xlsx"), sheet = sheets[[i]])
  top5_genes_cluster <- data.frame(genes = markers.clusters$gene[1:5],
                                   name = markers.clusters$cluster[1:5])
  top5_genes <- rbind(top5_genes,top5_genes_cluster)
}

immune <- ScaleData(immune, features=top5_genes$genes)
heatmap <- DoHeatmap(immune, features = top5_genes$genes, size=3, slot = "scale.data",disp.min = -2,disp.max = 2) + 
  NoLegend() + theme(plot.margin= unit(c(1, 2, 1, 1), "cm"), axis.text.y = element_text(face = "italic"),
                     panel.background = element_rect(fill = "transparent"),
                     plot.background = element_rect(fill = "transparent"))
ggsave("immune_markers_heatmap.png",p=heatmap,path = img.path, width = 6, height = 3, limitsize = FALSE, bg='transparent')


# fig 1G vascular
# violin plots
vascular_clusters<-  c("Infl_Endo","Endo","Lymph_Endo","SMC","Pericytes")
genes<- c("PECAM1", "CLDN5", "TIE1", "PROX1", "SELE", "IL6", "ACTA2","MYH11")

violin<-VlnPlot(vascular, features = genes, pt.size = 0, ncol = 4) & 
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

for (i in 1:length(violin)) {
  violin[[i]] <- violin[[i]] + 
    theme(axis.title = element_blank(),
          plot.title = element_text(size = 24, face = 'bold.italic'), 
          text = element_text(size = 24, face='plain'),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent")) + 
    geom_boxplot(width=0.1, fill="white",position=position_dodge(1))
}

ggsave(paste0("vascular_genes_violin_new.png"), p=violin, path = img.path, width = 15, height = 5, limitsize = FALSE, bg = "transparent")

# proportion
id_count_table <- as.data.frame(table(vascular@meta.data$new_clustering, vascular@meta.data$type)) %>% filter(Var1 %in% vascular_clusters)
colnames(id_count_table) <- c("cluster","type","frequency")
id_count_table_nerve <- id_count_table %>% filter(type=="nerve") %>% mutate(perc = frequency/sum(frequency))
id_count_table_neuroma <- id_count_table %>% filter(type=="neuroma") %>% mutate(perc = frequency/sum(frequency))
id_count_table <- bind_rows(id_count_table_nerve,id_count_table_neuroma)

histogram <- ggplot(id_count_table, aes(x = type, y = perc, fill = cluster)) +
  geom_bar(stat = "identity", color = "black") + 
  theme_classic() + 
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, size=14, hjust = 1),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent")) +
  scale_x_discrete(labels = c('nerve'='Nerve', 'neuroma'='Neuroma')) + 
  labs(fill='Cluster')

ggsave("vascular_clusters_histogram.png", plot = histogram, width = 2.5, height = 5.5, path = img.path, bg = "transparent")

# heatmap
sheets <- excel_sheets(here('outs', 'reanalysis', 'vascular', 'tables',"vascular_markers_clusters_seuratV4.xlsx"))
top5_genes <- data.frame()
for (i in 1:length(sheets)) {
  markers.clusters <- read_xlsx(here('outs', 'reanalysis', 'vascular', 'tables',"vascular_markers_clusters_seuratV4.xlsx"), sheet = sheets[[i]])
  top5_genes_cluster <- data.frame(genes = markers.clusters$gene[1:5],
                                   name = markers.clusters$cluster[1:5])
  top5_genes <- rbind(top5_genes,top5_genes_cluster)
}

vascular <- ScaleData(vascular, features=top5_genes$genes)
heatmap <- DoHeatmap(vascular, features = top5_genes$genes, size=3, slot = "scale.data",disp.min = -2,disp.max = 2) + 
  NoLegend() + theme(plot.margin= unit(c(1, 2, 1, 1), "cm"), axis.text.y = element_text(face = "italic"),
                     panel.background = element_rect(fill = "transparent"),
                     plot.background = element_rect(fill = "transparent"))
ggsave("vascular_markers_heatmap.png",p=heatmap,path = img.path, width = 6, height = 5, limitsize = FALSE, bg='transparent')


