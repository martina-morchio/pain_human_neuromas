####### snRNAseq annotation updated ###########
library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
library(enrichR)
library(writexl)
library(ggrepel)
library(qs)
library(tibble)
library(tidyr)
library(here)
library(data.table)
options(future.globals.maxSize = 1e10)
options(Seurat.object.assay.version = "v5")

here::i_am('pain_human_neuromas.Rproj')

img.path<- here('outs', 'reanalysis', 'all', 'images')
tbl.path<- here('outs', 'reanalysis', 'all', 'tables')

paths<- c(img.path,tbl.path)
for (path in paths) {
  if (!dir.exists(path)) {dir.create(path, recursive = T)}
}

# making images for visualisation of snRNAseq data
# load annotated R object
nerves.comb <- readRDS(here('data','nerves_annotated.rds'))
immune<- qread(here('outs', 'reanalysis', 'immune', 'immune.qs'))
fibro<- qread(here('outs', 'reanalysis', 'fibro', 'fibro.qs'))
schwann<- qread(here('outs', 'reanalysis', 'schwann', 'schwann.qs'))
vascular<- qread(here('outs', 'reanalysis', 'vascular', 'vascular.qs'))
others <- nerves.comb %>% subset(idents=c('MenF_1', 'MenF_2', 'SGC_1', 'SGC_2', 'SGC_3', 'SGC_4', 'Myo', 'Oligo', 'Astro'))

# get all annotation
immune.meta <- immune@meta.data %>% select(new_clustering)
fibro.meta <- fibro@meta.data %>% select(new_clustering)
schwann.meta <- schwann@meta.data %>% select(new_clustering)
vascular.meta <- vascular@meta.data %>% select(new_clustering)
others.meta <- others@meta.data %>% select(annotation) %>% rename('new_clustering'='annotation')

new.meta<- rbind(immune.meta, fibro.meta, schwann.meta, vascular.meta, others.meta)

# combine with old
nerves.comb@meta.data$new_clustering <- new.meta[rownames(nerves.comb@meta.data), "new_clustering"]

# use new annotation
Idents(nerves.comb) <- "new_clustering"

# reorder idents
order<- c("mSC_1","mSC_2", "nmSC", "damage_mSC", "damage_nmSC", "Astro", "Oligo", 'EndoF', 'PeriF_1', 
          'PeriF_2', 'PFF', 'MenF_1', 'MenF_2', "Endo", "Infl_Endo","Lymph_Endo","SMC","Pericytes",
          "Macro", "Lympho", "SGC_1", "SGC_2", "SGC_3", "SGC_4", "Myo")
levels(nerves.comb)<- order
nerves.comb<- nerves.comb %>% NormalizeData()

# make images
# UMAPs
DefaultAssay(nerves.comb) <- "integrated"
p1 <- DimPlot(nerves.comb, reduction = "umap", group.by = "type", order = c("neuroma","nerve")) + NoAxes() + ggtitle("")
p2 <- DimPlot(nerves.comb, reduction = "umap", label = TRUE, repel = TRUE, label.box = TRUE, label.size = 3) + NoAxes() + NoLegend()
umap<-cowplot::plot_grid(plotlist=list(p1,p2),nrow = 1)
ggsave("all_UMAP.png", p=umap, path=img.path, width = 12, height = 6, limitsize = FALSE)
p3 <- DimPlot(nerves.comb, reduction = "umap", split.by = "type") + NoAxes()
ggsave("splitbytype_UMAP.pdf", p=p3, path=img.path, width = 12, height = 6, limitsize = FALSE)

# QC by cluster
QC<-VlnPlot(nerves.comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "new_clustering", split.plot = TRUE)
ggsave("QC.pdf",plot=QC,width = 20, height = 5, path=img.path)

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

ggsave("nerves_clusters_histogram.png", plot=histogram, width= 5, height=7, path=img.path)

# save object
qsave(nerves.comb, here('outs', 'reanalysis', 'all', 'nerves.comb.qs'))
nerves.comb<- qread(here('outs', 'reanalysis', 'all', 'nerves.comb.qs'))

# print summary stats
meta <- nerves.comb@meta.data %>% rownames_to_column('cellID')
meta_by_sample <- meta %>% group_by(id, new_clustering) %>% summarise(n = n(), .groups = 'drop') %>% 
  pivot_wider(names_from = new_clustering, values_from = n, values_fill = 0)
meta_by_type <- meta %>% group_by(type, new_clustering) %>% summarise(n = n(), .groups = 'drop') %>% 
  pivot_wider(names_from = new_clustering, values_from = n, values_fill = 0)

# save
write_xlsx(list(sample = meta_by_sample, type = meta_by_type), path = here(tbl.path,"cell_type_distrib.xlsx"))
fwrite(meta, file= here(tbl.path, "snRNAseq_meta.csv"))
