########## Heatmap spatial clusters identified with Giotto ##########
library(dittoSeq)
library(scater)
library(patchwork)
library(cowplot)
library(viridis)
library(SingleCellExperiment)

working_directory<-"/Users/martina/Documents/PhD/Main_projects/Visium/analysis/giotto"
setwd(working_directory)
images<-"images/heatmap/"

# Read sce experiment generated with 07_sce_object_from_giotto
spe<- readRDS("giotto_integrated_visium_sce_norm.rds")

# select colour for each cluster
ann_cols <-c("Bcells"= '#B01E3A',
             "Endo"= '#ff9a36',
             "Fibro"= '#CCB1F1',
             "Macro"= 'lightpink',
             "Endo"= '#ff9a36',
             "Myo1"= '#D4D915',
             "Myo2"= '#E6C122',
             "Myo3"= '#AC8F14',
             "Myo4"= 'chocolate4',
             "Peri"= '#B95FBB',
             "SC1"= '#28CECA',
             "SC2"= '#25aff5',
             "SC3"= '#A4DFF2',
             "SC4"= '#03FCA5',
             "SC5"= '#2236E1',
             "SC6"= '#31C53F')

# subset to remove clusters 16 and 17
clusters<-names(ann_cols)
spe <- subset(spe, ,annotation %in% clusters)

## subsampling 4000 cells
set.seed(220818)
cur_cells <- sample(seq_len(ncol(spe)), 8000)

# select top genes for each cluster
marker_genes<- c("COL1A1", "SFRP2", "FBLN1", "COL1A2", "SFRP4", "AQP1", "CCL14", "IL6", 
"TM4SF1", "SELE", "HBA2", "MPZ", "HBA1", "PMP22", "MBP", "MB", 
"TNNT1", "TCAP", "CKM", "TTN", "PTGDS", "CLDN1", "SLC2A1", "IGFBP6", "PRX", "S100B", 
"TNNI1", "APOD", "ACTC1", "THBS4", "MYLPF", "CA3", "IGKC", "IGHG2", "IGHG1", "LYZ", "MMP9", 
"SPP1", "LAPTM5", "CHIT1")

order<- c("SC1","SC2","SC3","SC4","SC5","SC6","Fibro","Peri","Endo","Bcells","Macro","Myo1","Myo2","Myo3","Myo4")
levels(colData(spe)$annotation) <- order

## make heatmap with Ditto heatmap
# aggregate cell types
celltype_mean <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),  
                                      ids = spe$annotation, 
                                      statistics = "mean",
                                      use.assay.type = "counts", 
                                      subset.row = marker_genes)

aggr<- dittoHeatmap(celltype_mean,
             assay = "counts", 
             cluster_cols = FALSE, 
             scale = "none",
             heatmap.colors = viridis(100),
             annot.by = c("annotation"),
             annot.colors = ann_cols)

ggsave(plot=aggr,filename="heatmap_celltype_aggr.pdf",path=images,width=6,height=6)