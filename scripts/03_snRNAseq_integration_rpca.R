####### Integrative analysis using Seurat v5 ###########
library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
options(future.globals.maxSize = 1e10)
options(Seurat.object.assay.version = "v5")

# integration of single nuclei from neuromas and trigeminal nerve root samples using Seurat's rpca method
working_directory <- "/mnt/parscratch/users/mdq19mm/nucseq/combined_runs/integration"
setwd(working_directory)
filefolder <- "/mnt/parscratch/users/mdq19mm/nucseq/combined_runs/integration/"
images <- "images/"
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

##### Load data #####
# Load filtered snRNAseq samples
n1 <- readRDS(paste0(filefolder,"n1.rds"))
n2 <- readRDS(paste0(filefolder,"n2.rds"))
tg1 <- readRDS(paste0(filefolder,"tg1.rds"))
tg2 <- readRDS(paste0(filefolder,"tg2.rds"))

# Add metadata for each sample
n1$id <-  "n1"
n1$type <-  "neuroma"
n1$sex <-  "F"
n1$pain.status <-  "NP"

n2$id <-  "n2"
n2$type <-  "neuroma"
n2$sex <-  "F"
n2$pain.status <-  "P"

tg1$id <-  "tg1"
tg1$orig.ident <- "trigeminal_nerve_root"
tg1$type <-  "nerve"
tg1$sex <-  "F"
tg1$pain.status <-  "NP"

tg2$id <-  "tg2"
tg2$orig.ident <- "trigeminal_nerve_root"
tg2$type <-  "nerve"
tg2$sex <-  "M"
tg2$pain.status <-  "NP"

# merge samples
nerves.comb <- merge(n1,c(n2,tg1,tg2), add.cell.ids = c("n1","n2","tg1","tg2"), merge.data=TRUE)

# split the dataset into a list of two seurat objects (nerve and neuroma)
nerves.list <- SplitObject(nerves.comb, split.by = "type")

# SCT normalization of each sample type
nerves.list.SCT <- lapply(X = nerves.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = nerves.list.SCT, nfeatures = 3000)
nerves.list.SCT <- PrepSCTIntegration(object.list = nerves.list.SCT, anchor.features = features)
nerves.list.SCT <- lapply(X = nerves.list.SCT, FUN = RunPCA, features = features)

# reciprocal PCA integration with SCT normalization
nerves.anchors <- FindIntegrationAnchors(object.list = nerves.list.SCT, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "rpca")
nerves.comb.sct <- IntegrateData(anchorset = nerves.anchors, normalization.method = "SCT", dims = 1:30)

# Run the standard workflow for visualization and clustering
nerves.comb.sct <- ScaleData(nerves.comb.sct, verbose = FALSE)
nerves.comb.sct <- RunPCA(nerves.comb.sct, npcs = 30, verbose = FALSE)
nerves.comb.sct <- RunUMAP(nerves.comb.sct, reduction = "pca", dims = 1:30)
nerves.comb.sct <- FindNeighbors(nerves.comb.sct, reduction = "pca", dims = 1:30)
nerves.comb.sct <- FindClusters(nerves.comb.sct, resolution = 0.5)

saveRDS(nerves.comb.sct,"nerves_comb_SCT.rds")
nerves.comb.sct <- readRDS("nerves_comb_SCT.rds")

######### Visualization #########

p1 <- DimPlot(nerves.comb.sct, reduction = "umap", group.by = "type")
p2 <- DimPlot(nerves.comb.sct, reduction = "umap", label = TRUE, repel = TRUE)
umap<-cowplot::plot_grid(plotlist=list(p1,p2),nrow = 1)
p3 <- DimPlot(nerves.comb.sct, reduction = "umap", split.by = "type")
umap<-cowplot::plot_grid(plotlist=list(umap,p3),nrow = 2, labels = c('A', 'B'))
ggsave("rpca_UMAP.pdf", p=umap, path = images, width = 12, height = 10, limitsize = FALSE)

# QC by cluster
QC<-VlnPlot(nerves.comb.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "integrated_snn_res.0.5")
ggsave("QC.pdf",plot=QC,width = 15, height = 5, path=images)

# Distribution across samples
id_count_table <- as.data.frame(table(nerves.comb.sct@meta.data$integrated_snn_res.0.5, nerves.comb.sct@meta.data$type))
colnames(id_count_table) <- c("cluster","id","frequency")
histogram<- ggplot(id_count_table, aes(x=id, y=frequency, fill=cluster))+geom_bar(position="fill",stat="identity") +theme_classic()
ggsave("nerves_clusters_histogram.pdf", plot=histogram, width= 3, height=10, path=images)

# calculate propoertion of nuclei in each cluster
table(nerves.comb.sct$id, nerves.comb.sct$integrated_snn_res.0.5) %>% 
  write.csv(paste(working_directory,images,"cluster_distribution.csv",sep="/"),quote=F)

####### Annotating clusters after running 04_finding_markers.R script ######
# based on top DE in each cluster and on expression of known markers, the clusters are annotated with code below
# Cluster 26 is removed due to the high levels of neuronal genes indicating ambient RNA
# this is due to the experimental design, the neuroma and trigeminal nerve root samples were multiplexed with DRG samples, which led to this type of ambient RNA contamination
nerves.comb.sct<- subset(nerves.comb.sct, idents="26", invert=TRUE)

nerves.comb.sct <- RenameIdents(nerves.comb.sct,
                                "0"="mSC_1",
                                "2"="mSC_2",
                                "13"="mSC_3",
                                "1"="nmSC",
                                "3"="EndoF",
                                "5"="PeriF_1",
                                "7"="PeriF_2",
                                "8"="PeriF_3",
                                "9"="PeriF_4",
                                "6"="PFF",
                                "14"="MenF_1",
                                "19"="MenF_2",
                                "4"="Endo_1",
                                "16"="Endo_2",
                                "17"="Endo_3",
                                "18"="Endo_4",
                                "10"="Mural",
                                "11"="Macro",
                                "21"="Lympho",
                                "12"="SGC_1",
                                "23"="SGC_2",
                                "20"="SGC_3",
                                "15"="SGC4",
                                "22"="Myo",
                                "24"="Astro",
                                "25"="Oligo")

levels(nerves.comb) <- c("mMSC_1","mSC_2","mSC_3","nmSC","Astro","Oligo","EndoF","PeriF_1","PeriF_2","PeriF_3",
                         "PeriF_4","MenF_1","MenF_2","PFF","Endo_1","Endo_2","Endo_3","Endo_4","Mural","Macro",
                         "Lympho","SGC_1","SGC_2","SGC_3","SGC_4","Myo")
nerves.comb$annotation <- Idents(nerves.comb)
nerves.comb.sct$annotation <- Idents(nerves.comb.sct)

# further filtering of clusters
# small number of nuclei classified as meningeal, astrocytes or oligos from the neuroma samples are removed
nerves.comb.sct <- subset(nerves.comb.sct, subset = type =="neuroma" & 
                              (annotation=="MenF_1" | annotation=="MenF_2" | annotation=="Oligo" | annotation=="Astro"), invert=TRUE)
# small number of nuclei classified as salivary gland cells or myocytes from the nerve root samples are removed
nerves.comb.sct <- subset(nerves.comb.sct, subset = type =="nerve" & 
                              (annotation=="SGC_1" | annotation=="SGC_2" | annotation=="SGC_3" | annotation=="SGC_4" | annotation=="Myo"), invert=TRUE)

DefaultAssay(nerves.comb.sct) <- "RNA"
saveRDS(nerves.comb.sct,"nerves_annotated.rds")

######### Subset subtypes #########
schwann <- subset(nerves.comb.sct,idents=c("mSC_1","mSC_2", "mSC_3", "nmSC"))
fibro <- subset(nerves.comb.sct, idents=c("EndoF","PeriF_1","PeriF_2","PeriF_3","PeriF_4","PFF"))
vascular <- subset(nerves.comb.sct, idents=c("Endo_1","Endo_2","Endo_3","Endo_4","Mural"))
immune <- subset(nerves.comb.sct, idents=c("Macro","Lympho"))
other <- subset(nerves.comb.sct, idents=c("SGC_1","SGC_2","SGC_3","SGC_4","Myo","Oligo","Astro"))

saveRDS(schwann,"schwann.rds")
saveRDS(fibro,"fibro.rds")
saveRDS(vascular,"vascular.rds")
saveRDS(immune,"immune.rds")
saveRDS(other,"other.rds")

######## Separate TN and Neuroma objects #########
working_directory <- "/mnt/parscratch/users/mdq19mm/nucseq/combined_runs/integration"
setwd(working_directory)

nerves.comb.sct <- readRDS("nerves_annotated.rds")
Idents(nerves.comb.sct) <- "type"
nerves <- subset(nerves.comb.sct, ident="nerve")
neuromas <- subset(nerves.comb.sct, ident="neuroma")

saveRDS(nerves, "nerves.rds")
saveRDS(neuromas, "neuromas.rds")

############ Get matrix for cell type deconvolution ##########
sc_path <- '/mnt/parscratch/users/mdq19mm/nucseq/combined_runs/integration/nerves_annotated.rds' #use whole dataset
giotto_path <- '/mnt/parscratch/users/mdq19mm/giotto/combined_runs/'
neuromas_sc<- readRDS(sc_path)
Idents(neuromas_sc) <- "annotation"

sc_expression_norm <- GetAssayData(neuromas_sc, assay = "RNA", slot = "data") # normalized data matrix
meta <- neuromas_sc@meta.data

write.table(sc_expression_norm, file=paste0(giotto_path,'/sc_expression_norm.tsv'), quote=FALSE, sep='\t')
write.table(meta, file=paste0(giotto_path,'/sc_meta.txt'), quote=FALSE, sep='\t')
