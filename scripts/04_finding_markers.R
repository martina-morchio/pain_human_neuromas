######## Finding markers for annotation ##########
library(writexl)
library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
options(future.globals.maxSize = 1e10)
options(Seurat.object.assay.version = "v5")

# Identifying top DE genes in each cluster for annotation

# Set working environment
filefolder<-"/mnt/parscratch/users/mdq19mm/nucseq/combined_runs/integration/"
setwd(filefolder)
images <- "images/"
markers_folder <- "images/markers"

# Load integrated data and normalise
nerves.comb.sct <- readRDS("nerves_comb_SCT.rds")
Idents(nerves.comb.sct) <- "integrated_snn_res.0.5"
DefaultAssay(nerves.comb.sct) <- "RNA"
nerves.comb.sct <- NormalizeData(nerves.comb.sct)
nerves.comb.sct <- ScaleData(nerves.comb.sct, verbose = FALSE)

# Find markers
nerves.comb.sct.markers <- FindAllMarkers(nerves.comb.sct, only.pos = TRUE, logfc.threshold = 0)

signif.markers <- nerves.comb.sct.markers %>%
  filter(p_val_adj < 0.05) 

nerves.comb.sct.markers.clusters<- list()

for (i in 1:length(levels(nerves.comb.sct))) {
  nerves.comb.sct.markers.clusters[[i]] <- signif.markers %>% 
    filter(cluster==levels(nerves.comb.sct)[[i]]) %>% arrange(desc(avg_log2FC))
}

write_xlsx(nerves.comb.sct.markers.clusters, path=paste0(images,"nerves_markers_clusters_seuratV4.xlsx"))

########## Annotation with known markers ##########
# Generate dotplots with known markers for each cell type to confirm annotation
# Markers are obtained from l;iterature, view supplementary materials for more information

markers <- list(
  fibroblasts = c("DCN", "MFAP5", "GSN", "VIM","COL1A1","FGFR1","FN1"),
  epineurial.fibroblasts = c("SFRP2", "DPT", "PCOLCE2", "ADAMTS5", "PI16", "SFRP4", "PRRX1", "COMP", 
                             "LY6C1", "ZFHX4", "POSTN", "GDF10", "BMP4", "BMPER", "FGF10", "DPP4", 
                             "PTGES", "CXCL1", "TSPAN11", "PDGFRB"),
  endoneurial.fibroblasts = c("OSR2","ABCA8","ABCA9","ABCA10", "PLXDC1", "ALDH1A1","COL15A1"),
  differentiating.fibroblasts = c("DLK1", "CILP", "PLAGL1", "PTN"),
  injury.responsive.fibroblasts = c("PRRX1", "PDGFRA"),
  perineurial = c("CLDN1", "SLC2A1", "PTCH1", "LMO7","ITGB4", "KLF5"),
  
  endothelial = c("EGFL7", "PECAM1", "TIE1", "EMCN", "CDH5", "VWF", "CLDN5", "ECSCR"),
  epineurial.endothelial = c("SOX17", "SPOCK2", "RGCC"),
  endoneurial.endothelial = c("LRG1", "ICAM1"),
  lymphatic.endothelial = c("LYVE1", "MMRN1", "FLT4", "PROX1"),
  meningeal.dural= c("FXYD5", "DKK2", "ALPL", "CRABP2", "COL1A1"),
  meningeal.arachnoid= c("COL1A1","OGN","ALDH1A2","PTGDS"),
  meningeal.pial= c("COL1A1", "S100A6","LAMA1","NGRF"),
  
  proliferating = c("MKI67", "TOP2A", "FOXM1"),
  vascular.smooth.muscle = c("DES", "TPM2", "MYH11", "ACTA2", "MYLK", "MYOM1", "MYOCD"),
  pericytes = c("RGS5", "KCNJ8", "PDGFRB","OCLN","CLDN12","JAM1","TJP1","TJP2"),
  mesenchymal = c("ACTA2", "PDGFRA", "PDGFRB", "CSPG4", "THY1", "CD34", "ADAM12", "TWIST1"),
  endoneurial.mesenchymal = c("ETV1", "OSR2"),
  
  schwann = c("SOX10", "PLP1", "ERBB3", "NCAM1", "S100B"),
  non.myelinating = c("L1CAM", "NRXN1","NCAM1"),
  myelinating = c("MBP", "MPZ", "EGR2", "NCMAP"),
  
  immune = c("PTPRC", "CD52"),
  macrophages = c("AIF1", "CD68", "MRC1", "ADGRE1", "SIGLEC1",  "ITGAM"),
  epineurial.macrophages = c( "CLEC10A"),
  monocytes = c("CXCR1", "CSF1R", "CD300A", "CLEC4A"),
  neutrophils = c("S100A8", "S100A9", "CXCR2", "CXCL2"),
  mast.cells = c("KIT","CPA3","HPGDS","VWA5A","MS4A2","GATA2"),
  t.cells = c("CD3G", "CXCR6", "TRAC", "CD3E", "SKAP1", "THEMIS", "IL7R"),
  nk.cells = c("NKG7", "KLRK1", "NCR1"),
  b.cells = c("BANK1", "CBFA2T3", "TAOK3", "MS4A1", "CD19", "CD79A"),
  
  myocytes=c("ITGA3","MYH14","MYL1","TNNT1","TNNT3","TNNI1"),
  
  oligodendrocytes = c("OLIG2","OLIG1","MOG","CNP","PLP1"),
  astrocytes = c("GFAP"),
  salivary.gland.cells=c("MUC5B","BHLHA15","AQP5","KRT19","SLPI","KRT7","KRT14"),
  adipocytes = c("AQP7","ADIPOQ","GPAM","PLIN1")
)

# make dotplot of markers for each cell type

for (i in 1:length(markers)) {
  dotplot<-DotPlot(nerves.comb.sct,features=markers[[i]]) + labs(title=names(markers)[[i]]) +
    theme(axis.text.x = element_text(angle = 90), axis.title = element_blank())
  filename <- paste0(names(markers)[[i]],"_markers_")
  height <- ceiling(length(markers[[i]])/4)*5
  width<-length(markers[[i]])/8 + 4
  ggsave(paste0(filename,"dotplot.pdf"),p=dotplot,path = markers_folder, width = width, height = 5, limitsize = FALSE)
}
