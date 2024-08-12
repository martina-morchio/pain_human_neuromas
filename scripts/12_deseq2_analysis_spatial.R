######## DE Analysis on local machine #############
library(dplyr)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(tibble)
library(PCAtools)
library(ggrepel)

# This script is to run DE analysis on aggregated pseudo-bulk dataset from Visium and get images for visualisation

working.directory<-"/Users/martina/Documents/PhD/Main_projects/Visium/analysis/DEseq2"
setwd(working.directory)
images<- "images/"
DE_tables <- "DE_tables/"

# get aggregated matrix from 10_pseudobulk_preprocessing.R
counts<-read.table("spatial_aggr.txt")
metadata<-read.table("spatial_metadata.txt")
metadata <- metadata[order(row.names(metadata)),]

# create DEseq Dataset and collapse technical replicates
dds <- DESeqDataSetFromMatrix(counts, 
                              colData = metadata, 
                              design = ~ sex +pain.status)
dds<- collapseReplicates(dds, groupby = dds$sample)

####### Run differential expression analysis with DEseq2 #######
dds <- DESeq(dds)
vst <- assay(vst(dds))

# PCA plot
p <- pca(vst, metadata = colData(dds), removeVar = 0.1)
biplot(p, labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
ggsave("biplot_collapsed.pdf", path = images, width=8, height=8)

# Plot dispersion estimates
plotDispEsts(dds)

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, 
               name = "pain.status_P_vs_NP",
               alpha = 0.05)

# Shrink the log2 fold changes to be more appropriate using the apeglm method 
# should cite Zhu, A., Ibrahim, J.G., Love, M.I. (2018) https://doi.org/10.1093/bioinformatics/bty895 when using this method

res <- lfcShrink(dds, 
                 coef = "pain.status_P_vs_NP",
                 res=res,
                 type = "apeglm")

# Turn the DESeq2 results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl 

# Write all results to file
write.csv(res_tbl,
          paste0(working.directory,"/",DE_tables,"P_vs_NP_all_genes.csv"),
          quote = FALSE,
          row.names = FALSE)

#table of significant genes only
padj_cutoff <- 0.05 # Set thresholds
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Write significant results to file
write.csv(sig_res,
          paste0(working.directory,"/",DE_tables,"P_vs_NP_signif_genes_.csv"),
          quote = FALSE,
          row.names = FALSE)

#### Make Volcano plot ####
log2fc_cutoff<-0.5
res_table_thres <- res_tbl[!is.na(res_tbl$padj), ] %>% 
  mutate(threshold = padj < padj_cutoff)
res_table_thres$genelabels <- ifelse(res_table_thres$gene %in% sig_res$gene, T, F)

volplot <- ggplot(res_table_thres) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold), size=1) +
  ggtitle("Volcano plot of painful versus non painful human neuromas") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj)),
                  label = ifelse(res_table_thres$genelabels == TRUE, as.character(res_table_thres$gene),""), 
                  max.overlaps = 70, force=30) + 
  scale_y_continuous(limits = c(0, 8)) +
  scale_color_manual(values = c("grey60", "red3"), name= NULL, labels = c("Not significant","P_adj < 0.01")) + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff),
             linetype = "dashed") +
  xlim(-5,+5) +  theme_bw() +  theme(plot.title = element_text(face = "bold", size = 12, hjust=0.5))

ggsave("volcanoplot_nerve_fascicles.pdf", plot= volplot,path = images, width=6,height=5)

##### Images from single cell data ##### 
library(ggplot2)
library(tidyr)
library(forcats)

# validate DE genes in snRNAseq data
# use neuroma subset, generated in 03_snRNAseq_integration_rpca.R
images<- "./images/single_cell"
neuromas_sc<- readRDS("/Volumes/shared/boissonade/User/mdq19mm/snRNAseq/combined_runs/neuromas.rds")
DefaultAssay(neuromas_sc) <-"RNA"

## Use aggregated expression from each cluster
genes<- res_tbl %>% filter(padj<0.01, abs(log2FoldChange)>1.3) %>% arrange(desc(log2FoldChange)) %>% select(gene,log2FoldChange,padj)
aggr<- AggregateExpression(neuromas_sc, assays="RNA", features=genes$gene, group.by="id", normalization.method = "LogNormalize")
aggr_df <- as.data.frame(aggr$RNA) %>% rownames_to_column("gene")
aggr_log2fc<- inner_join(aggr_df, genes)

df<-aggr_df %>% gather(key="sample",value="aggregated.expression", -gene)
df$gene <- factor(df$gene)
df$sample<-factor(df$sample)
df$sample<-gsub("n1","NP",df$sample)
df$sample<-gsub("n2","P",df$sample)
ggplot(df, aes(color=sample)) + geom_point(aes(x=fct_inorder(gene), y=aggregated.expression, color=sample, fill=sample)) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Aggregated expression of genes upregulated un painful neuromas in the snRNAseq dataset")
  
ggsave("aggr_expression_sc_deg.pdf",path=images,width=5,height=4)

## Dotplots
genes_oi <- aggr_log2fc # Filtering genes of interest to validate
Idents(neuromas_sc) <- "annotation"

genes_oi <- genes %>% filter(padj<0.01)
neuromas_df<- subset(neuromas_sc,idents=c("MenF_1","MenF_2","Astro","Endo_4","Oligo","Myo_1","Myo_2","SGC_1","SGC_2","SGC_3"), invert=TRUE)

p1<-DotPlot(neuromas_df, features = c(genes_oi$gene), dot.scale = 8) +
  RotatedAxis() + theme(axis.text.y = element_text(size=15),axis.text.x = element_text(size=15))
ggsave("genes_oi_dotplot.pdf",p=p1,path = "./images/single_cell", width = 6, height = 5, limitsize = FALSE)

##################### Images from spatial data ##################### 
filedir <- "/Users/martina/Documents/PhD/Main_projects/Visium/all_data/"
sample.ids<- list.dirs(filedir,full.names = FALSE, recursive = FALSE)
spatial_all<- readRDS("/Volumes/shared/boissonade/User/mdq19mm/Visium_2023/analysis/spatial_combined.rds")
names(spatial_all) <-sample.ids
images <- "images/spatial/"
selected_samples<-c("LN1_D1","LN15_D1","LN2_C1","LN8_D1")
genes_oi_2 <- c("HLA-A","JMJD1C","NLRC5","CXCL2","TNFAIP3","ADGRG1","NOTCH1","KLF4","ANXA2","MMP19","NID2")

# Using only a few representative sections
spatial_all<- spatial_all[selected_samples]

# SCT transform so they're comparable
spatial_all <- lapply(spatial_all,SCTransform, assay = "Spatial", verbose = FALSE)

# Merge in one object
spatial.comb<-merge(spatial_all[[1]], list(spatial_all[[2]], spatial_all[[3]], spatial_all[[4]]))
DefaultAssay(spatial.comb) <- "SCT"
spatial.comb<- SCTransform(spatial.comb, assay = "Spatial")

# Get max expression value
counts<- GetAssayData(object = spatial.comb, assay = "SCT", slot = "data")

for (i in 11:length(genes_oi_2)) {
  
  n<-max(counts[genes_oi_2[i],])
  
  plot<-SpatialFeaturePlot(spatial.comb, features = genes_oi_2[i], alpha=0.5, pt.size.factor=1.2,image.alpha=0, stroke=0, crop = FALSE)
  
  for (j in 1:length(plot)) {
    plot[[j]]<- plot[[j]] + scale_fill_continuous(type="viridis",limits = c(0,n))
  }
  
  ggsave(paste0(genes_oi_2[i],"_transparent_spatial_plot.pdf"), p=plot, path =images, width = 20, height = 5, limitsize = FALSE)
}

