######## Cell type distribution in spatial dataset ########

#### Get cell type proportions from giotto object (on HPC) ####
library(Giotto)
library(dplyr)
# load giotto object
testcombo = loadGiotto(path_to_folder = "/mnt/parscratch/users/mdq19mm/giotto/results/objects/testcombo",reconnect_giottoImage = TRUE)

meta<- pDataDT(testcombo) #get metadata
# establish cluster names
cluster_ann<- c("1"="Fibro",
                "2"="Endo",
                "3"="SC1",
                "4"="Myo1",
                "5"="Peri",
                "6"="SC2",
                "7"="Myo2",
                "8"="SC3",
                "9"="Myo3",
                "10"="SC4",
                "11"="Myo4",
                "12"="SC5",
                "13"="Bcells",
                "14"="SC6",
                "15"="Macro")
meta<- meta %>% mutate(annotation= recode(leiden_harmony, !!!cluster_ann)) #convert cluster numbers to names

df <- table(meta$list_ID, meta$annotation) # make dataframe with barcode and annotation
# save csv
write.csv(df, "/mnt/parscratch/users/mdq19mm/giotto/results/harmony_cells_proportions/cell_type_distribution.csv", quote=F)  

######## Differential abundance analysis #########
# as described in Amezquita et al 2019
library(edgeR)
library(tidyverse)

working_directory <- "/Users/martina/Documents/PhD/Main_projects/Visium/analysis/giotto"
setwd(working_directory)
images <- "images/cell_type_distribution/"

# get cell type annotation obtained from above and tidy data frame
abundances <- read.csv("cell_type_distribution.csv") %>% t() %>% as.data.frame()
colnames(abundances) <- abundances[1,]
abundances<- abundances[-1,] %>% mutate_at(vars(1:27),as.numeric)

#adding metadata
samples<-c("LN1","LN1","LN1","LN1","LN12","LN12","LN13","LN13","LN15","LN15","LN15","LN15",
           "LN2","LN2","LN2","LN2","LN2","LN7","LN7","LN7","LN7","LN7","LN8","LN8","LN8","LN8","LN8")
status<-c("P","P","P","P","P","P","P","P","P","P","P","P",
          "NP","NP","NP","NP","NP","NP","NP","NP","NP","NP","NP","NP","NP","NP","NP")
merged<- data.frame(replicate=colnames(abundances), sample=samples, status=status) %>% column_to_rownames(var="replicate")

# create edgeR object and follow differential abundance protocol in Amezquita 2019
y.ab<- DGEList(abundances,samples=merged)
design <- model.matrix(~factor(status), y.ab$samples)

y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex=1)
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
plotQLDisp(fit.ab, cex=1)

res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
DA<-data.frame(topTags(res, n=17))

write.csv(DA, paste0(images,"DA_analysis_res.csv"), quote=F)

#### Barplots of cell-type distribution #####
abundances <- read.csv("cell_type_distribution.csv")
abundances$sample <- c("LN1","LN1","LN1","LN1","LN12","LN12","LN13","LN13","LN15","LN15","LN15","LN15",
               "LN2","LN2","LN2","LN2","LN2","LN7","LN7","LN7","LN7","LN7","LN8","LN8","LN8","LN8","LN8")
abundances <- abundances %>% dplyr::select(-X) %>% group_by(sample) %>% summarize_all(sum, na.rm=TRUE) %>% 
  as.data.frame() %>% dplyr::select(sample,SC1,SC2,SC3,SC4,SC5,SC6,Fibro,Peri,Endo,Bcells,Macro,Myo1,Myo2,Myo3,Myo4)

abundances<- abundances %>% rowwise() %>% mutate(total=sum(c_across(2:16))) %>% 
  ungroup() %>% mutate(across(2:16, ~./total))

histo_data<- data.frame()
groupings<- c("SC", "SC", "SC", "SC", "SC", "SC", "Fibro", "Fibro", "Endo", "Immune", "Immune", "Myo", "Myo", "Myo", "Myo")
for (i in 1:length(abundances$sample)) {
  histo_data_new<- data.frame(sample=as.factor(abundances$sample[i]),
                              celltype=as.factor(colnames(abundances)[2:16]),
                              groups=groupings,
                              freq=as.numeric(abundances[i,2:16]))
  histo_data <- rbind(histo_data,histo_data_new)
}

histo_data<- group_by(histo_data,groups)
histo_data$celltype <- fct_relevel(histo_data$celltype, as.character(histo_data$celltype[1:15]))

cell_colours <- c("Fibro"= '#CCB1F1',
                  "Endo"= '#ff9a36',
                  "SC1"= '#28CECA',
                  "Myo1"= '#D4D915',
                  "Peri"= '#B95FBB',
                  "SC2"= '#25aff5',
                  "Myo2"= '#E6C122',
                  "SC3"= '#A4DFF2',
                  "Myo3"= '#AC8F14',
                  "SC4"= '#03FCA5',
                  "Myo4"= 'chocolate4',
                  "SC5"= '#2236E1',
                  "Bcells"= '#B01E3A',
                  "SC6"= '#31C53F',
                  "Macro"= 'lightpink')


histogram<- ggplot(histo_data, aes(x=sample, fill=celltype, y=freq))+
  geom_col(position="stack") +theme_classic() +
  scale_fill_manual(values = cell_colours) + facet_wrap(~ groups, nrow = 1) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  labs(y="Frequency") 

ggsave("cell_type_proportion.pdf",plot = histogram,
       path = images, width = 14,height=6)

######### boxplots of change in all clusters#########
sig_cells <- DA %>% filter(FDR<0.01)
celltypes<- rownames(sig_cells)

for (i in 1:length(celltypes)) {
  box_data <- histo_data %>% filter(celltype==celltypes[i])
  box_data$status <- c("P","P","P","P","NP","NP","NP")
  pval<- sig_cells %>% filter(rownames(sig_cells)==celltypes[i]) %>% pull(FDR) %>% round(digits=8)
  box_data %>%
    ggplot( aes(x=status, y=freq, fill=status)) +
    geom_boxplot() +
    scale_fill_manual(values=c("P"="#E84B1E","NP"="#2A4C9C")) +
    geom_point(color="black", size=1,position=position_jitter(height=0, width = 1e-01)) +
    theme_classic() +
    theme(legend.position="none", plot.title = element_text(size=9),plot.subtitle = element_text(size=7), vjust=0) + 
    ggtitle(celltypes[i], subtitle = paste0("FDR=",pval)) +xlab("")
  ggsave(paste0(celltypes[i],"_abundance_boxplot.pdf"), path=images, width=2,heigh=3)
}

####### GSEA analysis of DA cell types ##########
library(topGO)
library(GO.db)
library(biomaRt)
library(Rgraphviz)
library(dplyr)

images <- "images/enrichment/"
# get de genes
marker_genes<- read.csv("scran_markers_clusters.csv")
cluster_ann<- c("1"="Fibro",
                "2"="Endo",
                "3"="SC1",
                "4"="Myo1",
                "5"="Peri",
                "6"="SC2",
                "7"="Myo2",
                "8"="SC3",
                "9"="Myo3",
                "10"="SC4",
                "11"="Myo4",
                "12"="SC5",
                "13"="Bcells",
                "14"="SC6",
                "15"="Macro")

# annotate cluster with cell-type
marker_genes<- marker_genes %>% mutate(celltype= recode(cluster, !!!cluster_ann)) %>% filter(!is.na(celltype))

## prepare objects for GO analysis
# get background genes
bg <- read.csv("background_genes_GSEA.csv")

# create GO db for genes to be used using biomaRt - please note that this takes a while
ensembl <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
go_bg_ids <- getBM(attributes=c("go_id","external_gene_name", "namespace_1003"), 
                   filters="external_gene_name", 
                   values=bg, 
                   mart=ensembl)

# build the gene 2 GO annotation list (needed to create topGO object)
gene_2_GO=unstack(go_bg_ids[,c(1,2)])

# GO analysis for each cluster

for (i in 1:length(cluster_ann)) {
  # get cluster name
  cell<- as.character(cluster_ann[i])
  res_tbl <- marker_genes %>% filter(celltype==cell)
  # make named factor showing which genes are of interest
  geneList=res_tbl$FDR
  names(geneList)= res_tbl$feats
  geneList<- geneList[!is.na(geneList)]
  
  # make function to select genes with FDR< 0.001
  topDiffGenes <- function(allScore) {
    return(allScore < 0.001)
  }
  # Make topGO data object
  sampleGOdata_BP<- new("topGOdata", ontology="BP",
                        allGenes=geneList, geneSel=topDiffGenes, 
                        nodeSize=5,
                        annot=annFUN.gene2GO, gene2GO=gene_2_GO)
  
  # Run test
  resultFisher <- runTest(sampleGOdata_BP,algorithm="classic",statistic="fisher")
  
  ## Analysis of results of biological process
  allGO=usedGO(sampleGOdata_BP)
  allRes <- GenTable(sampleGOdata_BP, fisher=resultFisher, 
                     orderBy = "fisher", topNodes = length(allGO), numChar=1000)
  # Add significant genes to table
  sel.terms <- sample(usedGO(sampleGOdata_BP), 10)
  ann.genes <- genesInTerm(sampleGOdata_BP)
  
  for (j in 1:length(ann.genes)) {
    ann.genes[[j]] <- ann.genes[[j]][ann.genes[[j]] %in% sigGenes(sampleGOdata_BP)]
  }
  
  allRes$genes<- NA
  ann.genes<-replace(ann.genes,ann.genes=="","NA")
  for (i in 1:length(ann.genes)) {
    GOid<-names(ann.genes)[i]
    genes<-paste0(ann.genes[i])
    genes<- gsub("c\\(","",genes)
    genes<- gsub("\\)","",genes)
    genes<- gsub('"',"",genes)
    allRes[which(allRes$GO.ID==GOid),"genes"]<-ifelse(length(genes)==0,"NA", genes)
  }
  
  #performing FDR correction on our p values
  p.adj=round(p.adjust(allRes$fisher,method="fdr"),digits = 4)
  
  # create the file with all the statistics from GO analysis
  all_res_final=cbind(allRes,p.adj)
  all_res_final=all_res_final[order(all_res_final$p.adj),]
  
  # add enrichment score
  all_res_final <- all_res_final %>% mutate("enrichment"=Significant/Expected) %>% filter(p.adj<0.01) %>% arrange(desc(enrichment))
  
  # save table
  write.table(all_res_final, file=paste0("images/enrichment/",cell,"_BP_fisher_test.txt"), sep="\t",row.names = F)
  
  ##### Barplot
  data<- all_res_final %>% head(n=15) %>% dplyr::select(GO.ID,Term,Annotated,Significant,fisher,p.adj,enrichment) %>% 
    mutate(across(everything(), ~replace(., .== 0e+00,2.220446e-16)))
  barplot<- ggplot(data, aes(x=reorder(Term, enrichment), y=enrichment, fill=p.adj)) + geom_bar(stat = "identity") + 
    coord_flip() +theme_classic() + theme(axis.title.y=element_blank())
  ggsave(paste0(cell,"_barplot_enrichment.pdf"),plot=barplot,path=images, width=6, height=3)
  
}


