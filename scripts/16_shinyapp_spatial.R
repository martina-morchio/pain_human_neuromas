######### Creating spatial shiny app #########
library(Seurat)
library(ShinyCell2)
library(here)
library(qs)
library(shiny)
here::i_am('pain_human_neuromas.Rproj')

# making shinyapp
citation<- list(
  author = 'Morchio, M., Sankaranarayanan, I., Tavares-Ferreira, D., et al.',
  title = 'Investigation of cellular and molecular changes linked with neuropathic pain in healthy and injured human trigeminal nerves.',
  journal = 'bioRxiv', 
  doi = 'https://doi.org/10.1101/2024.10.05.616798')

# config for snRNAseq data
sn<- qread(here('data', 'shinyapp','snRNAseq', 'nerves_annotated.qs'))
config_snrna <- createConfig(sn)
config_snrna = modDefault(config_snrna, "cell_type", "nFeature_RNA")

# let's create config for two of the spatial objects
sp <- qread(here('outs', 'spatial', 'spatial_annotated.qs'))
DefaultAssay(sp) <- 'Spatial'
samples.images <- c("LN2_C1","LN8_D1","LN1_A1","LN12_B1") 
ln<- list()
config_spatial <- list()

for (sample in samples.images) {
  cells_keep <- Cells(sp)[sp$list_ID == sample]
  ln[[sample]] <- subset(sp, cells = cells_keep) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims=1:30)
  config_spatial[[sample]] <- createConfig(ln[[sample]])
  config_spatial[[sample]] = modDefault(config_spatial[[sample]], "cell_type", "nr_feats")
}

# let's make shiny files
makeShinyFiles(sn, config_snrna, shiny.prefix="sn", , default.gene1 = 'PTGDS', default.gene2= 'MBP',
               shiny.dir="neuroma_atlas/")

for (sample in samples.images) {
  makeShinyFiles(ln[[sample]], scConf = config_spatial[[sample]], dimred.to.use = "umap", shiny.prefix = sample,
                 shiny.dir = "neuroma_atlas/", default.gene1 = "MBP", 
                 default.gene2 = "PTGDS", default.multigene = NA, default.dimred = NA)
}

# make shiny app code
makeShinyCodes(shiny.title="Atlas of human trigeminal nerves and neuromas ",shiny.prefix=c("sn",samples.images),
               shiny.headers = c("snRNAseq Nerves vs Neuromas", samples.images),
               shiny.dir="neuroma_atlas/", shiny.footnotes = citation)
