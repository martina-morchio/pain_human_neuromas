############ Prepare snRNAseq data for cell type deconvolution ##########
library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
library(Giotto)
options(future.globals.maxSize = 1e10)
options(Seurat.object.assay.version = "v5")

# this script is to generate giotto object from single nuclei dataset
# in order to run cell-type deconvolution with PAGE

sc_path <- here('outs','reanalysis', 'all', 'nerves.comb.qs') #use whole dataset
giotto_path <- here('outs','reanalysis', 'giotto')
if (!dir.exists(giotto_path)) {dir.create(giotto_path)}

neuromas_sc<- qread(sc_path)
Idents(neuromas_sc) <- "new_clustering"

sc_expression_norm <- GetAssayData(neuromas_sc, assay = "RNA", slot = "data") # normalized data matrix
meta <- neuromas_sc@meta.data

# save matrices
write.table(sc_expression_norm, file=paste0(giotto_path,'/sc_expression_norm.tsv'), quote=FALSE, sep='\t')
write.table(meta, file=paste0(giotto_path,'/sc_meta.txt'), quote=FALSE, sep='\t')

# read sc data
meta<-read.table(here(giotto_path, "sc_meta.txt"))
sc_expression_norm<-read.table(here(giotto_path, "sc_expression_norm.tsv"),header = T,row.names = 1)

# create Giotto object with single cells
my_python_path = "C:/Users/mmorchio/AppData/Local/r-miniconda/envs/giotto_env/python.exe" # set to NULL to use previously installed giotto environment
results_folder = here('outs','reanalysis', 'giotto')
if (!dir.exists(results_folder)) {dir.create(results_folder)}

instrs = createGiottoInstructions(python_path = my_python_path,
                                  save_dir = results_folder)
neuromas_sc <- createGiottoObject(raw_exprs = sc_expression_norm,instructions = instrs)
neuromas_sc <- addCellMetadata(neuromas_sc,
                               new_metadata = meta)
neuromas_sc <- normalizeGiotto(gobject = neuromas_sc, scalefactor = 6000, verbose = T)

# Save object
saveGiotto(gobject = neuromas_sc, method='qs', dir=results_folder, foldername = 'giotto_object')
