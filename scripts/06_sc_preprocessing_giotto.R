####### Processing of single cell data for cell-type deconvolution with Giotto ############
library(Giotto)

# this script is to generate giotto object from single nuclei dataset, in order to run cell-type deconvolution with PAGE
# a Giotto object is created from count matrix and metadata generated in 03_snRNAseq_integration_rpca.R

# read sc data
sc_meta<-read.table("/mnt/parscratch/users/mdq19mm/giotto/sc_meta.txt")
sc_data<-read.table("/mnt/parscratch/users/mdq19mm/giotto/sc_expression_norm.tsv",header = T,row.names = 1)

# create Giotto object with single cells
instrs = createGiottoInstructions(python_path = my_python_path)
neuromas_sc <- createGiottoObject(raw_exprs = sc_data,instructions = instrs)
neuromas_sc <- addCellMetadata(neuromas_sc,
                               new_metadata = sc_meta)
neuromas_sc <- normalizeGiotto(gobject = neuromas_sc, scalefactor = 6000, verbose = T)

# Save object
saveGiotto(gobject = neuromas_sc,
           dir = "/mnt/parscratch/users/mdq19mm/giotto/combined_runs/objects",
           foldername = 'nerves_sc')