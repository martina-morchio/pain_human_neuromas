######### Analysis of spatial transcriptomics with Giotto ########

## Useful tutorials:
# One section: https://giottosuite.readthedocs.io/en/latest/subsections/datasets/mouse_visium_brain.html
# Integrated sections: https://giottosuite.readthedocs.io/en/latest/subsections/datasets/visium_prostate_integration.html

## If installation needed, make sure you install most up-to-date version of Giotto suite
# if(!"devtools" %in% installed.packages()) {
#   install.packages("devtools")
# }
# 
# devtools::install_github("drieslab/Giotto@suite")

## Loading packages
library(Giotto)
library(ggplot2)
library(scran)
library(harmony)
library(cowplot)
library(dplyr)
library(stringr)

########### Basic Giotto workflow ##########

## Start Giotto environment
# create instructions for your python path, how to view your plots and
# parameters to save your plot if wanted
my_python_path = "/users/mdq19mm/.local/share/r-miniconda/envs/giotto_env/bin/python" # set to NULL to use previously installed giotto environment
results_folder = '/mnt/parscratch/users/mdq19mm/giotto/results/integration_test'
instrs = createGiottoInstructions(python_path = my_python_path,
                                  show_plot = FALSE,
                                  return_plot = TRUE,
                                  save_plot = TRUE,
                                  save_dir = results_folder,
                                  plot_format = 'pdf',
                                  dpi = 200,
                                  height = 9,
                                  width = 9)

## Create Giotto objects with expression data, location data and instructions
# List of paths 
filedir <- "/mnt/parscratch/users/mdq19mm/visium_data" # directory with subdirectory for each sample
sample.ids<- list.dirs(filedir,full.names = FALSE, recursive = FALSE) # get list of sample ids from directory names
sample.ids<- sample.ids[-c(1,4,5,9,22)] # removing LN1_C1
sample.paths<- paste(filedir,sample.ids,"outs", sep="/") # get path to each sample's data

samples.images <- c("LN2_C1","LN7_D1","LN8_D1","LN1_A1","LN12_B1","LN15_D1") # Representative sections from each sample to display

# make list of giotto object for each sample from outs directory after running Spaceranger
giotto_list<- list()
for (i in 1:length(sample.ids)) {
  giotto_list[i] <- createGiottoVisiumObject(visium_dir=sample.paths[i], 
                                             png_name='tissue_hires_image.png', 
                                             expr_data= 'filter', 
                                             gene_column_index = 2, 
                                             instructions = instrs)
}
names(giotto_list) <- sample.ids

## show associated images and metadata with giotto object
showGiottoImageNames(giotto_list$LN1_A1) # "image" is the default name
pDataDT(giotto_list$LN1_A1) ## check metadata
showGiottoImageNames(giotto_list$LN2_C1) # "image" is the default name
pDataDT(giotto_list$LN2_C1) ## check metadata

## Join giotto objects
testcombo = joinGiottoObjects(gobject_list = giotto_list,
                              gobject_names = sample.ids,
                              join_method = 'shift', x_padding = 1000)

# join info is stored in this slot
testcombo@join_info

# check joined Giotto object
fDataDT(testcombo)
pDataDT(testcombo)
showGiottoImageNames(testcombo)
showGiottoSpatLocs(testcombo)
showGiottoExpression(testcombo)

## Process object
# subset on in-tissue spots
metadata = pDataDT(testcombo)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
testcombo = subsetGiotto(testcombo, cell_ids = in_tissue_barcodes)

## filter
testcombo <- filterGiotto(gobject = testcombo,
                          expression_threshold = 1,
                          feat_det_in_min_cells = 2,
                          min_det_feats_per_cell = 100,
                          expression_values = c('raw'),
                          verbose = T)
## normalize
testcombo <- normalizeGiotto(gobject = testcombo, scalefactor = 6000, verbose=T)

## add gene & cell statistics
testcombo <- addStatistics(gobject = testcombo)

# to view statistics
head(fDataDT(testcombo))
head(pDataDT(testcombo))

## visualize
spatPlot2D(gobject = testcombo, show_image = T, point_alpha = 0.7,
           cell_color = 'nr_feats', color_as_factor = F, save_param = list(save_name = '2_b_spatplot_image'))

## Calculate highly variable features
testcombo <- calculateHVF(gobject = testcombo, save_plot = TRUE)
gene_metadata = fDataDT(testcombo)
# Use HVF expressed in > 3% of spots with a mean expression value in spots where feat is detected > 0.4
featgenes = gene_metadata[hvf == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$feat_ID

## run PCA on expression values (default)
testcombo <- runPCA(gobject = testcombo,
                 feats_to_use = featgenes)

## Visualise PCA
screePlot(testcombo, ncp = 30)
dimPlot2D(gobject = testcombo,dim_reduction_to_use = "pca", save_param = list(save_name = "3a_screeplot"))

########## Harmony integration ##########
testcombo = runGiottoHarmony(testcombo, vars_use = 'list_ID', do_pca = F)

## UMAP dimension reduction
testcombo = runUMAP(testcombo, dim_reduction_name = 'harmony', dim_reduction_to_use = 'harmony', name = 'umap_harmony')

## sNN network (default)
testcombo <- createNearestNetwork(gobject = testcombo,
                                  dim_reduction_to_use = 'harmony', dim_reduction_name = 'harmony', name = 'NN.harmony',
                                  dimensions_to_use = 1:10, k = 15)

## Leiden clustering
testcombo <- doLeidenCluster(gobject = testcombo,
                             network_name = 'NN.harmony', resolution = 0.4, n_iterations = 1000, name = 'leiden_harmony')

## Get images
dir.create(paste0(results_folder,'/plots'))
plots <- paste0(results_folder,'/plots')

cell_colours <- c("1"= '#CCB1F1',
                  "2"= '#ff9a36',
                  "3"= '#28CECA',
                  "4"= '#D4D915',
                  "5"= '#B95FBB',
                  "6"= '#25aff5',
                  "7"= '#E6C122',
                  "8"= '#A4DFF2',
                  "9"= '#AC8F14',
                  "10"= '#03FCA5',
                  "11"= 'chocolate4',
                  "12"= '#2236E1',
                  "13"= '#B01E3A',
                  "14"= '#31C53F',
                  "15"= 'lightpink',
                  "16"= '#aeadb3',
                  "17"= 'ivory')

umap<- plotUMAP(gobject = testcombo,
                cell_color = 'leiden_harmony', dim_reduction_name = "umap_harmony",point_size = 1, cell_color_code= cell_colours,
                point_shape='no_border', save_param = list(save_name = 'umap'))

spatial <- spatPlot2D(testcombo, cell_color = 'leiden_harmony', cell_color_code = cell_colours, image_name = samples.images, group_by = "list_ID",
                      point_shape='no_border', point_size=1.5, coord_fix_ratio = 0.75, save_param = list(save_name = 'spatial'))

ggsave("umap_harmony.pdf", plot=umap,path=plots,width=10, height=10)
ggsave("spatial_harmony.pdf", plot=spatial,path=plots,width=35, height=28,limitsize = F)

umap<-umap + theme_void() +theme(legend.position="none") +ggtitle('')
spatial<-spatial + theme_void()+ggtitle('')
plot<-cowplot::plot_grid(plotlist=list(umap,spatial),nrow = 1,rel_widths = c(1,6))
ggsave("umap_spatial_harmony.pdf", plot=plot,path=plots,width=35, height=5)

## Find markers with Scran
scran_markers_subclusters = findMarkers_one_vs_all(gobject = testcombo,
                                                   method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_harmony')
topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats

write.csv(scran_markers_subclusters, file = paste0(results_folder,"/scran_markers_clusters.csv"),row.names = F)

####### PAGE enrichment ########
### load giotto object with sc data from 06_sc_preprocessing_giotto.R
neuromas = loadGiotto(path_to_folder = "/mnt/parscratch/users/mdq19mm/giotto/combined_runs/objects/nerves_sc")

# Create PAGE matrix
# PAGE matrix should be a binary matrix with each row represent a gene marker and each column represent a cell type
# markers_scran is generated from single cell analysis ()
markers_scran = findMarkers_one_vs_all(gobject=neuromas,
                                       method="scran",
                                       expression_values="normalized",
                                       cluster_column='annotation',
                                       min_feats=3)

top_markers <- markers_scran[, head(.SD, 10), by="cluster"]
celltypes<-levels(factor(markers_scran$cluster))
sign_list<-list()

for (i in 1:length(celltypes)){
  sign_list[[i]]<-top_markers[which(top_markers$cluster == celltypes[i]),]$feats
}

# A count matrix and a vector for all cell labels will be needed
sc_expression_norm = getExpression(neuromas,
                                   values = "normalized",
                                   output = "matrix")

annotation_feats = pDataDT(neuromas)$annotation

######## PAGE Enrichment #######
PAGE_matrix = makeSignMatrixPAGE(sign_names = celltypes,
                                 sign_list = sign_list)

testcombo = runPAGEEnrich(gobject = testcombo,
                       sign_matrix = PAGE_matrix,
                       min_overlap_genes = 2)

cell_types_subset = colnames(PAGE_matrix)

# Plot PAGE enrichment result
testcombo@instructions$save_dir <- "/mnt/parscratch/users/mdq19mm/giotto/combined_runs/PAGE_enrichment/"
samples.images <- c("LN2_C1","LN7_D1","LN8_D1","LN1_A1","LN12_B1","LN15_D1") # Representative sections from each sample to display

#subset object
spatial_locs<-getSpatialLocations(testcombo,output="data.table") 

for (i in 1:length(samples.images)){
  spat_locs <- spatial_locs %>% filter(str_starts(cell_ID,samples.images[i]))
  subset <- subsetGiottoLocs(testcombo, 
                             x_max = max(spat_locs$sdimx), x_min = min(spat_locs$sdimx),
                             y_max = max(spat_locs$sdimy), y_min = min(spat_locs$sdimy))
  spatCellPlot2D(gobject = subset, 
                 spat_enr_names = 'PAGE',
                 point_shape='no_border',point_size=1, 
                 cell_annotation_values = cell_types_subset[1:17], cow_n_col = 2, coord_fix_ratio = 1,
                 save_param = list(save_name = paste0("PAGE_plot_",samples.images[i]), base_width=10, base_height=45, limitsize=FALSE))
  
}

# Reload object
testcombo = loadGiotto(path_to_folder = "/mnt/parscratch/users/mdq19mm/giotto/results/objects/testcombo",reconnect_giottoImage = TRUE)

# Clusters were annotated as follows from 1 to 17
cluster_ann<- c("Fibro",
                "Endo",
                "SC1",
                "Myo1",
                "Peri",
                "SC2",
                "Myo2",
                "SC3",
                "Myo3",
                "SC4",
                "Myo4",
                "SC5",
                "Bcells",
                "SC6",
                "Macro",
                "16",
                "17")

