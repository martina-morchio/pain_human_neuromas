########### Making SingleCellExperiment object from Giotto ############
library(Giotto)
library(SingleCellExperiment)
library(dplyr)

# obtaining metadata and count matrix from integrated Visium slides to be used for further analysis (heatmap, cellchat)
testcombo = loadGiotto(path_to_folder = "/mnt/parscratch/users/mdq19mm/giotto/results/objects/testcombo",reconnect_giottoImage = TRUE)

fDataDT(testcombo)
pDataDT(testcombo) %>% class()

counts_raw<- getExpression(testcombo, values="raw", output="matrix") 
counts_normalised<- getExpression(testcombo, values="normalized", output="matrix") 
counts_scaled<- getExpression(testcombo, values="scaled", output="matrix") 
meta<- pDataDT(testcombo)
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
names(cluster_ann)<- c(1:17)
meta <- meta %>% mutate(annotation=recode(leiden_harmony, !!!cluster_ann))

#Adding patient metadata
samples <- c("LN1","LN1","LN1","LN1","LN12","LN12","LN13","LN13","LN15","LN15","LN15","LN15",
             "LN2","LN2","LN2","LN2","LN2","LN7","LN7","LN7","LN7","LN7","LN8","LN8","LN8","LN8","LN8")
status<- c("P","P","P","P","P","P","P","P","P","P","P","P",
           "NP","NP","NP","NP","NP","NP","NP","NP","NP","NP","NP","NP","NP","NP","NP")
reps<- meta$list_ID %>% unique()
info<- data.frame(list_ID=reps, sample=samples,status=status)
meta<- left_join(meta,info, by="list_ID")

# save metadata table
write.csv(meta, "/mnt/parscratch/users/mdq19mm/giotto/results/metadata_giotto_harmony_integrated.csv")

# creating sce object with normalised data
sce_norm <- SingleCellExperiment(list(counts=counts_normalised),
                            colData=meta)

saveRDS(sce_norm, "/mnt/parscratch/users/mdq19mm/giotto/results/giotto_integrated_visium_sce_norm.rds")
