# Investigation of the cellular and molecular changes linked with neuropathic pain in healthy and injured human trigeminal nerves

## Authors
Martina Morchio<sup>1</sup>, Ishwarya Sankaranarayanan<sup>2</sup>, Diana Tavares-Ferreira<sup>2</sup>, Natalie Wong<sup>1</sup>, Simon Atkins<sup>1</sup>, Emanuele Sher<sup>3</sup>, Theodore J Price<sup>2</sup>, Daniel W Lambert<sup>1</sup>, Fiona M Boissonade<sup>1</sup>*


## Affiliations
<sup>1</sup>Neuroscience Institute and School of Clinical Dentistry, University of Sheffield, Sheffield, UK

<sup>2</sup>Department of Neuroscience and Center for Advanced Pain Studies, University of Texas at Dallas, Richardson, TX, USA.

<sup>3</sup>Lilly Research Centre, Eli Lilly and Company, Surrey, UK

Email: f.boissonade@sheffield.ac.uk



## Abstract
Injuries to the trigeminal nerve, responsible for sensory innervation to the face, may occur during routine dental procedures, resulting in the formation of a neuroma accompanied by loss of sensation and/or symptoms of pain. In order to gain insight into the molecular mechanisms underpinning the sensory changes, single nuclei RNA sequencing and spatial transcriptomics were employed to profile the transcriptional landscape at single cell resolution of human trigeminal nerves and neuromas. Cellular and transcriptional changes were identified that correlated with the presence of pain, including an expansion of endothelial cells with a pro-inflammatory phenotype and over-expression of HLA-A, CXCL2 and CXCL8. Interactome analysis highlighted signalling changes linked with the presence of pain. HLA-A protein expression was confirmed in neuromas and positively correlated with symptoms of pain. The atlas generated represents a valuable resource for pain research, highlighting the role of inflammation, endothelial cell dysfunction and chemokine signalling in neuropathic pain.

This repository contains qanalysis code for the study above.

## What's included in this directory
- [01_cellbender.sh](https://github.com/martina-boop/pain_human_neuromas/blob/main/scripts/01_cellbender.sh): Ambient RNA removal from snRNAseq data
- [02_doubletfinder_and_filtering.R](https://github.com/martina-boop/pain_human_neuromas/blob/main/scripts/02_doubletfinder_and_filtering.R): Filtering of low quality nuclei and doublets
- [03_snRNAseq_integration_rpca.R](https://github.com/martina-boop/pain_human_neuromas/blob/main/scripts/03_snRNAseq_integration_rpca.R): Workflow for integration of neuromas and trigeminal nerve roots nuclei using the rpca method
- [04_finding_markers.R](https://github.com/martina-boop/pain_human_neuromas/blob/main/scripts/04_finding_markers.R): Calculation of differentially expressed genes in each cluster and validation of the expression of known marker genes in the expected cell types
- [05_images_snRNAseq.R](https://github.com/martina-boop/pain_human_neuromas/blob/main/scripts/05_images_snRNAseq.R): snRNAseq data visualisation (figure 1)
- [06_sc_preprocessing_giotto.R](https://github.com/martina-boop/pain_human_neuromas/blob/main/scripts/06_sc_preprocessing_giotto.R): Code to process snRNAseq data for cell-type deconvolution with PAGE in the Giotto workflow
- [07_visium_sections_integration.R](https://github.com/martina-boop/pain_human_neuromas/blob/main/scripts/07_visium_sections_integration.R): Integration of Visium sections using the standard Giotto workflow (figure 3)
- [08_sce_object_from_giotto.R](https://github.com/martina-boop/pain_human_neuromas/blob/main/scripts/08_sce_object_from_giotto.R): Code to obtain a SingleCellExperiment object from a Giotto object for downstream data analysis
- [09_making_heatmap_spatial.R](https://github.com/martina-boop/pain_human_neuromas/blob/main/scripts/09_making_heatmap_spatial.R): Heatmap of spatial data (figure 3)
- [10_differential_abundance.R](https://github.com/martina-boop/pain_human_neuromas/blob/main/scripts/10_differential_abundance.R): Differential abundance analysis of clusters identified with spatial transcriptomics using EdgeR in painful and non-painful samples (figure 3)
- [11_pseudobulk_preprocessing.R](https://github.com/martina-boop/pain_human_neuromas/blob/main/scripts/11_pseudobulk_preprocessing.R): Pre-processing of spatial transcriptomics data for pseudo-bulk differential expression analysis
- [12_deseq2_analysis_spatial.R](https://github.com/martina-boop/pain_human_neuromas/blob/main/scripts/12_deseq2_analysis_spatial.R): Differential gene expression analysis of barcodes in nerve fascicles from painful and non-painful samples with DEseq2 (figure 4)
- [13_cellchat_objects_computation.R](https://github.com/martina-boop/pain_human_neuromas/blob/main/scripts/13_cellchat_objects_computation.R): Computation of cellchat objects and inference of communication network in painful and non-painful samples with cellchat v2
- [14_cellchat_images.R](https://github.com/martina-boop/pain_human_neuromas/blob/main/scripts/14_cellchat_images.R): Visualisation of cell-cell communication with cellchat v2 (figure 5)

## Data availability
Data will be available upon article acceptance.

## Related Packages
- [Seurat](https://satijalab.org/seurat/)
- [Giotto](https://github.com/drieslab/Giotto)
- [DEseq2](https://github.com/thelovelab/DESeq2)
- [Differential cluster abundance analysis](https://bioconductor.org/books/3.13/OSCA.multisample/differential-abundance.html)
- [EdgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
- [CellChat v2](https://github.com/jinworks/CellChat)

