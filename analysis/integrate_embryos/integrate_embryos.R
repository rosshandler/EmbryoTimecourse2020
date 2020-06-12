library(irlba)
library(scran)
library(Matrix)
library(batchelor)
library(SingleCellExperiment)

########################
########################
####
#### Integrate embryos and generate batch corrected PCA
####
########################

path2atlas_ext <-  "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/"
path2atlas     <- "/hps/research1/marioni/ivan/EmbryoTimeCourse/atlas/data/"
path2integ     <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/Integrated_Atlas/"
cluster <- "ebi"

path2atlas_ext <- "/mnt/b/iimaz/embryos/"
path2atlas     <- "/mnt/b/iimaz/embryos/"
path2integ     <- "/mnt/b/iimaz/embryos/"
cluster <- "ema"

source("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/core_scripts/core_functions.R")

### Integrate sce and generate several files (this is very intense memory wise)
sceIntegrateEmbryos(path2sce_atlas = path2atlas, path2sce_extension = path2atlas_ext, path2out = path2integ)

### Read Original Atlas
load(file = paste0(path2atlas, "embryo_data.RData"))

sce_atlas  <- sce
meta_atlas <- meta
logcounts(sce_atlas) <- Matrix(logcounts(sce_atlas), sparse=TRUE) ## change permantely original atlas logcounts to sparse? 

load(file = paste0(path2atlas_ext, "embryo_extension_data_hvgs_pca.RData"))
  
sce_ext  <- sce
meta_ext <- meta

rm(sce);rm(meta)

integrated_corrected_pcs <- embryoBatchCorrection(
  atlas_sce      = sce_atlas,
  atlas_meta     = meta_atlas,
  extension_sce  = sce_ext,
  extension_meta = meta_ext,
  normalised = TRUE, computer = cluster)

saveRDS(integrated_corrected_pcs,
  file = paste0(path2integ, "integrated_corrected_pcs.rds"))


integrated_corrected_pcs <- embryoBatchCorrection(
  atlas_sce      = sce_atlas,
  atlas_meta     = meta_atlas,
  extension_sce  = sce_ext,
  extension_meta = meta_ext,
  normalised = TRUE, computer = cluster,
  integrated_pca_overwrite = TRUE, mnn_order = "reverse")

saveRDS(integrated_corrected_pcs,
  file = paste0(path2integ, "integrated_reverse_corrected_pcs.rds"))

