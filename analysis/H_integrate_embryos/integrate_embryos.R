## Integrate embryos and generate batch corrected PCA

library(irlba)
library(scran)
library(Matrix)
library(batchelor)
library(SingleCellExperiment)

path2atlas     <- "/hps/research1/marioni/ivan/EmbryoTimeCourse/atlas/data/"
path2atlas_ext <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/"
path2integ     <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/Integrated_Atlas/"
cluster <- "ebi"

source("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/core_scripts/core_functions.R")

## Integrate sce and generate several files (this is memory intense!)
sceIntegrateEmbryosStep1(
  path2sce_atlas = paste0(path2atlas, "embryo_data.RData"), 
  path2sce_extension = paste0(path2atlas_ext, "embryo_extension_data_hvgs_pca.RData"),
  path2out = path2integ
)

sceIntegrateEmbryosStep2(path2working_dir = path2integ)

sceIntegrateEmbryosStep3counts(path2working_dir = path2integ)

sceIntegrateEmbryosStep3logcounts(path2working_dir = path2integ)

sceIntegrateEmbryosStep4(path2working_dir = path2integ)

## Batch corrected PCA
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
  normalised = TRUE, nPCs = 75, computer = cluster)

saveRDS(integrated_corrected_pcs,
  file = paste0(path2integ, "integrated_corrected_pcs.rds"))


integrated_corrected_pcs <- embryoBatchCorrection(
  atlas_sce      = sce_atlas,
  atlas_meta     = meta_atlas,
  extension_sce  = sce_ext,
  extension_meta = meta_ext,
  normalised = TRUE, nPCs = 75, computer = cluster,
  integrated_pca_overwrite = TRUE, mnn_order = "reverse")

saveRDS(integrated_corrected_pcs,
  file = paste0(path2integ, "integrated_reverse_corrected_pcs.rds"))


## Load original atlas metadata extended by mixed gastrulation stage mapping (QC passed)
meta_atlas <- readRDS("metadata_mixgast_extended.rds")

## Integrate metadata
integrated_meta <- embryoPreMetaData(atlas_meta = meta_atlas, extension_meta = meta_ext)
saveRDS(integrated_meta, file = "integrated_premeta_v0.rds")

integrated_meta_e95  <- integrated_meta[integrated_meta$stage == "E9.5",]
integrated_meta_e925 <- integrated_meta[integrated_meta$stage == "E9.25",]

## Collapsing E9.5 and E9.25
integrated_meta$stage        <- gsub("E9.5", "E9.25", integrated_meta$stage)
integrated_meta$stage.mapped <- gsub("E9.5", "E9.25", integrated_meta$stage.mapped)

saveRDS(integrated_meta, file="integrated_premeta_stage_collapsed.rds")

meta1 <- readRDS("integrated_premeta_stage_collapsed.rds")
meta2 <- readRDS("integrated_premeta_v0.rds")

meta2$stage.collapsed        <- meta1$stage
meta2$stage.mapped.collapsed <- meta1$stage.mapped

names_order <- c("cell", "sample", "stage", "stage.mapped", "stage.collapsed", "stage.mapped.collapsed",
  "somite_count", "tube_name", "celltype", "doub.density")

saveRDS(meta2[,names_order], file="integrated_premeta_v1.rds")


