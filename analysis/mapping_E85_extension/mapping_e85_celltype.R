library(irlba)
library(scran)
library(scater)
library(Matrix)
library(batchelor)
library(SingleCellExperiment)

path2atlas      <- "/hps/research1/marioni/ivan/EmbryoTimeCourse/atlas/data/"
path2atlas_ext  <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/"
path2integ      <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/Integrated_Atlas/"

## Read Original Atlas
load(file = paste0(path2atlas, "embryo_data.RData"))

sce_atlas  <- sce
meta_atlas <- readRDS(file = paste0(path2atlas, "metadata_mixgast_extended.rds"))
logcounts(sce_atlas) <- Matrix(logcounts(sce_atlas), sparse=TRUE) ## change permantely original atlas logcounts to sparse? 

load(file = paste0(path2atlas_ext, "embryo_extension_data_hvgs_pca.RData"))
  
sce_ext  <- sce
meta_ext <- meta

rm(sce);rm(meta)

res <- mapWrap_extended(
  atlas_sce   = sce_atlas,
  atlas_meta  = meta_atlas,
  map_sce     = sce_ext,
  map_meta    = meta_ext,
  mapstage_x  = "E8.5", 
  map2stage_x = "E8.5", 
  normalised  = TRUE,
  integrated_pca_overwrite = TRUE,
  mixed_gastrulation_map   = TRUE)

saveRDS(res, file = paste0(path2integ, "maping_withinE85.rds"))


## Update metadata

meta <- readRDS(paste0(path2integ, "integrated_premeta_v1.rds"))

e85index <- match(gsub("map_", "ext_", res$mapping$cell), meta$cell)

celltype.extended <-  meta$celltype
celltype.extended[e85index] <- res$mapping$celltype.mapped
meta <- cbind(meta, celltype.extended)

names_order <- c("cell", "sample", "stage", "stage.mapped", "stage.collapsed", "stage.mapped.collapsed",
"somite_count", "tube_name", "celltype", "celltype.extended", "doub.density")

meta <- meta[, names_order]

saveRDS(meta, paste0(path2integ, "integrated_premeta_e85mapping.rds"))
