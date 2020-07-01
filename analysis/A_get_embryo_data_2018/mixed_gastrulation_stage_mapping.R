library(irlba)
library(scran)
library(scater)
library(Matrix)
library(batchelor)
library(BiocNeighbors)
library(SingleCellExperiment)

path2atlas  <- "/hps/research1/marioni/ivan/EmbryoTimeCourse/atlas/data/"

source(paste0(path2code, "embryo_marionilab_extended.r"))

load(paste0(path2atlas, "embryo_data.RData"))

dim(sce)
#[1]  29452 116312

sce_mg  <- sce[, meta$stage == "mixed_gastrulation"]
meta_mg <- meta[meta$stage  == "mixed_gastrulation", ]

sce_atlas  <- sce[, (!meta$stage == "mixed_gastrulation")]
meta_atlas <- meta[!(meta$stage  == "mixed_gastrulation"), ]

dim(sce_mg)
#[1] 29452  7455

## Mapping mixed_gastrulation cells to embryo stages in the atlas time course

res <- mapWrap_extended(atlas_sce = sce_atlas, atlas_meta = meta_atlas,
  map_sce = sce_mg, map_meta = meta_mg, k = 30)

saveRDS(res, file = paste0(path2atlas, "mapping_mixgast.rds"))

meta_extended <- meta
stage.mapped  <- meta$stage
stage.mapped[meta$stage == "mixed_gastrulation"] <- as.character(res$mapping$stage.mapped)
meta_extended <- cbind(meta[,c("cell", "barcode", "sample", "stage")], stage.mapped,
  meta[,5:ncol(meta)])

saveRDS(meta_extended, file = paste0(path2atlas, "metadata_mixgast_extended.rds"))
