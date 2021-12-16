# ROUTINELY USED PACKAGES
#library(topGO)
library(irlba)
library(scran)
library(limma)
library(scater)
library(Matrix)
library(biomaRt)
library(cowplot)
library(ggplot2)
library(batchelor)
library(data.table)
library(BiocParallel)
library(SingleCellExperiment) 

# Set number of cores
setParallel <- function(ncores = 1){
  library(BiocParallel)
  mcparam = SnowParam(workers = ncores)
  register(mcparam)
}

# Colour transparency
ac <- function(x, alpha=1, ...) { y <- adjustcolor(x, alpha.f=alpha, ...); names(y) <- names(x); return(y)}

load_data2020 = function(normalise = TRUE, remove_doublets = FALSE, remove_stripped = FALSE, load_corrected = FALSE){
  
  if(load_corrected & (!remove_doublets | !remove_stripped)){
    message("Using corrected PCs, also removing doublets + stripped now.")
    remove_doublets = TRUE
    remove_stripped = TRUE
  }
    
  path2data  <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/"
  counts     <- readRDS(paste0(path2data, "raw_counts_qc.rds"))
  genes      <- readLines(paste0(path2data, "genes.tsv"))
# genes      <- read.table(paste0(path2data, "genes.tsv"), stringsAsFactors = F)
  meta       <- read.table(paste0(path2data, "meta.tab"), header = TRUE, sep = "\t",
    stringsAsFactors = FALSE, comment.char = "$")

  rownames(counts) <- genes #ensembl  
# rownames(counts) <- genes[,1] #ensembl
  colnames(counts) <- meta$cell

  sce <- SingleCellExperiment(assays = list("counts" = counts))
  
  if(normalise){
    sfs <- read.table(paste0(path2data, "sizefactors.tab"), stringsAsFactors = F)[,1]
    sizeFactors(sce) <- sfs
    sce <- logNormCounts(sce)
  }
  
  if(remove_doublets){
    sce  <- logNormCounts(sce[,!meta$doublet])
    meta <- meta[!meta$doublet,]
  }
  
  if(remove_stripped){
    sce  <- logNormCounts((sce[,!meta$stripped])
    meta <- meta[!meta$stripped, ]
  }
  
  if(load_corrected){
    corrected <- readRDS(paste0(path2data, "corrected_pcas.rds"))
    assign("corrected", corrected, envir = .GlobalEnv)  
  }
    
  
  assign("genes", genes, envir = .GlobalEnv)
  assign("meta", meta, envir = .GlobalEnv)
  assign("sce", sce, envir = .GlobalEnv)
  
  
  invisible(0)
}


