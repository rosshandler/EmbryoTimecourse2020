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

# Biomart reference
mouse_ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

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



sceIntegrateEmbryosStep1 <- function(path2sce_atlas, path2sce_extension, path2out){

  load(path2atlas)
  atlas_sce  <- sce
  atlas_meta <- meta

  load(path2extension)
  ext_sce    <- sce
  ext_meta   <- meta
  
  #mouse_ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
 
  shared_genes  <- intersect(rownames(atlas_sce), rownames(ext_sce))

  atlas_sce  <- atlas_sce[shared_genes, ]
  ext_sce    <- ext_sce[shared_genes, ]

  #prevent duplicate cell names
  colnames(ext_sce) <- paste0("ext_", colnames(ext_sce))
  
  writeLines(shared_genes, con = paste0(path2out, "shared_genes.txt"))
  
  gene_map <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol"),
    filters = "ensembl_gene_id", values = shared_genes, mart = mouse_ensembl)
  
  write.table(gene_map, file = paste0(path2out,  "shared_genesyms.txt"), row.names = FALSE)

  counts    <- cbind(counts(atlas_sce), counts(ext_sce))
  logcounts <- cbind(Matrix(logcounts(atlas_sce), sparse=TRUE), logcounts(ext_sce))
 
  saveRDS(counts, paste0(path2out, "integrated_counts_sparse.rds"))
  saveRDS(logcounts, paste0(path2out, "integrated_logcounts_sparse.rds"))

}
