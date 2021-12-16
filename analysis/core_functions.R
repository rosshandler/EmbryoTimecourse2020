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

sceIntegrateEmbryosStep2 <- function(path2working_dir){

  gene_map <- read.delim(paste0(path2working_dir, "shared_genesyms.txt"), header = TRUE, sep=" ")

  counts    <- readRDS(paste0(path2working_dir, "integrated_counts_sparse.rds"))
  counts    <- counts[rowSums(counts) != 0,]
  genes     <- intersect(unique(gene_map$ensembl_gene_id), rownames(counts))
  counts    <- counts[genes,]

  logcounts <- readRDS(paste0(path2working_dir, "integrated_logcounts_sparse.rds"))
  logcounts <- logcounts[rowSums(logcounts) != 0,]
  genes     <- intersect(unique(gene_map$ensembl_gene_id), rownames(logcounts))
  logcounts <- logcounts[genes,]

  writeLines(colnames(counts), con = paste0(path2working_dir, "integrated_barcodes.txt"))
  writeLines(genes, con = paste0(path2working_dir, "integrated_genes.txt"))

  sce_sparse <- SingleCellExperiment(assays = list(counts = counts,
    logcounts = logcounts))
  saveRDS(sce_sparse, file = paste0(path2working_dir, "integrated_sce_sparse.rds"))

  writeMM(counts, file = paste0(path2working_dir, "integrated_raw_counts_qc.mtx"))
  writeMM(t(counts), file = paste0(path2working_dir, "integrated_raw_counts_qc_transposed.mtx"))

}

sceIntegrateEmbryosStep3counts <- function(path2working_dir){

  counts    <- counts(readRDS(paste0(path2working_dir, "integrated_sce_sparse.rds")))

  c1 <- as.data.table(as.matrix(counts[,1:50000]), keep.rownames=TRUE)
  c2 <- as.data.table(as.matrix(counts[,50001:100000]), keep.rownames=FALSE)
  c3 <- as.data.table(as.matrix(counts[,100001:150000]), keep.rownames=FALSE)
  c4 <- as.data.table(as.matrix(counts[,150001:200000]), keep.rownames=FALSE)
  c5 <- as.data.table(as.matrix(counts[,200001:250000]), keep.rownames=FALSE)
  c6 <- as.data.table(as.matrix(counts[,250001:300000]), keep.rownames=FALSE)
  c7 <- as.data.table(as.matrix(counts[,300001:350000]), keep.rownames=FALSE)
  c8 <- as.data.table(as.matrix(counts[,350001:400000]), keep.rownames=FALSE)
  c9 <- as.data.table(as.matrix(counts[,400001:ncol(counts)]), keep.rownames=FALSE)
  dt <- cbind(c1,c2,c3,c4,c5,c6,c7,c8,c9)
  colnames(dt)[1] <- "ensembl_id"
  fwrite(dt, file = paste0(path2working_dir, "integrated_counts_datatable.tab"),
    sep = "\t", row.names = FALSE, col.names = TRUE)
  fwrite(t(dt), file = paste0(path2working_dir, "integrated_counts_datatable_transposed.tab"),
    sep = "\t", row.names = FALSE, col.names = TRUE)
}

sceIntegrateEmbryosStep3logcounts <- function(path2working_dir){

  logcounts <- logcounts(readRDS(paste0(path2working_dir, "integrated_sce_sparse.rds")))

  c1 <- as.data.table(as.matrix(logcounts[,1:50000]), keep.rownames=TRUE)
  c2 <- as.data.table(as.matrix(logcounts[,50001:100000]), keep.rownames=FALSE)
  c3 <- as.data.table(as.matrix(logcounts[,100001:150000]), keep.rownames=FALSE)
  c4 <- as.data.table(as.matrix(logcounts[,150001:200000]), keep.rownames=FALSE)
  c5 <- as.data.table(as.matrix(logcounts[,200001:250000]), keep.rownames=FALSE)
  c6 <- as.data.table(as.matrix(logcounts[,250001:300000]), keep.rownames=FALSE)
  c7 <- as.data.table(as.matrix(logcounts[,300001:350000]), keep.rownames=FALSE)
  c8 <- as.data.table(as.matrix(logcounts[,350001:400000]), keep.rownames=FALSE)
  c9 <- as.data.table(as.matrix(logcounts[,400001:ncol(logcounts)]), keep.rownames=FALSE)
  dt <- cbind(c1,c2,c3,c4,c5,c6,c7,c8,c9)
  colnames(dt)[1] <- "ensembl_id"
  fwrite(dt, file = paste0(path2working_dir, "integrated_logcounts_datatable.tab"),
    sep = "\t", row.names = FALSE, col.names = TRUE)

}

sceIntegrateEmbryosStep4 <- function(path2working_dir){

  counts    <- fread(file = paste0(path2working_dir, "integrated_counts_datatable.tab"))
  logcounts <- fread(file = paste0(path2working_dir, "integrated_logcounts_datatable.tab"))

  sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts[,-1]),
   logcounts = as.matrix(logcounts[,-1])))
  rownames(sce) <- counts$ensembl_id

  saveRDS(sce, file =paste0(path2working_dir, "integrated_sce.rds"))
}

sceIntegrateEmbryosADATA <- function(path2working_dir){

  logcounts <- readRDS(paste0(path2working_dir, "integrated_logcounts_sparse.rds"))

  c1 <- as.data.table(as.matrix(logcounts[,1:50000]), keep.rownames=TRUE)
  c2 <- as.data.table(as.matrix(logcounts[,50001:100000]), keep.rownames=FALSE)
  c3 <- as.data.table(as.matrix(logcounts[,100001:150000]), keep.rownames=FALSE)
  c4 <- as.data.table(as.matrix(logcounts[,150001:200000]), keep.rownames=FALSE)
  c5 <- as.data.table(as.matrix(logcounts[,200001:250000]), keep.rownames=FALSE)
  c6 <- as.data.table(as.matrix(logcounts[,250001:300000]), keep.rownames=FALSE)
  c7 <- as.data.table(as.matrix(logcounts[,300001:350000]), keep.rownames=FALSE)
  c8 <- as.data.table(as.matrix(logcounts[,350001:400000]), keep.rownames=FALSE)
  c9 <- as.data.table(as.matrix(logcounts[,400001:ncol(logcounts)]), keep.rownames=FALSE)
  dt <- cbind(c1,c2,c3,c4,c5,c6,c7,c8,c9)
  colnames(dt)[1] <- "ensembl_id"
  fwrite(dt, file = paste0(path2working_dir, "integrated_logcounts_datatable_adata.tab"),
    sep = "\t", row.names = FALSE, col.names = TRUE)

}  

