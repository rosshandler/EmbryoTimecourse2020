library(Matrix)
library(scran)
library(scater)
library(igraph)
library(BiocParallel)

parallel = FALSE
#set it up for scran to be properly parallelised
if (parallel == TRUE){
  library(BiocParallel)
  ncores = 10
  mcparam = SnowParam(workers = ncores)
  register(mcparam)
}else{
  library(BiocParallel)
  ncores = 1
  mcparam = SnowParam(workers = ncores)
  register(mcparam)
}

path2data  <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/"

load_data2020 = function(normalise = TRUE, remove_doublets = FALSE, remove_stripped = FALSE, load_corrected = FALSE){
  
  if(load_corrected & (!remove_doublets | !remove_stripped)){
    message("Using corrected PCs, also removing doublets + stripped now.")
    remove_doublets = TRUE
    remove_stripped = TRUE
  }
  
  require(scran)
  require(scater)
  require(SingleCellExperiment)
  require(Matrix)
  
  path2data  <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/"
  counts     <- readRDS(paste0(path2data,"raw_counts_qc.rds"))
  genes      <- read.table(paste0(path2data, "genes.tsv"), stringsAsFactors = F)
  meta       <- read.table(paste0(path2data, "meta.tab"), header = TRUE, sep = "\t",
    stringsAsFactors = FALSE, comment.char = "$")
  
  rownames(counts) = genes[,1] #ensembl
  colnames(counts) = meta$cell

  sce = SingleCellExperiment(assays = list("counts" = counts))
  
  if(normalise){
    sfs = read.table(paste0(path2data, "sizefactors.tab"), stringsAsFactors = F)[,1]
    sizeFactors(sce) = sfs
    sce = scater::normalize(sce)
  }
  
  if(remove_doublets){
    sce = scater::normalize(sce[,!meta$doublet])
    meta = meta[!meta$doublet,]
  }
  
  if(remove_stripped){
    sce = scater::normalize(sce[,!meta$stripped])
    meta = meta[!meta$stripped, ]
  }
  
  if(load_corrected){
    corrected = readRDS(paste0(path2data, "corrected_pcas.rds"))
    assign("corrected", corrected, envir = .GlobalEnv)  
  }
    
  
  assign("genes", genes, envir = .GlobalEnv)
  assign("meta", meta, envir = .GlobalEnv)
  assign("sce", sce, envir = .GlobalEnv)
  
  
  invisible(0)
}

load_data2020(normalise = FALSE)

lib.sizes <- Matrix::colSums(counts(sce))
sce       <- sce[calcAverage(sce)>0.1,]

clusts = as.numeric(quickCluster(sce, method = "igraph", min.size = 100, BPPARAM = mcparam))

path2data  <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/"
#saveRDS(clusts, file = paste0(path2data, "quickCluster.rds"))
#clusts <- readRDS(paste0(path2data, "quickCluster.rds"))

min.clust <- min(table(clusts))/2
new_sizes <- c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce <- computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

ggplot(data = data.frame(X = lib.sizes, Y = sizeFactors(sce)),
              mapping = aes(x = X, y = Y)) +
  geom_point() +
  scale_x_log10(breaks = c(5000, 10000, 50000, 100000), labels = c("5,000", "10,000", "50,000", "100,000") ) +
  scale_y_log10(breaks = c(0.2, 1, 5)) +
  labs(x = "Number of UMIs", y = "Size Factor")

ggsave(paste0(path2data, "QC/sizefactors.pdf")

write.table(sizeFactors(sce), quote = F, col.names = F, row.names = F,
  file = paste0(path2data, "sizefactors.tab"))
