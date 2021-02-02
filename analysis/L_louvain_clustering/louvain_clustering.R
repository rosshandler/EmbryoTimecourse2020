library(Matrix)
library(igraph)

library(scran)
#set it up for scran to be properly parallelised
library(BiocParallel)
ncores  = 16
mcparam = SnowParam(workers = ncores)
register(mcparam)

path2integ <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/Integrated_Atlas/"

corrected <- readRDS(paste0(path2integ, "integrated_corrected_pcs.rds"))
meta      <- readRDS(paste0(path2integ, "integrated_premeta_e85wot_somites.rds"))

graph <- buildSNNGraph(corrected, BPPARAM = mcparam, d = NA, transposed = TRUE)

set.seed(42)
clustered    <- cluster_louvain(graph)
meta$cluster <- as.vector(membership(clustered))
clust.sizes  <- sapply(1:length(unique(meta$cluster)), function(x) sum(meta$cluster == x))

sub_clusters <- lapply(1:length(unique(meta$cluster)), function(x){
  #possibly this needs recalculation of HVGs
  pc_sub <- corrected[meta$cluster == x,]
  
  graph  <- buildSNNGraph(pc_sub, d = NA, BPPARAM = mcparam, transposed = TRUE)
  set.seed(42)
  clusts <- cluster_louvain(graph)
  return(as.numeric(membership(clusts)))
})

# Save to metadata
meta$cluster.sub <- NA
for( i in 1:length(unique(meta$cluster))){
  meta$cluster.sub[meta$cluster == i] = sub_clusters[[i]]
}

saveRDS(meta, paste0(path2integ, "integrated_meta.rds"))
