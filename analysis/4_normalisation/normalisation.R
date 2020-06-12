library(Matrix)
library(scran)
library(scater)
library(igraph)
library(BiocParallel)

setParallel(ncores = 1)

source("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/core_scripts/core_functions.R")

load_data2020(normalise = FALSE)

path2data  <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/"

lib.sizes <- Matrix::colSums(counts(sce))
sce       <- sce[calcAverage(sce)>0.1,]

clusts <- as.numeric(quickCluster(sce, method = "igraph", min.size = 100, BPPARAM = mcparam))

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
