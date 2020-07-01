library(irlba)
library(Rtsne)
library(Matrix)
library(igraph)
library(scater)
library(ggplot2)
library(biomaRt)
library(matrixStats)

setParallel(ncores = 1)

source("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/core_scripts/core_functions.R")

load_data2020()

path2data  <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/"

#Doublets
##Computing doublet scores

sub_sces <- lapply(unique(meta$sample), function(x) return(scater::normalize(sce[, meta$sample == x])))
hvg.list <- lapply(sub_sces, getHVGs)
names(hvg.list) <- unique(meta$sample)

set.seed(42)
scores_hvgs <- lapply(1:length(sub_sces), function(i) doubletCells(sub_sces[[i]],
                 approximate = TRUE, subset.row = rownames(sub_sces[[i]]) %in% hvg.list[[i]]))

scores_hvgs <- do.call(c, scores_hvgs)
scores      <- scores_hvgs

#saveRDS(scores,readRDS(paste0(path2data, "doublet_scores.rds"))
#scores <- readRDS(readRDS(paste0(path2data, "doublet_scores.rds"))

## Scoring of cells

sampsize <- melt(table(meta$sample))
order    <- order(sampsize$value, decreasing = FALSE)

ggplot(data.frame(score = log2(scores+1), sample = meta$sample, stage = meta$stage),
  aes (x = factor(sample, levels = sampsize$Var1[order]), y = score, fill = stage)) +
    geom_boxplot(col = "black") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(y = "log2(doublet score + 1)", x = "Sample") +
    scale_fill_manual(values = stage_colours_extension, labels = stage_labels, name = "Stage")


ggsave(readRDS(paste0(path2data, "QC/doubletScores.pdf"))

tsnes  <- lapply(unique(meta$sample), function(x){
  sub  <- scater::normalize(sce[, meta$sample == x])
  hvgs <- getHVGs(sub)
  pca  <- prcomp_irlba(t(logcounts(sub[hvgs,])), n = 50)
  tsne <- Rtsne(pca$x, pca = FALSE)
  return(tsne$Y)
})
names(tsnes) <- unique(meta$sample)
colmax       <- max(log2(scores + 1))

score_plots <- lapply(unique(meta$sample), function(x){
  
  scr <- log2(scores[meta$sample == x] + 1)
  ord <- order(scr)
  p <- ggplot(as.data.frame(tsnes[[as.character(x)]])[ord,],
         aes(x = V1, y= V2, col = scr[ord])) +
    geom_point(size = 0.6) +
    scale_color_gradient2(name = "log2(score+1)", mid = "cornflowerblue", low = "gray75", high = "black", midpoint = max(colmax)/2) +
    # scale_color_viridis() +
    ggtitle(paste0("Sample ", x)) +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(colour = "black", linetype = 1, size = 0.5))
  
  return(p)
})

ggsave(readRDS(paste0(path2data, "QC/doubletScores_tsne_samplewise.pdf"))

##Clustering doublets

clusters <- lapply(unique(meta$sample), function(x){
  sub_sce    <- scater::normalize(sce[,meta$sample == x])
  sub_meta   <- meta[meta$sample == x,]
  graph      <- buildSNNGraph(sub_sce, pc.approx = TRUE)
  clusters   <- cluster_louvain(graph)
  vec        <- as.numeric(membership(clusters))
  names(vec) <- sub_meta$cell
  return(vec)
})
names(clusters) <- unique(meta$sample)

clusters_sub <- lapply(unique(meta$sample), function(x){
  sub_sce  <- scater::normalize(sce[,meta$sample == x])
  sub_meta <- meta[meta$sample == x,]
  clusts   <- clusters[[as.character(x)]]
  sub_clusts <- lapply(unique(clusts), function(y){
    sub_sub_sce  <- scater::normalize(sub_sce[,clusts == y])
    sub_sub_meta <- sub_meta[clusts == y,]
    graph        <- buildSNNGraph(sub_sub_sce, pc.approx = TRUE, d = min(c(ncol(sub_sub_sce)-1, 50)))
    clusters     <- cluster_louvain(graph)
    vec <- as.numeric(membership(clusters))
    vec <- paste0(y, ".", vec)
    names(vec) <- sub_sub_meta$cell
    return(vec)
  })
  clusts <- do.call(c, sub_clusts)
  clusts <- clusts[match(meta$cell[meta$sample == x], names(clusts))]
  return(clusts)
})
