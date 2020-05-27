
library(irlba)
library(Rtsne)
library(Matrix)
library(igraph)
library(scater)
library(ggplot2)
library(biomaRt)
library(matrixStats)

parallel = FALSE
#set it up for scran to be properly parallelised
if (parallel == TRUE){
  library(BiocParallel)
  ncores  = 10
  mcparam = SnowParam(workers = ncores)
  register(mcparam)
}else{
  library(BiocParallel)
  ncores  = 1
  mcparam = SnowParam(workers = ncores)
  register(mcparam)
}

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
                                                                  approximate = TRUE,
                                                                  subset.row = rownames(sub_sces[[i]]) %in% hvg.list[[i]]))
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


clusters          <- clusters_sub
meta$doub.density <- scores

maxscore <- max(log2(scores+1))

plots <- lapply(unique(meta$sample), function(x){
  p <- ggplot(data = data.frame(clusts = clusters[meta$sample == x],
        scores = log2(scores[meta$sample == x] + 1)),
        mapping = aes(x = factor(clusts), y=scores, fill = factor(clusts))) +
        geom_boxplot() +
        theme(axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = "none") +
        lims(y = c(0, maxscore)) +
        scale_fill_Publication() +
        ggtitle(paste0("Sample, ", x, "(n=", sampsize$value[match(x, sampsize$Var1)], ")"))
  
  return(p)
})

plot_grid(plotlist = plots, ncol = 3)


pdf <- aggregate(meta$doub.density, list(meta$sample, clusters), median)
names(pdf) <- c("sample", "cluster", "median.score")

mad_upper <- function(x){
  x <- x-median(x)
  return(mad(x[x>0], center = 0))
}

tests <- lapply(unique(meta$sample), function(x){
  sub <- pdf[pdf$sample == x,]
  scores_sample <- meta$doub.density[meta$sample == x]
  sub$p.value   <- pnorm(sub$median.score, mean = median(sub$median.score), sd = mad_upper(sub$median.score), lower.tail = FALSE)
  return(sub)
})
pdf <- do.call(rbind, tests)
pdf$fdr     <- p.adjust(pdf$p.value, method = "fdr")
pdf$n.cells <- sapply(1:nrow(pdf), function(row){
  sum(clusters == pdf$cluster[row] & meta$sample == pdf$sample[row])
})

pdf$frac.cells <- sapply(1:nrow(pdf), function(row){
  pdf$n.cells[row]/sum(pdf$n.cells[pdf$sample == pdf$sample[row]])
})


ggplot(pdf, aes(x = sample, y = log2(median.score + 1), col = factor(cluster))) +
  geom_point(data = pdf[pdf$fdr < 0.1,], size = 2) +
  geom_jitter(data = pdf[pdf$fdr >= 0.1,], col = "darkgrey", size = 0.4, width = 0.2, height = 0) +
  scale_colour_Publication() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = 1:max(pdf$sample)) +
  labs(x = "Sample", y = "log2(cluster median score + 1)")

ggsave(paste0(path2data, "QC/medianDoubletScores_clustersamplewise.pdf"), height = 10, width = 10)

## Testing doublet calls

doub.call <- paste(meta$sample, clusters) %in% paste(pdf$sample, pdf$cluster)[pdf$fdr < 0.1]
libs      <- Matrix::colSums(counts(sce))

ggplot(mapping = aes(x = doub.call, y = libs)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("FALSE" = "Singlets", "TRUE" = "Doublets")) +
  theme(axis.title.x = element_blank()) +
  labs(y = "#UMIs") +
  scale_y_log10()

ggsave(paste0(path2data, "QC/UMIsBySingletsDoublets.pdf"))

fd <- sapply(unique(meta$sample), function(x){
  sum(doub.call[meta$sample == x])/sum(meta$sample == x)
})

nc <- sapply(unique(meta$sample), function(x){
  sum(meta$sample == x)
})

frac_doubs <- ggplot(mapping = aes(x = nc, y = fd, fill = meta$stage[match(unique(meta$sample), meta$sample)])) +
  geom_point(shape = 21, size = 4, alpha = 1) +
  scale_fill_manual(values = stage_colours_extension, labels = stage_labels, name = "Stage") +
  labs(x = "Number of cells", y = "Fraction of called doublets")

frac_doubs

ggsave(paste0(path2data, "QC/DoubletCellsBySample.pdf"))

#Xist
xist_on <- counts(sce)[match("ENSMUSG00000086503", as.character(genes)),]>0

mouse_ensembl <- useMart("ensembl")
mouse_ensembl <- useDataset("mmusculus_gene_ensembl", mart = mouse_ensembl)

gene_map  <- getBM(attributes = c("ensembl_gene_id", "chromosome_name"), filters = "ensembl_gene_id", values = genes, mart = mouse_ensembl)
sex_genes <- c(gene_map[gene_map[,2] == "Y",1])
sex_genes <- sex_genes[sex_genes != "ENSMUSG00000096768"] #this gene has an X chromosome paralogue Erdr1 - exclude

y_on <- Matrix::colSums(counts(sce)[match(sex_genes, genes),] > 0 )>0
sex  <- xist_on & y_on

pdf$sex.frac <- sapply(1:nrow(pdf), function(row){
  select <- clusters == pdf$cluster[row] & meta$sample == pdf$sample[row]
  return(sum(sex[select])/sum(select))
})

pdf$stage <- meta$stage[match(pdf$sample, meta$sample)]
pdf$call  <- pdf$fdr < 0.1

xist_plot <- ggplot(pdf, aes(y = sex.frac, x = stage, fill = call)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Paired", labels = c("FALSE" = "Singlet", "TRUE" = "Doublet"), name = "") +
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(labels = stage_labels) +
  labs(y = "Fraction of cells coexpressing Xist + Y-CHR")

xist_plot

ggsave(paste0(path2data, "QC/mixed_sex_doublets.pdf"))

## Refining doublet calls

hvgs <- getHVGs(sce)
pca  <- prcomp_irlba(t(logcounts(sce)[hvgs,]), n = 50)

order_df <- meta[!duplicated(meta$sample), c("stage", "sample")]
order_df$ncells <- sapply(order_df$sample, function(x) sum(meta$sample == x))
order_df$stage  <- factor(order_df$stage,
                        levels = rev(c("E9.5", 
                                   "E9.25", 
                                   "E9.0", 
                                   "E8.75", 
                                   "E8.5")))

order_df <- order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
order_df$stage <- as.character(order_df$stage)

all_correct <- doBatchCorrect(counts = logcounts(sce)[rownames(sce) %in% hvgs,],
                             timepoints = meta$stage, 
                             samples = meta$sample, 
                             timepoint_order = order_df$stage, 
                             sample_order = order_df$sample, 
                             npc = 50)


tsne  <- Rtsne(all_correct, pca = FALSE)$Y
graph <- buildSNNGraph(all_correct, d = NA, transposed = TRUE)

set.seed(42)
clusts <- as.numeric(membership(cluster_louvain(graph)))
names(clusts) <- meta$cell

tab <- table(clusts, doub.call)
tab <- as.data.frame(sweep(tab, 1, rowSums(tab), "/"))
tab <- tab[as.logical(tab$doub.call),c(1,3)]
tab$p   <- pnorm(tab$Freq, mean = median(tab$Freq), sd = mad_upper(tab$Freq), lower.tail = FALSE)
tab$fdr <- p.adjust(tab$p, method = "fdr")

#Fraction of called doublets in all-data subclusters. Subclusters coloured red were considered to be composed of doublets."}
ggplot(mapping = aes(x = factor(rownames(tab), levels = rownames(tab)), y = tab[,2], fill = tab$fdr < 0.1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "darkgrey")) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank()) +
  labs(y = "Fraction doublets in cluster")

#ggsave(paste0(path2data, "QC/fraction_subclusters_with_doublets.pdf"))
ggsave(paste0(path2data, "QC/fraction_of_doublets_subclusters.pdf"))

cluster_calls <- clusts %in% tab$clusts[tab$fdr < 0.1]
cell_calls    <- doub.call
state <- rep("Singlet", nrow(meta))
state[cluster_calls] <- "Cluster-doublet"
state[cell_calls]    <- "Sample-doublet"
#ord <- order(factor(state, levels = c("Singlet", "Cluster-doublet", "Sample-doublet")))
ord <- sample(length(state), length(state))

ggsave(paste0(path2data, "QC/tsne_doublets.pdf"))

#Stripped nuclei

mt.counts   <- counts(sce)[which(genes %in% gene_map$ensembl_gene_id[gene_map$chromosome_name=="MT"]), ]
mt.fraction <- Matrix::colSums(mt.counts)/Matrix::colSums(counts(sce))
libsizes    <- Matrix::colSums(counts(sce))

meds <- data.frame(mt = sapply(unique(clusts), function(x) median(mt.fraction[clusts == x])),
                  lib = sapply(unique(clusts), function(x) median(libsizes[clusts == x])),
                  cluster = unique(clusts))

meds$stripped = meds$mt < 0.005

ggplot(meds, aes (x = mt, y = lib, col = stripped)) +
  geom_point() +
  labs(x = "Cluster median mitochondrial fraction",
       y = "Cluster median library size") +
  scale_x_log10() + 
  scale_y_log10() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))
meta$stripped = clusts %in% meds$cluster[meds$stripped]

ggsave(paste0(path2data, "QC/Stripped_nuclei.pdf"))

#Save the results
write.table(meta, file = paste0(path2data, "meta.tab"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
