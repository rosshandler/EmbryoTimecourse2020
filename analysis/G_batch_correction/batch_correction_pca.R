library(Rtsne)
library(irlba)
library(scran)
library(Matrix)
library(ggplot2)
library(batchelor)
library(easyGgplot2)
library(SingleCellExperiment)

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

########################
########################
#### Extended atlas pre-processing:
####
#### Load QC filtered SingleCellExperiment and metadata
#### Compute Highly Variable Genes and PCA
#### Compute batch corrected PCA (only atlas extension)
#### Compute tSNE to inspect batch correction
########################

source("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/core_scripts/core_functions.R")

load_data2020(remove_doublets = TRUE, remove_stripped = TRUE)

path2data  <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/"

chrY.genes.file <- "/hps/research1/marioni/ivan/EmbryoTimeCourse/atlas/data/ygenes.tab"

# HVGs, PCA and Batch correction

trend  <- scran::trendVar(sce, use.spikes = FALSE, loess.args = list(span = 0.05))
decomp <- scran::decomposeVar(sce, fit = trend)
decomp <- decomp[decomp$mean > 1e-3,]

xist <- "ENSMUSG00000086503"
ychr <- read.table(chrY.genes.file, stringsAsFactors = FALSE)[,1]

decomp <- decomp[!rownames(decomp) %in% c(xist, ychr, other),]
decomp$FDR <- p.adjust(decomp$p.value, method = "fdr")

hvgs_colour <- rep("black",nrow(decomp))
#hvgs_colour[decomp$p.value < 0.05] <- "yellow3"
hvgs_colour[decomp$FDR < 0.05] <- "yellow3"

pdf("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/hvgs_trend.pdf")
pl.index <- order(decomp$p.value, decreasing=TRUE)
plot(trend$mean[pl.index], trend$var[pl.index], col=hvgs_colour[pl.index], pch=19, cex=.75, 
     xlab="Mean log expression", ylab="Variance of log expression", bty="n", ylim=c(0,3))
#x <- sort(trend$mean)
#lines(x, trend$trend(x), col="dodgerblue", lwd=2)
curve(trend$trend(x), col="red", lwd=2, add=TRUE)
grid()
dev.off()

dec <- modelGeneVar(sce)
hvgs_colour <- rep("black",nrow(dec))
hvgs_colour[dec$FDR < 0.05] <- "yellow3"

pdf("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/hvgs_trend_2020.pdf")
pl.index <- order(dec$FDR, decreasing=TRUE)
plot(dec$mean[pl.index], dec$total[pl.index], col=hvgs_colour[pl.index], pch=19, cex=.75, 
     xlab="Mean log expression", ylab="Variance of log expression", bty="n", ylim=c(0,3))
curve(metadata(dec)$trend(x), add=TRUE, col="red")
grid()
dev.off()



nPCs <- 50
hvgs     <- getHVGs(sce, computer = "ebi")
base_pca <- prcomp_irlba(t(logcounts(sce)[rownames(sce) %in% hvgs,]), n = nPCs)$x

#Get order: oldest to youngest; most cells to least cells
order_df <- meta[!duplicated(meta$sample), c("stage", "sample")]
order_df$ncells <- sapply(order_df$sample, function(x) sum(meta$sample == x))
order_df$stage  <- factor(order_df$stage, 
                        levels = rev(c("E9.5", 
                                   "E9.25", 
                                   "E9.0", 
                                   "E8.75", 
                                   "E8.5")))
order_df       <- order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
order_df$stage <- as.character(order_df$stage)

corrected_pca <- doBatchCorrect(counts = logcounts(sce)[rownames(sce) %in% hvgs,], 
                             timepoints = meta$stage, 
                             samples = meta$sample, 
                             timepoint_order = order_df$stage, 
                             sample_order = order_df$sample, 
                             npc = nPCs,
                             BPPARAM = mcparam)
                             
saveRDS(corrected_pca, file = paste0(path2data, "corrected_pcas.rds"))

save(sce, meta, hvgs, base_pca, corrected_pca, file=paste0(path2data, "embryo_extension_data_hvgs_pca.RData"))

### tSNE
#corrected_pca <- readRDS(paste0(path2data, "corrected_pcas.rds"))

base_tsne      <- Rtsne(base_pca, pca = FALSE)$Y
corrected_tsne <- Rtsne(corrected_pca, pca = FALSE)$Y

save(base_tsne, corrected_tsne, file=paste0(path2data, "embryo_extension_tsne.RData"))

base_tsne <- as.data.frame(base_tsne)
base_tsne$sample <- meta$sample
base_tsne$sample <- meta$stage
base_tsne$state  <- "Uncorrected"
corrected_tsne        <- as.data.frame(corrected_tsne)
corrected_tsne$sample <- meta$sample
corrected_tsne$stage  <- corrected_tsne$stage
corrected_tsne$state  <- "Corrected"

bc_tsne <- rbind(base_tsne, corrected_tsne)
reorder <- sample(nrow(bc_tsne), nrow(bc_tsne))

ggplot(bc_tsne[reorder,], aes(x = V1, y = V2, col = factor(sample))) +
  geom_point(size = 0.4) +
  scale_colour_Publication() +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  facet_wrap(~state, nrow = 2)
  
ggsave(paste0(path2data, "QC/batchcorrection_tsne.pdf"), height=15, width=10)

ggplot(bc_tsne[reorder,], aes(x = V1, y = V2, col = factor(stage))) +
  geom_point(size = 0.4) +
  scale_color_manual(values=stage_colours_extension,name = "Embryo stage", labels = names(stage_colours_extension)) +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  facet_wrap(~state, nrow = 2)
