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

# HVGs, PCA and Batch correction 
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
