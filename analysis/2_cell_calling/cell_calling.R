##
## Cell calling

library(DropletUtils)
library(Matrix)

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

sample_sheet_file <- paste0(path2data, "embryo_extension_sample_metadata.txt")
sample_sheet      <- read.table(sample_sheet_file, header = TRUE)

sample_index     <- paste0("SIGA", sample_sheet$index)
sequencing_round <- gsub("-", "_", sample_sheet$sequencing_round)

sequencing_rounds <- unique(sequencing_round)

paths <- paste0("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/cellranger_output/",
  sequencing_round, "/", sample_index)

mtx_loc  <- paste0(paths, "/outs/matrix_unswapped.mtx")
bc_loc   <- paste0(paths, "/outs/barcodes_unswapped.tsv")

matrices <- bplapply(mtx_loc, readMM)
bcs      <- bplapply(bc_loc, function(x) read.table(x, header = FALSE, stringsAsFactors = FALSE)[,1])

#correct barcode sample number
for(i in 1:length(bcs)){
  bcs[[i]] <- paste0(bcs[[i]], "-", i)
}

set.seed(42)
#do call
outs <- lapply(matrices, emptyDrops, niters = 20000, ignore = 4999, BPPARAM = mcparam, lower = 100, retain = Inf)

#identify cells
sigs <- lapply(outs, function(x) x$FDR <= 0.01 & !is.na(x$FDR))
#subset the cells
cells    <- lapply(1:length(matrices), function(i) matrices[[i]][, sigs[[i]]])
barcodes <- lapply(1:length(bcs), function(i) bcs[[i]][sigs[[i]]])
#append
counts   <- do.call(cbind, cells)
barcodes <- do.call(c, barcodes)

writeMM(counts, file = paste0(path2data, "raw_counts.mtx"))
write.table(barcodes, file = paste0(path2data,"barcodes.tsv",
  col.names = FALSE, row.names = FALSE, quote = FALSE)
holder = file.copy(from = "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/cellranger_output/SLX_18551/SIGAA11/outs/genes_unswapped.tsv", 
                   to = paste0(path2data, "genes.tsv"),
                   overwrite = TRUE)

sample     <- 1:68
exp_design <- sample_sheet

exp_design$time_point[grep("^9$",exp_design$time_point)] <-
  paste0(exp_design$time_point[grep("^E9$",exp_design$time_point)],"9.0")
            
exp_design$time_point <- paste0("E",exp_design$time_point)

exp_design <- cbind(sample, exp_design)
colnames(exp_design)[7] <- "stage"
colnames(exp_design)[8] <- "embryo_tube"

tube_name  <- strsplit(as.character(exp_design$embryo_tube), split="_")
tube_name  <- unlist(lapply(tube_name, function(x) x[2]))
exp_design <- cbind(exp_design[,1:7],tube_name,exp_design[,8:10])

write.table(exp_design, file = paste0(path2data, "embryo_extension_experimental_design.txt"), quote = FALSE, row.names = FALSE, sep="\t")
                 
summary_df       <- data.frame(sample = 1:length(matrices), value = sapply(cells, ncol))
summary_df$stage <- exp_design$stage[match(summary_df$sample, exp_design$sample)]

tab <- acast(summary_df, sample ~ stage, fill = 0)
tab <- rbind(tab, colSums(tab))
rownames(tab)[nrow(tab)] <- "Total"

print(summary_df)

