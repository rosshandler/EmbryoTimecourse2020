##
## Unswap molecule barcodes

library(Matrix)
library(DropletUtils)

sample_sheet_file <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/embryo_extension_sample_metadata.txt"
sample_sheet      <- read.table(sample_sheet_file, header = TRUE)

sample_index     <- paste0("SIGA", sample_sheet$index)
sequencing_round <- gsub("-", "_", sample_sheet$sequencing_round)

sequencing_rounds <- unique(sequencing_round)

paths <- paste0("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/cellranger_output/",
  sequencing_round, "/", sample_index)

index    <- grep(sequencing_rounds[1],paths)
mol_loc  <- paste0(paths[index], "/outs/molecule_info.h5")
out_loc  <- paste0(paths[index], "/outs/matrix_unswapped.mtx")
bc_loc   <- paste0(paths[index], "/outs/barcodes_unswapped.tsv")
gene_loc <- paste0(paths[index], "/outs/genes_unswapped.tsv")

unswapped <- swappedDrops(mol_loc, get.swapped = TRUE)

ratios <- sapply(1:length(unswapped$cleaned), function(i){
  sum(unswapped$swapped[[i]])/(sum(unswapped$cleaned[[i]]) + sum(unswapped$swapped[[i]]))
})

ratios_SLX_18551 <- ratios

saveRDS(ratios_SLX_18551, paste0("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/","swapping_ratios_SLX_18551.rds"))

for(i in 1:length(mol_loc)){
  writeMM(unswapped$cleaned[[i]], file = out_loc[i])
  write.table(colnames(unswapped$cleaned[[i]]), file = bc_loc[i], col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(rownames(unswapped$cleaned[[i]]), file = gene_loc[i], col.names = FALSE, row.names = FALSE, quote = FALSE)
}

index    <- grep(sequencing_rounds[2],paths)
mol_loc  <- paste0(paths[index], "/outs/molecule_info.h5")
out_loc  <- paste0(paths[index], "/outs/matrix_unswapped.mtx")
bc_loc   <- paste0(paths[index], "/outs/barcodes_unswapped.tsv")
gene_loc <- paste0(paths[index], "/outs/genes_unswapped.tsv")

unswapped <- swappedDrops(mol_loc, get.swapped = TRUE)

ratios = sapply(1:length(unswapped$cleaned), function(i){
  sum(unswapped$swapped[[i]])/(sum(unswapped$cleaned[[i]]) + sum(unswapped$swapped[[i]]))
})

ratios_SLX_18597 <- ratios
saveRDS(ratios_SLX_18597, paste0("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/","swapping_ratios_SLX_18597.rds"))

for(i in 1:length(mol_loc)){
  writeMM(unswapped$cleaned[[i]], file = out_loc[i])
  write.table(colnames(unswapped$cleaned[[i]]), file = bc_loc[i], col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(rownames(unswapped$cleaned[[i]]), file = gene_loc[i], col.names = FALSE, row.names = FALSE, quote = FALSE)
}

index    <- grep(sequencing_rounds[3],paths)
mol_loc  <- paste0(paths[index], "/outs/molecule_info.h5")
out_loc  <- paste0(paths[index], "/outs/matrix_unswapped.mtx")
bc_loc   <- paste0(paths[index], "/outs/barcodes_unswapped.tsv")
gene_loc <- paste0(paths[index], "/outs/genes_unswapped.tsv")

unswapped <- swappedDrops(mol_loc, get.swapped = TRUE)

ratios = sapply(1:length(unswapped$cleaned), function(i){
  sum(unswapped$swapped[[i]])/(sum(unswapped$cleaned[[i]]) + sum(unswapped$swapped[[i]]))
})

ratios_SLX_18943 <- ratios
saveRDS(ratios_SLX_18943, paste0("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/","swapping_ratios_SLX_18943.rds"))

for(i in 1:length(mol_loc)){
  writeMM(unswapped$cleaned[[i]], file = out_loc[i])
  write.table(colnames(unswapped$cleaned[[i]]), file = bc_loc[i], col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(rownames(unswapped$cleaned[[i]]), file = gene_loc[i], col.names = FALSE, row.names = FALSE, quote = FALSE)
}

ratios_SLX_18551 <- readRDS(paste0("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/","swapping_ratios_SLX_18551.rds"))

ratios_SLX_18597 <- readRDS(paste0("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/","swapping_ratios_SLX_18597.rds"))

ratios_SLX_18943 <- readRDS(paste0("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/","swapping_ratios_SLX_18943.rds"))

pdf(paste0("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/", "swappping_ratio_hist.pdf"))
hist(c(ratios_SLX_18551, ratios_SLX_18597, ratios_SLX_18943)*100, col = "gray", xlab="Swapping molecules percentage", ylab="Numner of samples", main="")
dev.off()
