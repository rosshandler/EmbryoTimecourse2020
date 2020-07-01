path2trajs  <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/W-OT/E85_trajs/"
path2integ  <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/Integrated_Atlas/"

trajs <- read.table(paste0(path2trajs, "e85_trajectory.txt"), header = TRUE)

premeta_e85annot    <- readRDS(paste0(path2integ, "integrated_premeta_e85mapping.rds"))
celltype.descendant <- as.character(premeta_e85annot$celltype.extended)

meta_tmp <- premeta_e85annot[
  premeta_e85annot$stage == "E8.75" |
  premeta_e85annot$stage == "E9.0"  |
  premeta_e85annot$stage == "E9.25" |
  premeta_e85annot$stage == "E9.5", ]

index <- match(meta_tmp$cell, as.character(trajs$id))
index <- index[!is.na(index)]
trajs <- trajs[index, ]

descendant <- apply(trajs[,-1], 1, function(x) {
  colnames(trajs[,-1])[which(x == max(x))]
  })

descendant  <- gsub("\\.", " ", descendant)
descendant  <- gsub("Def  endoderm", "Def. endoderm", descendant)
descendant  <- gsub("Forebrain Midbrain Hindbrain", "Forebrain/Midbrain/Hindbrain", descendant)

cell <- as.character(trajs$id)

descendant <- data.frame(cbind(cell, descendant))

index <- match(descendant$cell, premeta_e85annot$cell)
celltype.descendant[index] <- as.character(descendant$descendant)

meta <- data.frame(cbind(premeta_e85annot, celltype.descendant))

names_order <- c("cell", "sample", "stage", "stage.mapped", "stage.collapsed", "stage.mapped.collapsed",
"somite_count", "tube_name", "celltype", "celltype.extended", "celltype.descendant", "doub.density")

meta <- meta[, names_order]

saveRDS(meta, paste0(path2integ, "integrated_premeta_e85wot.rds"))
