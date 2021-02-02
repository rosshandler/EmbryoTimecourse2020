path2trajs  <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/W-OT/E85_trajs/"
path2integ  <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/Integrated_Atlas/"
path2trajs_somites <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/W-OT/somite_trajs/

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

#####
### Extension of somites using published subcluster annotations (Carolina, Jonny, et al.)
###

meta  <- readRDS(paste0(path2integ, "integrated_premeta_e85wot.rds"))

trajs <- read.table(paste0(path2trajs_somites, "somites_e85_trajectory.txt")), header = TRUE)

meta_tmp <- meta[
  meta$stage == "E8.5"  |
  meta$stage == "E8.75" |
  meta$stage == "E9.0"  |
  meta$stage == "E9.25" |
  meta$stage == "E9.5", ]

index1  <- grep("Paraxial mesoderm", meta_tmp$celltype.descendant)
index2  <- grep("Somitic mesoderm", meta_tmp$celltype.descendant)

meta_tmp <- meta_tmp[c(index1, index2), ]

index <- match(meta_tmp$cell, as.character(trajs$id))
index <- index[!is.na(index)]
traj_somites <- trajs[index, ]
rownames(traj_somites) <- traj_somites$id

#set.seed(42) if equally probable cells for more than one trajectory
#descendant <- unlist(apply(traj_somites[,-1], 1, function(x) {
#  colnames(traj_somites[,-1])[sample(which(x == max(x)))[1]]
#}))

descendant <- unlist(apply(traj_somites[,-1], 1, function(x) {
  colnames(traj_somites[,-1])[which(x == max(x))]
}))

cell <- as.character(traj_somites$id)

descendant <- data.frame(cbind(cell, descendant))

celltype.descendant.somites <- as.character(meta$celltype.descendant)
index <- match(descendant$cell, meta$cell)

descendant$descendant <- gsub("Head_mesoderm", "Cranial mesoderm",as.character(descendant$descendant))
descendant$descendant <- gsub("Anterior.most_somites", "Anterior Somitic Tissues",as.character(descendant$descendant))
descendant$descendant <- gsub("Posterior.most_somites", "Posterior Somitic Tissues",as.character(descendant$descendant))
descendant$descendant <- gsub("Presomitic_mesoderm", "Presomitic mesoderm",as.character(descendant$descendant))
                  
celltype.descendant.somites[index] <- as.character(descendant$descendant)

meta <- data.frame(cbind(meta[,1:(ncol(meta)-1)], celltype.descendant.somites, doub.density=meta[,ncol(meta)]))

saveRDS(meta, paste0(path2integ, "integrated_premeta_e85wot_somites.rds"))
