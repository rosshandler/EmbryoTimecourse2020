###############
# Tube_Names correction
#
# After chatting with some corrections were made to tube_names 

path2integ     <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/Integrated_Atlas/"

path2integ     <- "/mnt/b/iimaz/embryos/"

meta <- readRDS(paste0(path2integ, "integrated_meta.rds"))

#Corrected tube names: anterior, medial, posterior (instead of head, torso and tail)

tube_name_corrected <- as.character(meta$tube_name)
tube_name_corrected <- gsub("anterior", "Anterior",tube_name_corrected)
tube_name_corrected <- gsub("posterior", "Posterior",tube_name_corrected)
tube_name_corrected <- gsub("head", "Anterior_section",tube_name_corrected)
tube_name_corrected <- gsub("torso", "Medial_section",tube_name_corrected)
tube_name_corrected <- gsub("tail", "Posterior_section",tube_name_corrected)

meta <- cbind(meta[,1:8],tube_name_corrected,meta[,9:ncol(meta)])

saveRDS(meta, paste0(path2integ, "integrated_meta.rds"))
