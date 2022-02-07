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

#Also prepare experimental design for submission
df <- read.table("embryo_extension_experimental_design.txt", header=TRUE)

tube_name_corrected <- as.character(df$tube)
tube_name_corrected <- gsub("anterior", "Anterior",tube_name_corrected)
tube_name_corrected <- gsub("posterior", "Posterior",tube_name_corrected)
tube_name_corrected <- gsub("head", "Anterior section",tube_name_corrected)
tube_name_corrected <- gsub("torso", "Medial section",tube_name_corrected)
tube_name_corrected <- gsub("tail", "Posterior section",tube_name_corrected)

embryo_tube_name_corrected <- as.character(df$embryo_tube)
embryo_tube_name_corrected <- gsub("_"," ",gsub("anterior", "Anterior",embryo_tube_name_corrected))
embryo_tube_name_corrected <- gsub("_"," ",gsub("posterior", "Posterior",embryo_tube_name_corrected))
embryo_tube_name_corrected <- gsub("_"," ",gsub("head", "Anterior section",embryo_tube_name_corrected))
embryo_tube_name_corrected <- gsub("_"," ",gsub("torso", "Medial section",embryo_tube_name_corrected))
embryo_tube_name_corrected <- gsub("_"," ",gsub("tail", "Posterior section",embryo_tube_name_corrected))

organism <- rep("Mus musculus", nrow(df))
strain   <- rep("C57BL/6J", nrow(df))
genotype <- rep("wild type genotype", nrow(df))
sex <- rep("undetermined", nrow(df))
material_type <- rep("organism part", nrow(df))

df <- cbind(df, embryo_tube_name_corrected, tube_name_corrected)

df_sub <- cbind(df[,"sample"], organism, strain, paste0("embryonic day ", df[,"stage"]), as.character(df[,"embryo_tube_name_corrected"]), as.character(df[,"tube_name_corrected"]), genotype, sex, material_type)
colnames(df_sub) <- c("source name", "organism", "strain", "developmental stage", "individual", "organism part",	"genotype",	"sex", "material type")
df_sub <- as.data.frame(df_sub)

write.table(df_sub, file = "embryo_extension_experimental_design_submission.txt", quote=FALSE, row.names=FALSE, sep="\t")
