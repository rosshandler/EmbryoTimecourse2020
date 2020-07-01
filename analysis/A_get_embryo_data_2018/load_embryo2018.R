celltype_colours <- c("Epiblast" = "#635547",
                     "Primitive Streak" = "#DABE99",
                     "Caudal epiblast" = "#9e6762",
                     
                     "PGC" = "#FACB12",
                     
                     "Anterior Primitive Streak" = "#c19f70",
                     "Notochord" = "#0F4A9C",
                     "Def. endoderm" = "#F397C0",
                     "Gut" = "#EF5A9D",
                     
                     "Nascent mesoderm" = "#C594BF",
                     "Mixed mesoderm" = "#DFCDE4",
                     "Intermediate mesoderm" = "#139992",
                     "Caudal Mesoderm" = "#3F84AA",
                     "Paraxial mesoderm" = "#8DB5CE",
                     "Somitic mesoderm" = "#005579",
                     "Pharyngeal mesoderm" = "#C9EBFB",
                     "Cardiomyocytes" = "#B51D8D",
                     "Allantois" = "#532C8A",
                     "ExE mesoderm" = "#8870ad",
                     "Mesenchyme" = "#cc7818",
                     
                     "Haematoendothelial progenitors" = "#FBBE92",
                     "Endothelium" = "#ff891c",
                     "Blood progenitors 1" = "#f9decf",
                     "Blood progenitors 2" = "#c9a997",
                     "Erythroid1" = "#C72228",
                     "Erythroid2" = "#f79083",
                     "Erythroid3" = "#EF4E22",
                     
                     "NMP" = "#8EC792",
                     
                     "Rostral neurectoderm" = "#65A83E",
                     "Caudal neurectoderm" = "#354E23",
                     "Neural crest" = "#C3C388",
                     "Forebrain/Midbrain/Hindbrain" = "#647a4f",
                     "Spinal cord" = "#CDE088",
                     
                     "Surface ectoderm" = "#f7f79e",
                     
                     "Visceral endoderm" = "#F6BFCB",
                     "ExE endoderm" = "#7F6874",
                     "ExE ectoderm" = "#989898",
                     "Parietal endoderm" = "#1A1A1A"
                     
)

haem_colours <- c(
  "Mes1"= "#c4a6b2",
  "Mes2"= "#ca728c",
  
  "Cardiomyocytes" =  "#B51D8D",  
  
  "BP1" = "#6460c5",
  "BP2" = "#96b8e4",
  "Haem3"= "#02f9ff",
  "BP3" = "#07499f",
  "BP4" = "#036ef9",
  
  "Haem1"= "#bb22a7",
  "Haem2" = "#f695e9",
  "Haem4" = "#4c4a81",
  
  "EC1"= "#006737",
  
  "EC2" = "#5eba46",
  "EC3" = "#818068",
  "EC4"="#d6de22",
  "EC5"="#5f7238",
  "EC6"="#2ab572",
  "EC7"="#000000",
  "EC8"="#a0cb3b",
  
  "Ery1"="#f67a58",
  "Ery2" ="#a26724",
  "Ery3"="#cdaf7f",
  "Ery4"= "#625218",
  
  "My" = "#c62127",
  "Mk"= "#f6931d"
)

stage_colours <- c("E6.5" = "#D53E4F",
                  "E6.75" = "#F46D43",
                  "E7.0" = "#FDAE61",
                  "E7.25" = "#FEE08B",
                  "E7.5" = "#FFFFBF",
                  "E7.75" = "#E6F598",
                  "E8.0" = "#ABDDA4",
                  "E8.25" = "#66C2A5",
                  "E8.5" = "#3288BD",
                  "mixed_gastrulation" = "#A9A9A9")

stage_labels <- names(stage_colours)
names(stage_labels) <- names(stage_colours)
stage_labels[10]    <- "Mixed"

load_data <- function(normalise = TRUE, remove_doublets = FALSE, remove_stripped = FALSE, load_corrected = FALSE, path2atlas=NULL){
  
  if(is.null(path2atlas)){
    path2atlas <- "/hps/research1/marioni/ivan/EmbryoTimeCourse/atlas/data/"

  }
  
  if(load_corrected & (!remove_doublets | !remove_stripped)){
    message("Using corrected PCs, also removing doublets + stripped now.")
    remove_doublets <- TRUE
    remove_stripped <- TRUE
  }
  

  #counts  <- readRDS(paste0(path2atlas, "raw_counts.rds"))  
  counts  <- Matrix::readMM(paste0(path2atlas, "raw_counts.mtx"))
  genes   <- read.table(paste0(path2atlas, "genes.tsv"), stringsAsFactors = F)
  meta    <- read.table(paste0(path2atlas, "meta.tab"), header = TRUE, 
             sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
  
  rownames(counts) <- genes[,1] #ensembl
  colnames(counts) <- meta$cell

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list("counts" = counts))
  
  if(normalise){
    sfs <- read.table(paste0(path2atlas, "sizefactors.tab"), stringsAsFactors = F)[,1]
    SingleCellExperiment::sizeFactors(sce) <- sfs
    sce <- scater::normalize(sce)
  }
  
  if(remove_doublets){
    sce  <- scater::normalize(sce[,!meta$doublet])
    meta <- meta[!meta$doublet,]
  }
  
  if(remove_stripped){
    sce  <- scater::normalize(sce[,!meta$stripped])
    meta <- meta[!meta$stripped, ]
  }
  
  if(load_corrected){
    corrected <- readRDS(paste0(path2atlas, "corrected_pcas.rds"))
    assign("corrected", corrected, envir = .GlobalEnv)
    
  }
  
  assign("genes", genes, envir = .GlobalEnv)
  assign("meta", meta, envir = .GlobalEnv)
  assign("sce", sce, envir = .GlobalEnv)
  
  
  invisible(0)
}

load_data(normalise = TRUE, remove_doublets = TRUE, remove_stripped = TRUE, load_corrected = TRUE, path2atlas = NULL)

save.image(file = '/hps/research1/marioni/ivan/EmbryoTimeCourse/atlas/data/embryo_data.RData')
