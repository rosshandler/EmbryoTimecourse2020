library(Matrix)
library(biomaRt)
library(ggplot2)
library(scales)
library(viridis)

path2data  <- "/hps/research1/marioni/ivan/EmbryoTimeCourse2020/data/"
counts     <- readMM(paste0(path2data, "raw_counts.mtx"))
genes      <- readLines(paste0(path2data, "genes.tsv"))
barcodes   <- readLines(paste0(path2data, "barcodes.tsv"))
exp_design <- read.table(paste0(path2data, "embryo_extension_experimental_design.txt"),
  sep = "\t", header = TRUE)

lib.sizes    <- colSums(counts)
ngenes       <- colSums(counts > 0)
split_bc     <- strsplit(as.character(barcodes), "-", fixed = TRUE)
samples      <- sapply(split_bc, function(x) x[2])
sample_names <- as.character(exp_design$sample_name[match(samples, exp_design$sample)])
sequencing_round  <- gsub("-", "_", exp_design$sequencing_round)
sequencing_rounds <- unique(sequencing_round)
sample_index      <- paste0("SIGA", exp_design$index)

mouse_ensembl <- useMart("ensembl")
mouse_ensembl <- useDataset("mmusculus_gene_ensembl", mart = mouse_ensembl)

paths <- paste0("/hps/research1/marioni/ivan/EmbryoTimeCourse2020/cellranger_output/",
  sequencing_round, "/", sample_index,"/outs/metrics_summary.csv")

names(paths) <- as.character(unique(exp_design$sample_name))

#Comparing sequencing depth

## UMI counts

batch = sapply(samples, function(x){
  if(x %in% 1:16){
    return(1)
  } else if(x %in% 17:32){
    return(2)
  } else{
    return(3)
  }
})

plot_df <- data.frame(lib = lib.sizes, sample = samples, batch = batch)

# batches
plot_df$batch <- gsub("^1$","SLX_18551", plot_df$batch)
plot_df$batch <- gsub("^2$","SLX_18597", plot_df$batch)
plot_df$batch <- gsub("^3$","SLX_18943", plot_df$batch)

ggplot(plot_df, aes (x = factor(batch), y = lib)) +
  geom_boxplot() +
  theme_bw() + 
labs(x = "Batch", y = "Number of UMIs") +
scale_y_log10(breaks = c(1000, 5000, 10000, 50000, 100000),
    labels = c("1,000", "5,000", "10,000", "50,000", "100,000"))

ggsave(paste0(path2data, "QC/UMISbyBatch.pdf"))
#ggsave(paste0(path2data, "QC/UMISbyBatch2018.pdf"))

## Sequencing saturation

pdf <- data.frame(
  sample   = as.character(exp_design$sample_name),
  med_lib  = sapply(as.character(exp_design$sample_name),function(x) median(lib.sizes[sample_names == x])),
  mean_lib = sapply(as.character(exp_design$sample_name), function(x) mean(lib.sizes[sample_names == x])),
  stage    = sapply(as.character(exp_design$sample_name),
    function(x) exp_design$stage[match(x, exp_design$sample_name)]),
  saturation = sapply(as.character(unique(exp_design$sample_name)), function(x){
    tab  = read.table(paths[x], header = TRUE, sep = ",")
    stat = tab[1, "Sequencing.Saturation"]
    return(as.numeric(substr(stat, 1, 4)))}
  )
)
                 
colors <- c(brewer_pal(palette = "Spectral")(length(unique(pdf$stage))-1), "lightgrey")

stage_colours_extension <- c(colorRampPalette(c("red", "orange", "yellow", "green", "blue"),
  space = "Lab")(13),"#A9A9A9")

names(stage_colours_extension) <- c(
  "E6.5",
  "E6.75",
  "E7.0",
  "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5",
  "E8.75",
  "E9.0",
  "E9.25",
  "E9.5",
  "mixed_gastrulation")

ggplot(pdf, aes(x = med_lib, y = saturation, label = 1:nrow(pdf), fill = stage)) +
  geom_label(alpha = 0.8) +
  labs(x = "Median library size", y = "Sequencing Saturation (%)") +
  theme_bw() +
  scale_fill_manual(values = colors)

ggsave(paste0(path2data, "QC/sequencing_saturation.pdf"))
ggsave(paste0(path2data, "QC/sequencing_saturation2018.pdf"))

ggplot(pdf, aes(x = med_lib, y = saturation, label = 1:nrow(pdf), fill = stage)) +
  geom_label(alpha = 0.8) +
  labs(x = "Median library size", y = "Sequencing Saturation (%)") +
  theme_bw() +
  scale_fill_manual(values = stage_colours_extension)

ggsave(paste0(path2data, "QC/sequencing_saturation_colours_ext.pdf"))
#ggsave(paste0(path2data, "QC/sequencing_saturation2018_colours_ext.pdf"))

#"Cells dropped according to UMI threshold"
plot_change <- function(barcodes, logical_keep){

  sample_names <- as.character(exp_design$sample_name[match(samples, exp_design$sample)])

  pdf <- data.frame(Sample  = as.character(exp_design$sample_name),
                   Total    = sapply(as.character(exp_design$sample_name), function(x) sum(sample_names == x)),
                   Dropped  = sapply(as.character(exp_design$sample_name),
                                    function(x) sum(!logical_keep[sample_names == x])),
                   Retained = sapply(as.character(exp_design$sample_name),
                                     function(x) sum(logical_keep[sample_names == x])))

  p <- ggplot(data = pdf) +
          geom_bar(mapping = aes(y = Total, 
                                 x = factor(Sample, levels = as.character(exp_design$sample_name))), 
                   fill = "darkgrey",
                   stat = "identity") +
          geom_bar(mapping = aes(y = Retained, 
                                 x = factor(Sample, levels = as.character(exp_design$sample_name))), 
                   fill = "coral",
                   stat = "identity") +
          geom_segment(mapping = aes(y.   = Total, 
                                     yend = Retained, 
                                     x    = factor(Sample, levels = as.character(exp_design$sample_name)),  
                                     xend = factor(Sample, levels = as.character(exp_design$sample_name))),
                       arrow = arrow(length = unit(0.1, "inches"))) +
          theme_bw() +
          labs(y = "Number of cells", x= "Sample") + theme(axis.text.x=element_text(angle = -90, hjust = 0))
                             
  pdf <- rbind(pdf, data.frame(Sample  = "Total",
                              Total    = length(barcodes), 
                              Dropped  = sum(!logical_keep), 
                              Retained = sum(logical_keep)))
    
  return(list(plot = p, df = pdf))
}


# Cell complexity thresholding

qplot(lib.sizes, ngenes, col = ifelse(ngenes < 1000, "drop", "keep")) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() + 
  theme(text = element_text(size=20),legend.position = "none")  +
  labs(x = "UMI count", y = "Number of expressed genes") +
  scale_color_manual(values = c("drop" = "grey50", "keep" = "black"), name = "")


ggsave(paste0(path2data, "QC/cell_complexity.pdf"))
#ggsave(paste0(path2data, "QC/cell_complexity2018.pdf"))

gene_drop = plot_change(barcodes, logical_keep = ngenes > 1000 )
print(gene_drop$plot)

write.table(gene_drop$df, file = paste0(path2data, "QC/cellsBysample.txt",
  quote =FALSE, row.names =FALSE, sep="\t")

ggsave(paste0(path2data, "QC/cellsBysample.pdf"))
#ggsave(paste0(path2data, "QC/cellsBysample2018.pdf"))

#####
##### Filter 1
counts    <- counts[, ngenes > 1000]
barcodes  <- barcodes[ngenes > 1000]
ngenes    <- ngenes[ngenes > 1000]
lib.sizes <- colSums(counts)

# Mitochondrial gene expression
sym_map   <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id", values = genes, mart = mouse_ensembl)
write.table(sym_map, file = paste0(path2data, "genes_symsmap.tab"), quote =FALSE, row.names =FALSE, sep="\t")

gene_map  <- getBM(attributes=c("ensembl_gene_id", "chromosome_name"),
  filters = "ensembl_gene_id", values = genes, mart = mouse_ensembl)
write.table(gene_map, file = paste0(path2data, "genes_chrsmap.tab"), quote =FALSE, row.names =FALSE, sep="\t")

mt.index    <- gene_map$chromosome_name == "MT"
mt.counts   <- counts[which(genes %in% gene_map$ensembl_gene_id[mt.index]), ]
mt.fraction <- colSums(mt.counts)/lib.sizes

mt.p   <- pnorm(mt.fraction, mean = median(mt.fraction), sd = mad(mt.fraction), lower.tail = FALSE)
mt.lim <- min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 0.05)])

#Threhdold
mt.lim
#[1] 0.04680187
            
qplot(lib.sizes, mt.fraction, col = ifelse(mt.fraction>mt.lim, "drop", "keep")) +
  scale_x_log10() +
  labs(x = "UMI count", y = "MT read fraction") +
  theme_minimal() + 
  theme(text = element_text(size=20),legend.position = "none")  +
  scale_color_manual(values = c("drop" = "grey50", "keep" = "black"), name = "")

ggsave(paste0(path2data, "QC/mtreadfraction.pdf"))

mt_drop = plot_change(barcodes, logical_keep = mt.fraction<mt.lim )
print(mt_drop$plot)

ggsave(paste0(path2data, "QC/mtreadfraction_cellsdropped.pdf"))

write.table(mt_drop$df, file = paste0(path2data, "QC/cellsBysample_mt.txt"),
  quote =FALSE, row.names =FALSE, sep="\t")

#####
##### Filter 2
counts      <- counts[, mt.fraction < mt.lim]
barcodes    <- barcodes[mt.fraction < mt.lim]
mt.fraction <- mt.fraction[mt.fraction < mt.lim]

#Number of detected genes

lib.sizes <- colSums(counts)
n.genes   <- colSums(counts>0)

split_bc  <- strsplit(as.character(barcodes), "-", fixed = T)
bcs       <- sapply(split_bc, function(x) x[1])
samples   <- sapply(split_bc, function(x) x[2])
sample_names <- as.character(exp_design$sample_name[match(samples, exp_design$sample)])

batch = sapply(samples, function(x){
  if(x %in% 1:16){
    return(1)
  } else if(x %in% 17:32){
    return(2)
  } else{
    return(3)
  }
})

plot_df = data.frame(lib = lib.sizes, genes = n.genes, sample = samples)
plot_df = plot_df[sample(nrow(plot_df), nrow(plot_df), replace = FALSE),]

p = ggplot(plot_df, aes(x = lib, y = genes)) +
  geom_point(alpha = 0.5) +
  scale_colour_Publication() +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  facet_wrap(~factor(sample, levels = unique(sample)[order(nchar(unique(as.character(sample))), unique(sample))]), ncol = 3) +
  labs(x = "Number of UMIs", y = "Number of detected genes")
suppressWarnings(plot(p))

#Inter-sample comparisons

batch = gsub("^1$","SLX_18551", batch)
batch = gsub("^2$","SLX_18597", batch)
batch = gsub("^3$","SLX_18943", batch)

p = ggplot(data.frame(lib = lib.sizes, sample = samples), aes (x = factor(sample, levels = unique(sample)), 
                                                           y = lib.sizes, fill = batch)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_log10() +
  labs(x = "Sample", y = "Number of UMIs")+
  scale_color_brewer(palette = "Set1", name = "Sequencing\nbatch")
suppressWarnings(plot(p))

ggsave(paste0(path2data, "QC/UMIs_distribution.pdf"))

p = ggplot(data.frame(mt = mt.fraction, sample = samples), aes (x = factor(sample, levels = unique(sample)), 
                                                           y = mt, fill = factor(batch))) +
  geom_boxplot() +
  theme_bw() +
  scale_y_log10() +
  labs(x = "Sample", y = "Mitochondrial gene fraction") +
  scale_color_brewer(palette = "Set1", name = "Sequencing\nbatch")
suppressWarnings(plot(p))

ggsave(paste0(path2data, "QC/UMIs_distribution_mt.pdf"))


pdf = data.frame(mt = mt.fraction, sample = samples, lib.size = lib.sizes,
  stage = exp_design$stage[match(samples, exp_design$sample)])

p = ggplot(pdf, aes(x = factor(stage), y= lib.size)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_log10() +
  labs(x = "Stage", y = "Number of UMIs")
suppressWarnings(plot(p))

ggsave(paste0(path2data,"QC/UMISbyStage.pdf"))

p = ggplot(pdf, aes(x = factor(stage), y=mt.fraction)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_log10() +
  labs(x = "Stage", y = "Mitchondrial gene fraction")
suppressWarnings(plot(p))

ggsave(paste0(path2data,"QC/UMISbyStage_mt.pdf"))

#Make metadata

meta <- data.frame(cell = paste0("cell_", 1:ncol(counts)),
  barcode      = bcs,
  sample       = as.numeric(samples),
  sample_name  = sample_names,
  sample_index = exp_design$sample_index[match(samples, as.character(exp_design$sample))],
  stage        = exp_design$stage[match(samples, as.character(exp_design$sample))],
  tube_name    = exp_design$tube_name[match(samples, as.character(exp_design$sample))],
  embryo_tube  = exp_design$embryo_tube[match(samples, as.character(exp_design$sample))],
  batch        = exp_design$sequencing_round[match(samples, as.character(exp_design$sample))],
  somite_count = exp_design$somite_count[match(samples, as.character(exp_design$sample))]
)

#Write summary

tab <- as.matrix(table(meta$sample, meta$stage))
new_rownames  <- c("Total", rownames(tab))

tab <- rbind(sapply(colnames(tab), function(x) sum(meta$stage == x)), tab)
rownames(tab) <- new_rownames

write.table(paste0(path2data,"QC/cellsBySampleAndStage.txt"), row.names = TRUE, quote = FALSE)

exp <- exp_design
exp$obs   <- sapply(exp$sample, function(x) sum(meta$sample == x))
exp$cells <- as.numeric(exp$estimated_captured_60perc)

#Expected vs. retained cells

ggplot(exp, aes(x = cells, y = obs, col = embryo_tube)) +
  geom_point() +
  #scale_colour_Publication(name = "Sample") +
  scale_x_log10(breaks = c(500, 1000, 5000)) +
  scale_y_log10(breaks = c(500, 1000, 5000)) +
  labs(x = "Predicted number of cells", y = "Observed number of cells") +
  geom_abline(slope = 1, intercept = 0)

ggsave(paste0(path2data,"QC/UMIS_expectedVSobservedEmbryoTube.pdf"))

ggplot(exp, aes(x = cells, y = obs, col = tube_name)) +
  geom_point() +
  #scale_colour_Publication(name = "Sample") +
  scale_x_log10(breaks = c(500, 1000, 5000)) +
  scale_y_log10(breaks = c(500, 1000, 5000)) +
  labs(x = "Predicted number of cells", y = "Observed number of cells") +
  geom_abline(slope = 1, intercept = 0)

ggsave(paste0(path2data, "QC/UMIS_expectedVSobservedTubeName.pdf"))

## Write output

writeMM(counts, file = paste0(path2data, "raw_counts_qc.mtx"))

saveRDS(as(counts, "dgCMatrix"), file = paste0(path2data, "raw_counts_qc.rds"))

write.table(barcodes, file = paste0(path2data, "barcodes_qc.tsv"),
  row.names = F, col.names = F, quote = F, sep = "\t")

write.table(meta, file = paste0(path2data, "meta.tab"), row.names = F, col.names = T, quote = F, sep = "\t")

