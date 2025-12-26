#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
  library(readr)
  library(plyr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(scCustomize)
  library(cowplot)
  library(optparse) # For parsing command line arguments
})

# --- 1. Argument Parsing ---
option_list = list(
  make_option(c("--run_name"), type="character", default=NULL, help="Nextflow run name (e.g., A1)"),
  make_option(c("--matrix_h5"), type="character", default=NULL, help="Path to the filtered_feature_bc_matrix.h5 file"),
  make_option(c("--rds"), type="character", default=NULL, help="Path to a Seurat object .rds file"),
  make_option(c("--demux_csv"), type="character", default=NULL, help="Path to the demuxalot posterior probabilities CSV"),
  make_option(c("--sample_names"), type="character", default=NULL, help="Comma-separated list of patient SM IDs (e.g., 'P1','P2')"),
  make_option(c("--sexes"), type="character", default=NULL, help="Comma-separated list of patient sexes (e.g., 'Male','Female')")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$run_name) || is.null(opt$demux_csv)) {
  stop("Run name and demux CSV path must be provided.", call.=FALSE)
}

if (is.null(opt$matrix_h5) && is.null(opt$rds)) {
  stop("Either --matrix_h5 or --rds must be provided.", call.=FALSE)
}

# Clean and format arguments
experiment1 <- opt$run_name
demuxalot_path <- opt$demux_csv
sample_names_vec <- gsub("'", "", strsplit(opt$sample_names, ",")[[1]])
sexes_vec <- gsub("'", "", strsplit(opt$sexes, ",")[[1]])

# Create temporary sex dataframe
sex <- data.frame(
  patient = sample_names_vec,
  sex = sexes_vec
)

# --- 2. Load Data ---
cat(paste0("Processing run: ", experiment1, "\n"))
cat(paste0("Loading Demuxalot from: ", demuxalot_path, "\n"))

if (!file.exists(demuxalot_path)) {
  stop(paste("Error: Demuxalot CSV not found at", demuxalot_path))
}
demuxalot <- read.csv(demuxalot_path)

# Load Seurat object
if (!is.null(opt$rds)) {
  cat(paste0("Loading Seurat object from RDS: ", opt$rds, "\n"))
  obj <- readRDS(opt$rds)
} else if (!is.null(opt$matrix_h5)) {
  cat(paste0("Loading matrix from H5: ", opt$matrix_h5, "\n"))
  if (!file.exists(opt$matrix_h5)) stop("H5 file not found!")
  h5 <- Read10X_h5(opt$matrix_h5)
  obj <- CreateSeuratObject(h5$`Gene Expression`)
} else {
  stop("No valid input provided.")
}

meta <- obj@meta.data
meta$cell_id <- rownames(meta)

# --- 3. Demuxalot Assignment ---
demuxalot$assignment <- NA
# for(i in 1:nrow(demuxalot)){
#   idx <- which.max(c(demuxalot[i, 2:ncol(demuxalot)]))[[1]] 
#   demuxalot$assignment[i] <- colnames(demuxalot)[1 + idx]
# }

# Find the column index with the max value for each row (excluding BARCODE)
max_indices <- apply(demuxalot[, 2:ncol(demuxalot)], 1, which.max)
# Map the index back to the column name (1 + index for 1-based indexing)
demuxalot$assignment <- colnames(demuxalot)[1 + max_indices]

demuxalot$assignment <- sub("X","",demuxalot$assignment)
demuxalot$assignment <- gsub(".","-",demuxalot$assignment, fixed = T)

cat("Assignment Table:\n")
print(table(demuxalot$assignment))

# --- 4. Merge Demultiplexing ---
demuxalot_1 <- subset(demuxalot, select = c("BARCODE", "assignment"))
order <- rownames(meta)
meta <- merge(meta, demuxalot_1, by.x = "cell_id", by.y = "BARCODE", all.x = TRUE)
meta <- meta %>% slice(match(order, cell_id))
rownames(meta) <- meta$cell_id
obj@meta.data <- meta

# --- 5. PCA ---
rownames(demuxalot) <- demuxalot$BARCODE
clusters_pca <- subset(demuxalot, select = -c(BARCODE, assignment))
clusters_pca <- log2(abs(clusters_pca) + 1)
clusters_pca <- prcomp(clusters_pca, scale = TRUE, center = TRUE)

PCi <- data.frame(clusters_pca$x)
PCi$cell_id <- rownames(PCi)
PCi <- merge(PCi, meta, by = "cell_id", all.x = TRUE)
PCi <- merge(PCi, sex, by.x = "assignment", by.y = "patient", all.x = TRUE)

# --- 6. PCA Visualization ---
pdf(file = paste0(experiment1, "_PCA_plots.pdf"), width = 10, height = 12)

p1 <- ggplot(PCi, aes(x = PC1, y = PC2, col = sex)) +
  geom_point(size = 1, alpha = 0.7) +
  theme_classic() +
  scale_color_manual(values = c("#F8A19FFF", "#325A9BFF", "#27e627d4")) +
  ggtitle(experiment1)

p2 <- ggplot(PCi, aes(x = PC1, y = PC2, col = assignment)) +
  geom_point(size = 1, alpha = 0.5) +
  theme_classic() +
  scale_color_manual(values = DiscretePalette_scCustomize(8, "glasbey", shuffle_pal = TRUE)) +
  ggtitle(experiment1)

print(plot_grid(p1, p2, nrow = 1))
dev.off()

# --- 7. Sex Gene Expression ---
pdf(file = paste0(experiment1, "_SexGene_plots.pdf"), width = 12, height = 18)
obj <- NormalizeData(obj)
obj <- ScaleData(obj)
genes <- c("XIST", "DDX3Y", "RPS4Y1", "USP9Y", "UTY", "ZFY")
gene <- FetchData(obj, vars = genes)
gene$cell_id <- rownames(gene)
names(gene) <- gsub("-", "_", names(gene))
genes <- gsub("-", "_", genes)
PCi_plot <- merge(PCi, gene, by = "cell_id", all.y = TRUE)

create_gene_plots <- function(PCx, PCy, genes){
  plots <- list()
  for(i in 1:length(genes)){
    plots[[i]] <- ggplot(PCi_plot, aes_string(x = PCx, y = PCy)) +
      geom_point(aes_string(col = genes[i]), alpha = 0.7, size = 1) +
      theme_classic() +
      scale_colour_gradient2(low = "orange", mid = "grey", high = "blue", guide = "colourbar")
  }
  return(plots)
}

print(plot_grid(plotlist = create_gene_plots("PC1", "PC2", genes), nrow = 3))
print(plot_grid(plotlist = create_gene_plots("PC1", "PC3", genes), nrow = 3))
print(plot_grid(plotlist = create_gene_plots("PC2", "PC3", genes), nrow = 3))
dev.off()
