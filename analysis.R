# Load required libraries
library(limma)
library(edgeR)
library(tidyverse)
library(openxlsx)
library(pheatmap)
library(data.table)
library(org.Hs.eg.db)
library(Homo.sapiens)
library(RColorBrewer)
library(ggrepel)
library(hues)
library(Biobase)
library(stringr)
library(ggfortify)
library(qs)
library(DT)
library(dplyr)
library(NOISeq)
library(AnnotationDbi)
library(org.Mm.eg.db)
source("HRK_funcs.R")
select <- dplyr::select
annotation_obj <- get(report_params$annotation_db, envir = asNamespace(report_params$annotation_db)) 

# Ensure output directory exists
if (!dir.exists(report_params$out_dir)) dir.create(report_params$out_dir)

#Step 1: Read in the tab-separated files
files.list <- list.files(path = report_params$rsem_dir, pattern = 'genes.results$', full.names = T)

#Gets just the base sample name
orig_names <- gsub(".genes.results","",basename(files.list))

list_of_data <- lapply(files.list, function(file) {
  read.table(file, header = TRUE, sep = "\t", row.names = 1)
})

names(list_of_data) <- orig_names

DGE <- readDGE(files = files.list, columns = c(1, 5))

# Rename samples
orig_names <- str_sub(string = basename(colnames(DGE)), end = -7)
orig_names <- str_sub(string = orig_names, start = 1, end = 18)
rownames(DGE$samples) <- orig_names
colnames(DGE$counts) <- orig_names
DGE$samples$SampleName <- orig_names

# Load sample keys
sample.keys <- read_csv(report_params$sample_data)

# Merge sample metadata
DGE$samples <- merge(DGE$samples, sample.keys, by="SampleName")
rownames(DGE$samples) <- orig_names

# Annotate genes
ensemblid <- rownames(DGE)

# Get gene annotations
genes <- AnnotationDbi::select(annotation_obj, 
                               keys = ensemblid, 
                               columns = c("ENSEMBL", "SYMBOL"), 
                               keytype = "ENSEMBL")

# Remove duplicated Ensembl IDs (keeping the first occurrence)
genes <- genes[!duplicated(genes$ENSEMBL), ]

# Ensure we retain all IDs and row order in DGE
matched_genes <- genes[match(ensemblid, genes$ENSEMBL), ]

# Assign the matched gene annotations as a new column in DGE
DGE$genes <- matched_genes

# Define treatments
treatment.all <- as.factor(DGE$samples[[report_params$group_var]])

# Filter genes based on expression level (edgeR)
keep.exprs <- filterByExpr(DGE, group = treatment.all)
DGE.edgeRfilt <- DGE[keep.exprs, , keep.lib.sizes = FALSE] 

# Normalize using TMM method (edgeR)
DGE.edgeRfilt <- calcNormFactors(DGE.edgeRfilt, method = 'TMM')

# Filter genes based on NOISeq
DGE.NOIseqfilt <- DGE
DGE.NOIseqfilt$counts <- filtered.data(
  DGE$counts, 
  factor = DGE$samples[[report_params$group_var]], 
  norm = FALSE, 
  depth = NULL, 
  method = 1, 
  cv.cutoff = 100, 
  cpm = 1, 
  p.adj = "fdr"
)

DGE.NOIseqfilt$samples$lib.size <- apply(DGE.NOIseqfilt$counts, 2, sum)
DGE.NOIseqfilt <- calcNormFactors(DGE.NOIseqfilt, method = 'TMM')

# Create design matriDGE
design.mat <- model.matrix(~ 0 + DGE$samples[[report_params$group_var]])
rownames(design.mat) <- DGE$samples$SampleName
colnames(design.mat) <- make.names(levels(as.factor(DGE$samples[[report_params$group_var]])))

# Perform voom transformation
# Create voom plot and save voom object
png(filename = paste0(report_params$out_dir, "/voom_plot.png"))  # Open PNG device
v <- voom(DGE.NOIseqfilt, design.mat, plot = TRUE)
dev.off()  # Close PNG device

# Fit linear model
vfit <- lmFit(v, design.mat)

# Extract unique contrast strings and split multiple contrasts in each entry
contrast_strings <- unique(unlist(strsplit(DGE.NOIseqfilt$samples$Contrast, ";")))

# Create a named list of contrasts dynamically
contrasts_list <- setNames(
  lapply(contrast_strings, function(contrast) {
    groups <- unlist(strsplit(contrast, "-"))  # Split into individual groups
    makeContrasts(contrasts = paste0(groups[1], " - ", groups[2]), 
                  levels = colnames(coef(vfit)))  # Use levels from vfit
  }),
  paste0(gsub("-", "_vs_", contrast_strings))  # Name informatively
)

# Print the list to check
print(contrasts_list)


# Compute contrasts and apply eBayes
efit_list <- lapply(contrasts_list, function(contr) {
  vfit_contr <- contrasts.fit(vfit, contr)
  eBayes(vfit_contr)
})

# Save differential expression results

library(clusterProfiler)

# Create result data frames dynamically
efit_results_list <- setNames(
  lapply(seq_along(efit_list), function(i) {
    create_efit_results_df(efit_list[[i]])
  }),
  paste0("efit_", names(efit_list), "_results_df") # Informative names
)

# Function to extract top DE genes and convert to Entrez IDs
top_DE_entrezIDs <- function(df, direction = "up") {
  if (direction == "up") {
    filtered_genes <- df %>% filter(logFC > 0.58 & adj.P.value < 0.05)
  } else if (direction == "down") {
    filtered_genes <- df %>% filter(logFC < -0.58 & adj.P.value < 0.05)
  } else {
    stop("Invalid direction. Choose either 'up' or 'down'.")
  }
  
  # Return empty vector if no significant genes
  if (nrow(filtered_genes) == 0) {
    return(character(0))  # Empty list instead of failure
  }
  
  entrez_ids <- mapIds(
    annotation_obj,
    keys = filtered_genes$ensembleID,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  return(entrez_ids)
}

# Apply function dynamically for each contrast and store results
top_DE_entrezIDs_list <- setNames(
  lapply(efit_results_list, function(df) {
    list(
      up = top_DE_entrezIDs(df, "up"),
      down = top_DE_entrezIDs(df, "down")
    )
  }),
  paste0("top_DE_entrezIDs_", names(efit_results_list))
)

universe_entrez <- mapIds(
  annotation_obj,
  keys = rownames(DGE.NOIseqfilt),  # Ensembl universe
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

dge_list_edgeR <- DGE.edgeRfilt
dge_list_raw <- DGE
dge_list_NOISeq <- DGE.NOIseqfilt
efit_list <- efit_list
efit_results_dfs <- efit_results_list
entrez_ids_list <- top_DE_entrezIDs_list


