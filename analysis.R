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

# Extract unique contrast names
contrast_strings <- unique(unlist(strsplit(DGE.NOIseqfilt$samples$Contrast, ";")))

# Generate contrast list
contrasts_list <- setNames(
  lapply(contrast_strings, function(contrast) {
    groups <- trimws(unlist(strsplit(contrast, "-")))  # Remove whitespace
    
    print(paste("Processing contrast:", paste(groups, collapse = " - ")))  # Debugging step
    
    # Directly pass the contrast formula (no eval/parse needed!)
    makeContrasts(contrasts = paste(groups[1], "-", groups[2]), 
                  levels = colnames(vfit$design))
  }),
  paste0(gsub("-", "_vs_", contrast_strings))  # Clean names for list elements
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
# Create a flat list of DE gene lists with contrast-specific names
top_DE_entrezIDs_list <- do.call(c, lapply(c("up", "down"), function(direction) {  
  unlist(lapply(names(efit_results_list), function(contrast_name) {
    df <- efit_results_list[[contrast_name]]
    contrast_name_clean <- gsub("^efit_|_results_df$", "", contrast_name)
    
    # Extract genes based on the current direction
    gene_list <- top_DE_entrezIDs(df, direction)
    list(gene_list)  # Store in a list
  }), recursive = FALSE)
}))

# Assign names in the correct order
names(top_DE_entrezIDs_list) <- unlist(lapply(c("up", "down"), function(direction) {  
  lapply(names(efit_results_list), function(contrast_name) {
    contrast_name_clean <- gsub("^efit_|_results_df$", "", contrast_name)
    paste0(contrast_name_clean, ".DE.", direction)  # Create meaningful names
  })
}))


# Create a flat list of non-DE gene lists
non_DE_entrezIDs_list <- setNames(
  lapply(names(efit_results_list), function(contrast_name) {
    df <- efit_results_list[[contrast_name]]
    non_DE_entrezIDs(df)  # Should return a character vector
  }),
  gsub("^efit_|_results_df$", "", names(efit_results_list)) %>%
    paste0(".nonDE")
)

#Extract and format gene lists specified in report_params$gene_lists

# Identify all gene lists available in the session
all_entrez_lists <- ls(pattern = "entrez")
entrez_lists <- all_entrez_lists[sapply(all_entrez_lists, function(x) is.list(get(x)))]

# Extract user-specified gene lists from report_params
specified_gene_lists <- eval(parse(text = report_params$gene_lists))

# Create a structured named list for gene sets, keeping groups separate
gene_lists <- list()

for (pattern_vec in specified_gene_lists) {
  pattern <- pattern_vec[[1]]  # Extract pattern from list
  regex_pattern <- gsub("\\*", ".*", pattern)  # Convert wildcard to regex
  
  # Find matching objects
  matching_objects <- entrez_lists[grepl(regex_pattern, entrez_lists)]
  
  # Store matching objects in a **sub-list** within gene_lists
  gene_lists[[pattern]] <- list()  
  
  if (length(matching_objects) > 0) {
    for (list_name in matching_objects) {
      contrast_names <- names(get(list_name))
      
      for (contrast in contrast_names) {
        gene_lists[[pattern]][[contrast]] <- get(list_name)[[contrast]]
      }
    }
  }
}

# Save universe for enrichment funcs
universe_entrez <- mapIds(
  annotation_obj,
  keys = rownames(DGE.NOIseqfilt),  # Ensembl universe
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Rename some objs to work with exisiting .Rmd code
dge_list_edgeR <- DGE.edgeRfilt
dge_list_raw <- DGE
dge_list_NOISeq <- DGE.NOIseqfilt
efit_list <- efit_list
efit_results_dfs <- efit_results_list
entrez_ids_list <- top_DE_entrezIDs_list


