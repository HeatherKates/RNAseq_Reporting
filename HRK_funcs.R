library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(plotly)
library(heatmaply)
library(org.Mm.eg.db)

generate_enrichment_plot <- function(gene_lists, de_results_df, universe_entrez, ont_category, significance_threshold = 0.05, top_n = 10,annotation_db) {
  annotation_obj <- get(annotation_db, envir = asNamespace(annotation_db)) 
  
  # Ensure gene lists are named and define contrast order
  if (is.null(names(gene_lists))) stop("Each gene list must be named!")
  contrast_order <- names(gene_lists)
  
  # Convert DE results dataframe to a lookup table for mapping Entrez IDs
  de_results_df <- de_results_df %>%
    mutate(ENTREZID = mapIds(annotation_obj, keys = ensembleID, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first"))
  
  # Run GO enrichment using compareCluster()
  formula_res <- compareCluster(
    Entrez ~ Contrast,
    data = bind_rows(lapply(contrast_order, function(contrast) {
      data.frame(
        Entrez = gene_lists[[contrast]],
        Contrast = contrast
      )
    })),
    fun = "enrichGO",
    universe = na.omit(universe_entrez),
    OrgDb = annotation_obj,
    keyType = "ENTREZID",
    ont = ont_category,
    pvalueCutoff = 0.05
  )
  
  # Ensure clusters are ordered correctly
  formula_res@compareClusterResult$Cluster <- factor(
    formula_res@compareClusterResult$Cluster,
    levels = contrast_order
  )
  
  # **Filter only by significance threshold FIRST (to preserve all significant results)**
  filtered_results <- subset(
    formula_res@compareClusterResult,
    p.adjust <= significance_threshold
  )
  
  # Convert Entrez IDs to Gene Symbols and select top 20 by p-value **within the correct cluster & direction**
    filtered_results$GeneSymbols <- sapply(filtered_results$geneID, function(gene_list) {
      # Split "/" separated Entrez IDs into a vector
      entrez_ids <- unlist(strsplit(gene_list, "/"))
      
      # Convert to Gene Symbols
      gene_symbols <- mapIds(annotation_obj, 
                             keys = entrez_ids, 
                             column = "SYMBOL", 
                             keytype = "ENTREZID", 
                             multiVals = "first")
      
      # Remove NAs and format output
      gene_symbols <- na.omit(gene_symbols)
      
      if (length(gene_symbols) > 0) {
        paste(gene_symbols, collapse = "<br>")  # Return as formatted string
      } else {
        NA_character_  # Return NA if empty
      }
    })
    
  
  # **Save full filtered results for downloading**
  download_go_results <- filtered_results %>%
    select(Cluster, Description, p.adjust, GeneSymbols, everything())  # Reorder columns
  
  # **Identify top `n` GO terms across all clusters for plotting**
  top_GO_terms <- filtered_results %>%
    group_by(Cluster) %>%
    slice_min(order_by = p.adjust, n = top_n) %>%
    pull(Description) %>%
    unique()
  
  # **Filter for plotting (only keeping top GO terms)**
  formula_res@compareClusterResult <- filtered_results %>%
    filter(Description %in% top_GO_terms)
  
  # Convert GeneRatio to numeric
  formula_res@compareClusterResult <- formula_res@compareClusterResult %>%
    mutate(GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
  
  # **Restore GO Term Ordering Using Hierarchical Clustering**
 
  reorder_GO_terms <- function(df) {
    term_matrix <- table(df$Description, df$Cluster)  
    
    if (nrow(term_matrix) > 1 && length(unique(df$Description)) > 1) {
      term_dist <- dist(term_matrix, method = "binary")  
      term_hclust <- hclust(term_dist, method = "ward.D2")  
      term_order <- rownames(term_matrix)[term_hclust$order]
    } else {
      term_order <- unique(df$Description)  # Ensure a valid factor level
    }
    
    # Safely convert to a factor with the correct order
    df$Description <- factor(df$Description, levels = term_order)  
    return(df)
  }
  
  formula_res@compareClusterResult <- reorder_GO_terms(formula_res@compareClusterResult)
  
  # **Define GeneRatio bins**
  bin_breaks <- c(0, 0.01, 0.05, 0.10, max(formula_res@compareClusterResult$GeneRatio, na.rm = TRUE) + 0.01)  
  bin_labels <- c("≤0.01", "0.01 - 0.05", "0.05 - 0.10", "≥0.10")
  
  formula_res@compareClusterResult <- formula_res@compareClusterResult %>%
    mutate(GeneRatioCategory = cut(GeneRatio, breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE))
  
  # **Define reference dot sizes**
  size_mapping <- c("≤0.01" = 2, "0.01 - 0.05" = 4, "0.05 - 0.10" = 6, "≥0.10" = 8)
  
  # **Ensure all necessary variables are characters/numeric**
  formula_res@compareClusterResult <- formula_res@compareClusterResult %>%
    mutate(
      p.adjust = as.numeric(as.character(p.adjust)),
      GeneRatio = as.numeric(as.character(GeneRatio)),
      tooltip_text = paste(
        "Cluster: ", as.character(Cluster), "<br>",  
        "GO Term: ", as.character(Description), "<br>",  
        "p.adjust: ", signif(p.adjust, 3), "<br>",
        "GeneRatio: ", signif(GeneRatio, 3), "<br>",
        "Top Genes:<br>", as.character(GeneSymbols)  
      )
    )
  
  # **Create ggplot object**
  p <- ggplot(formula_res@compareClusterResult, aes(
    x = Cluster, 
    y = Description, 
    size = GeneRatioCategory,  
    color = p.adjust,
    text = tooltip_text  
  )) +
    geom_point(alpha = 0.8) +  
    scale_size_manual(name = "Gene Ratio", values = size_mapping) +  
    scale_color_gradient(low = "red", high = "blue", 
                         limits = c(min(formula_res@compareClusterResult$p.adjust, na.rm = TRUE), 
                                    max(formula_res@compareClusterResult$p.adjust, na.rm = TRUE)),
                         name = "p.adjust") +  
    guides(color = guide_colorbar(title = "p.adjust")) +  
    ggtitle(paste("GO Term Enrichment (", ont_category, ")", sep = "")) +
    xlab("DE Gene list (p.adj < 0.05, abs(log2FC) > 0.58)") +
    ylab("GO Term") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # **Convert to interactive plot (ensure tooltips work)**
  interactive_plot <- ggplotly(p, tooltip = "text") %>%  
    layout(legend = list(title = list(text = "Gene Ratio")))
  
  # **Static High-Resolution Plot (for manuscript)**
  static_plot <- dotplot(formula_res, showCategory = top_n) +
    ggtitle(paste("GO Term Enrichment (", ont_category, ")", sep = "")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # **Return interactive plot, static plot, and full GO results**
  return(list(
    interactive_plot = interactive_plot, 
    static_plot = static_plot, 
    go_results = download_go_results
  ))
}

# Run the function
#GO_BP_results <- generate_enrichment_plot(
#  gene_lists = gene_lists, 
#  de_results_df = bind_rows(efit.1.results_df, efit.2.results_df, efit.3.results_df), 
#  universe_entrez = universe_entrez,
#  ont_category = "BP"
#)
# Run the function
#GO_MF_results <- generate_enrichment_plot(
#  gene_lists = gene_lists, 
#  de_results_df = bind_rows(efit.1.results_df, efit.2.results_df, efit.3.results_df), 
#  universe_entrez = universe_entrez,
#  ont_category = "MF"
#)

# Save high-resolution PNG for manuscript
#ggsave("GO_Enrichment_BP.png", plots$static_plot, width = 10, height = 12, dpi = 300,bg = "white")

generate_kegg_enrichment_plot <- function(gene_lists, de_results_df, universe_entrez, significance_threshold = 0.05, top_n = 10,annotation_db) {
  annotation_obj <- get(annotation_db, envir = asNamespace(annotation_db)) 
  
  # Ensure gene lists are named and define contrast order
  if (is.null(names(gene_lists))) stop("Each gene list must be named!")
  contrast_order <- names(gene_lists)
  
  # Convert DE results dataframe to a lookup table for mapping Entrez IDs
  de_results_df <- de_results_df %>%
    mutate(ENTREZID = mapIds(annotation_obj, keys = ensembleID, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first"))
  
  # Run KEGG enrichment using compareCluster()
  kegg_res <- compareCluster(
    Entrez ~ Contrast,
    data = bind_rows(lapply(contrast_order, function(contrast) {
      data.frame(
        Entrez = gene_lists[[contrast]],
        Contrast = contrast
      )
    })),
    fun = "enrichKEGG",
    universe = na.omit(universe_entrez),
    organism = report_params$organism,
    keyType = "kegg",
    pvalueCutoff = 0.05
  )
  
  # Ensure clusters are ordered correctly
  kegg_res@compareClusterResult$Cluster <- factor(
    kegg_res@compareClusterResult$Cluster,
    levels = contrast_order
  )
  
  # **Filter only by significance threshold FIRST (to preserve all significant results)**
  filtered_results <- subset(
    kegg_res@compareClusterResult,
    p.adjust <= significance_threshold
  )
  
  # Convert Entrez IDs to Gene Symbols and select top 20 by p-value **within the correct cluster & direction**
  # Convert Entrez IDs to Gene Symbols and select top 20 by p-value **within the correct cluster & direction**
  filtered_results$GeneSymbols <- sapply(filtered_results$geneID, function(gene_list) {
    # Split "/" separated Entrez IDs into a vector
    entrez_ids <- unlist(strsplit(gene_list, "/"))
    
    # Convert to Gene Symbols
    gene_symbols <- mapIds(annotation_obj, 
                           keys = entrez_ids, 
                           column = "SYMBOL", 
                           keytype = "ENTREZID", 
                           multiVals = "first")
    
    # Remove NAs and format output
    gene_symbols <- na.omit(gene_symbols)
    
    if (length(gene_symbols) > 0) {
      paste(gene_symbols, collapse = "<br>")  # Return as formatted string
    } else {
      NA_character_  # Return NA if empty
    }
  })
  
  
  # **Save full filtered results for downloading**
  download_kegg_results <- filtered_results %>%
    select(Cluster, Description, p.adjust, GeneSymbols, everything())  # Reorder columns
  
  # **Identify top `n` KEGG pathways across all clusters for plotting**
  top_KEGG_terms <- filtered_results %>%
    group_by(Cluster) %>%
    slice_min(order_by = p.adjust, n = top_n) %>%
    pull(Description) %>%
    unique()
  
  # **Filter for plotting (only keeping top KEGG pathways)**
  kegg_res@compareClusterResult <- filtered_results %>%
    filter(Description %in% top_KEGG_terms)
  
  # Convert GeneRatio to numeric
  kegg_res@compareClusterResult <- kegg_res@compareClusterResult %>%
    mutate(GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
  
  # **Restore KEGG Pathway Ordering Using Hierarchical Clustering**
  reorder_KEGG_terms <- function(df) {
    term_matrix <- table(df$Description, df$Cluster)  
    
    if (nrow(term_matrix) > 1) {
      term_dist <- dist(term_matrix, method = "binary")  
      term_hclust <- hclust(term_dist, method = "ward.D2")  
      term_order <- rownames(term_matrix)[term_hclust$order]
    } else {
      term_order <- df$Description  
    }
    
    df$Description <- factor(df$Description, levels = term_order)  
    return(df)
  }
  
  kegg_res@compareClusterResult <- reorder_KEGG_terms(kegg_res@compareClusterResult)
  
  # **Define GeneRatio bins**
  bin_breaks <- c(0, 0.01, 0.05, 0.10, max(kegg_res@compareClusterResult$GeneRatio, na.rm = TRUE) + 0.01)  
  bin_labels <- c("≤0.01", "0.01 - 0.05", "0.05 - 0.10", "≥0.10")
  
  kegg_res@compareClusterResult <- kegg_res@compareClusterResult %>%
    mutate(GeneRatioCategory = cut(GeneRatio, breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE))
  
  # **Define reference dot sizes**
  size_mapping <- c("≤0.01" = 2, "0.01 - 0.05" = 4, "0.05 - 0.10" = 6, "≥0.10" = 8)
  
  # **Ensure all necessary variables are characters/numeric**
  kegg_res@compareClusterResult <- kegg_res@compareClusterResult %>%
    mutate(
      p.adjust = as.numeric(as.character(p.adjust)),
      GeneRatio = as.numeric(as.character(GeneRatio)),
      tooltip_text = paste(
        "Cluster: ", as.character(Cluster), "<br>",  
        "KEGG Pathway: ", as.character(Description), "<br>",  
        "p.adjust: ", signif(p.adjust, 3), "<br>",
        "GeneRatio: ", signif(GeneRatio, 3), "<br>",
        "Top Genes:<br>", as.character(GeneSymbols)  
      )
    )
  
  # **Create ggplot object**
  p <- ggplot(kegg_res@compareClusterResult, aes(
    x = Cluster, 
    y = Description, 
    size = GeneRatioCategory,  
    color = p.adjust,
    text = tooltip_text  
  )) +
    geom_point(alpha = 0.8) +  
    scale_size_manual(name = "Gene Ratio", values = size_mapping) +  
    scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +  
    guides(color = guide_colorbar(title = "p.adjust")) +  
    ggtitle("KEGG Pathway Enrichment") +
    xlab("DE Gene list (p.adj < 0.05, abs(log2FC) > 0.58)") +
    ylab("KEGG Pathway") +
    theme_minimal()+
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) 
  
  interactive_plot <- ggplotly(p, tooltip = "text") %>%  
    layout(legend = list(
      title = list(text = "Gene Ratio"),  # Ensure Gene Ratio is labeled
      colorbar = list(title = "p.adjust") # Ensure p.adjust is labeled for color
    ))
  
  
  
  # **Static KEGG Plot**
  static_plot <- dotplot(kegg_res, showCategory = top_n) +
    ggtitle("KEGG Pathway Enrichment") +
    theme_minimal()
  
  return(list(interactive_plot = interactive_plot, static_plot = static_plot, kegg_results = download_kegg_results))
}

# Run the function
#kegg_results <- generate_kegg_enrichment_plot(
#  gene_lists = gene_lists, 
#  de_results_df = bind_rows(efit.1.results_df, efit.2.results_df, efit.3.results_df), 
#  universe_entrez = universe_entrez
#)

#kegg_results$interactive_plot

generate_volcano_plot <- function(efit_results_df, contrast_name) {
  library(plotly)
  
  volcano.res <- efit_results_df
  volcano.res$group <- "NotSignificant"
  volcano.res[which(volcano.res$adj.P.value < 0.05 & volcano.res['logFC'] > 0.58), "group"] <- "Upregulated"
  volcano.res[which(volcano.res$adj.P.value < 0.05 & volcano.res['logFC'] < -0.58), "group"] <- "Downregulated"
  
  # Apply small offset to avoid log(1) = 0 issues
  volcano.res$neglog10p <- -log10(volcano.res$adj.P.value + 1e-10)  # Ensure no log(1)
  # Check if all y-values are small (<1) and set range accordingly
  if (max(volcano.res$neglog10p, na.rm = TRUE) < 1) {
    y_axis_range <- c(0, 1)  # Force range to 0-1
  } else {
    y_axis_range <- c(0, max(volcano.res$neglog10p, na.rm = TRUE) + 0.5)  # Dynamic range otherwise
  }
  plot_ly(
    data = volcano.res,
    x = volcano.res$logFC, 
    y = volcano.res$neglog10p,  # Use modified neg log10 p-values
    hoverinfo = "text", 
    text = paste(
      "Gene: ", volcano.res$SYMBOL, 
      "<br>Ensemble: ",volcano.res$ensembleID,
      "<br>logFC: ", round(volcano.res$logFC, 2),
      "<br>adj.P.val: ", round(volcano.res$adj.P.value, 2)
    ),
    mode = "markers", 
    color = volcano.res$group,
    colors = c(NotSignificant = "gray", Upregulated = "red", Downregulated = "blue")
  ) %>% 
    layout(
      title = paste(contrast_name, "(Adj. P-value vs. Fold-Change)"),
      xaxis = list(title = 'logFC'), 
      yaxis = list(
        title = '-log10(adj.PVal)',
        exponentformat = "none",  # Disable scientific notation
        separatethousands = TRUE,  # Keep numbers readable
        tickmode = "linear",  # Ensure even tick spacing
        range = y_axis_range 
      ),
      margin = list(l = 75, t = 150)
    )
  
  
}

generate_heatmap <- function(efit_results_df, lcpm_matrix, dge_list_NOISeq, title, num_genes = 50, fontsize_row = 12) {
  
  # Extract contrast name dynamically
  contrast_name <- unique(efit_results_df$contrast)[1]  # Ensure it's a single contrast
  contrast_groups <- unlist(strsplit(contrast_name, " - "))  # Extract two groups being compared
  
  # Filter samples that belong to these groups
  selected_samples <- dge_list_NOISeq$samples %>%
    filter(!!sym(report_params$group_var) %in% contrast_groups) %>%  # Select samples in the contrast
    select(SampleName, all_of(report_params$group_var))
  
  # Step 1: Select top 50 DE genes
  top50_sigOE_genes <- efit_results_df %>%
    arrange(adj.P.value) %>%
    slice_head(n = num_genes) %>%
    select(ensembleID, SYMBOL)
  
  # Step 2: Subset the lcpm matrix for the selected samples
  lcpm_filtered <- lcpm_matrix[, colnames(lcpm_matrix) %in% selected_samples$SampleName]
  
  # Step 3: Ensure genes are properly filtered
  lcpm_filtered <- lcpm_filtered[rownames(lcpm_filtered) %in% top50_sigOE_genes$ensembleID, ]
  
  # Step 4: Transpose `lcpm_filtered`
  lcpm_filtered <- as.data.frame(t(lcpm_filtered))
  
  # Step 5: Ensure sample metadata matches
  sidecols <- selected_samples %>% select(SampleName, all_of(report_params$group_var))
  rownames(sidecols) <- sidecols$SampleName
  
  # Step 6: Reorder `lcpm_filtered` based on sorted `sidecols`
  sidecols <- sidecols %>% arrange(all_of(report_params$group_var))  
  lcpm_filtered <- t(lcpm_filtered)
  
  # Ensure order matches
  matching_indices <- match(rownames(sidecols), colnames(lcpm_filtered))  
  lcpm_filtered <- lcpm_filtered[, matching_indices, drop = FALSE]  
  rownames(lcpm_filtered) <- top50_sigOE_genes$SYMBOL  
  
  # Step 7: Define dynamic colors for the annotation
  group_levels <- unique(selected_samples[[report_params$group_var]])  
  color_palette <- brewer.pal(min(length(group_levels), 8), "Set2")  
  col_side_palette <- colorRampPalette(color_palette)  # Create a color palette function
  
  # Step 8: Create heatmap
  hm <- heatmaply(
    plot_method = "plotly",
    lcpm_filtered,
    show_dendrogram = FALSE,
    showticklabels = c(FALSE, TRUE),
    scale = "row",
    colors = c("darkblue", "white", "darkred"),
    col_side_colors = sidecols %>% select(all_of(report_params$group_var)),  
    xlab = "Samples",
    ylab = "Genes",
    main = title,
    margins = c(0, 0, 60, 0),
    fontsize_row = fontsize_row,
    cex.main = 1,
    col_side_palette = col_side_palette,  # Corrected palette function
    labRow = rownames(lcpm_filtered),
    side_color_colorbar_len = 0.2,   
    side_color_width = 0.05,
    height=700
  )
  
  return(hm)
}


#source("rnaseq_functions.R")
library(gtable)

# Function to generate MDS plots for multiple DGELists
generate_mds_plots <- function(dge_list, titles, top = 1000, group = "treatment") {
  if (length(dge_list) != length(titles)) {
    stop("Error: Length of DGEList and titles must be the same")
  }
  
  plots <- lapply(seq_along(dge_list), function(i) {
    dge <- dge_list[[i]]
    title <- titles[i]
    
    plot_MDS_2(lcpm = cpm(dge, log = TRUE), dge = dge, top = top, group = group) +
      theme(legend.position = "none") + 
      geom_point(size = 3) +
      ggtitle(title)
  })
  
  return(plots)
}

# Function to extract and adjust legend
get_legend <- function(p) {
  g <- ggplotGrob(
    p + theme(
      legend.key.width = unit(2, "cm"),  # Increase width of legend keys
      legend.spacing.x = unit(0.5, "cm"), # Increase spacing between legend items
      legend.text = element_text(size = 12),  # Adjust text size
      legend.title = element_text(size = 14, face = "bold") # Adjust title size
    )
  )
  legend <- gtable_filter(g, "guide-box")
  return(legend)
}

library(ggplot2)
library(base64enc)
library(htmltools)

download_button_png <- function(plot_object, output_name = "plot", width = 12, height = 6, dpi = 300) {
  # Step 1: Save the ggplot object as a temporary PNG file
  temp_file <- tempfile(fileext = ".png")
  ggsave(temp_file, plot = plot_object, width = width, height = height, dpi = dpi)
  
  # Step 2: Convert the image to Base64 encoding
  encoded_img <- base64enc::dataURI(file = temp_file, mime = "image/png")
  
  # Step 3: Create a styled HTML download button
  download_button <- HTML(paste0(
    '<a href="', encoded_img, '" download="', output_name, '.png" ',
    'class="btn btn-primary" style="padding:10px; font-size:16px; text-decoration:none; ',
    'color:white; background-color:#007BFF; border-radius:5px;">',
    '📥 Download ', output_name, '</a>'
  ))
  
  # Return the HTML download button
  return(download_button)
}

library(DT)
library(dplyr)

display_de_result_table <- function(efit, contrast_name = "Contrast") {
  # Step 1: Extract relevant columns from the efit object
  efit_results_df <- data.frame(
    ensembleID = rownames(efit$coefficients),  # Extract rownames as ensemble IDs
    logFC = efit$coefficients[, 1],           # Log fold change
    AveExpr = efit$Amean,                      # Average expression
    t = efit$t[, 1],                           # t-statistic
    P.value = efit$p.value[, 1],               # Raw p-value
    adj.P.value = p.adjust(efit$p.value[, 1], method = "BH"), # Adjusted p-value (BH)
    B = efit$lods[, 1],                         # Log-odds of differential expression
    contrast = contrast_name
  )
  
  # Step 2: Merge with gene symbols from efit$genes
  if (!is.null(efit$genes)) {
    efit_results_df <- merge(efit_results_df, efit$genes, by.x = "ensembleID", by.y = "ENSEMBL", all.x = TRUE)
  }
  
  # Step 3: Relocate SYMBOL and ensembleID for better display
  efit_results_df <- efit_results_df %>%
    relocate(SYMBOL) %>%
    relocate(ensembleID, .after = contrast)
  
  # Step 4: Create an interactive datatable
  datatable(
    efit_results_df %>% arrange(adj.P.value),  # Sort by adjusted p-value
    filter = "top",
    extensions = c('Buttons','Scroller'),
    caption = paste("Differential expression results for", contrast_name),
    options = list(
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel'),
      #pageLength = 10,
      autoWidth = TRUE,
      deferRender = TRUE,
      scrollY = 200,
      scroller = TRUE,
      order = list(list(which(names(efit_results_df) == "adj.P.value") - 1, "asc")),  # Sort by adj.P.value column
      columnDefs = list(list(className = 'dt-center', targets = "_all"))
    )
  ) %>%
    formatRound(columns = c("logFC", "AveExpr", "t", "P.value", "adj.P.value", "B"), digits = 3) %>%
    formatStyle(
      'logFC',
      background = styleInterval(c(-0.58, 0.58), c('lightblue', 'white', 'lightpink'))
    )
}

library(DT)
library(limma)
library(dplyr)
library(tidyr)

display_de_summary <- function(efit, lfc = 0.58, pval = 0.05, adjust_method = "BH",contrast) {
  # Step 1: Get summary of DE results (produces a table object)
  de_summary <- summary(decideTests(efit, lfc = lfc, adjust.method = adjust_method, p.value = pval))
  
  # Step 2: Convert table to data frame with explicit variable names
  de_summary_df <- as.data.frame(de_summary, responseName = "Number of Genes")
  
  # Step 3: Rename columns explicitly
  colnames(de_summary_df) <- c("Direction", "Contrast", "Number of Genes")
  
  # Step 4: Remove leading "X" and also any "X" after a "-" if followed by a digit
  de_summary_df$Contrast <- gsub("^X| X", "", de_summary_df$Contrast, perl = TRUE)
  
  de_summary_df <- de_summary_df %>% relocate(Contrast)
  
  # Step 5: Create an interactive datatable
  datatable(
    de_summary_df,
    rownames = FALSE,  # Hide row names
    caption = paste("Summary of Differential Expression Results (abs(lfc) > 0.58, adj.p.val < 0.05) for contrast",contrast),
    options = list(
      dom = 'Bfrtip',
      buttons = c('copy'),
      pageLength = 10,
      autoWidth = TRUE
    )
  )
}

library(dplyr)

create_efit_results_df <- function(efit) {
  # Step 1: Extract contrast name dynamically
  contrast_name <- gsub("^X| X", "", colnames(efit$contrasts))
  
  # Step 2: Construct results data frame
  efit_results_df <- data.frame(
    ensembleID = rownames(efit$coefficients),  # Extract rownames as ensemble IDs
    logFC = efit$coefficients[, 1],           # Log fold change
    AveExpr = efit$Amean,                     # Average expression
    t = efit$t[, 1],                          # t-statistic
    P.value = efit$p.value[, 1],              # Raw p-value
    adj.P.value = p.adjust(efit$p.value[, 1], method = "BH"), # Adjusted p-value (BH)
    B = efit$lods[, 1],                       # Log-odds of differential expression
    contrast = contrast_name                  # Use extracted contrast name
  )
  
  # Step 3: Merge with gene symbols from efit$genes
  if (!is.null(efit$genes)) {
    efit_results_df <- merge(efit_results_df, efit$genes, by.x = "ensembleID", by.y = "ENSEMBL", all.x = TRUE)
  }
  
  # Step 4: Reorder columns for readability
  efit_results_df <- efit_results_df %>% relocate(SYMBOL, .after = ensembleID)
  
  return(efit_results_df)
}

library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)

# Function to add Entrez IDs to a single efit results df
add_entrez_ids <- function(results_df) {
  # Extract Ensembl IDs
  ensembl_ids <- results_df$ensembleID
  
  # Convert to Entrez ID
  entrez_ids <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")  # Use "first" to handle multiple mappings
  
  # Add Entrez IDs to DataFrame
  results_df <- results_df %>%
    mutate(EntrezID = entrez_ids)
  
  return(results_df)
}

library(plotly)
library(RColorBrewer)
library(car)

plot_pca <- function(dge, title, grp_var=report_params$group_var, show_legend = TRUE, combine_plots = FALSE) {
  # Extract log-transformed CPM values
  PCA_DATA <- t(cpm(dge, log = TRUE))  # Transpose to have samples as rows
  numsonly <- as.data.frame(PCA_DATA)
  numsonly <- as.data.frame(lapply(numsonly, as.numeric))  # Ensure numeric values
  
  # Perform PCA
  pca_res <- prcomp(numsonly, center = TRUE, scale. = FALSE, rank. = 2)
  scores <- pca_res$x
  
  # Extract group information (batch or treatment)
  group_labels <- dge$samples[[grp_var]]  
  colors <- RColorBrewer::brewer.pal(length(unique(group_labels)), "Set1")
  
  # Initialize plotly object
  fig <- plot_ly()
  
  # Function to plot ellipses
  plot_ellipse_function <- function(target) {
    target_data <- as.data.frame(scores[group_labels == target, 1:2])
    target_color <- colors[which(unique(group_labels) == target)]
    
    ellipse_coords <- car::dataEllipse(target_data$PC1, target_data$PC2, 
                                       levels = 0.68, plot.points = FALSE, add = TRUE, draw = FALSE)
    
    fig <<- fig %>%
      add_polygons(x = ellipse_coords[,1], y = ellipse_coords[,2],
                   line = list(color = target_color, dash = "dot"),
                   fillcolor = target_color, opacity = 0.3,
                   showlegend = FALSE, hoverinfo = "skip")
  }
  
  # Function to plot points
  plot_points_function <- function(target) {
    target_data <- as.data.frame(scores[group_labels == target, 1:2])
    rownames(target_data) <- rownames(PCA_DATA)[group_labels == target]
    target_color <- colors[which(unique(group_labels) == target)]
    
    text_data <- as.character(rownames(target_data))
    PC1 <- as.numeric(target_data$PC1)
    PC2 <- as.numeric(target_data$PC2)
    
    fig <<- fig %>%
      add_trace(data = target_data, x = ~PC1, y = ~PC2, 
                type = "scatter", mode = "markers", name = target,
                marker = list(color = target_color),
                text = ~text_data, hoverinfo = "text", showlegend = show_legend)
  }
  
  # Add ellipses and points
  invisible(lapply(unique(group_labels), plot_ellipse_function))
  invisible(lapply(unique(group_labels), plot_points_function))
  
  # Add title
  fig <- fig %>% layout(title = title, xaxis = list(title = "PC1"), yaxis = list(title = "PC2"))
  
  # **If not combining, return single plot**
  if (!combine_plots) return(fig)
  
  # **Generate PCA plots for each filtering step**
  dge_list <- dge_list
  pca_plots <- mapply(plot_pca, dge_list, grp_var = grp_var, show_legend = c(TRUE, FALSE, FALSE), title = "", SIMPLIFY = FALSE)
  
  # **Combine the three plots in a single row**
  combined_plot <- subplot(pca_plots[[1]], pca_plots[[2]], pca_plots[[3]], 
                           nrows = 1, shareY = FALSE)
  
  # **Add annotations**
  annotations = list( 
    list(x = 0.15, y = 1.0, text = "PCA Before Filtering", xref = "paper", yref = "paper", 
         xanchor = "center", yanchor = "bottom", showarrow = FALSE),  
    list(x = 0.5, y = 1, text = "PCA After Filtering: NOISeq", xref = "paper", yref = "paper", 
         xanchor = "center", yanchor = "bottom", showarrow = FALSE),
    list(x = 0.85, y = 1, text = "PCA After Filtering: edgeR", xref = "paper", yref = "paper", 
         xanchor = "center", yanchor = "bottom", showarrow = FALSE)
  )
  
  # **Apply annotations and return combined plot**
  combined_plot <- combined_plot %>% layout(annotations = annotations)
  return(combined_plot)
}
# Define function to subset samples by contrast and run PCA
plot_pca_by_contrast <- function(dge_list, contrast_name, group_var) {
  # Extract groups from contrast name
  groups <- unlist(strsplit(gsub("efit_|_results_df", "", contrast_name), "_vs_"))
  
  # Subset samples in the DGEList based on the groups in the contrast
  filtered_dge_list <- lapply(dge_list, function(dge) {
    sample_subset <- dge$samples[dge$samples[[report_params$group_var]] %in% groups, ]
    dge_filtered <- dge
    dge_filtered$counts <- dge$counts[, colnames(dge$counts) %in% sample_subset$SampleName]
    dge_filtered$samples <- sample_subset
    return(dge_filtered)
  })
  
  # Generate PCA plots for raw, NOISeq-filtered, and edgeR-filtered data
  pca_plots <- mapply(plot_pca, filtered_dge_list, 
                      title = paste0("PCA for ", gsub("efit_", "", contrast_name)),  
                      grp_var = report_params$group_var, 
                      show_legend = c(TRUE, FALSE, FALSE), 
                      SIMPLIFY = FALSE)
  
  # Combine plots into one row
  combined_plot <- subplot(pca_plots[[1]], pca_plots[[2]], pca_plots[[3]], 
                           nrows = 1, shareY = FALSE)
  
  # Add annotations
  annotations = list( 
    list( 
      x = 0.15, y = 1.0,  
      text = "PCA Before Filtering", xref = "paper", yref = "paper", 
      xanchor = "center", yanchor = "bottom", showarrow = FALSE 
    ),  
    list( 
      x = 0.5, y = 1,  
      text = "PCA After Filtering: NOISeq", xref = "paper", yref = "paper",  
      xanchor = "center", yanchor = "bottom", showarrow = FALSE
    ),
    list( 
      x = 0.85, y = 1,  
      text = "PCA After Filtering: edgeR", xref = "paper", yref = "paper",  
      xanchor = "center", yanchor = "bottom", showarrow = FALSE
    )
  )
  
  # Apply annotations
  combined_plot <- combined_plot %>% layout(
    annotations = annotations,
    title = paste0("PCA for ", gsub("efit_", "", contrast_name))
  )
  
  return(combined_plot)
}

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

# Function to extract top DE genes and convert to Entrez IDs
non_DE_entrezIDs <- function(df) {
    filtered_genes <- df %>% filter(adj.P.value > 0.05|abs(logFC) < 0.58)
  
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



