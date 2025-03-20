# RNAseq Reporting Pipeline

This repository provides a streamlined workflow for generating reports for RNA-seq differential expression analysis. The pipeline is designed to automate the process of statistical analysis and visualization, ensuring reproducibility and consistency.

## Overview

The analysis is controlled via a parameter file (`Report_params.txt`), and the main report generation script is `RNAseq_report.Rmd`. Custom functions for plotting and additional processing steps are stored in `HRK_funcs.R`. The report is rendered in RMarkdown and generates a detailed HTML summary of the RNA-seq results.

## Files

- **\`RNAseq_report.Rmd\`**: Main RMarkdown script that runs the full analysis and generates an HTML report.
- **\`analysis.R\`**: Core analysis script that processes RNA-seq data, performs differential expression analysis, and applies normalization methods.
- **\`HRK_funcs.R\`**: Custom R functions for visualization and other auxiliary tasks.
- **\`Report_params.txt\`**: Configuration file where all user-defined variables are specified.

## Usage

### Running the Pipeline

Ensure that **R** is loaded on your system. Then, run the following command to execute the report:

\`\`\`sh
module load R
Rscript -e 'rmarkdown::render("RNAseq_report.Rmd", 
                              output_file = "CX.novogene_brain.Report.html", 
                              params = list(param_file="My_Report_Params.txt"))'
\`\`\`

This will generate an HTML report summarizing the RNA-seq analysis.

### Customization

- Modify **\`Report_params.txt\`** to set experiment-specific parameters.
- Additional test files and explanations for parameters will be added soon.

## Dependencies

This pipeline requires R with the following packages:

- \`rmarkdown\`
- \`ggplot2\`
- \`plotly\`
- \`edgeR\`
- \`limma\`
- \`NOISeq\`
- \`clusterProfiler\`
- \`org.Hs.eg.db\` (for human datasets)
- \`org.Mm.eg.db\` (for mouse datasets)

To install missing dependencies, run:

\`\`\`r
install.packages(c("rmarkdown", "ggplot2", "plotly"))
BiocManager::install(c("edgeR", "limma", "NOISeq", "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db"))
\`\`\`

---

### Notes
- The analysis is **fully automated** once \`Report_params.txt\` is configured.
- The HTML report provides detailed PCA plots, heatmaps, volcano plots, and enrichment analyses.

More details and example data will be added soon!

