# Rat Gut Microbiome Analysis (HFHS Diet, Wistar vs Sprague Dawley)

R scripts for the longitudinal analysis of gut microbiota composition and function in Wistar and Sprague Dawley rats fed an obesogenic high-fat, high-sucrose (HFHS) diet. This code was developed to evaluate whether and how host genetic background modulates the gut microbiome's response to diet, and can be adapted by others studying diet-microbiome interactions across rat strains or similar rodent model systems.

## What this code does

- Processes gut microbiome abundance, taxonomy, and sample metadata
- Computes alpha and beta diversity metrics across strains and time points
- Tests differential taxon abundance across strain and time using MaAsLin2 (samples from different time points were modeled as independent observations; rat ID was not included as a random effect)
- Tests differential abundance of predicted functional pathways using DESeq2
- Generates figures summarizing microbiome composition and diversity over time

## Requirements

- **R version:** 4.4.2 (should also work on other recent R 4.x versions)
- **Operating system:** Linux (tested on the CSC computing environment); should also run on macOS/Windows with the same R packages installed
- **R packages:**
  - phyloseq
  - ape
  - microbiome
  - vegan
  - ggplot2
  - dplyr
  - tidyr
  - tidyverse
  - lme4
  - lmerTest
  - emmeans
  - ggpubr
  - effectsize
  - patchwork
  - DESeq2
  - pheatmap
  - igraph
  - ggraph
  - data.table

## Installation

1. Clone or download this repository:
   ```bash
   git clone https://github.com/ljforeror/rat_microbiome_analysis.git
   ```
2. Open R or RStudio and install the required packages:
   ```r
   install.packages(c("vegan", "ggplot2", "dplyr", "tidyr", "tidyverse", "lme4",
                       "lmerTest", "emmeans", "ggpubr", "effectsize", "patchwork",
                       "pheatmap", "igraph", "ggraph", "data.table"))

   # phyloseq, ape, microbiome, and DESeq2 are distributed via Bioconductor:
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install(c("phyloseq", "ape", "microbiome", "DESeq2"))
   ```


## Testing

To confirm the code runs correctly on your system, run the test script in `tests/`, which executes the main analysis on the example data in `example_data/` and compares the result against the pre-computed expected output.

## License

This project is released under the MIT License. See the Copyright (c) 2026 Lady Johanna Forero Rodriguez et, al., Faculty of Medicine, University of Helsinki and Helsinki University Hospital (LICENSE) file for details.

## Generative AI disclosure

Portions of this code were developed with the assistance of a generative AI tool (Claude, Anthropic), in accordance with the ASM Generative AI Policy.

## Support

This repository is actively maintained. Issues and pull requests are welcome without prior approval. The authors commit to maintaining and supporting this resource for at least three years following publication.

## Citation

If you use this code, please cite [xxxxxxx].
