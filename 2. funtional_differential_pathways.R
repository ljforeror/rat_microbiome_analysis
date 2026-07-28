############################################################
############ PICRUSt2 – DESeq2 STRAIN PIPELINE ############
############ MEMORY-OPTIMIZED VERSION #####################
############################################################

library(data.table)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggraph)
library(ggrepel)
library(igraph)
library(patchwork)
library(pheatmap)
library(stringr)
library(tidyr)
library(tidyverse)

############################################################
############ CREATE RESULTS DIRECTORY
############################################################

if(!dir.exists("Results")){
  dir.create("Results", recursive = TRUE)
}

############################################################
# COLORS
############################################################

strain_colors <- c("#7B61FF", "#FF4FA3")


names(strain_colors) <- c("SD","Wi")

heat_colors <- colorRampPalette(
  c("#7A3CFF", "#F9FCFF", "#00CFC1")
)(100)

############################################################
# LOAD PICRUSt UNSTRATIFIED MATRIX
############################################################

picrust_path <- "/your_data_folder/path_abun_unstrat_with_desc.tsv.gz"

pathways_raw <- read_tsv(picrust_path, show_col_types = FALSE)

pathway_annotation <- pathways_raw %>%
  select(pathway, description)

pathways_mat <- pathways_raw %>%
  select(-description) %>%
  column_to_rownames("pathway")

count_matrix <- as.matrix(pathways_mat)
storage.mode(count_matrix) <- "numeric"

picrust_df<- read_tsv("/your_data_folder/path_abun_contrib.tsv.gz")

picrust_clases <- read_tsv(
  "/your_data_folder/your_metacyc_classes.tsv",
  show_col_types = FALSE
)
picrust_clases <- picrust_clases %>%
  dplyr::rename(pathway = feature)

picrust_clases <- picrust_clases %>%
  as.data.frame() %>%
  mutate(
    pathway = as.character(pathway),
    class1 = as.character(class1),
    class2 = as.character(class2)
  )
############################################################
# LOAD METADATA
############################################################

metadata <- read.csv(
  "/your_data_folder/Metadata.csv"
)
############################################################
# PIPELINE PICRUST - LONGITUDINAL + NETWORK + FINAL FIGURE
############################################################

dir.create("Results/picrust_results", recursive = TRUE, showWarnings = FALSE)
out <- "Results/picrust_results/"

############################################################
# COLORs
############################################################

strain_colors <- c("SD" = "#FF6EC7", "Wi" = "#6EC6FF")
pathway_color <- "#FFD1DC"
edge_color    <- "#BBBBBB"

############################################################
#  METADATA
############################################################

metadata <- metadata_df1 %>%
  dplyr::filter(Diet == "High Fat/Sucrose") %>%
  dplyr::filter(!is.na(Body_weight), !is.na(Time_point), !is.na(Strain)) %>%
  dplyr::mutate(Sample = make.names(Sample))

colnames(count_matrix) <- make.names(colnames(count_matrix))

common_samples <- intersect(metadata$Sample, colnames(count_matrix))
metadata       <- metadata[metadata$Sample %in% common_samples, ]
count_matrix   <- count_matrix[, common_samples]
metadata       <- metadata[match(colnames(count_matrix), metadata$Sample), ]

stopifnot(all(colnames(count_matrix) == metadata$Sample))

metadata$Strain <- factor(metadata$Strain, levels = c("SD", "Wi"))

metadata$Time_point <- factor(metadata$Time_point,
                              levels = c("0 weeks", "12 weeks", "24 weeks"),
                              labels = c("0weeks", "12weeks", "24weeks"))

metadata$Body_weight_scaled <- scale(metadata$Body_weight)

strain_ref  <- levels(metadata$Strain)[1]  # SD
strain_test <- levels(metadata$Strain)[2]  # Wi

set.seed(1991)

############################################################
# PICRUST 
############################################################

picrust_df <- picrust_df %>%
  dplyr::rename(pathway = `function`)

############################################################
# PATHWAY MAP
############################################################

pathway_map <- picrust_clases %>%
  dplyr::select(pathway, pathway_name) %>%
  dplyr::distinct()

############################################################
# FUNCTiON HELPER
############################################################

safe_deseq_df <- function(res) {
  data.frame(
    pathway        = rownames(res),
    baseMean       = as.numeric(res$baseMean),
    log2FoldChange = as.numeric(res$log2FoldChange),
    lfcSE          = as.numeric(res$lfcSE),
    stat           = as.numeric(res$stat),
    pvalue         = as.numeric(res$pvalue),
    padj           = as.numeric(res$padj),
    stringsAsFactors = FALSE
  )
}

############################################################
#  DESEQ2 LONGITUDINAL
############################################################

dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData   = metadata,
  design    = ~ Time_point + Strain + Body_weight_scaled + Strain:Time_point
)
dds <- DESeq(dds)

############################################################
#  Interaction terms
############################################################

interaction_terms <- grep("Time_point.*Strain|Strain.*Time_point",
                          resultsNames(dds), value = TRUE)

res_list <- lapply(interaction_terms, function(n) {
  df <- safe_deseq_df(results(dds, name = n))
  df$coef_name <- n
  df <- df[!is.na(df$padj), ]
  return(df)
})

res_interaction_df <- dplyr::bind_rows(res_list)

res_interaction_df <- res_interaction_df %>%
  dplyr::mutate(
    Time_point = gsub("Time_point|\\.Strain.*|Strain.*", "", coef_name),
    Time_point = gsub("^\\.", "", Time_point),
    Direction  = ifelse(log2FoldChange > 0, strain_test, strain_ref)
  )

############################################################
#  SIGNIFICANtS
############################################################

res_interaction_sig <- res_interaction_df[res_interaction_df$padj < 0.05, ]

write.csv(res_interaction_sig,
          paste0(out, "Table_Longitudinal_Significant.csv"),
          row.names = FALSE)

############################################################
#  TOP PATHWAYS (top 5 per timepoint and direction)
############################################################

top_pathways_df <- res_interaction_sig %>%
  dplyr::group_by(Time_point, Direction) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()

top_pathways <- unique(top_pathways_df$pathway)

############################################################
#  HEATMAP — FIGURE A
############################################################

vsd <- varianceStabilizingTransformation(dds[top_pathways, ], blind = FALSE)
mat <- assay(vsd)
mat <- t(scale(t(mat)))
mat[is.na(mat)] <- 0

long_df <- as.data.frame(mat) %>%
  tibble::rownames_to_column("pathway") %>%
  dplyr::left_join(pathway_map, by = "pathway") %>%
  dplyr::mutate(pathway_name = ifelse(is.na(pathway_name), pathway, pathway_name)) %>%
  tidyr::pivot_longer(-c(pathway, pathway_name),
                      names_to  = "Sample",
                      values_to = "value") %>%
  dplyr::mutate(Sample = as.character(Sample)) %>%
  dplyr::left_join(metadata %>% dplyr::mutate(Sample = as.character(Sample)),
                   by = "Sample")

if (any(is.na(long_df$Time_point))) stop("❌ Mismatch in Sample names")

p_long_heat <- ggplot(long_df, aes(Time_point, pathway_name, fill = value)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~Strain) +
  scale_fill_gradient2(low = "#B388FF", mid = "white", high = "#FF80AB", midpoint = 0) +
  theme_minimal(base_size = 16) +
  theme(panel.grid  = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text  = element_text(face = "bold")) +
  labs(title = "A  Longitudinal interaction pathways", fill = "Z-score")

ggsave(paste0(out, "Figure_A_Heatmap.png"), p_long_heat, width = 12, height = 8, dpi = 300)

############################################################
#  FUNCTIONAL CLASSES — FIGURE B
############################################################

class_df <- res_interaction_sig %>%
  dplyr::left_join(picrust_clases, by = "pathway") %>%
  dplyr::filter(!is.na(class1)) %>%
  dplyr::mutate(
    Time_point = gsub("Time_point|\\.Strain.*|Strain.*", "", coef_name),
    Time_point = gsub("^\\.", "", Time_point),
    Strain     = ifelse(log2FoldChange > 0, strain_test, strain_ref),
    Time_point = factor(Time_point, levels = c("0weeks", "12weeks", "24weeks"))
  ) %>%
  dplyr::group_by(Time_point, Strain, class1) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

p_class_time <- ggplot(class_df, aes(reorder(class1, n), n, fill = Strain)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  coord_flip() +
  facet_wrap(~Time_point) +
  scale_fill_manual(values = strain_colors) +
  theme_minimal(base_size = 16) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold")) +
  labs(title = "B  Functional classes over time", x = NULL, y = "N pathways")

ggsave(paste0(out, "Figure_B_Classes_Time.png"), p_class_time, width = 12, height = 8, dpi = 300)

############################################################
#  NETWORK — FIGURE C
############################################################

net_df <- picrust_df %>%
  dplyr::filter(pathway %in% top_pathways,
                taxon   %in% taxa_sig_features) %>%
  dplyr::left_join(pathway_map, by = "pathway") %>%
  dplyr::mutate(pathway_name = ifelse(is.na(pathway_name), pathway, pathway_name)) %>%
  dplyr::left_join(tax_map, by = c("taxon" = "feature")) %>%
  dplyr::mutate(DisplayName = ifelse(is.na(Genus), taxon, Genus)) %>%
  dplyr::left_join(metadata %>% dplyr::select(Sample, Time_point, Strain),
                   by = c("sample" = "Sample")) %>%
  dplyr::filter(!is.na(Time_point)) %>%
  dplyr::mutate(Significance = ifelse(pathway %in% res_interaction_sig$pathway,
                                      "Significant", "Not significant"))

top_taxa <- net_df %>%
  dplyr::group_by(Time_point, Strain, DisplayName) %>%
  dplyr::summarise(total = sum(taxon_function_abun), .groups = "drop") %>%
  dplyr::group_by(Time_point, Strain) %>%
  dplyr::slice_max(total, n = 3)

net_df <- net_df %>%
  dplyr::semi_join(top_taxa, by = c("Time_point", "Strain", "DisplayName"))

net_list <- lapply(levels(metadata$Time_point), function(tp) {
  
  df_tp <- net_df %>% dplyr::filter(Time_point == tp)
  if (nrow(df_tp) == 0) return(NULL)
  
  plots_strain <- lapply(levels(metadata$Strain), function(st) {
    
    df_st <- df_tp %>% dplyr::filter(Strain == st)
    if (nrow(df_st) == 0) return(NULL)
    
    g <- igraph::graph_from_data_frame(
      df_st %>% dplyr::select(DisplayName, pathway_name, taxon_function_abun),
      directed = FALSE
    )
    
    V(g)$type <- ifelse(V(g)$name %in% df_st$DisplayName, "Bacteria", "Pathway")
    
    pathway_sig <- df_st %>%
      dplyr::select(pathway_name, Significance) %>%
      dplyr::distinct()
    
    V(g)$color <- ifelse(
      V(g)$type == "Bacteria",
      strain_colors[st],
      ifelse(
        pathway_sig$Significance[match(V(g)$name, pathway_sig$pathway_name)] == "Significant",
        "#FF4FA3",
        "#D3D3D3"
      )
    )
    
    V(g)$size <- ifelse(V(g)$type == "Bacteria", 8, 5)
    
    ggraph(g, layout = "kk") +
      geom_edge_link(aes(width = taxon_function_abun),
                     color = edge_color, alpha = 0.4) +
      geom_node_point(aes(size = size), color = V(g)$color) +
      geom_node_text(
        aes(label = name),
        repel         = TRUE,
        size          = 5,
        fontface      = "bold",
        box.padding   = 0.5,
        point.padding = 0.3
      ) +
      scale_size(range = c(5, 12)) +
      theme_void(base_size = 16) +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      ) +
      labs(title = paste(tp, "-", st))
  })
  
  plots_strain <- Filter(Negate(is.null), plots_strain)
  patchwork::wrap_plots(plots_strain, nrow = 1)
})

net_list <- Filter(Negate(is.null), net_list)
p_network_final <- patchwork::wrap_plots(net_list, ncol = 1)

ggsave(paste0(out, "Figure_C_Network.png"),
       p_network_final, width = 24, height = 20, dpi = 400)

############################################################
#  MAIN FIGURE (A / B / C) — LARGER
############################################################

Figure_MAIN <- p_long_heat /
  p_class_time /
  patchwork::wrap_elements(p_network_final) +
  patchwork::plot_layout(heights = c(1, 1, 3)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 20, face = "bold")
    )
  )

ggsave(paste0(out, "Figure_MAIN.png"),
       Figure_MAIN,
       width  = 28,
       height = 36,
       dpi    = 300)

ggsave(paste0(out, "Figure_MAIN.pdf"),
       Figure_MAIN,
       width  = 28,
       height = 36)

cat("✔ Pipeline complete — results in Results/picrust_results/\n")
############################################################
#  MAIN FIGURE AB (A and B only)
############################################################

Figure_MAIN_AB <- p_long_heat /
  p_class_time +
  patchwork::plot_layout(heights = c(1, 1)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 20, face = "bold")
    )
  )

ggsave(paste0(out, "Figure_MAIN_AB.png"),
       Figure_MAIN_AB,
       width  = 16,
       height = 20,
       dpi    = 300)

ggsave(paste0(out, "Figure_MAIN_AB.pdf"),
       Figure_MAIN_AB,
       width  = 16,
       height = 20)

cat("✔ Figure AB saved in Results/picrust_results/\n")

############################################################
#  FIGURE C — NETWORK SEPARATE
############################################################

ggsave(paste0(out, "Figure_C_Network_FINAL.png"),
       p_network_final,
       width  = 28,
       height = 24,
       dpi    = 300)

ggsave(paste0(out, "Figure_C_Network_FINAL.pdf"),
       p_network_final,
       width  = 28,
       height = 24)

cat("✔ Figure C saved in Results/picrust_results/\n")

############################################################
# SANKEY — Bacteria → Functional class → Strain
############################################################

library(ggalluvial)

sankey_df <- picrust_df %>%
  dplyr::filter(pathway %in% res_interaction_sig$pathway,
                taxon   %in% taxa_sig_features) %>%
  dplyr::left_join(tax_map, by = c("taxon" = "feature")) %>%
  dplyr::mutate(DisplayName = ifelse(is.na(Genus), taxon, Genus)) %>%
  dplyr::left_join(picrust_clases, by = "pathway") %>%
  dplyr::filter(!is.na(class1)) %>%
  dplyr::left_join(metadata %>% dplyr::select(Sample, Strain),
                   by = c("sample" = "Sample")) %>%
  dplyr::filter(!is.na(Strain)) %>%
  dplyr::group_by(DisplayName, class1, Strain) %>%
  dplyr::summarise(contribution = sum(taxon_function_abun, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::semi_join(
    dplyr::group_by(., DisplayName) %>%
      dplyr::summarise(total = sum(contribution)) %>%
      dplyr::slice_max(total, n = 8),
    by = "DisplayName"
  )

p_sankey <- ggplot(sankey_df,
                   aes(axis1 = DisplayName,
                       axis2 = class1,
                       axis3 = Strain,
                       y     = contribution)) +
  geom_alluvium(aes(fill = Strain), alpha = 0.7, width = 0.3) +
  geom_stratum(fill = "grey95", color = "grey50", width = 0.3) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size     = 4,
            fontface = "bold",
            color    = "#333333") +
  scale_x_discrete(
    limits = c("Bacteria", "Functional class", "Strain"),
    expand = c(0.1, 0.1)
  ) +
  scale_fill_manual(values = strain_colors) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid      = element_blank(),
    axis.text.y     = element_blank(),
    axis.title.y    = element_blank(),
    legend.position = "bottom",
    plot.title      = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(title = "Bacteria → Functional class → Strain",
       fill  = "Strain")

ggsave(paste0(out, "Figure_D_Sankey.png"),
       p_sankey, width = 16, height = 12, dpi = 300)

ggsave(paste0(out, "Figure_D_Sankey.pdf"),
       p_sankey, width = 16, height = 12)

cat("✔ Sankey saved\n")

############################################################
# SUPPLEMENTARY TABLES
############################################################
supp_table1 <- res_interaction_sig %>%
  dplyr::left_join(picrust_clases, by = "pathway") %>%
  dplyr::mutate(
    pathway_name = ifelse(is.na(pathway_name), pathway, pathway_name),
    Direction    = ifelse(log2FoldChange > 0, strain_test, strain_ref),
    Time_point   = gsub("Time_point|\\.Strain.*|Strain.*", "", coef_name),
    Time_point   = gsub("^\\.", "", Time_point)
  ) %>%
  dplyr::select(
    pathway,
    pathway_name,
    class1,
    class2,
    Time_point,
    Direction,
    log2FoldChange,
    pvalue,
    padj
  ) %>%
  dplyr::arrange(padj)

write.csv(supp_table1,
          paste0(out, "Supplementary_Table_Significant_Pathways.csv"),
          row.names = FALSE)

supp_table2 <- supp_table1 %>%
  dplyr::filter(!is.na(class1)) %>%
  dplyr::group_by(class1, class2, Direction, Time_point) %>%
  dplyr::summarise(
    n_pathways  = dplyr::n(),
    mean_log2FC = round(mean(log2FoldChange), 3),
    median_padj = round(median(padj), 4),
    .groups     = "drop"
  ) %>%
  dplyr::arrange(class1, Time_point, Direction)

write.csv(supp_table2,
          paste0(out, "Supplementary_Table_Functional_Classes.csv"),
          row.names = FALSE)

cat("✔ Supplementary tables saved\n")
cat(paste0("  → Significant pathways: ", nrow(supp_table1), "\n"))
cat(paste0("  → Functional classes:      ", nrow(supp_table2), "\n"))
############################################################
# SUPPLEMENTARY TABLES — ALL MODEL TERMS
############################################################

# view all available terms
resultsNames(dds)

# ── function to extract any term ────────────────

extract_term <- function(dds, term_name) {
  safe_deseq_df(results(dds, name = term_name)) %>%
    dplyr::filter(!is.na(padj), padj < 0.05) %>%
    dplyr::left_join(picrust_clases, by = "pathway") %>%
    dplyr::mutate(
      pathway_name = ifelse(is.na(pathway_name), pathway, pathway_name),
      term         = term_name
    ) %>%
    dplyr::select(
      term,
      pathway,
      pathway_name,
      class1,
      class2,
      log2FoldChange,
      pvalue,
      padj
    ) %>%
    dplyr::arrange(padj)
}

# ── extract each term ──────────────────────────────────

# Strain
supp_strain <- extract_term(dds, grep("Strain", resultsNames(dds), value = TRUE) %>%
                              grep("Time", ., value = TRUE, invert = TRUE) %>%
                              .[1])

# Body_weight
supp_Body_weight <- extract_term(dds, "Body_weight_scaled")

# Interactions (all Strain:Time terms)
interaction_terms <- grep("Time_point.*Strain|Strain.*Time_point",
                          resultsNames(dds), value = TRUE)

supp_interaction <- dplyr::bind_rows(
  lapply(interaction_terms, function(t) extract_term(dds, t))
)

# ── general summary table ─────────────────────────────────

supp_all <- dplyr::bind_rows(
  supp_strain      %>% dplyr::mutate(Variable = "Strain"),
  supp_Body_weight      %>% dplyr::mutate(Variable = "Body_weight"),
  supp_interaction %>% dplyr::mutate(Variable = "Strain:Timepoint interaction")
)

# ── count summary ────────────────────────────────────

supp_summary <- supp_all %>%
  dplyr::group_by(Variable, class1) %>%
  dplyr::summarise(n_pathways = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(Variable, dplyr::desc(n_pathways))

# ── save ───────────────────────────────────────────────

write.csv(supp_strain,
          paste0(out, "Supplementary_Table_Strain_Significant.csv"),
          row.names = FALSE)

write.csv(supp_Body_weight,
          paste0(out, "Supplementary_Table_Body_weight_Significant.csv"),
          row.names = FALSE)

write.csv(supp_interaction,
          paste0(out, "Supplementary_Table_Interaction_Significant.csv"),
          row.names = FALSE)

write.csv(supp_all,
          paste0(out, "Supplementary_Table_ALL_Terms.csv"),
          row.names = FALSE)

write.csv(supp_summary,
          paste0(out, "Supplementary_Table_Summary_by_Variable.csv"),
          row.names = FALSE)

# ── console report ────────────────────────────────────

cat("✔ Supplementary tables saved\n")
cat(paste0("  → Strain:               ", nrow(supp_strain),      " pathways\n"))
cat(paste0("  → Body_weight:               ", nrow(supp_Body_weight),      " pathways\n"))
cat(paste0("  → Strain:Time interact: ", nrow(supp_interaction), " pathways\n"))
cat(paste0("  → Unique total:         ",
           length(unique(supp_all$pathway)),                      " pathways\n"))

############################################################
# SUMMARY HEATMAP BY STRAIN — DESeq2 + MaAsLin annotation
############################################################
############################################################
# SUMMARY HEATMAP BY STRAIN — DESeq2 + MaAsLin annotation
############################################################

library(pheatmap)

# ── 1. PREPARE TAXA LABEL (MaAsLin style) ───────────────

tax_map_label <- tax_map %>%
  dplyr::mutate(Taxa_label = dplyr::case_when(
    !is.na(Family) & Family != "" & !is.na(Genus) & Genus != "" ~
      paste0(Family, " | ", Genus),
    !is.na(Family) & Family != "" ~
      paste0(Family, " (", feature, ")"),
    TRUE ~ feature
  ))

maaslin_df <- maaslin_sig %>%
  dplyr::left_join(tax_map_label, by = "feature") %>%
  dplyr::mutate(coef = as.numeric(coef))

# ── 2. PREPARE DESeq2 MATRIX ─────────────────────────────

heat_df <- supp_all %>%
  dplyr::filter(!is.na(class1)) %>%
  dplyr::mutate(
    term_label = dplyr::case_when(
      Variable == "Strain"                       ~ "Strain",
      Variable == "Body_weight"                       ~ "Body_weight",
      Variable == "Strain:Timepoint interaction" ~
        gsub("^\\.|Time_point|\\.Strain.*|Strain.*", "", term),
      TRUE ~ term
    )
  )

heat_mat <- heat_df %>%
  dplyr::select(pathway_name, term_label, log2FoldChange) %>%
  dplyr::distinct(pathway_name, term_label, .keep_all = TRUE) %>%
  tidyr::pivot_wider(
    names_from  = term_label,
    values_from = log2FoldChange,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("pathway_name")

# z-score per row
heat_mat_z <- t(scale(t(as.matrix(heat_mat))))
heat_mat_z[is.na(heat_mat_z)] <- 0

# ── 3. ROW ANNOTATION — functional class ──────────────────

row_annotation <- heat_df %>%
  dplyr::select(pathway_name, class1) %>%
  dplyr::distinct(pathway_name, .keep_all = TRUE) %>%
  tibble::column_to_rownames("pathway_name")

row_annotation <- row_annotation[rownames(heat_mat_z), , drop = FALSE]

# ── 4. ROW ANNOTATION — MaAsLin ─────────────────────────

maaslin_driven_pathways <- picrust_df %>%
  dplyr::filter(taxon %in% taxa_sig_features) %>%
  dplyr::left_join(picrust_clases, by = "pathway") %>%
  dplyr::filter(!is.na(pathway_name)) %>%
  dplyr::pull(pathway_name) %>%
  unique()

# remove old column if it exists
row_annotation$MaAsLin_driven <- NULL

# add with correct name and labels in English
row_annotation$MaAsLin2 <- ifelse(
  rownames(row_annotation) %in% maaslin_driven_pathways,
  "Significant", "Not significant"
)

# ── 5. COLORS ────────────────────────────────────────────

n_classes  <- length(unique(row_annotation$class1))
class_cols <- setNames(
  colorRampPalette(c("#FF6EC7", "#6EC6FF", "#FFD700",
                     "#90EE90", "#FF8C00", "#DA70D6"))(n_classes),
  unique(row_annotation$class1)
)

ann_colors <- list(
  class1   = class_cols,
  MaAsLin2 = c("Significant" = "#7B61FF", "Not significant" = "grey90")
)

# ── 6. PLOT ───────────────────────────────────────────────

n_paths <- nrow(heat_mat_z)

for(ext in c("png", "pdf")) {
  
  if(ext == "png") {
    png(paste0(out, "Figure_Heatmap_DESeq2_MaAsLin.png"),
        width = 14, height = max(10, 6 + n_paths * 0.2),
        units = "in", res = 300)
  } else {
    pdf(paste0(out, "Figure_Heatmap_DESeq2_MaAsLin.pdf"),
        width = 14, height = max(10, 6 + n_paths * 0.2))
  }
  
  pheatmap(
    heat_mat_z,
    annotation_row    = row_annotation,
    annotation_colors = ann_colors,
    color             = colorRampPalette(c("#6EC6FF", "white", "#FF6EC7"))(100),
    cluster_cols      = FALSE,
    cluster_rows      = TRUE,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = 8,
    fontsize_col      = 12,
    fontsize          = 11,
    border_color      = NA,
    main              = "Integrated Analysis of Functional Pathways: DESeq2 Terms and MaAsLin2 Significant Taxa"
  )
  
  dev.off()
}

cat("✔ Heatmap saved\n")
cat(paste0("  → Total pathways:      ", n_paths, "\n"))
cat(paste0("  → MaAsLin2 significant: ",
           sum(row_annotation$MaAsLin2 == "Significant"), "\n"))

############################################################
#  HEATMAP DESEQ2 + MAASLIN2
############################################################

library(pheatmap)

# ── PREPARE MATRIX ───────────────────────────────────────

heat_df_pub <- supp_all %>%
  dplyr::filter(!is.na(class1)) %>%
  dplyr::mutate(
    term_label = dplyr::case_when(
      Variable == "Strain"                       ~ "Strain",
      Variable == "Body_weight"                       ~ "Body_weight",
      Variable == "Strain:Timepoint interaction" ~
        gsub("^\\.|Time_point|\\.Strain.*|Strain.*", "", term),
      TRUE ~ term
    )
  )

heat_mat_pub <- heat_df_pub %>%
  dplyr::select(pathway_name, term_label, log2FoldChange) %>%
  dplyr::distinct(pathway_name, term_label, .keep_all = TRUE) %>%
  tidyr::pivot_wider(
    names_from  = term_label,
    values_from = log2FoldChange,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("pathway_name")

# z-score per row
heat_mat_z_pub <- t(scale(t(as.matrix(heat_mat_pub))))
heat_mat_z_pub[is.na(heat_mat_z_pub)] <- 0

# ── ROW ANNOTATION ───────────────────────────────────────

row_ann_pub <- heat_df_pub %>%
  dplyr::select(pathway_name, class1) %>%
  dplyr::distinct(pathway_name, .keep_all = TRUE) %>%
  tibble::column_to_rownames("pathway_name")

row_ann_pub <- row_ann_pub[rownames(heat_mat_z_pub), , drop = FALSE]

# MaAsLin2
maaslin_driven_pub <- picrust_df %>%
  dplyr::filter(taxon %in% taxa_sig_features) %>%
  dplyr::left_join(picrust_clases, by = "pathway") %>%
  dplyr::filter(!is.na(pathway_name)) %>%
  dplyr::pull(pathway_name) %>%
  unique()

row_ann_pub$MaAsLin2 <- ifelse(
  rownames(row_ann_pub) %in% maaslin_driven_pub,
  "Significant", "Not significant"
)

# ── COLORS ───────────────────────────────────────────────

n_classes_pub  <- length(unique(row_ann_pub$class1))
class_cols_pub <- setNames(
  colorRampPalette(c("#FF6EC7", "#6EC6FF", "#FFD700",
                     "#90EE90", "#FF8C00", "#DA70D6"))(n_classes_pub),
  unique(row_ann_pub$class1)
)

ann_colors_pub <- list(
  class1   = class_cols_pub,
  MaAsLin2 = c("Significant" = "#7B61FF", "Not significant" = "grey90")
)

n_paths_pub <- nrow(heat_mat_z_pub)
# remove MaAsLin2 from annotation
row_ann_pub$MaAsLin2 <- NULL
ann_colors_pub$MaAsLin2 <- NULL

# ── PLOT ──────────────────────────────────────────────────

for(ext in c("png", "pdf")) {
  
  if(ext == "png") {
    png(paste0(out, "Figure_1_Heatmap_DESeq2_MaAsLin2.png"),
        width  = 14,
        height = max(10, 4 + n_paths_pub * 0.15),
        units  = "in", res = 300)
  } else {
    pdf(paste0(out, "Figure_1_Heatmap_DESeq2_MaAsLin2.pdf"),
        width  = 14,
        height = max(10, 4 + n_paths_pub * 0.15))
  }
  
  pheatmap(
    heat_mat_z_pub,
    annotation_row    = row_ann_pub,
    annotation_colors = ann_colors_pub,
    color             = colorRampPalette(c("#6EC6FF", "white", "#FF6EC7"))(100),
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = 7,
    fontsize_col      = 12,
    fontsize          = 11,
    border_color      = NA,
    main              = "Integrated Analysis of Functional Pathways: DESeq2 Terms and MaAsLin2 Significant Taxa"
  )
  
  dev.off()
}

cat("✔ Heatmap DESeq2 + MaAsLin2 saved\n")
cat(paste0("  -> Pathways:         ", n_paths_pub, "\n"))
cat(paste0("  -> MaAsLin2 sig:     ",
           sum(row_ann_pub$MaAsLin2 == "Significant"), "\n"))

############################################################
#  SANKEY — Bacteria -> Functional class -> Strain
############################################################

library(ggalluvial)

sankey_df <- picrust_df %>%
  dplyr::filter(pathway %in% res_interaction_sig$pathway,
                taxon   %in% taxa_sig_features) %>%
  dplyr::left_join(tax_map, by = c("taxon" = "feature")) %>%
  dplyr::mutate(
    DisplayName = dplyr::case_when(
      !is.na(Family) & Family != "" & !is.na(Genus) & Genus != "" ~
        paste0(Family, " | ", Genus),
      !is.na(Family) & Family != "" ~ Family,
      TRUE ~ taxon
    )
  ) %>%
  dplyr::left_join(picrust_clases, by = "pathway") %>%
  dplyr::filter(!is.na(class1)) %>%
  dplyr::left_join(metadata %>% dplyr::select(Sample, Strain),
                   by = c("sample" = "Sample")) %>%
  dplyr::filter(!is.na(Strain)) %>%
  dplyr::group_by(DisplayName, class1, Strain) %>%
  dplyr::summarise(contribution = sum(taxon_function_abun, na.rm = TRUE),
                   .groups = "drop")

# top 8 bacteria
top8_bacteria <- sankey_df %>%
  dplyr::group_by(DisplayName) %>%
  dplyr::summarise(total = sum(contribution), .groups = "drop") %>%
  dplyr::slice_max(total, n = 8) %>%
  dplyr::pull(DisplayName)

sankey_df <- sankey_df %>%
  dplyr::filter(DisplayName %in% top8_bacteria)

p_sankey <- ggplot(sankey_df,
                   aes(axis1 = DisplayName,
                       axis2 = class1,
                       axis3 = Strain,
                       y     = contribution)) +
  geom_alluvium(aes(fill = Strain), alpha = 0.7, width = 0.3) +
  geom_stratum(fill = "grey95", color = "grey50", width = 0.3) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size     = 4,
            fontface = "bold",
            color    = "#333333") +
  scale_x_discrete(
    limits = c("Bacteria", "Functional class", "Strain"),
    expand = c(0.1, 0.1)
  ) +
  scale_fill_manual(values = strain_colors) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid      = element_blank(),
    axis.text.y     = element_blank(),
    axis.title.y    = element_blank(),
    legend.position = "bottom",
    plot.title      = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(title = "Bacteria - Functional Class - Strain",
       fill  = "Strain")

ggsave(paste0(out, "Figure_2_Sankey.png"),
       p_sankey, width = 16, height = 12, dpi = 300)
ggsave(paste0(out, "Figure_2_Sankey.pdf"),
       p_sankey, width = 16, height = 12)

cat("✔ Sankey saved\n")

############################################################
#  SUPPLEMENTARY TABLE — DETAILED CONTRIBUTIONS
############################################################

supp_contributions <- picrust_df %>%
  dplyr::filter(pathway %in% res_interaction_sig$pathway,
                taxon   %in% taxa_sig_features) %>%
  dplyr::left_join(tax_map, by = c("taxon" = "feature")) %>%
  dplyr::mutate(
    DisplayName = dplyr::case_when(
      !is.na(Family) & Family != "" & !is.na(Genus) & Genus != "" ~
        paste0(Family, " | ", Genus),
      !is.na(Family) & Family != "" ~ Family,
      TRUE ~ taxon
    )
  ) %>%
  dplyr::left_join(picrust_clases, by = "pathway") %>%
  dplyr::left_join(metadata %>% dplyr::select(Sample, Strain, Time_point),
                   by = c("sample" = "Sample")) %>%
  dplyr::filter(!is.na(Strain)) %>%
  dplyr::group_by(pathway, pathway_name, class1, class2,
                  DisplayName, Strain, Time_point) %>%
  dplyr::summarise(
    mean_contribution = mean(taxon_function_abun, na.rm = TRUE),
    total_contribution = sum(taxon_function_abun, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    res_interaction_sig %>%
      dplyr::select(pathway, log2FoldChange, padj, coef_name) %>%
      dplyr::mutate(
        Interaction_timepoint = dplyr::case_when(
          grepl("24weeks", coef_name) ~ "24 weeks",
          grepl("12weeks", coef_name) ~ "12 weeks",
          TRUE ~ "Unknown"
        )
      ) %>%
      dplyr::select(pathway, log2FoldChange, padj, Interaction_timepoint) %>%
      dplyr::distinct(pathway, .keep_all = TRUE),
    by = "pathway"
  ) %>%
  dplyr::arrange(padj, pathway_name, Strain, Time_point)

write.csv(supp_contributions,
          paste0(out, "Supplementary_Table_PiCRUST_Contributions.csv"),
          row.names = FALSE)

cat("✔ Supplementary table saved\n")
cat(paste0("  -> Rows:     ", nrow(supp_contributions), "\n"))
cat(paste0("  -> Pathways:  ", length(unique(supp_contributions$pathway)), "\n"))
cat(paste0("  -> Bacteria: ", length(unique(supp_contributions$DisplayName)), "\n"))

############################################################
# COMBINED FIGURE — Heatmap + Sankey + Class barplot
############################################################

library(pheatmap)
library(ggalluvial)
library(patchwork)
library(grid)
library(gridExtra)

# ── A. TOP 20 HEATMAP BY PADJ ───────────────────────────

top20_paths <- res_interaction_sig %>%
  dplyr::arrange(padj) %>%
  dplyr::distinct(pathway, .keep_all = TRUE) %>%
  dplyr::slice(1:20) %>%
  dplyr::pull(pathway)

heat_top20 <- heat_mat_z_pub[
  rownames(heat_mat_z_pub) %in%
    (picrust_clases$pathway_name[picrust_clases$pathway %in% top20_paths]),
  , drop = FALSE
]

# if they don't match by pathway_name, use directly
if(nrow(heat_top20) == 0) {
  top20_names <- supp_all %>%
    dplyr::filter(pathway %in% top20_paths) %>%
    dplyr::pull(pathway_name) %>%
    unique()
  heat_top20 <- heat_mat_z_pub[
    rownames(heat_mat_z_pub) %in% top20_names, , drop = FALSE
  ]
}

row_ann_top20 <- row_ann_pub[rownames(heat_top20), , drop = FALSE]

# ── B. SANKEY ─────────────────────────────────────────────

# you already have p_sankey from the previous block

# ── C. FUNCTIONAL CLASSES BARPLOT BY STRAIN ─────────────

class_strain_df <- res_interaction_sig %>%
  dplyr::left_join(picrust_clases %>%
                     dplyr::select(pathway, class1),
                   by = "pathway") %>%
  dplyr::filter(!is.na(class1)) %>%
  dplyr::mutate(
    Strain = ifelse(log2FoldChange > 0, strain_test, strain_ref)
  ) %>%
  dplyr::group_by(class1, Strain) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

p_class_strain <- ggplot(class_strain_df,
                         aes(x = reorder(class1, n),
                             y = n,
                             fill = Strain)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  coord_flip() +
  scale_fill_manual(values = strain_colors) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(title    = "C  |  Functional classes by Strain",
       x        = NULL,
       y        = "Number of significant pathways",
       fill     = "Strain enriched")

# ── SAVE HEATMAP AS GROB ─────────────────────────────

heatmap_grob <- function() {
  pheatmap(
    heat_top20,
    annotation_row    = row_ann_top20,
    annotation_colors = ann_colors_pub,
    color             = colorRampPalette(c("#6EC6FF", "white", "#FF6EC7"))(100),
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = 16,
    fontsize_col      = 16,
    fontsize          = 16,
    border_color      = NA,
    main              = "A  |  Top 10 Pathways by Strain:Timepoint Interaction",
    silent            = TRUE
  )$gtable
}

# ── SANKEY WITH LARGE TEXT AND ADJUSTED BOXES ───────────

p_sankey <- ggplot(sankey_df,
                   aes(axis1 = DisplayName,
                       axis2 = class1,
                       axis3 = Strain,
                       y     = contribution)) +
  geom_alluvium(aes(fill = Strain), alpha = 0.7, width = 0.4) +
  geom_stratum(fill = "grey95", color = "grey50", width = 0.4) +
  geom_text(stat    = "stratum",
            aes(label = after_stat(stratum)),
            size     = 6,
            fontface = "bold",
            color    = "#333333",
            check_overlap = FALSE) +
  scale_x_discrete(
    limits = c("Bacteria", "Functional class", "Strain"),
    expand = c(0.15, 0.15)
  ) +
  scale_fill_manual(values = strain_colors) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid      = element_blank(),
    axis.text.y     = element_blank(),
    axis.title.y    = element_blank(),
    axis.text.x     = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    legend.text     = element_text(size = 14),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 16)
  ) +
  labs(title = "B  |  Bacteria - Functional Class - Strain",
       fill  = "Strain")

# ── COMBINED FIGURE ──────────────────────────────────────

png(paste0(out, "Figure_MAIN_PiCRUST.png"),
    width = 26, height = 24, units = "in", res = 300)

grid.arrange(
  heatmap_grob(),
  p_sankey,
  nrow    = 2,
  heights = c(1.2, 1)
)

dev.off()

pdf(paste0(out, "Figure_MAIN_PiCRUST.pdf"),
    width = 26, height = 24)

grid.arrange(
  heatmap_grob(),
  p_sankey,
  nrow    = 2,
  heights = c(1.2, 1)
)

dev.off()

cat("✔ Combined PiCRUST figure saved with large text\n")


############################################################
# PALETTES
############################################################

strain_colors <- c("SD" = "#FF6EC7", "Wi" = "#6EC6FF")
heat_colors   <- colorRampPalette(c("#A78BFF", "white", "#FFAA80"))(100)

class_palette <- c(
  "#FF6EC7", "#6EC6FF", "#A78BFF", "#FFAA80",
  "#FF9EE5", "#FFE566", "#5EEAD4", "#FF7C7C",
  "#FFC0F0", "#7EB8FF", "#FFD4A8", "#C4B5FD"
)

############################################################
# CORRECTED NAMES
############################################################

pathway_name_fixes <- c(
  "PWY-7294" = "D-xylose degradation IV",
  "PWY-7376" = "Cob(II)yrinate a,c-diamide biosynthesis II (late cobalt incorporation)",
  "P381-PWY" = "Adenosylcobalamin biosynthesis II (aerobic)"
)

clean_pathway_name <- function(x) {
  for (id in names(pathway_name_fixes)) {
    x <- gsub(id, pathway_name_fixes[id], x, fixed = TRUE)
  }
  x <- sapply(x, function(nm) {
    paste0(toupper(substr(nm, 1, 1)), substr(nm, 2, nchar(nm)))
  })
  return(x)
}

############################################################
# PREPARE PANEL A DATA
# Mean Z-score by strain and time point
# Top 20 pathways ordered by padj
############################################################

top_paths_ids <- res_interaction_sig %>%
  dplyr::arrange(padj) %>%
  dplyr::distinct(pathway, .keep_all = TRUE) %>%
  dplyr::slice(1:10) %>%
  dplyr::pull(pathway)

top_pathways_ordered <- res_interaction_sig %>%
  dplyr::arrange(padj) %>%
  dplyr::distinct(pathway, .keep_all = TRUE) %>%
  dplyr::slice(1:10) %>%
  dplyr::left_join(pathway_map, by = "pathway") %>%
  dplyr::mutate(
    pathway_name = ifelse(is.na(pathway_name), pathway, pathway_name),
    pathway_name = clean_pathway_name(pathway_name)
  ) %>%
  dplyr::pull(pathway_name)

# VSD
vsd_top <- varianceStabilizingTransformation(
  dds[top_paths_ids, ], blind = FALSE
)
mat_top <- assay(vsd_top)
mat_top <- t(scale(t(mat_top)))
mat_top[is.na(mat_top)] <- 0

# Long format
long_df_top <- as.data.frame(mat_top) %>%
  tibble::rownames_to_column("pathway") %>%
  dplyr::left_join(pathway_map, by = "pathway") %>%
  dplyr::mutate(
    pathway_name = ifelse(is.na(pathway_name), pathway, pathway_name),
    pathway_name = clean_pathway_name(pathway_name)
  ) %>%
  tidyr::pivot_longer(-c(pathway, pathway_name),
                      names_to  = "Sample",
                      values_to = "zscore") %>%
  dplyr::mutate(Sample = as.character(Sample)) %>%
  dplyr::left_join(
    metadata %>%
      dplyr::mutate(Sample = as.character(Sample)) %>%
      dplyr::select(Sample, Strain, Time_point),
    by = "Sample"
  ) %>%
  dplyr::filter(!is.na(Strain), !is.na(Time_point))

# Mean by group
long_df_mean <- long_df_top %>%
  dplyr::group_by(pathway_name, Strain, Time_point) %>%
  dplyr::summarise(mean_z = mean(zscore, na.rm = TRUE),
                   .groups = "drop")

# Labels with padj
padj_labels <- res_interaction_sig %>%
  dplyr::arrange(padj) %>%
  dplyr::distinct(pathway, .keep_all = TRUE) %>%
  dplyr::slice(1:20) %>%
  dplyr::left_join(pathway_map, by = "pathway") %>%
  dplyr::mutate(
    pathway_name = ifelse(is.na(pathway_name), pathway, pathway_name),
    pathway_name = clean_pathway_name(pathway_name),
    padj_label   = paste0(pathway_name,
                          " (q=", formatC(padj, format = "e", digits = 1), ")")
  ) %>%
  dplyr::select(pathway_name, padj_label)

long_df_mean <- long_df_mean %>%
  dplyr::left_join(padj_labels, by = "pathway_name")

long_df_mean$padj_label <- factor(
  long_df_mean$padj_label,
  levels = rev(unique(padj_labels$padj_label))
)

############################################################
# PANEL A — Heatmap mean Z-score SD vs Wi by time point
############################################################

p_5A <- ggplot(long_df_mean,
               aes(x = Time_point, y = padj_label, fill = mean_z)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(mean_z, 1)),
            size = 4, color = "grey20") +
  facet_wrap(~Strain, nrow = 1) +
  scale_fill_gradient2(
    low      = "#A78BFF",
    mid      = "white",
    high     = "#FFAA80",
    midpoint = 0,
    name     = "Mean\nZ-score"
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1,
                                    size = 13, face = "bold"),
    axis.text.y      = element_text(size = 11),
    strip.text       = element_text(face = "bold", size = 15),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title       = element_text(face = "bold", size = 16, hjust = 0),
    plot.subtitle    = element_text(size = 11, color = "grey40", hjust = 0),
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 13),
    legend.position  = "right"
  ) +
  labs(
    title    = "A  |  Top 10 Strain × Time interaction pathways (mean Z-score)",
    subtitle = "Ordered by adjusted p-value (most significant at top) | SD and Wi shown separately",
    x        = "Time point",
    y        = NULL
  )

############################################################
# PANEL B — Functional classes barplot by time and strain
############################################################

class_df <- res_interaction_sig %>%
  dplyr::left_join(picrust_clases %>%
                     dplyr::select(pathway, class1),
                   by = "pathway") %>%
  dplyr::filter(!is.na(class1)) %>%
  dplyr::mutate(
    Time_point = gsub("Time_point|\\.Strain.*|Strain.*", "", coef_name),
    Time_point = gsub("^\\.", "", Time_point),
    Time_point = factor(Time_point,
                        levels = c("0weeks", "12weeks", "24weeks"),
                        labels = c("0 weeks", "12 weeks", "24 weeks")),
    Strain     = ifelse(log2FoldChange > 0, strain_test, strain_ref),
    class1     = clean_pathway_name(class1)
  ) %>%
  dplyr::group_by(Time_point, Strain, class1) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

p_5B <- ggplot(class_df,
               aes(x = reorder(class1, n), y = n, fill = Strain)) +
  geom_bar(stat = "identity", position = "dodge",
           color = "black", linewidth = 0.4) +
  coord_flip() +
  facet_wrap(~Time_point, nrow = 1) +
  scale_fill_manual(values = strain_colors) +
  theme_classic(base_size = 15) +
  theme(
    strip.text       = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title       = element_text(face = "bold", hjust = 0, size = 16),
    axis.text.y      = element_text(size = 12),
    axis.text.x      = element_text(size = 12),
    axis.title.x     = element_text(size = 14),
    legend.text      = element_text(size = 13),
    legend.title     = element_text(size = 14),
    legend.position  = "right"
  ) +
  labs(
    title = "B  |  Functional pathway classes across time points",
    x     = NULL,
    y     = "Number of significant pathways (q < 0.05)",
    fill  = "Enriched in"
  )

############################################################
# PANEL C — Sankey bacteria → class → strain
############################################################

sankey_df <- picrust_df %>%
  dplyr::filter(pathway %in% res_interaction_sig$pathway,
                taxon   %in% taxa_sig_features) %>%
  dplyr::left_join(tax_map, by = c("taxon" = "feature")) %>%
  dplyr::mutate(
    DisplayName = dplyr::case_when(
      !is.na(Family) & Family != "" & !is.na(Genus) & Genus != "" ~
        paste0(Family, " | ", Genus),
      !is.na(Family) & Family != "" ~ Family,
      TRUE ~ taxon
    )
  ) %>%
  dplyr::left_join(picrust_clases, by = "pathway") %>%
  dplyr::filter(!is.na(class1)) %>%
  dplyr::mutate(class1 = clean_pathway_name(class1)) %>%
  dplyr::left_join(metadata %>% dplyr::select(Sample, Strain),
                   by = c("sample" = "Sample")) %>%
  dplyr::filter(!is.na(Strain)) %>%
  dplyr::group_by(DisplayName, class1, Strain) %>%
  dplyr::summarise(contribution = sum(taxon_function_abun, na.rm = TRUE),
                   .groups = "drop")

top8_bacteria <- sankey_df %>%
  dplyr::group_by(DisplayName) %>%
  dplyr::summarise(total = sum(contribution), .groups = "drop") %>%
  dplyr::slice_max(total, n = 8) %>%
  dplyr::pull(DisplayName)

sankey_df <- sankey_df %>%
  dplyr::filter(DisplayName %in% top8_bacteria)

p_5C <- ggplot(sankey_df,
               aes(axis1 = DisplayName,
                   axis2 = class1,
                   axis3 = Strain,
                   y     = contribution)) +
  geom_alluvium(aes(fill = Strain), alpha = 0.7, width = 0.4) +
  geom_stratum(fill = "grey95", color = "grey50", width = 0.4) +
  geom_text(stat         = "stratum",
            aes(label    = after_stat(stratum)),
            size         = 5,
            fontface     = "bold",
            color        = "#333333",
            check_overlap = FALSE) +
  scale_x_discrete(
    limits = c("Bacteria", "Functional class", "Strain"),
    expand = c(0.15, 0.15)
  ) +
  scale_fill_manual(values = strain_colors) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid      = element_blank(),
    axis.text.y     = element_blank(),
    axis.title.y    = element_blank(),
    axis.text.x     = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.text     = element_text(size = 13),
    legend.title    = element_text(size = 14),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 16)
  ) +
  labs(
    title = "C  |  Microbial taxa contributions to functional classes by strain",
    fill  = "Strain"
  )

#############################################################
# FIXING THE FIGURE BEFORE
#############################################################
############################################################
# REQUIRED LIBRARIES
############################################################
library(ggplot2)
library(ggalluvial)
library(ggfittext)
library(magick)
library(dplyr)

############################################################
# PANEL A — Heatmap mean Z-score SD vs Wi by time point
############################################################

p_5A <- ggplot(long_df_mean,
               aes(x = Time_point, y = padj_label, fill = mean_z)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(mean_z, 1)),
            size = 4, color = "grey20") +
  facet_wrap(~Strain, nrow = 1) +
  scale_fill_gradient2(
    low      = "#A78BFF",
    mid      = "white",
    high     = "#FFAA80",
    midpoint = 0,
    name     = "Mean\nZ-score"
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1,
                                    size = 13, face = "bold"),
    axis.text.y      = element_text(size = 11),
    strip.text       = element_text(face = "bold", size = 15),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title       = element_text(face = "bold", size = 16, hjust = 0),
    plot.subtitle    = element_text(size = 11, color = "grey40", hjust = 0),
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 13),
    legend.position  = "right"
  ) +
  labs(
    title    = "A  |  Top 10 Strain × Time interaction pathways (mean Z-score)",
    subtitle = "Ordered by adjusted p-value (most significant at top) | SD and Wi shown separately",
    x        = "Time point",
    y        = NULL
  )

############################################################
# PANEL B — Functional classes barplot by time and strain
############################################################

class_df <- res_interaction_sig %>%
  dplyr::left_join(picrust_clases %>%
                     dplyr::select(pathway, class1),
                   by = "pathway") %>%
  dplyr::filter(!is.na(class1)) %>%
  dplyr::mutate(
    Time_point = gsub("Time_point|\\.Strain.*|Strain.*", "", coef_name),
    Time_point = gsub("^\\.", "", Time_point),
    Time_point = factor(Time_point,
                        levels = c("0weeks", "12weeks", "24weeks"),
                        labels = c("0 weeks", "12 weeks", "24 weeks")),
    Strain     = ifelse(log2FoldChange > 0, strain_test, strain_ref),
    class1     = clean_pathway_name(class1)
  ) %>%
  dplyr::group_by(Time_point, Strain, class1) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

p_5B <- ggplot(class_df,
               aes(x = reorder(class1, n), y = n, fill = Strain)) +
  geom_bar(stat = "identity", position = "dodge",
           color = "black", linewidth = 0.4) +
  coord_flip() +
  facet_wrap(~Time_point, nrow = 1) +
  scale_fill_manual(values = strain_colors) +
  theme_classic(base_size = 15) +
  theme(
    strip.text       = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title       = element_text(face = "bold", hjust = 0, size = 16),
    axis.text.y      = element_text(size = 12),
    axis.text.x      = element_text(size = 12),
    axis.title.x     = element_text(size = 14),
    legend.text      = element_text(size = 13),
    legend.title     = element_text(size = 14),
    legend.position  = "right"
  ) +
  labs(
    title = "B  |  Functional pathway classes across time points",
    x     = NULL,
    y     = "Number of significant pathways (q < 0.05)",
    fill  = "Enriched in"
  )

############################################################
# PANEL C — Sankey bacteria → class → strain (CORRECTED)
############################################################

sankey_df <- picrust_df %>%
  dplyr::filter(pathway %in% res_interaction_sig$pathway,
                taxon   %in% taxa_sig_features) %>%
  dplyr::left_join(tax_map, by = c("taxon" = "feature")) %>%
  dplyr::mutate(
    DisplayName = dplyr::case_when(
      !is.na(Family) & Family != "" & !is.na(Genus) & Genus != "" ~
        paste0(Family, "\n", Genus),
      !is.na(Family) & Family != "" ~ Family,
      TRUE ~ taxon
    )
  ) %>%
  dplyr::left_join(picrust_clases, by = "pathway") %>%
  dplyr::filter(!is.na(class1)) %>%
  dplyr::mutate(class1 = clean_pathway_name(class1)) %>%
  dplyr::left_join(metadata %>% dplyr::select(Sample, Strain),
                   by = c("sample" = "Sample")) %>%
  dplyr::filter(!is.na(Strain)) %>%
  dplyr::group_by(DisplayName, class1, Strain) %>%
  dplyr::summarise(contribution = sum(taxon_function_abun, na.rm = TRUE),
                   .groups = "drop")

top8_bacteria <- sankey_df %>%
  dplyr::group_by(DisplayName) %>%
  dplyr::summarise(total = sum(contribution), .groups = "drop") %>%
  dplyr::slice_max(total, n = 8) %>%
  dplyr::pull(DisplayName)

sankey_df <- sankey_df %>%
  dplyr::filter(DisplayName %in% top8_bacteria)


p_5C <- ggplot(sankey_df,
               aes(axis1 = DisplayName,
                   axis2 = class1,
                   axis3 = Strain,
                   y     = contribution)) +
  geom_alluvium(aes(fill = Strain), alpha = 0.7, width = 0.4) +
  geom_stratum(fill = "grey95", color = "grey50", width = 0.4) +
  ggfittext::geom_fit_text(
    stat          = "stratum",
    aes(label     = after_stat(stratum)),
    width         = 0.4,
    min.size      = 1,
    size          = 10,
    fontface      = "bold",
    color         = "#333333",
    place         = "center",
    reflow        = TRUE,
    grow          = FALSE
  ) +
  scale_x_discrete(
    limits = c("Bacteria", "Functional class", "Strain"),
    expand = c(0.2, 0.2)
  ) +
  scale_fill_manual(values = strain_colors) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid      = element_blank(),
    axis.text.y     = element_blank(),
    axis.title.y    = element_blank(),
    axis.text.x     = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.text     = element_text(size = 13),
    legend.title    = element_text(size = 14),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 15),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    title = "C  |  Microbial taxa contributions to functional classes by strain",
    fill  = "Strain"
  )
p_5C

############################################################
# SAVE FIGURE 5 — HIGH RESOLUTION (patchwork, no magick)
############################################################
p_5C <- p_5C +
  scale_x_discrete(
    limits = c("Bacteria", "Functional class", "Strain"),
    expand = c(0.05, 0.05)
  ) +
  theme(
    aspect.ratio = NULL,
    plot.margin  = margin(5, 5, 5, 5)
  )
library(patchwork)

# Increase font size proportional to the combined figure
p_5A <- p_5A + theme(text = element_text(size = rel(1.15)))
p_5B <- p_5B + theme(text = element_text(size = rel(1.15)))
p_5C <- p_5C + theme(text = element_text(size = rel(1.15)))

# Combine the three panels vertically
Figure_5 <- p_5A / p_5B / p_5C +
  plot_layout(heights = c(1.5, 0.7, 1.3))

gc()

# Combined figure — PNG and PDF
ggsave(paste0(out, "Figure5_Functional_Overview.png"),
       Figure_5, width = 15, height = 18, dpi = 400)

ggsave(paste0(out, "Figure5_Functional_Overview.pdf"),
       Figure_5, width = 15, height = 18)

# Individual PDFs per panel
ggsave(paste0(out, "Figure5A_heatmap.pdf"), p_5A, width = 16, height = 9)
ggsave(paste0(out, "Figure5B_barplot.pdf"), p_5B, width = 16, height = 7)
ggsave(paste0(out, "Figure5C_sankey.pdf"),  p_5C, width = 22, height = 8)

cat("✔ Figure 5 saved successfully\n")
cat("  Combined PNG: Figure5_Functional_Overview.png\n")
cat("  Combined PDF: Figure5_Functional_Overview.pdf\n")
cat("  Individual PDFs: Figure5A_heatmap.pdf / Figure5B_barplot.pdf / Figure5C_sankey.pdf\n")


############################################################
# SUPPLEMENTARY FIGURE S1 — FULL LONG HEATMAP
############################################################

clases_supp     <- unique(row_ann_pub$class1)
class_cols_supp <- setNames(
  class_palette[seq_along(clases_supp)],
  clases_supp
)

row_ann_supp <- data.frame(
  Class = row_ann_pub$class1,
  row.names = rownames(row_ann_pub)
)

ann_colors_supp <- list(Class = class_cols_supp)
alto_supp       <- max(16, 4 + n_paths_pub * 0.20)

for (ext in c("png", "pdf")) {
  
  if (ext == "png") {
    png(paste0(out, "FigureS1_Heatmap_Complete.png"),
        width  = 20,
        height = alto_supp,
        units  = "in",
        res    = 300)
  } else {
    pdf(paste0(out, "FigureS1_Heatmap_Complete.pdf"),
        width  = 20,
        height = alto_supp)
  }
  
  pheatmap(
    heat_mat_z_pub,
    annotation_row    = row_ann_supp,
    annotation_colors = ann_colors_supp,
    color             = heat_colors,
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = 9,
    fontsize_col      = 14,
    fontsize          = 13,
    border_color      = NA,
    main              = "Supplementary S1  |  All significant pathways — Strain, Strain × Time and Body_weight effects (DESeq2, q < 0.05)"
  )
  
  dev.off()
}

cat("✔ FigureS1 saved\n")

############################################################
# SUPPLEMENTARY FIGURE S2 — NETWORK
############################################################

theme_network <- theme_void(base_size = 18) +
  theme(
    plot.title      = element_text(size = 18, face = "bold",
                                   hjust = 0.5, margin = margin(b = 8)),
    plot.background = element_rect(fill = "white", color = NA),
    legend.text     = element_text(size = 14),
    legend.title    = element_text(size = 15, face = "bold"),
    legend.position = "bottom"
  )

net_list_final <- lapply(levels(metadata$Time_point), function(tp) {
  
  df_tp <- net_df %>% dplyr::filter(Time_point == tp)
  if (nrow(df_tp) == 0) return(NULL)
  
  plots_strain <- lapply(levels(metadata$Strain), function(st) {
    
    df_st <- df_tp %>% dplyr::filter(Strain == st)
    if (nrow(df_st) == 0) return(NULL)
    
    g <- igraph::graph_from_data_frame(
      df_st %>% dplyr::select(DisplayName, pathway_name, taxon_function_abun),
      directed = FALSE
    )
    
    V(g)$type <- ifelse(V(g)$name %in% df_st$DisplayName,
                        "Bacteria", "Pathway")
    
    pathway_sig <- df_st %>%
      dplyr::select(pathway_name, Significance) %>%
      dplyr::distinct()
    
    V(g)$color <- ifelse(
      V(g)$type == "Bacteria",
      strain_colors[st],
      ifelse(
        pathway_sig$Significance[
          match(V(g)$name, pathway_sig$pathway_name)
        ] == "Significant",
        "#FFAA80",
        "#E0E0E0"
      )
    )
    
    V(g)$size <- ifelse(V(g)$type == "Bacteria", 10, 7)
    
    ggraph(g, layout = "kk") +
      geom_edge_link(aes(width = taxon_function_abun),
                     color = "#BBBBBB", alpha = 0.5) +
      geom_node_point(aes(size = size), color = V(g)$color) +
      geom_node_text(
        aes(label = name),
        repel         = TRUE,
        size          = 5.5,
        fontface      = "bold",
        box.padding   = 0.6,
        point.padding = 0.4,
        color         = "#222222"
      ) +
      scale_size(range = c(6, 14)) +
      theme_network +
      labs(title = paste0(gsub("_", " ", tp), " — ", st))
  })
  
  plots_strain <- Filter(Negate(is.null), plots_strain)
  patchwork::wrap_plots(plots_strain, nrow = 1) +
    patchwork::plot_annotation(
      title = paste0("Time point: ", gsub("_", " ", tp)),
      theme = theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
      )
    )
})

net_list_final <- Filter(Negate(is.null), net_list_final)

p_network_final_improved <- patchwork::wrap_plots(
  net_list_final, ncol = 1
) +
  patchwork::plot_annotation(
    title    = "Supplementary S2  |  Taxa–pathway network by time point and strain",
    subtitle = "Nodes: bacteria (strain color) and pathways (orange = significant q < 0.05)\nEdge width proportional to taxon functional contribution",
    theme    = theme(
      plot.title    = element_text(size = 22, face = "bold",    hjust = 0.5),
      plot.subtitle = element_text(size = 16, color = "grey40", hjust = 0.5)
    )
  )

ggsave(paste0(out, "FigureS2_Network_Taxa_Pathways.png"),
       p_network_final_improved,
       width  = 28,
       height = 10 * length(net_list_final),
       dpi    = 300)

ggsave(paste0(out, "FigureS2_Network_Taxa_Pathways.pdf"),
       p_network_final_improved,
       width  = 28,
       height = 10 * length(net_list_final))

cat("✔ FigureS2 saved\n")
cat("\n=== FINAL SUMMARY ===\n")
cat("  Figure5:   Heatmap Z-score (A) + Class barplot (B) + Sankey (C)\n")
cat("  FigureS1:  Complete heatmap, all pathways\n")
cat("  FigureS2:  Network taxa-pathways\n")
