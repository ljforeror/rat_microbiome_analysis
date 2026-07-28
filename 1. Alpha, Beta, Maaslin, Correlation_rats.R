################################################################################
######################## COMPLETE MICROBIOME ANALYSIS ##########################
################################################################################

################################################################################
############################ 1. LIBRARIES ######################################
################################################################################

.libPaths(c("/projappl/project_2007408/project_rpackages_Microbiome16srat", .libPaths()))

library(phyloseq)
library(ape)
library(microbiome)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggpubr)
library(effectsize)
library(patchwork)

################################################################################
############################ 2. DIRECTORIES ####################################
################################################################################

base_path <- "/scratch/project_2007408/Microbiome_Mice/Microbiome_Endotarget_mice/I"
analysis_path <- paste0(base_path,"/R_analysis_1")

Figures  <- paste0(analysis_path,"/New_Results/")
Tablesen <- paste0(analysis_path,"/New_Results/")

dir.create(Figures, showWarnings = FALSE, recursive = TRUE)

################################################################################
############################ 3. LOAD METADATA ##################################
################################################################################
# ── Clinical variables ───────────────────────────────────────────────
datos_variables <- read.table(
  "/scratch/project_2007408/Microbiome_Mice/Microbiome_Endotarget_mice/I/R_analysis_1/df_completo_4.csv",
  header           = TRUE,
  sep              = ",",
  stringsAsFactors = FALSE
)

metadata <- read.csv(paste0(analysis_path,"/Metadata.csv"), row.names = 1)

metadata$Time_point <- factor(metadata$Time_point,
                              levels = c("Baseline","Midpoint","Endpoint"),
                              labels = c("0 weeks","12 weeks","24 weeks"))

metadata$Surgery.status <- factor(metadata$Surgery.status,
                                  levels = c("sham","OA groove"),
                                  labels = c("Sham","Groove"))

metadata$Strain <- as.factor(metadata$Strain)
metadata$Rat.ID <- as.factor(metadata$Rat.ID)

################################################################################
############################ 4. LOAD PHYLOSEQ ##################################
################################################################################

taxa <- read.csv(paste0(base_path,"/Results/_taxa_final_1_mod.csv"), row.names = 1)
tax_mat <- as.matrix(taxa)

seq <- read.csv(paste0(base_path,"/Results/_taxa_asv_final_1.csv"), row.names = 1)
otu_mat <- as.matrix(t(seq))

tree <- read.tree(paste0(base_path,"/Results/GTR.phy"))

taxa1 <- read.csv(paste0(base_path,"/Results/_taxa_final_1.csv"), row.names = 1)
name_dict <- setNames(rownames(taxa1), taxa1$seq)
tree$tip.label <- name_dict[tree$tip.label]

OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
SAMPLES <- sample_data(metadata)

ps <- phyloseq(OTU, TAX, SAMPLES, phy_tree(tree))

################################################################################
############################ 5. FILTERING ######################################

############################################################
# PIPELINE COMPLETE— MICROBIOME ANALYSIS
# Taxonomic composition, Alpha/Beta diversity, MaAsLin2
############################################################

library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(tibble)
library(Maaslin2)
library(lme4)
library(lmerTest)
library(rstatix)
library(ggpubr)
library(car)

dir.create("Results", showWarnings = FALSE, recursive = TRUE)
Figures  <- "Results/"
Tablesen <- "Results/"

############################################################
# COLORS
############################################################

strain_colors <- c("SD" = "#7B61FF", "Wi" = "#FF4FA3")
surgery_colors <- c("Sham" = "#7B61FF", "Groove" = "#FF4FA3")
coef_colors    <- c("Positive" = "#D73027", "Negative" = "#4575B4")
gradient_low   <- "#00CFE8"
gradient_mid   <- "white"
gradient_high  <- "#007BFF"

############################################################
# 🦠FILTERING
############################################################

ps <- subset_taxa(ps, !is.na(Phylum) & Phylum != "" & Phylum != "uncharacterized")

prev      <- apply(otu_table(ps), 1, function(x) sum(x > 0))
keep_taxa <- names(prev[prev >= 0.05 * nsamples(ps)])
ps        <- prune_taxa(keep_taxa, ps)
ps        <- subset_samples(ps, Diet == "High Fat/Sucrose")

############################################################
# 🦠 RAREFACTION
############################################################

set.seed(123)
ps_rarefied <- rarefy_even_depth(ps,
                                 sample.size = min(sample_sums(ps)),
                                 rngseed     = 123,
                                 replace     = FALSE,
                                 verbose     = FALSE)

ps_rarefied <- subset_samples(ps_rarefied, Diet == "High Fat/Sucrose")

############################################################
# 🦠 TAXONOMIC COMPOSITION
############################################################

# Agglomerate
ps_genus  <- tax_glom(ps, taxrank = "Genus")
ps_family <- tax_glom(ps, taxrank = "Family")
ps_phylum <- tax_glom(ps, taxrank = "Phylum")

cat("ASVs after filtering:    ", ntaxa(ps),       "\n")
cat("Genera after glom:       ", ntaxa(ps_genus),  "\n")
cat("Families after glom:     ", ntaxa(ps_family), "\n")
cat("Phyla after glom:        ", ntaxa(ps_phylum), "\n")

# Relative abundance
ps_genus_rel  <- transform_sample_counts(ps_genus,  function(x) x / sum(x))
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))

genus_df  <- psmelt(ps_genus_rel)
phylum_df <- psmelt(ps_phylum_rel)

# Top taxa
top_phyla  <- phylum_df %>%
  group_by(Phylum) %>%
  summarise(total = mean(Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 7) %>%
  pull(Phylum)

top_genera <- genus_df %>%
  group_by(Genus) %>%
  summarise(total = mean(Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

phylum_df$Phylum <- ifelse(phylum_df$Phylum %in% top_phyla, phylum_df$Phylum, "Other")
genus_df$Genus   <- ifelse(genus_df$Genus   %in% top_genera, genus_df$Genus,  "Other")

# Summaries
phylum_summary <- phylum_df %>%
  group_by(Strain, Phylum) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop")

genus_summary_tp <- genus_df %>%
  filter(Genus %in% top_genera) %>%
  group_by(Strain, Time_point, Genus) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  mutate(percentage = round(mean_abundance * 100, 2))

write.csv(phylum_summary,  paste0(Tablesen, "Table_Phylum_RelAbundance.csv"),   row.names = FALSE)
write.csv(genus_summary_tp, paste0(Tablesen, "Table_Genus_Top10_pct.csv"),       row.names = FALSE)

############################################################
# 🦠 FIGURE 1 — TAXONOMIC COMPOSITION (A + B)
############################################################
phylum_colors <- c(
  "Firmicutes"                  = "#E63946",  # rojo
  "Bacteroidetes"               = "#1D7874",  # verde azulado oscuro
  "Actinobacteria"               = "#F4A300",  # naranja/mostaza
  "Proteobacteria"               = "#3A6EA5",  # azul
  "Verrucomicrobia"              = "#8E44AD",  # púrpura
  "Tenericutes"                  = "#2A9D8F",  # verde esmeralda
  "Candidatus Saccharibacteria"  = "#D4A5A5",  # rosa arcilla
  "Other"                        = "#B0B0B0"   # gris neutro
)
phylum_colors <- c(
    "Firmicutes"                  = "#F4A6A6",  # rosa-rojo pastel
    "Bacteroidetes"               = "#8BE6DF",  # turquesa pastel (el tuyo original)
    "Actinobacteria"               = "#FFD8A8",  # naranja pastel
    "Proteobacteria"               = "#A8C6FA",  # azul pastel
    "Verrucomicrobia"              = "#C9A7EB",  # púrpura pastel
    "Tenericutes"                  = "#B4E6B0",  # verde pastel
    "Candidatus Saccharibacteria"  = "#F9E2AE",  # amarillo pastel (el tuyo original)
    "Other"                        = "#DDDDDD"   # gris neutro
)
  
genus_base_colors <- c(
    "#F4A6A6", "#8BE6DF", "#FFD8A8", "#A8C6FA", "#C9A7EB",
    "#B4E6B0", "#F9E2AE", "#F7B6D2", "#A0D9D9", "#E6C3A0",
    "#B9A7FF", "#9FD8B0", "#F4C7C3", "#AEE0F5", "#D9C48A"
)
  
genus_palette <- setNames(
    genus_base_colors[1:length(unique(genus_df$Genus))],
    unique(genus_df$Genus)
)
  
genus_palette["Other"] <- "#DDDDDD"



# Reordenar los niveles para que "Other" quede al final (abajo del stack)
phylum_df <- phylum_df %>%
  mutate(Phylum = factor(Phylum, levels = c(setdiff(names(phylum_colors), "Other"), "Other")))

genus_df <- genus_df %>%
  mutate(Genus = factor(Genus, levels = c(setdiff(unique(Genus), "Other"), "Other")))

p_phylum <- ggplot(phylum_df, aes(Time_point, Abundance, fill = Phylum)) +
  stat_summary(fun = mean, geom = "bar", position = position_stack(reverse = TRUE), color = "black", linewidth = 0.2) +
  #stat_summary(fun = mean, geom = "bar", position = "stack", color = "black", linewidth = 0.2) +
  facet_wrap(~Strain) +
  scale_fill_manual(values = phylum_colors) +
  scale_y_continuous(labels = percent_format()) +
  theme_classic(base_size = 20) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 22, face = "bold"),
        legend.title = element_blank()) +
  labs(title = "A  |  Phylum-level composition",
       y = "Mean relative abundance (%)", x = "Time point")


p_genus <- ggplot(genus_df, aes(Time_point, Abundance, fill = Genus)) +
  stat_summary(fun = mean, geom = "bar", position = position_stack(reverse = TRUE), color = "black", linewidth = 0.2) +
  #stat_summary(fun = mean, geom = "bar", position = "stack", color = "black", linewidth = 0.2) +
  facet_wrap(~Strain) +
  scale_fill_manual(values = genus_palette) +
  scale_y_continuous(labels = percent_format()) +
  theme_classic(base_size = 20) +
  theme(strip.text = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 22, face = "bold"),
        legend.title = element_blank()) +
  labs(title = "B  |  Top 10 genera",
       y = "Mean relative abundance (%)", x = "Time point")

figure1_composition <- p_phylum / p_genus +
  plot_annotation(
    title = "Strain-specific gut microbiome composition across time",
    theme = theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5))
  )

ggsave(paste0(Figures, "Figure_1_Taxonomic_Composition.png"),
       figure1_composition, width = 18, height = 17, dpi = 300)
ggsave(paste0(Figures, "Figure_1_Taxonomic_Composition.pdf"),
       figure1_composition, width = 18, height = 17)

cat("✔ Figure 1 (Taxonomic composition) done\n")

# --- Pooled phylum plot ---
phylum_pooled <- phylum_df %>%
  group_by(Strain, Phylum) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

p_phylum_pooled <- ggplot(phylum_pooled, aes(x = Strain, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), color = "black", linewidth = 0.2) +
  scale_fill_manual(values = phylum_colors) +
  scale_y_continuous(labels = percent_format()) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(size = 22, face = "bold"),
        legend.title = element_blank()) +
  labs(title = "Pooled phylum-level composition", y = "Mean relative abundance (%)", x = "Strain")

# --- Pooled genus plot ---
genus_pooled <- genus_df %>%
  group_by(Strain, Genus) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

p_genus_pooled <- ggplot(genus_pooled, aes(x = Strain, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), color = "black", linewidth = 0.2) +
  scale_fill_manual(values = genus_palette) +
  scale_y_continuous(labels = percent_format()) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(size = 22, face = "bold"),
        legend.title = element_blank()) +
  labs(title = "Pooled top 10 genera", y = "Mean relative abundance (%)", x = "Strain")

figure_pooled <- p_phylum_pooled / p_genus_pooled +
  plot_annotation(
    title = "Pooled gut microbiome composition (all time points combined)",
    theme = theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5))
  )

ggsave(paste0(Figures, "Figure_1B_Pooled_Composition.png"),
       figure_pooled, width = 14, height = 14, dpi = 300)
ggsave(paste0(Figures, "Figure_1B_Pooled_Composition.pdf"),
       figure_pooled, width = 14, height = 14)
############################################################
# 🦠 ALPHA DIVERSITY
############################################################
alpha_df <- estimate_richness(ps_rarefied,
                              measures = c("Observed", "Chao1", "Shannon", "Simpson")) %>%
  mutate(
    Time_point     = sample_data(ps_rarefied)$Time_point,
    Strain         = sample_data(ps_rarefied)$Strain,
    Surgery.status = sample_data(ps_rarefied)$Surgery.status,
    Rat.ID         = sample_data(ps_rarefied)$Rat.ID
  )
alpha_df$Strain     <- factor(alpha_df$Strain, levels = c("SD", "Wi"))
alpha_df$Time_point <- factor(alpha_df$Time_point)

alpha_df <- alpha_df %>%
  mutate(Pielou_evenness = Shannon / log(Observed))

# ── Measures: Shannon principal + resto suplementarios ──────────────────────
alpha_measures_all  <- c("Shannon", "Observed", "Chao1", "Simpson", "Pielou_evenness")
alpha_measures_supp <- c("Observed", "Chao1", "Simpson", "Pielou_evenness")   # para figuras suplementarias

# ── Mixed models (todos los índices) ────────────────────────────────────────
alpha_mixed_results <- dplyr::bind_rows(lapply(alpha_measures_all, function(measure) {
  model     <- lmer(as.formula(paste0(measure, " ~ Time_point * Strain + (1|Rat.ID)")),
                    data = alpha_df)
  anova_res           <- as.data.frame(anova(model))
  anova_res$Term      <- rownames(anova_res)
  anova_res$Measure   <- measure
  anova_res
}))

# FDR sobre p-valores del modelo mixto
alpha_mixed_results$p.FDR <- p.adjust(alpha_mixed_results[["Pr(>F)"]], method = "BH")

write.csv(alpha_mixed_results,
          paste0(Tablesen, "TableS2_AlphaDiversity_MixedModels_FDR.csv"),
          row.names = FALSE)

# ── Modelo Shannon principal ─────────────────────────────────────────────────
model_shannon <- lmer(Shannon ~ Time_point * Strain + (1|Rat.ID), data = alpha_df)
anova_shannon <- anova(model_shannon)
interaction_p <- anova_shannon["Time_point:Strain", "Pr(>F)"]

# ── Tabla S1: Wilcoxon no pareado por tiempo (SD vs Wi en cada time-point) ──
wilcox_by_time <- dplyr::bind_rows(lapply(alpha_measures_all, function(measure) {
  dplyr::bind_rows(lapply(levels(alpha_df$Time_point), function(tp) {
    sub  <- alpha_df %>% filter(Time_point == tp)
    grp1 <- sub[[measure]][sub$Strain == "SD"]
    grp2 <- sub[[measure]][sub$Strain == "Wi"]
    wt   <- wilcox.test(grp1, grp2, exact = FALSE)
    data.frame(
      Measure    = measure,
      Time_point = tp,
      n_SD       = length(grp1),
      n_Wi       = length(grp2),
      W          = wt$statistic,
      p.value    = wt$p.value
    )
  }))
}))

# FDR dentro de cada índice
wilcox_by_time <- wilcox_by_time %>%
  group_by(Measure) %>%
  mutate(p.FDR = p.adjust(p.value, method = "BH")) %>%
  ungroup()

write.csv(wilcox_by_time,
          paste0(Tablesen, "TableS1_AlphaDiversity_Wilcoxon_byTime_FDR.csv"),
          row.names = FALSE)

# ── Surgery effect at 24 weeks (Tabla S3) ────────────────────────────────────
alpha_24 <- alpha_df %>% filter(Time_point == "24 weeks")

alpha_surgery_results <- dplyr::bind_rows(lapply(alpha_measures_all, function(measure) {
  model     <- lm(as.formula(paste0(measure, " ~ Strain * Surgery.status")),
                  data = alpha_24)
  anova_res           <- as.data.frame(anova(model))
  anova_res$Term      <- rownames(anova_res)
  anova_res$Measure   <- measure
  anova_res
}))

# FDR sobre los p-valores del modelo lineal (24 weeks)
alpha_surgery_results$p.FDR <- p.adjust(alpha_surgery_results[["Pr(>F)"]], method = "BH")

write.csv(alpha_surgery_results,
          paste0(Tablesen, "TableS3_AlphaDiversity_24weeks_Surgery_FDR.csv"),
          row.names = FALSE)

############################################################
# 🦠 ALPHA FIGURES
############################################################

# ── Etiquetas Wilcoxon para figura principal (Shannon) ──────────────────────
shannon_wilcox_labels <- wilcox_by_time %>%
  filter(Measure == "Shannon") %>%
  mutate(label = case_when(
    p.FDR < 0.001 ~ "***",
    p.FDR < 0.01  ~ "**",
    p.FDR < 0.05  ~ "*",
    TRUE          ~ paste0("p=", signif(p.FDR, 2))
  ))

# ── Figura principal: Shannon ────────────────────────────────────────────────
p_alpha_main <- ggplot(alpha_df, aes(Time_point, Shannon, fill = Strain)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA,
               color = "black", linewidth = 0.7) +
  geom_jitter(shape = 21, aes(fill = Strain), color = "black",
              stroke = 0.7, size = 3.4, width = 0.12, alpha = 0.95) +
  stat_compare_means(
    aes(group = Strain),
    method = "wilcox.test",
    label  = "p.format",
    size   = 4,
    paired = FALSE
  ) +
  scale_fill_manual(values = strain_colors) +
  theme_classic(base_size = 14) +   # reducido de 18 a 14
  theme(
    aspect.ratio      = 1.2,        # más alto que ancho
    legend.position   = c(0.88, 0.80),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.key.size   = unit(0.5, "cm"),
    legend.text       = element_text(size = 11),
    plot.title        = element_text(size = 13, hjust = 0),
    plot.subtitle     = element_text(size = 9,  hjust = 0, color = "grey40"),
    axis.text         = element_text(size = 11),
    axis.title        = element_text(size = 12),
    plot.margin       = margin(5, 5, 5, 5)
  ) +
  labs(
    title    = "A  |  Shannon diversity",
    subtitle = paste0("Mixed model interaction p = ", signif(interaction_p, 3),
                      " | Wilcoxon FDR-adjusted (BH)"),
    x = "Time",
    y = "Shannon diversity index"
  )

# ── Figuras suplementarias: Observed, Chao1, Simpson ────────────────────────
make_supp_alpha <- function(measure) {
  ggplot(alpha_df, aes(Time_point, .data[[measure]], fill = Strain)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, color = "black", linewidth = 0.6) +
    geom_jitter(shape = 21, aes(fill = Strain), color = "black",
                stroke = 0.6, size = 3, width = 0.12, alpha = 0.9) +
    stat_compare_means(
      aes(group = Strain),
      method = "wilcox.test",
      label  = "p.signif",
      paired = FALSE
    ) +
    scale_fill_manual(values = strain_colors) +
    theme_classic(base_size = 15) +
    labs(title = paste0("Supplementary – ", measure),
         x = "Time", y = measure)
}

p_observed <- make_supp_alpha("Observed")
p_chao1    <- make_supp_alpha("Chao1")
p_simpson  <- make_supp_alpha("Simpson")
p_Pielou_evenness  <- make_supp_alpha("Pielou_evenness")

# ── Save figures ──────────────────────────────────────────────────────────
ggsave(paste0(Figures, "Figure_Alpha_Shannon.png"),
       p_alpha_main, width = 8, height = 6, dpi = 300)

ggsave(paste0(Figures, "FigureS_Alpha_Observed_Chao1_Simpson.png"),
       p_observed / p_chao1 / p_simpson / p_Pielou_evenness,
       width = 8, height = 14, dpi = 300)
############################################################
############################################################
# 🦠 BETA DIVERSITY
############################################################
ps_clr         <- microbiome::transform(ps, "clr")
dist_aitchison <- phyloseq::distance(ps_clr, method = "euclidean")
meta           <- as(sample_data(ps_clr), "data.frame")
meta           <- meta[match(rownames(meta), rownames(as.matrix(dist_aitchison))), ]

# Orden explícito de time points
meta$Time_point <- factor(meta$Time_point,
                          levels = c("0 weeks", "12 weeks", "24 weeks"))

# ── PERMANOVA global ──────────────────────────────────────────────────────────
# Correr por separado sin interacción para obtener efectos principales
adonis_time_strain <- adonis2(dist_aitchison ~ Strain + Time_point,
                              data         = meta,
                              permutations = 9999,
                              strata       = meta$Rat.ID,
                              by           = "margin")

# Interacción por separado
adonis_int <- adonis2(dist_aitchison ~ Strain * Time_point,
                      data         = meta,
                      permutations = 9999,
                      strata       = meta$Rat.ID,
                      by           = "terms")

print(adonis_time_strain)
print(adonis_int)
# Construir tabla combinada correctamente
adonis_combined <- data.frame(
  Term     = c("Time_point", "Strain", "Time_point:Strain"),
  Df       = c(adonis_time_strain$Df[1],     # Time_point
               adonis_time_strain$Df[2],     # Strain
               adonis_int$Df[3]),            # Interacción
  SumOfSqs = c(adonis_time_strain$SumOfSqs[1],
               adonis_time_strain$SumOfSqs[2],
               adonis_int$SumOfSqs[3]),
  R2       = c(adonis_time_strain$R2[1],
               adonis_time_strain$R2[2],
               adonis_int$R2[3]),
  F        = c(adonis_time_strain$F[1],
               adonis_time_strain$F[2],
               adonis_int$F[3]),
  p.value  = c(adonis_time_strain$`Pr(>F)`[1],   # margin
               adonis_time_strain$`Pr(>F)`[2],   # margin
               adonis_int$`Pr(>F)`[3])            # terms
)
adonis_combined$p.FDR <- p.adjust(adonis_combined$p.value, method = "BH")

print(adonis_combined)

write.csv(adonis_combined,
          paste0(Tablesen, "TableS4_PERMANOVA_Global_FDR.csv"),
          row.names = FALSE)

# Extraer valores
R2_strain <- adonis_combined$R2[adonis_combined$Term == "Strain"]
F_strain  <- adonis_combined$F[adonis_combined$Term == "Strain"]
p_strain  <- adonis_combined$p.value[adonis_combined$Term == "Strain"]

R2_time   <- adonis_combined$R2[adonis_combined$Term == "Time_point"]
F_time    <- adonis_combined$F[adonis_combined$Term == "Time_point"]
p_time    <- adonis_combined$p.value[adonis_combined$Term == "Time_point"]

R2_int    <- adonis_combined$R2[adonis_combined$Term == "Time_point:Strain"]
F_int     <- adonis_combined$F[adonis_combined$Term == "Time_point:Strain"]
p_int     <- adonis_combined$p.value[adonis_combined$Term == "Time_point:Strain"]

# ── Post-hoc: PERMANOVA Strain dentro de cada time point ─────────────────────
tp_levels  <- levels(meta$Time_point)   # "Baseline", "12 weeks", "24 weeks"

results_tp <- dplyr::bind_rows(lapply(tp_levels, function(tp) {
  idx      <- meta$Time_point == tp
  meta_sub <- meta[idx, ]
  dist_sub <- as.dist(as.matrix(dist_aitchison)[idx, idx])
  ad       <- adonis2(dist_sub ~ Strain,
                      data         = meta_sub,
                      permutations = 9999,
                      by           = "margin")
  data.frame(
    Time_point = tp,
    Df         = ad$Df[1],
    SumOfSqs   = round(ad$SumOfSqs[1], 3),
    F          = round(ad$F[1], 3),
    R2         = round(ad$R2[1], 3),
    p.value    = ad$`Pr(>F)`[1]
  )
}))
results_tp$p.FDR <- p.adjust(results_tp$p.value, method = "BH")

write.csv(results_tp,
          paste0(Tablesen, "TableS4b_PERMANOVA_PostHoc_byTimepoint.csv"),
          row.names = FALSE)

# ── BetaDisper ────────────────────────────────────────────────────────────────
disp_strain      <- betadisper(dist_aitchison, meta$Strain)
disp_time        <- betadisper(dist_aitchison, meta$Time_point)
disp_interaction <- betadisper(dist_aitchison,
                               interaction(meta$Strain, meta$Time_point))

betadisper_results <- rbind(
  data.frame(as.data.frame(anova(disp_strain)),      Factor = "Strain"),
  data.frame(as.data.frame(anova(disp_time)),         Factor = "Time_point"),
  data.frame(as.data.frame(anova(disp_interaction)),  Factor = "Strain x Time_point")
)
write.csv(betadisper_results,
          paste0(Tablesen, "TableS5_BetaDisper_Homogeneity.csv"),
          row.names = TRUE)

# Extraer p-values betadisper para párrafo
p_disp_strain <- anova(disp_strain)$`Pr(>F)`[1]
p_disp_time   <- anova(disp_time)$`Pr(>F)`[1]
p_disp_int    <- anova(disp_interaction)$`Pr(>F)`[1]

# ── PERMANOVA Surgery 24 weeks ────────────────────────────────────────────────
meta_24 <- meta[meta$Time_point == "24 weeks", ]
dist_mat <- as.matrix(dist_aitchison)
dist_24  <- as.dist(dist_mat[meta$Time_point == "24 weeks",
                             meta$Time_point == "24 weeks"])

adonis_surgery_24    <- adonis2(dist_24 ~ Strain * Surgery.status,
                                data         = meta_24,
                                permutations = 9999,
                                by           = "margin")
adonis_surgery_df        <- as.data.frame(adonis_surgery_24)
adonis_surgery_df$Term   <- rownames(adonis_surgery_df)
adonis_surgery_df$p.FDR  <- p.adjust(adonis_surgery_df$`Pr(>F)`, method = "BH")

write.csv(adonis_surgery_df,
          paste0(Tablesen, "TableS6_PERMANOVA_24weeks_Surgery_FDR.csv"),
          row.names = FALSE)

R2_surg_int <- round(adonis_surgery_df$R2[adonis_surgery_df$Term == "Strain:Surgery.status"], 3)
p_surg_int  <- adonis_surgery_df$`Pr(>F)`[adonis_surgery_df$Term == "Strain:Surgery.status"]
fdr_surg_int <- adonis_surgery_df$p.FDR[adonis_surgery_df$Term == "Strain:Surgery.status"]

# ── PCoA ordination ───────────────────────────────────────────────────────────
ordination <- ordinate(ps_clr, method = "PCoA", distance = dist_aitchison)

eig        <- ordination$values$Eigenvalues
pct_var    <- round(100 * eig / sum(eig[eig > 0]), 1)
pcoa1_var  <- pct_var[1]
pcoa2_var  <- pct_var[2]
total_var2 <- round(pcoa1_var + pcoa2_var, 1)

# Scores para ggplot manual (más control que plot_ordination)
pcoa_scores <- as.data.frame(ordination$vectors[, 1:2])
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")
pcoa_scores$Strain        <- meta[rownames(pcoa_scores), "Strain"]
pcoa_scores$Time_point    <- meta[rownames(pcoa_scores), "Time_point"]
pcoa_scores$Surgery.status <- meta[rownames(pcoa_scores), "Surgery.status"]

############################################################
# 🦠 FIGURES BETA DIVERSITY
############################################################
# Forzar conversión total
pcoa_scores <- data.frame(
  PCoA1          = as.numeric(ordination$vectors[, 1]),
  PCoA2          = as.numeric(ordination$vectors[, 2]),
  Strain         = as.character(meta[rownames(ordination$vectors), "Strain"]),
  Time_point     = as.character(meta[rownames(ordination$vectors), "Time_point"]),
  Surgery.status = as.character(meta[rownames(ordination$vectors), "Surgery.status"]),
  row.names      = rownames(ordination$vectors),
  stringsAsFactors = FALSE
)

# Verificar que no hay listas
cat("Clases de columnas:\n")
print(sapply(pcoa_scores, class))

# Solución sin dplyr::slice
hulls_main <- do.call(rbind, lapply(unique(pcoa_scores$Strain), function(s) {
  df_s <- pcoa_scores[pcoa_scores$Strain == s, ]
  df_s[chull(df_s$PCoA1, df_s$PCoA2), ]
}))

cat("Hulls calculados:\n")
print(table(hulls_main$Strain))
# ── Figura principal: PCoA global (todos los TPs, color = Strain) ────────────

# Subtítulo con resultados PERMANOVA globales
subtitle_main <- paste0(
  "PERMANOVA — Strain: F=", F_strain, ", R²=", R2_strain, ", p=", p_strain,
  " | Time: F=", F_time,   ", R²=", R2_time,   ", p=", p_time,
  " | Interaction: F=", F_int, ", R²=", R2_int, ", p=", p_int
)

p_beta_main <- ggplot(pcoa_scores, aes(PCoA1, PCoA2,
                                       fill  = Strain,
                                       color = Strain)) +
  geom_polygon(data  = hulls_main,
               aes(group = Strain), alpha = 0.10, color = NA) +
  geom_point(shape = 21, color = "black",
             size = 4, stroke = 0.8, alpha = 0.95) +
  stat_ellipse(aes(group = Strain), type = "t",
               level = 0.95, linewidth = 1.2) +
  scale_fill_manual(values  = strain_colors) +
  scale_color_manual(values = strain_colors) +
  theme_classic(base_size = 18) +
  theme(aspect.ratio    = 1,
        plot.subtitle   = element_text(size = 9, color = "grey40")) +
  labs(title    = "B  |  Beta diversity — Aitchison distance (PCoA)",
       subtitle = subtitle_main,
       x = paste0("PCoA1 (", pcoa1_var, "%)"),
       y = paste0("PCoA2 (", pcoa2_var, "%)"))

ggsave(paste0(Figures, "Figure2b_Beta_PCoA_Main.png"),
       p_beta_main, width = 7, height = 6, dpi = 300)

# ── Figura Suplementaria S1: PCoA faceteado por time point (post-hoc) ─────────
label_fun <- function(x) {
  row <- results_tp[results_tp$Time_point == x, ]
  paste0(x,
         "\nF=", row$F, ", R²=", row$R2,
         "\np=", signif(row$p.value, 2),
         ", FDR=", signif(row$p.FDR, 2))
}

p_beta_facet <- ggplot(pcoa_scores,
                       aes(PCoA1, PCoA2, fill = Strain, color = Strain)) +
  geom_point(shape = 21, color = "black",
             size = 3.5, stroke = 0.7, alpha = 0.9) +
  stat_ellipse(aes(group = Strain), type = "t",
               level = 0.95, linewidth = 1) +
  scale_fill_manual(values  = strain_colors) +
  scale_color_manual(values = strain_colors) +
  facet_wrap(~Time_point,
             nrow      = 1,
             labeller  = labeller(Time_point = label_fun)) +
  theme_classic(base_size = 14) +
  theme(aspect.ratio  = 1,
        strip.text    = element_text(size = 11, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = NA)) +
  labs(title = "Supplementary Figure S1 | Beta diversity by time point (post-hoc PERMANOVA)",
       x = paste0("PCoA1 (", pcoa1_var, "%)"),
       y = paste0("PCoA2 (", pcoa2_var, "%)"))

ggsave(paste0(Figures, "FigureS1_Beta_PCoA_byTimepoint.png"),
       p_beta_facet, width = 13, height = 5, dpi = 300)

# ── Figura Suplementaria S2: PCoA 24w — Surgery status ───────────────────────
pcoa_24 <- pcoa_scores %>% filter(Time_point == "24 weeks")

surgery_colors <- c("Sham" = "#4DAF4A", "OA" = "#E41A1C")  # ajusta si tienes otros niveles

p_beta_surgery <- ggplot(pcoa_24,
                         aes(PCoA1, PCoA2,
                             fill  = Surgery.status,
                             shape = Strain)) +
  geom_point(size = 4, stroke = 0.9,
             color = "black", alpha = 0.95) +
  stat_ellipse(aes(group  = interaction(Strain, Surgery.status),
                   color  = Surgery.status),
               type = "t", level = 0.95,
               linewidth = 1, linetype = "dashed") +
  scale_fill_manual(values  = surgery_colors) +
  scale_color_manual(values = surgery_colors) +
  scale_shape_manual(values = c(SD = 21, Wi = 24)) +
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 1) +
  labs(
    title    = "Supplementary Figure S2 | Beta diversity at 24 weeks by Surgery status",
    subtitle = paste0("PERMANOVA Strain × Surgery: R²=", R2_surg_int,
                      ", p=", p_surg_int,
                      ", FDR=", round(fdr_surg_int, 3)),
    x = paste0("PCoA1 (", pcoa1_var, "%)"),
    y = paste0("PCoA2 (", pcoa2_var, "%)")
  )

ggsave(paste0(Figures, "FigureS2_Beta_PCoA_24w_Surgery.png"),
       p_beta_surgery, width = 7, height = 6, dpi = 300)

############################################################
# 🦠 EFFECT SIZES
############################################################

# ── Alpha: Cohen's d por time point (Shannon) ────────────────────────────────
alpha_effect <- alpha_df %>%
  group_by(Time_point) %>%
  group_modify(~ {
    d <- rstatix::cohens_d(data = .x, formula = Shannon ~ Strain, ci = TRUE)
    tibble(Cohens_d = d$effsize, CI_low = d$conf.low, CI_high = d$conf.high)
  }) %>%
  ungroup()

# ── Beta: R² de adonis_combined ──────────────────────────────────────────────
# Ya no existe adonis_main ni adonis_margin — usamos adonis_combined
R2_strain_beta <- adonis_combined$R2[adonis_combined$Term == "Strain"]
R2_time_beta   <- adonis_combined$R2[adonis_combined$Term == "Time_point"]
R2_int_beta    <- adonis_combined$R2[adonis_combined$Term == "Time_point:Strain"]
p_strain_beta  <- adonis_combined$p.value[adonis_combined$Term == "Strain"]
p_time_beta    <- adonis_combined$p.value[adonis_combined$Term == "Time_point"]
p_int_beta     <- adonis_combined$p.value[adonis_combined$Term == "Time_point:Strain"]

# ── Figura: Cohen's d por time point ─────────────────────────────────────────
p_alpha_effect <- ggplot(alpha_effect, aes(Time_point, Cohens_d)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             linewidth = 0.8, color = "grey40") +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.12, linewidth = 1, color = "#7A3CFF") +
  geom_point(size = 5, shape = 21, fill = "#7A3CFF",
             color = "black", stroke = 0.8) +
  theme_classic(base_size = 17) +
  labs(title = "A  |  Alpha diversity effect size (Shannon)",
       y     = "Cohen's d (Wi − SD)",
       x     = "Time point")

# ── Figura: R² PERMANOVA  ─────────────────────────────────────────
beta_plot_df <- data.frame(
  Effect = c("Strain", "Time point", "Time × Strain"),
  R2     = c(R2_strain_beta, R2_time_beta, R2_int_beta),
  p      = c(p_strain_beta,  p_time_beta,  p_int_beta)
)
beta_plot_df$label <- paste0("R²=", round(beta_plot_df$R2, 3),
                             "\np=", beta_plot_df$p)

p_beta_effect <- ggplot(beta_plot_df, aes(Effect, R2, fill = Effect)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.7, width = 0.6) +
  geom_text(aes(label = label), vjust = -0.3, size = 4.5) +
  scale_fill_manual(values = c(
    "Strain"        = "#EFA6C4",
    "Time point"    = "#7EC8E3",
    "Time × Strain" = "#00CFC1"
  )) +
  theme_classic(base_size = 17) +
  theme(legend.position = "none") +
  ylim(0, 0.30) +
  labs(title = "B  |  Variance explained in beta diversity",
       y     = "PERMANOVA R²",
       x     = NULL)

# ── Figura combined ──────────────────────────────────────────────────────────
p_effect_final <- p_alpha_effect / p_beta_effect +
  plot_annotation(tag_levels = "A")

ggsave(paste0(Figures, "Supplementary_EffectSizes.png"),
       p_effect_final, width = 8, height = 11, dpi = 300)

# ── Supplementary tables ─────────────────────────────────────────────────────
write.csv(
  alpha_effect %>%
    dplyr::rename(Effect_size = Cohens_d, CI_Lower = CI_low, CI_Upper = CI_high) %>%
    mutate(Analysis = "Alpha_Shannon_Cohens_d", Comparison = "Wi − SD"),
  paste0(Tablesen, "Supplementary_Alpha_EffectSizes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    Analysis = "Beta_PERMANOVA",
    Term     = c("Strain", "Time_point", "Time_point:Strain"),
    R2       = c(R2_strain_beta, R2_time_beta, R2_int_beta),
    p.value  = c(p_strain_beta,  p_time_beta,  p_int_beta)
  ),
  paste0(Tablesen, "Supplementary_Beta_EffectSizes.csv"),
  row.names = FALSE
)

############################################################
# 🦠 PREPARE METADATA + MAASLIN
############################################################

datos_variables <- read.table(
  "/scratch/project_2007408/Microbiome_Mice/Microbiome_Endotarget_mice/I/R_analysis_1/df_completo_4.csv",
  header = TRUE, sep = ",", stringsAsFactors = FALSE
)

if ("HOMA-IR" %in% colnames(datos_variables)) {
  colnames(datos_variables)[colnames(datos_variables) == "HOMA-IR"] <- "HOMA.IR"
}

numeric_vars <- c("LPS", "Zonulin", "IL6", "TNFa", "IL1b", "IL10",
                  "HOMA.IR", "Triglicerides", "Weight", "Total.OARSI.score")

for (v in numeric_vars) {
  if (v %in% colnames(datos_variables)) {
    datos_variables[[v]] <- as.numeric(gsub(",", ".", datos_variables[[v]]))
  }
}

datos_variables$Time_point <- dplyr::recode(
  as.character(datos_variables$Time_point),
  "Baseline" = "0 weeks", "Midpoint" = "12 weeks", "Endpoint" = "24 weeks"
)
# Renombrar la columna "weight" en datas_variables
datos_variables <- datos_variables %>%
  dplyr::rename(Body_weight = "Weight")  # ajusta "weight" al nombre EXACTO actual (revisa mayúsculas/minúsculas)

# Phyloseq genus level
physeq_tax <- tax_glom(ps_rarefied, "Genus")

ab_table <- as(otu_table(physeq_tax), "matrix")
if (taxa_are_rows(physeq_tax)) ab_table <- t(ab_table)
ab_df <- as.data.frame(ab_table)

metadata_df1 <- metadata_df <- data.frame(sample_data(physeq_tax)) %>%
  rownames_to_column("Sample") %>%
  left_join(
    datos_variables %>%
      select(Sample, LPS, Zonulin, IL6, TNFa, IL1b, IL10,
             HOMA.IR, Triglicerides, Body_weight, Total.OARSI.score),
    by = "Sample"
  )

common_samples <- intersect(rownames(ab_df), metadata_df$Sample)
ab_df          <- ab_df[common_samples, ]
metadata_df    <- metadata_df[match(common_samples, metadata_df$Sample), ]
rownames(metadata_df) <- metadata_df$Sample

metadata_df$Time_point <- factor(metadata_df$Time_point,
                                 levels = c("0 weeks", "12 weeks", "24 weeks"))
metadata_df$Rat.ID     <- factor(metadata_df$Rat.ID)
metadata_df$Strain     <- factor(metadata_df$Strain)

# Taxonomy map (UNA SOLA VEZ)
tax_table_df <- as.data.frame(tax_table(physeq_tax))
tax_table_df$feature <- rownames(tax_table_df)
tax_table_df$DisplayName <- with(tax_table_df,
                                 ifelse(!is.na(Genus)  & Genus  != "", Genus,
                                        ifelse(!is.na(Family) & Family != "", Family, feature)))

tax_map <- tax_table_df %>% select(feature, Family, Genus)

############################################################
# 🦠 MAASLIN2 MODELS
############################################################

RESULTS_DIR <- file.path("Results", "Maaslin_Models")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

metadata_df <- metadata_df %>%
  mutate(Strain     = as.factor(Strain),
         Time_point = as.factor(Time_point),
         Body_weight     = as.numeric(Body_weight)) %>%
  filter(!is.na(Body_weight))

taxa_table    <- ab_df
common_samples <- intersect(rownames(metadata_df), rownames(taxa_table))
metadata_df    <- metadata_df[common_samples, ]
taxa_table     <- taxa_table[common_samples, ]

# 🦠── Model without body mass ──────────────────────────────────────
Maaslin2(taxa_table, metadata_df,
         output          = file.path(RESULTS_DIR, "No_Body_weight"),
         fixed_effects   = c("Time_point", "Strain"),
         normalization   = "CLR",
         transform       = "LOG",
         min_prevalence  = 0.1,
         reference       = c("Time_point,0 weeks"),
         standardize     = TRUE)

# 🦠 ── Model with Body_mass ──────────────────────────────────────
Maaslin2(taxa_table, metadata_df,
         output          = file.path(RESULTS_DIR, "With_Body_weight"),
         fixed_effects   = c("Time_point", "Strain", "Body_weight"),
         normalization   = "CLR",
         transform       = "LOG",
         min_prevalence  = 0.1,
         reference       = c("Time_point,0 weeks"),
         standardize     = TRUE)

# 🦠── Model Zonulin ─────────────────────────────────────────
meta_z <- metadata_df %>% filter(!is.na(Zonulin))
ab_z   <- ab_df[rownames(meta_z), ]

Maaslin2(ab_z, meta_z,
         output          = file.path(RESULTS_DIR, "Zonulin"),
         fixed_effects   = c("Time_point", "Strain", "Zonulin"),
         normalization   = "CLR",
         transform       = "LOG",
         min_prevalence  = 0.1,
         reference       = c("Time_point,0 weeks"),
         standardize     = TRUE)

# 🦠── Model LPS ─────────────────────────────────────────────
meta_l <- metadata_df %>% filter(!is.na(LPS))
ab_l   <- ab_df[rownames(meta_l), ]

Maaslin2(ab_l, meta_l,
         output          = file.path(RESULTS_DIR, "LPS"),
         fixed_effects   = c("Time_point", "Strain", "LPS"),
         normalization   = "CLR",
         transform       = "LOG",
         min_prevalence  = 0.1,
         reference       = c("Time_point,0 weeks"),
         standardize     = TRUE)

# 🦠── Modelo OA 24w ──────────────────────────────────────────
meta_oa24 <- metadata_df %>% filter(Time_point == "24 weeks")
ab_oa24   <- ab_df[rownames(meta_oa24), ]

Maaslin2(ab_oa24, meta_oa24,
         output          = file.path(RESULTS_DIR, "OA_24weeks"),
         fixed_effects   = c("Strain", "Surgery.status"),
         normalization   = "CLR",
         transform       = "LOG",
         min_prevalence  = 0.1,
         reference       = c("Surgery.status,Sham"),
         standardize     = TRUE)

############################################################
# 🦠🦠🦠🦠🦠LOAD MAASLIN RESULTS
############################################################

res_no_Body_weight   <- read.delim(file.path(RESULTS_DIR, "No_Body_weight/all_results.tsv"))
res_with_Body_weight <- read.delim(file.path(RESULTS_DIR, "With_Body_weight/all_results.tsv"))
res_zonulin     <- read.delim(file.path(RESULTS_DIR, "Zonulin/all_results.tsv"))
res_lps         <- read.delim(file.path(RESULTS_DIR, "LPS/all_results.tsv"))
res_oa          <- read.delim(file.path(RESULTS_DIR, "OA_24weeks/all_results.tsv"))

sig_main    <- res_with_Body_weight %>% filter(qval < 0.05)
sig_zonulin <- res_zonulin     %>% filter(qval < 0.05)
sig_lps     <- res_lps         %>% filter(qval < 0.05)
sig_oa      <- res_oa          %>% filter(qval < 0.05)

# save tables
write.csv(sig_main,    file.path(RESULTS_DIR, "Significant_Main_Model.csv"),    row.names = FALSE)
write.csv(sig_zonulin, file.path(RESULTS_DIR, "Significant_Zonulin_Model.csv"), row.names = FALSE)
write.csv(sig_lps,     file.path(RESULTS_DIR, "Significant_LPS_Model.csv"),     row.names = FALSE)
write.csv(sig_oa,      file.path(RESULTS_DIR, "Significant_OA_Model.csv"),      row.names = FALSE)

write.csv(dplyr::bind_rows(
  sig_main    %>% mutate(Model = "Main"),
  sig_zonulin %>% mutate(Model = "Zonulin"),
  sig_lps     %>% mutate(Model = "LPS"),
  sig_oa      %>% mutate(Model = "OA")
), file.path(RESULTS_DIR, "Table_All_Significant_Models.csv"), row.names = FALSE)

############################################################
# 🦠 ADD TAXONOMY TO MAASLIN RESULTS
############################################################

add_taxonomy <- function(df) {
  df %>%
    left_join(tax_table_df %>% select(feature, DisplayName, Family, Genus),
              by = "feature") %>%
    mutate(DisplayName = ifelse(is.na(DisplayName), feature, DisplayName))
}

sig_main_t    <- add_taxonomy(sig_main)
sig_zonulin_t <- add_taxonomy(sig_zonulin)
sig_lps_t     <- add_taxonomy(sig_lps)
sig_oa_t      <- add_taxonomy(sig_oa)

maaslin_main_sig    <- sig_main_t
maaslin_zonulin_sig <- sig_zonulin_t
maaslin_lps_sig     <- sig_lps_t
maaslin_oa_sig      <- sig_oa_t

############################################################
# 🦠FIGURE — HEATMAP ALL MAASLIN (WITH Body_mass)
############################################################

df_W <- res_with_Body_weight %>%
  filter(qval < 0.05) %>%
  mutate(coef = as.numeric(coef)) %>%
  select(feature, metadata, coef) %>%
  left_join(tax_map, by = "feature") %>%
  mutate(Taxa_label = case_when(
    !is.na(Family) & Family != "" & !is.na(Genus) & Genus != "" ~
      paste0(Family, " | ", Genus, " (", feature, ")"),
    !is.na(Family) & Family != "" ~ paste0(Family, " (", feature, ")"),
    TRUE ~ feature
  ))

taxa_order_W <- df_W %>%
  group_by(Taxa_label) %>%
  summarise(mean_coef = mean(abs(coef), na.rm = TRUE)) %>%
  arrange(desc(mean_coef)) %>%
  pull(Taxa_label)

df_W$Taxa_label <- factor(df_W$Taxa_label, levels = rev(taxa_order_W))

n_taxa_W <- length(taxa_order_W)

p_C <- ggplot(df_W, aes(x = metadata, y = Taxa_label, fill = coef)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = gradient_low, mid = "white",
                       high = gradient_high, name = "Coefficient") +
  theme_minimal(base_size = 16) +
  theme(panel.grid   = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y  = element_text(size = 12)) +
  labs(title = "C  |  Significant taxa (model adjusted for Body_weight)",
       x = "Model term", fill = "Coefficient")

############################################################
# 🦠FIGURE 1 MAIN — A / B / C
############################################################

############################################################
# 🦠FIGURE 1 MAIN — A / B / C
############################################################
p_alpha_main <- p_alpha_main +
  labs(title = "A  |  Shannon diversity")
p_beta_facet <- p_beta_facet +
  labs(title = "B  |  Beta diversity (Aitchison distance)")
p_C <- p_C +
  labs(title = "C  |  Differential abundant taxa")

# 
# ── dont include leyend of A (B have it) ───────────────────────────────────────

# ── Tema base compartido ──────────────────────────────────────────────────────
tema_base <- theme_classic(base_size = 13) +
  theme(
    axis.text      = element_text(size = 11),
    axis.title     = element_text(size = 12),
    axis.title.y   = element_text(size = 12, margin = margin(r = 2)),
    plot.title     = element_text(size = 13, hjust = 0, face = "bold"),
    plot.subtitle  = element_text(size = 9,  hjust = 0, color = "grey40"),
    strip.text     = element_text(size = 9, face = "bold"),
    legend.text    = element_text(size = 9),
    legend.title   = element_text(size = 10),
    plot.margin    = margin(5, 5, 5, 5)
  )

# ── Apply to panel ──────────────────────────────────────────────────────
p_alpha_main <- p_alpha_main + tema_base +
  theme(
    legend.position  = "none",
    axis.title.y     = element_text(size = 12, margin = margin(r = 1, l = -10)),
    plot.margin      = margin(5, 5, 5, -10)  # margen izquierdo negativo
  )

p_beta_facet <- p_beta_facet + tema_base +
  theme(
    aspect.ratio      = NULL,
    legend.position   = c(0.97, 0.75),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.key.size   = unit(0.4, "cm"),
    plot.margin       = margin(5, 5, 5, 0)
  )

p_C <- p_C + tema_base +
  theme(
    legend.key.size = unit(0.4, "cm"),
    plot.margin     = margin(5, 5, 5, 5)
  )

# ── Ensamblar ─────────────────────────────────────────────────────────────────


# ──  patchwork ───────────────────────────────────────────────────
fila_sup <- (p_alpha_main | p_beta_facet) +
  plot_layout(widths = c(1, 1.8))

Figure_1 <- fila_sup / p_C +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title = "Gut microbiome dynamics across strains and time",
    theme = theme(
      plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
  )

ggsave(paste0(Figures, "Figure1_Main_2.png"),
       Figure_1,
       width  = 14,
       height = 10,
       dpi    = 300)
####### end of the figure


p_alpha_main <- p_alpha_main +
  labs(title = "A  |  Shannon diversity")

p_beta_facet <- p_beta_facet +
  labs(title = "B  |  Beta diversity (Aitchison distance)")

p_C <- p_C +
  labs(title = "C  |  Differential abundant taxa")

# fixe proportions — B more space
# ── 
tema_base <- theme_classic(base_size = 13) +
  theme(
    axis.text     = element_text(size = 11),
    axis.title    = element_text(size = 12),
    axis.title.y  = element_text(size = 12, margin = margin(r = 2)),
    plot.title    = element_text(size = 13, hjust = 0, face = "bold"),
    plot.subtitle = element_text(size = 9,  hjust = 0, color = "grey40"),
    plot.margin   = margin(5, 5, 5, 5)
  )

# ── Apply──────────────────────────────────────────────────────
p_alpha_main <- p_alpha_main + tema_base +
  theme(
    legend.position  = "none",
    axis.title.y     = element_text(size = 12, margin = margin(r = 1, l = -10)),
    plot.margin      = margin(5, 5, 5, -10)  # margen izquierdo negativo
  )
p_beta_facet <- p_beta_facet + tema_base +
  theme(
    aspect.ratio      = NULL,
    strip.text        = element_text(size = 9, face = "bold"),
    legend.position   = c(0.97, 0.75),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.key.size   = unit(0.4, "cm"),
    legend.text       = element_text(size = 9),
    plot.margin       = margin(5, 5, 5, 0)
  )

p_C <- p_C + tema_base

# ── Ensemble ─────────────────────────────────────────────────────────────────
# ── patchwork ───────────────────────────────────────────────────
fila_sup <- (p_alpha_main | p_beta_facet) +
  plot_layout(widths = c(1, 1.8))

Figure_1 <- fila_sup / p_C +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title = "Gut microbiome dynamics across strains and time",
    theme = theme(
      plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
  )

ggsave(paste0(Figures, "Figure1_Main.png"),
       Figure_1,
       width  = 14,
       height = 10,
       dpi    = 300)



# ── Short names y C ──────────────────────────────────────────────────────
p_C <- p_C +
  scale_y_discrete(labels = function(x) {
    x <- stringr::str_remove(x, ".*\\| ")    # quita "Familia | "
    stringr::str_remove(x, " \\(ASV.*\\)")   # quita "(ASV_XXX)"
  }) +
  tema_base +
  theme(
    axis.text.y  = element_text(size = 12),
    axis.text.x  = element_text(size = 12),
    axis.title   = element_text(size = 13),
    plot.title   = element_text(size = 14, hjust = 0, face = "bold"),
    legend.text  = element_text(size = 10),
    legend.title = element_text(size = 11)
  )

# ── Ensemble ─────────────────────────────────────────────────────────────────
fila_sup <- (p_alpha_main | p_beta_facet) +
  plot_layout(widths = c(1, 1.8))

Figure_1 <- fila_sup / p_C +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title = "Gut microbiome dynamics across strains and time",
    theme = theme(
      plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
  )

ggsave(paste0(Figures, "Figure1_Main.png"),
       Figure_1,
       width  = 14,
       height = 10,
       dpi    = 300)




#
# 🦠 VIF depends on the matrix of predictors not on the answer
################################################################
modelo_vif <- lm(Body_weight ~ Time_point + Strain, data = metadata_df)
vif(modelo_vif)
############################################################
#🦠 SUPPLEMENTARY — MAASLIN DOTPLOTS (A/B/C)
############################################################

make_dotplot <- function(df, title_text) {
  df <- df %>%
    mutate(coef     = as.numeric(coef),
           neglogq  = -log10(qval),
           Direction = ifelse(coef > 0, "Positive", "Negative"))
  df$DisplayName <- make.unique(as.character(df$DisplayName))
  df <- df[order(df$coef), ]
  df$DisplayName <- factor(df$DisplayName, levels = df$DisplayName)
  
  ggplot(df, aes(x = coef, y = DisplayName, size = neglogq, color = Direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(alpha = 0.9) +
    scale_color_manual(values = coef_colors) +
    theme_classic(base_size = 14) +
    labs(title = title_text, x = "Standardized coefficient",
         y = NULL, size = "-log10(q)")
}

p_dot_A <- make_dotplot(sig_main_t,    "A  |  Longitudinal model")
p_dot_B <- make_dotplot(sig_zonulin_t, "B  |  Zonulin associations")
p_dot_C <- make_dotplot(sig_lps_t,     "C  |  LPS associations")

Figure_Supp_Maaslin <- p_dot_A / p_dot_B / p_dot_C +
  plot_annotation(
    title = "Differential microbial genera identified by MaAsLin2",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

ggsave(paste0(Figures, "Supplementary_Maaslin_ABC.png"),
       Figure_Supp_Maaslin, width = 10, height = 18, dpi = 300)
ggsave(paste0(Figures, "Supplementary_Maaslin_ABC.pdf"),
       Figure_Supp_Maaslin, width = 10, height = 18)

############################################################
# 🦠SUPPLEMENTARY — OA SURGERY (A + B)
############################################################

ps_24  <- subset_samples(ps_clr, Time_point == "24 weeks")
ord_24 <- ordinate(ps_24, method = "PCoA", distance = "euclidean")

p_SA <- ggplot(alpha_24, aes(Surgery.status, Shannon, fill = Surgery.status)) +
  geom_boxplot(alpha = 0.75, color = "black", linewidth = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, shape = 21, size = 3.2, color = "black", stroke = 0.6) +
  scale_fill_manual(values = surgery_colors) +
  theme_classic(base_size = 16) +
  labs(title = "A  |  Alpha diversity by surgery (24 weeks)",
       x = "Surgery status", y = "Shannon diversity")

p_SB <- plot_ordination(ps_24, ord_24, color = "Surgery.status") +
  geom_point(shape = 21, aes(fill = Surgery.status), size = 4,
             stroke = 0.8, color = "black", alpha = 0.95) +
  stat_ellipse(aes(group = Surgery.status, fill = Surgery.status),
               geom = "polygon", alpha = 0.15, color = NA) +
  stat_ellipse(aes(group = Surgery.status, color = Surgery.status), linewidth = 1.2) +
  scale_color_manual(values = surgery_colors) +
  scale_fill_manual(values = surgery_colors) +
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 1) +
  labs(title = "B  |  Beta diversity by surgery (24 weeks)", x = "PCoA1", y = "PCoA2")

Figure_S2_OA <- p_SA | p_SB

ggsave(paste0(Figures, "FigureS2_OA_Alpha_Beta.png"),
       Figure_S2_OA, width = 14, height = 6, dpi = 300)

############################################################
# 🦠BIOMARKERS — ZONULIN + LPS OVER TIME
############################################################

bio_df <- metadata_df %>%
  select(Sample, Rat.ID, Time_point, Strain, Zonulin, LPS) %>%
  pivot_longer(cols = c(Zonulin, LPS), names_to = "Biomarker", values_to = "Value")

stat_tests_bio <- bio_df %>%
  group_by(Biomarker, Time_point) %>%
  wilcox_test(Value ~ Strain) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  left_join(bio_df %>% group_by(Biomarker, Time_point) %>%
              summarise(y.position = max(Value, na.rm = TRUE) * 1.05, .groups = "drop"),
            by = c("Biomarker", "Time_point"))

p_biomarkers <- ggplot(bio_df, aes(Time_point, Value, color = Strain, group = Strain)) +
  stat_summary(fun = mean, geom = "line",  linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  geom_jitter(width = 0.08, alpha = 0.5, size = 2) +
  stat_pvalue_manual(stat_tests_bio, label = "p.adj.signif", tip.length = 0.01) +
  facet_wrap(~Biomarker, scales = "free_y") +
  scale_color_manual(values = strain_colors) +
  theme_classic(base_size = 16) +
  labs(title = "Temporal dynamics of gut permeability markers",
       x = "Time point", y = "Concentration", color = "Strain")

ggsave(paste0(Figures, "Figure_Zonulin_LPS_Time_Strain.png"),
       p_biomarkers, width = 10, height = 6, dpi = 300)

############################################################
#🦠PICRUST PREPARATION
############################################################

maaslin_sig       <- sig_main
taxa_sig_features <- unique(maaslin_sig$feature)

write.csv(maaslin_sig,
          paste0(Figures, "Maaslin_Main_Significant_Taxa_for_PICRUSt.csv"),
          row.names = FALSE)

write.csv(sig_main_t %>%
            mutate(Model = "Time + Strain") %>%
            select(feature, DisplayName, metadata, coef, pval, qval, Model),
          paste0(Tablesen, "Table_S1_Maaslin_Main.csv"), row.names = FALSE)

cat("✔ Pipeline completo — listo para PiCRUST\n")
########
############################################################
# 🦠BOXPLOTS Y DOTPLOTS — MAASLIN2
############################################################

# función nombres
map_taxa_names <- function(features) {
  genus <- tax_table_df$DisplayName[match(features, tax_table_df$feature)]
  genus[is.na(genus)] <- features[is.na(genus)]
  make.unique(genus)
}

# función para obtener top 5 por variable
get_top <- function(df, var) {
  df %>%
    filter(metadata == var, qval < 0.05) %>%
    arrange(qval) %>%
    slice_head(n = 5) %>%
    pull(feature)
}

# función para preparar datos
prep_df <- function(features) {
  tax_names <- map_taxa_names(features)
  df <- taxa_table[, features, drop = FALSE] %>% as.data.frame()
  colnames(df) <- tax_names
  df$SampleID <- rownames(df)
  metadata_df$SampleID <- rownames(metadata_df)
  df <- merge(df, metadata_df, by = "SampleID")
  pivot_longer(df, cols = all_of(tax_names),
               names_to = "Taxa", values_to = "Abundance")
}

# top 5 por variable
top5_strain <- get_top(res_with_Body_weight, "Strain")
top5_time   <- get_top(res_with_Body_weight, "Time_point")
top5_Body_weight <- get_top(res_with_Body_weight, "Body_weight")

top5_strain <- intersect(top5_strain, colnames(taxa_table))
top5_time   <- intersect(top5_time,   colnames(taxa_table))
top5_Body_weight <- intersect(top5_Body_weight, colnames(taxa_table))

df_strain <- prep_df(top5_strain)
df_time   <- prep_df(top5_time)
df_Body_weight <- prep_df(top5_Body_weight)

# ── PLOTS ─────────────────────────────────────────────────

col_time <- c("#8DD3C7", "#80B1D3", "#BEBADA", "#FDB462", "#FB8072")

p_strain <- ggplot(df_strain, aes(Strain, Abundance, fill = Strain)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  facet_wrap(~Taxa, scales = "free_y") +
  scale_fill_manual(values = strain_colors) +
  theme_classic(base_size = 14) +
  labs(title = "A  |  Top 5 taxa by Strain", y = "CLR Abundance", x = NULL)

p_time <- ggplot(df_time, aes(Time_point, Abundance, fill = Time_point)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  facet_wrap(~Taxa, scales = "free_y") +
  scale_fill_manual(values = col_time) +
  theme_classic(base_size = 14) +
  labs(title = "B  |  Top 5 taxa by Time", y = "CLR Abundance", x = NULL)

p_Body_weight <- ggplot(df_Body_weight, aes(Body_weight, Abundance, color = Strain)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~Taxa, scales = "free_y") +
  scale_color_manual(values = strain_colors) +
  theme_classic(base_size = 14) +
  labs(title = "C  |  Top 5 taxa by Body_weight", y = "CLR Abundance")

# ── DOTPLOT LONGITUDINAL ──────────────────────────────────

top_long_features <- intersect(
  res_with_Body_weight %>%
    filter(qval < 0.05) %>%
    arrange(qval) %>%
    slice_head(n = 10) %>%
    pull(feature),
  colnames(taxa_table)
)

long_abundance <- taxa_table[, top_long_features, drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  tidyr::pivot_longer(-Sample,
                      names_to  = "feature",
                      values_to = "Abundance") %>%
  dplyr::left_join(metadata_df %>% dplyr::select(Sample, Time_point, Strain),
                   by = "Sample") %>%
  dplyr::left_join(tax_table_df %>% dplyr::select(feature, DisplayName),
                   by = "feature")

p_dot_long <- ggplot(long_abundance,
                     aes(Time_point, Abundance, color = Strain)) +
  geom_point(alpha = 0.6, size = 2,
             position = position_jitter(width = 0.15)) +
  facet_wrap(~DisplayName, scales = "free_y") +
  scale_color_manual(values = strain_colors) +
  theme_classic(base_size = 14) +
  labs(title = "Top 10 longitudinal taxa",
       y = "CLR abundance", x = "Time point")

# ── BOXPLOT TOP 5 GENERAL ─────────────────────────────────

top5_general <- res_with_Body_weight %>%
  filter(qval < 0.05, metadata == "Strain") %>%
  arrange(qval) %>%
  slice_head(n = 5) %>%
  pull(feature)

top5_general  <- intersect(top5_general, colnames(taxa_table))
tax_names_box <- map_taxa_names(top5_general)

df_box <- taxa_table[, top5_general, drop = FALSE] %>% as.data.frame()
colnames(df_box) <- tax_names_box
df_box$Sample <- rownames(df_box)  # ← usar Sample en vez de SampleID

df_box <- merge(df_box, 
                metadata_df %>% dplyr::select(Sample, Strain, Time_point),
                by = "Sample")  # ← join por Sample

df_box_long <- pivot_longer(df_box, cols = all_of(tax_names_box),
                            names_to = "Taxa", values_to = "Abundance")

p_box <- ggplot(df_box_long, aes(Strain, Abundance, fill = Strain)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  facet_wrap(~Taxa, scales = "free_y") +
  scale_fill_manual(values = strain_colors) +
  theme_classic(base_size = 14) +
  labs(title = "Top 5 taxa by Strain", y = "CLR Abundance", x = NULL)

# ── Save ───────────────────────────────────────────────

ggsave(paste0(Figures, "Fig_main.png"),
       p_strain / p_time / p_Body_weight,
       width = 12, height = 16, dpi = 300)

ggsave(paste0(Figures, "Fig_main.pdf"),
       p_strain / p_time / p_Body_weight,
       width = 12, height = 16)

ggsave(paste0(Figures, "Fig_boxplots.png"),
       p_box, width = 10, height = 6, dpi = 300)

ggsave(paste0(Figures, "Fig_boxplots.pdf"),
       p_box, width = 10, height = 6)

ggsave(paste0(Figures, "Dotplot_longitudinal.png"),
       p_dot_long, width = 12, height = 8, dpi = 300)

cat("✔ Boxplots y dotplots guardados\n")
####
# top 7 por cada variable
top7_strain <- res_with_Body_weight %>%
  filter(metadata == "Strain", qval < 0.05) %>%
  arrange(qval) %>%
  slice_head(n = 7) %>%
  pull(feature) %>%
  intersect(colnames(taxa_table))

top7_time <- res_with_Body_weight %>%
  filter(metadata == "Time_point", qval < 0.05) %>%
  arrange(qval) %>%
  slice_head(n = 7) %>%
  pull(feature) %>%
  intersect(colnames(taxa_table))

top7_Body_weight <- res_with_Body_weight %>%
  filter(metadata == "Body_weight", qval < 0.05) %>%
  arrange(qval) %>%
  slice_head(n = 7) %>%
  pull(feature) %>%
  intersect(colnames(taxa_table))

df_strain7 <- prep_df(top7_strain)
df_time7   <- prep_df(top7_time)
df_Body_weight7 <- prep_df(top7_Body_weight)

p_strain7 <- ggplot(df_strain7,
                    aes(Time_point, Abundance,
                        fill  = Strain,
                        color = Strain)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  facet_wrap(~Taxa, scales = "free_y", ncol = 4) +
  scale_fill_manual(values  = strain_colors) +
  scale_color_manual(values = strain_colors) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text  = element_text(face = "bold", size = 11)) +
  labs(title = "A  |  Top 7 Strain-associated taxa across time",
       y = "CLR Abundance", x = "Time point")

p_time7 <- ggplot(df_time7,
                  aes(Time_point, Abundance,
                      fill  = Strain,
                      color = Strain)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  facet_wrap(~Taxa, scales = "free_y", ncol = 4) +
  scale_fill_manual(values  = strain_colors) +
  scale_color_manual(values = strain_colors) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text  = element_text(face = "bold", size = 11)) +
  labs(title = "B  |  Top 7 Time-associated taxa across time",
       y = "CLR Abundance", x = "Time point")

p_Body_weight7 <- ggplot(df_Body_weight7,
                    aes(Time_point, Abundance,
                        fill  = Strain,
                        color = Strain)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  facet_wrap(~Taxa, scales = "free_y", ncol = 4) +
  scale_fill_manual(values  = strain_colors) +
  scale_color_manual(values = strain_colors) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text  = element_text(face = "bold", size = 11)) +
  labs(title = "C  |  Top 7 Body_weight-associated taxa across time",
       y = "CLR Abundance", x = "Time point")

# Save
ggsave(paste0(Figures, "Fig_Top7_Strain_time.png"),
       p_strain7, width = 16, height = 10, dpi = 300)

ggsave(paste0(Figures, "Fig_Top7_Time_time.png"),
       p_time7, width = 16, height = 10, dpi = 300)

ggsave(paste0(Figures, "Fig_Top7_Body_weight_time.png"),
       p_Body_weight7, width = 16, height = 10, dpi = 300)

# save combined
ggsave(paste0(Figures, "Fig_Top7_Combined.png"),
       p_strain7 / p_time7 / p_Body_weight7,
       width = 16, height = 28, dpi = 300)

ggsave(paste0(Figures, "Fig_Top7_Combined.pdf"),
       p_strain7 / p_time7 / p_Body_weight7,
       width = 16, height = 28)

cat("✔ Top 7 taxa por variable guardados\n")
cat(paste0("  -> Strain taxa:  ", length(top7_strain), "\n"))
cat(paste0("  -> Time taxa:    ", length(top7_time),   "\n"))
cat(paste0("  -> Body_weight taxa:  ", length(top7_Body_weight), "\n"))
######
# top 15 bacteria Strain
top15_strain <- get_top(res_with_Body_weight, "Strain")

#top 15
top15_strain <- res_with_Body_weight %>%
  filter(metadata == "Strain", qval < 0.05) %>%
  arrange(qval) %>%
  slice_head(n = 15) %>%
  pull(feature)

top15_strain <- intersect(top15_strain, colnames(taxa_table))

df_strain15 <- prep_df(top15_strain)

p_strain_time <- ggplot(df_strain15,
                        aes(Time_point, Abundance,
                            fill  = Strain,
                            color = Strain)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  facet_wrap(~Taxa, scales = "free_y", ncol = 5) +  # ← 5 columnas para 15 taxa
  scale_fill_manual(values  = strain_colors) +
  scale_color_manual(values = strain_colors) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text  = element_text(face = "bold", size = 12)
  ) +
  labs(title = "Strain-associated taxa across time",
       y     = "CLR Abundance",
       x     = "Time point")

ggsave(paste0(Figures, "Fig_Strain15_taxa_over_time.png"),
       p_strain_time, width = 18, height = 12, dpi = 300)

ggsave(paste0(Figures, "Fig_Strain15_taxa_over_time.pdf"),
       p_strain_time, width = 18, height = 12)

cat("✔ Top 15 strain taxa over time guardado\n")


res_with_Body_weight %>%
  filter(metadata == "Time_point", qval < 0.05) %>%
  nrow()
####


# count total and by variable
cat("=== RESUMEN MAASLIN2 (res_with_Body_weight) ===\n\n")

cat("Total significativos (qval < 0.05):", 
    nrow(res_with_Body_weight %>% filter(qval < 0.05)), "\n\n")

cat("Por variable:\n")
sig_W <- res_with_Body_weight[res_with_Body_weight$qval < 0.05 & !is.na(res_with_Body_weight$qval), ]

cat("=== RESUMEN MAASLIN2 (res_with_Body_weight) ===\n\n")

cat("Total significativos (qval < 0.05):", nrow(sig_W), "\n\n")

cat("Por variable:\n")
print(table(sig_W$metadata))

cat("\nTaxa unicos significativos:", length(unique(sig_W$feature)), "\n")

cat("\nTaxa unicos por variable:\n")
print(tapply(sig_W$feature, sig_W$metadata, function(x) length(unique(x))))
##########################################################


supp_taxa_table <- sig_W %>%
  as.data.frame() %>%
  dplyr::select(feature, metadata, coef, qval) %>%
  dplyr::left_join(
    tax_table_df %>% dplyr::select(feature, Family, Genus, DisplayName),
    by = "feature"
  ) %>%
  dplyr::mutate(
    Genus       = ifelse(is.na(Genus)  | Genus  == "", "Unknown", Genus),
    Family      = ifelse(is.na(Family) | Family == "", "Unknown", Family),
    DisplayName = ifelse(is.na(DisplayName) | DisplayName == "", feature, DisplayName),
    Direction   = ifelse(as.numeric(coef) > 0, "Enriched in Wi", "Enriched in SD"),
    ASV_label   = paste0(Family, " | ", Genus, " (", feature, ")")
  ) %>%
  dplyr::select(
    ASV         = feature,
    ASV_label,
    Family,
    Genus,
    DisplayName,
    Variable    = metadata,
    Direction,
    Coefficient = coef,
    qval
  ) %>%
  dplyr::arrange(Variable, qval)

# imprimir en consola
cat("=== TAXA SIGNIFICATIVOS POR VARIABLE ===\n\n")

for(var in unique(supp_taxa_table$Variable)) {
  cat("--- ", var, "---\n")
  df_var <- supp_taxa_table[supp_taxa_table$Variable == var, ]
  cat("N total:", nrow(df_var), "\n")
  print(df_var[, c("ASV", "Family", "Genus", "Direction", "qval")])
  cat("\n")
}

# save
write.csv(supp_taxa_table,
          paste0(Tablesen, "Table_Significant_Taxa_detail.csv"),
          row.names = FALSE)

cat("✔ Tabla guardada\n")

#######################################################################3
############################################################
# 🦠🦠🦠🦠🦠🦠🦠CORRELATIONS 24 weeks
# (with prevalence filter and exclusion of perfect rho )
############################################################



# name HOMA-IR
if ("HOMA-IR" %in% colnames(datos_variables)) {
  colnames(datos_variables)[colnames(datos_variables) == "HOMA-IR"] <- "HOMA.IR"
}

# all variables
all_clinical_vars <- c(
  "LPS", "G.CSF", "Eotaxin", "GM.CSF", "IL1a", "MIP1a", "IL4",
  "IL1b", "IL2", "IL6", "IL13", "IL10", "IL12", "IFNy", "IL5",
  "IL17a", "IL18", "MCP1", "IP10", "GRO.KC", "VEGF", "Fractalkine",
  "MIP2", "TNFa", "LIX", "Leptin", "RANTES",
  "Total.OARSI.score", "OARSI.Cartilage.degeneration",
  "OARSI.Synovial.inflammation", "OARSI.Subchondral.bone.score",
  "OARSI.osteophyte.size", "Total.osteophyte.diameter", "CISI.score",
  "Body_weight", "Triglicerides", "HOMA.IR", "Zonulin"
)

# To numéric
for (v in all_clinical_vars) {
  if (v %in% colnames(datos_variables)) {
    datos_variables[[v]] <- suppressWarnings(
      as.numeric(gsub(",", ".", as.character(datos_variables[[v]])))
    )
  }
}

#  time point labels
datos_variables$Time_point <- dplyr::recode(
  as.character(datos_variables$Time_point),
  "Baseline" = "0 weeks",
  "Midpoint" = "12 weeks",
  "Endpoint" = "24 weeks"
)

# ──metadata_df ──────────────────────────────────────────────────────
clinical_vars <- datos_variables %>%
  select(Sample, any_of(all_clinical_vars))

# avoid duplicates
cols_ya_presentes <- intersect(all_clinical_vars, colnames(metadata_df))
if (length(cols_ya_presentes) > 0) {
  clinical_vars <- clinical_vars %>%
    select(-any_of(cols_ya_presentes))
}

metadata_df <- metadata_df %>%
  left_join(clinical_vars, by = "Sample")

# CRÍTIC:  Sample as rownames
rownames(metadata_df) <- metadata_df$Sample

# Variables d
blood_vars_present <- intersect(all_clinical_vars, colnames(metadata_df))

cat("Variables clínicas disponibles (", length(blood_vars_present), "):\n")
print(blood_vars_present)

# ── Subset 24 weeks ──────────────────────────────────────────────────────
meta_24_cor <- metadata_df %>%
  filter(Time_point == "24 weeks")

# taxa_table
common_24   <- intersect(rownames(meta_24_cor), rownames(taxa_table))
meta_24_cor <- meta_24_cor[common_24, ]
ab_24_cor   <- taxa_table[common_24, , drop = FALSE]

cat("\nMuestras a 24 semanas alineadas:", nrow(meta_24_cor), "\n")
cat("Taxa antes de filtro prevalencia:", ncol(ab_24_cor), "\n")

# ── Filter prevalence (≥ 30% samples with > 0) ──────────────────
prevalencia_minima <- 0.30

taxa_prevalentes <- colnames(ab_24_cor)[
  apply(ab_24_cor, 2, function(x) mean(x > 0, na.rm = TRUE) >= prevalencia_minima)
]

cat("Taxa con prevalencia ≥ 30%:", length(taxa_prevalentes),
    "de", ncol(ab_24_cor), "\n")

ab_24_cor_filt <- ab_24_cor[, taxa_prevalentes, drop = FALSE]

# ── enough data ─────────────────────────
blood_vars_24 <- blood_vars_present[sapply(blood_vars_present, function(v) {
  vec <- as.numeric(meta_24_cor[[v]])
  sum(!is.na(vec)) >= 5
})]

cat("Variables clínicas con ≥5 obs:", length(blood_vars_24), "\n")

# ── spearmann correlation by strain ──────────────────────────────────────────
cor_results <- list()

for (s in unique(meta_24_cor$Strain)) {
  
  meta_s <- meta_24_cor[meta_24_cor$Strain == s, ]
  ab_s   <- ab_24_cor_filt[rownames(meta_s), , drop = FALSE]
  
  cat("Cepa:", s, "| n =", nrow(meta_s),
      "| taxa prevalentes:", ncol(ab_s), "\n")
  
  for (blood in blood_vars_24) {
    
    vec_blood <- as.numeric(meta_s[[blood]])
    
    if (sum(!is.na(vec_blood)) < 5)        next
    sd_blood <- sd(vec_blood, na.rm = TRUE)
    if (is.na(sd_blood) || sd_blood == 0)  next
    
    for (g in colnames(ab_s)) {
      
      vec_taxa <- as.numeric(ab_s[[g]])
      
      if (sum(!is.na(vec_taxa)) < 5)       next
      sd_taxa <- sd(vec_taxa, na.rm = TRUE)
      if (is.na(sd_taxa) || sd_taxa == 0)  next
      
      # Pares completos suficientes
      complete_pairs <- sum(!is.na(vec_blood) & !is.na(vec_taxa))
      if (complete_pairs < 5)              next
      
      test <- tryCatch(
        suppressWarnings(
          cor.test(vec_taxa, vec_blood,
                   method = "spearman",
                   use    = "complete.obs",
                   exact  = FALSE)
        ),
        error = function(e) NULL
      )
      
      if (is.null(test))                               next
      if (is.na(test$estimate) || is.na(test$p.value)) next
      
      cor_results[[length(cor_results) + 1]] <- data.frame(
        feature  = g,
        Strain   = s,
        Variable = blood,
        rho      = as.numeric(test$estimate),
        pval     = test$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
}

cat("\nTotal correlaciones calculadas:", length(cor_results), "\n")

if (length(cor_results) == 0) {
  stop("No se calcularon correlaciones. Revisa los diagnósticos arriba.")
}

cor_df <- dplyr::bind_rows(cor_results)

# ── Excluid perfect rho (artefacts) ──────────────────────────────
n_antes <- nrow(cor_df)
cor_df   <- cor_df %>% filter(abs(rho) < 0.99)
n_despues <- nrow(cor_df)

cat("Correlaciones removidas (|rho| ≥ 0.99):", n_antes - n_despues, "\n")
cat("Correlaciones restantes:               ", n_despues, "\n")

# ── FDR by strain ──────────────────────────────────────────────────────────────
cor_df <- cor_df %>%
  group_by(Strain) %>%
  mutate(p.FDR = p.adjust(pval, method = "BH")) %>%
  ungroup()

# ── Taxonomy ──────────────────────────────────────────────────────────
cor_df <- cor_df %>%
  left_join(
    tax_table_df %>%
      mutate(Taxa_label = case_when(
        !is.na(Family) & Family != "" & !is.na(Genus) & Genus != "" ~
          paste0(Family, " | ", Genus),
        !is.na(Family) & Family != "" ~ Family,
        TRUE ~ feature
      )) %>%
      select(feature, Taxa_label),
    by = "feature"
  ) %>%
  mutate(Taxa_label = ifelse(is.na(Taxa_label), feature, Taxa_label))

# ── Biological categories ─────────────────────────────────────────────────────
cor_df <- cor_df %>%
  mutate(Category = case_when(
    Variable %in% c("Body_weight", "Triglicerides", "HOMA.IR", "Leptin") ~
      "Metabolic",
    Variable %in% c("Zonulin", "LPS") ~
      "Gut barrier",
    Variable %in% c("IL1a", "IL1b", "IL2", "IL4", "IL5", "IL6", "IL10",
                    "IL12", "IL13", "IL17a", "IL18", "TNFa", "IFNy",
                    "MCP1", "MIP1a", "MIP2", "IP10", "GRO.KC", "RANTES",
                    "Fractalkine", "Eotaxin", "G.CSF", "GM.CSF",
                    "VEGF", "LIX") ~
      "Inflammation",
    Variable %in% c("Total.OARSI.score", "OARSI.Cartilage.degeneration",
                    "OARSI.Synovial.inflammation", "OARSI.Subchondral.bone.score",
                    "OARSI.osteophyte.size", "Total.osteophyte.diameter",
                    "CISI.score") ~
      "OA pathology",
    TRUE ~ "Other"
  ))

cat("\nResumen correlaciones:\n")
cat("  Total:             ", nrow(cor_df), "\n")
cat("  pval < 0.05:       ", nrow(cor_df %>% filter(pval  < 0.05)), "\n")
cat("  FDR  < 0.05:       ", nrow(cor_df %>% filter(p.FDR < 0.05)), "\n")

# ── Tablas ────────────────────────────────────────────────────────────────────
write.csv(
  cor_df %>%
    arrange(pval) %>%
    mutate(rho   = round(rho, 3),
           pval  = signif(pval,  3),
           p.FDR = signif(p.FDR, 3)),
  paste0(Tablesen, "TableS_Correlations_24weeks_Exploratory.csv"),
  row.names = FALSE
)

write.csv(
  cor_df %>%
    filter(p.FDR < 0.05) %>%
    arrange(Strain, p.FDR) %>%
    select(Strain, Taxa_label, Variable, Category, rho, pval, p.FDR) %>%
    mutate(rho   = round(rho, 3),
           pval  = signif(pval,  3),
           p.FDR = signif(p.FDR, 3)),
  paste0(Tablesen, "Table_Correlations_FDR_significant.csv"),
  row.names = FALSE
)

############################################################
#🦠FIGUReS
############################################################

# ── significant data ─────────────────────────────────────────────
cor_fdr_sig <- cor_df %>%
  filter(p.FDR < 0.05) %>%
  arrange(Strain, p.FDR) %>%
  mutate(sig_label = case_when(
    p.FDR < 0.001 ~ "***",
    p.FDR < 0.01  ~ "**",
    p.FDR < 0.05  ~ "*",
    TRUE          ~ ""
  ))

cat("\nCorrelations FDR < 0.05:", nrow(cor_fdr_sig), "\n")
cat("Por cepa:\n")
print(table(cor_fdr_sig$Strain))

# ── Figura principal: top 20 pval < 0.05 ─────────────────────────────────────
top_features_cor <- cor_df %>%
  filter(pval < 0.05) %>%
  group_by(feature) %>%
  summarise(pval_min = min(pval), .groups = "drop") %>%
  arrange(pval_min) %>%
  slice_head(n = 20) %>%
  pull(feature)

plot_cor_main <- cor_df %>%
  filter(feature %in% top_features_cor) %>%
  mutate(
    Taxa_label = factor(Taxa_label,
                        levels = rev(unique(Taxa_label[
                          order(abs(rho), decreasing = TRUE)]))),
    sig_label  = case_when(
      p.FDR < 0.05 ~ "**",
      pval  < 0.05 ~ "*",
      TRUE         ~ ""
    )
  )

p_cor_main <- ggplot(plot_cor_main,
                     aes(x = Variable, y = Taxa_label, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label),
            size = 4, color = "black", vjust = 0.8) +
  facet_grid(. ~ Strain) +
  scale_fill_gradient2(
    low    = gradient_low,
    mid    = gradient_mid,
    high   = gradient_high,
    name   = "Spearman ρ",
    limits = c(-1, 1)
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y      = element_text(size = 10),
    strip.text       = element_text(size = 13, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title       = element_text(size = 14, face = "bold", hjust = 0),
    plot.subtitle    = element_text(size = 10, color = "grey40"),
    legend.position  = "right"
  ) +
  labs(
    title    = "Exploratory microbiota–clinical correlations (24 weeks)",
    subtitle = "Top 20 taxa | * p<0.05  ** FDR<0.05 | Spearman | Prevalence ≥ 30%",
    x        = "Clinical variable",
    y        = NULL
  )

ggsave(paste0(Figures, "Figure_Correlations_24w_Main.png"),
       p_cor_main,
       width  = max(12, length(unique(plot_cor_main$Variable)) * 0.6),
       height = max(6,  length(unique(plot_cor_main$Taxa_label)) * 0.4),
       dpi    = 300)

# ── Figura FDR significativos ─────────────────────────────────────────────────
if (nrow(cor_fdr_sig) == 0) {
  
  cat("⚠️ No correlations with FDR < 0.05 after prevalence filter\n")
  
} else {
  
  taxa_order_fdr <- cor_fdr_sig %>%
    group_by(Taxa_label) %>%
    summarise(mean_abs = mean(abs(rho)), .groups = "drop") %>%
    arrange(desc(mean_abs)) %>%
    pull(Taxa_label)
  
  cor_fdr_sig$Taxa_label <- factor(cor_fdr_sig$Taxa_label,
                                   levels = rev(taxa_order_fdr))
  
  p_fdr_main <- ggplot(cor_fdr_sig,
                       aes(x = Variable, y = Taxa_label, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sig_label),
              size = 4.5, color = "black", vjust = 0.8) +
    facet_grid(. ~ Strain) +
    scale_fill_gradient2(
      low    = gradient_low,
      mid    = gradient_mid,
      high   = gradient_high,
      name   = "Spearman ρ",
      limits = c(-1, 1)
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y      = element_text(size = 11),
      strip.text       = element_text(size = 13, face = "bold"),
      strip.background = element_rect(fill = "grey95", color = NA),
      plot.title       = element_text(size = 14, face = "bold", hjust = 0),
      plot.subtitle    = element_text(size = 10, color = "grey40"),
      legend.position  = "right"
    ) +
    labs(
      title    = "FDR-significant microbiota–clinical correlations (24 weeks)",
      subtitle = "FDR < 0.05 | Spearman | Prevalence ≥ 30% | |ρ| < 0.99",
      x        = "Clinical variable",
      y        = NULL
    )
  
  ggsave(paste0(Figures, "Figure_Correlations_FDR_significant.png"),
         p_fdr_main,
         width  = max(10, length(unique(cor_fdr_sig$Variable)) * 0.8),
         height = max(5,  length(unique(cor_fdr_sig$Taxa_label)) * 0.5),
         dpi    = 300)
  
  # ── Suplementary by biological category ──────────────────────────
  categorias <- c("Metabolic", "Gut barrier", "Inflammation", "OA pathology")
  e  plots_por_categoria <- lapply(categorias, function(cat) {
    
    df_cat <- cor_df %>%
      filter(pval < 0.05, Category == cat) %>%
      mutate(sig_label = case_when(
        p.FDR < 0.05 ~ "**",
        pval  < 0.05 ~ "*",
        TRUE         ~ ""
      ))
    
    if (nrow(df_cat) == 0) return(NULL)
    
    taxa_order_cat <- df_cat %>%
      group_by(Taxa_label) %>%
      summarise(mean_abs = mean(abs(rho)), .groups = "drop") %>%
      arrange(desc(mean_abs)) %>%
      pull(Taxa_label)
    
    df_cat$Taxa_label <- factor(df_cat$Taxa_label,
                                levels = rev(taxa_order_cat))
    
    ggplot(df_cat, aes(x = Variable, y = Taxa_label, fill = rho)) +
      geom_tile(color = "white", linewidth = 0.4) +
      geom_text(aes(label = sig_label),
                size = 4, color = "black", vjust = 0.8) +
      facet_grid(. ~ Strain) +
      scale_fill_gradient2(
        low    = gradient_low,
        mid    = gradient_mid,
        high   = gradient_high,
        name   = "Spearman ρ",
        limits = c(-1, 1)
      ) +
      theme_classic(base_size = 12) +
      theme(
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y      = element_text(size = 10),
        strip.text       = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = NA),
        plot.title       = element_text(size = 13, face = "bold", hjust = 0),
        legend.position  = "right"
      ) +
      labs(
        title = paste0(cat, " variables"),
        x     = "Clinical variable",
        y     = NULL
      )
  })
  
  plots_por_categoria <- Filter(Negate(is.null), plots_por_categoria)
  
  if (length(plots_por_categoria) > 0) {
    
    Figure_Cor_byCategory <- wrap_plots(plots_por_categoria, ncol = 1) +
      plot_annotation(
        title    = "Microbiota–clinical correlations by biological category (24 weeks)",
        subtitle = "* p<0.05  ** FDR<0.05 | Spearman | Prevalence ≥ 30%",
        theme    = theme(
          plot.title    = element_text(size = 14, face = "bold",    hjust = 0.5),
          plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0.5)
        )
      )
    
    ggsave(paste0(Figures, "FigureS_Correlations_24w_byCategory.png"),
           Figure_Cor_byCategory,
           width  = 16,
           height = 6 * length(plots_por_categoria),
           dpi    = 300)
  }
}

cat("✔ Correlaciones exploratorias completadas\n")
cat("  Tabla completa:       TableS_Correlations_24weeks_Exploratory.csv\n")
cat("  Tabla FDR sig:        Table_Correlations_FDR_significant.csv\n")
cat("  Figura principal:     Figure_Correlations_24w_Main.png\n")
cat("  Figura FDR:           Figure_Correlations_FDR_significant.png\n")
cat("  Figura por categoría: FigureS_Correlations_24w_byCategory.png\n")

# ── 🦠biological category ──────────────────────────────
categorias <- c("Metabolic", "Gut barrier", "Inflammation", "OA pathology")

plots_por_categoria <- lapply(categorias, function(cat) {
  
  df_cat <- cor_df %>%
    filter(pval < 0.05, Category == cat) %>%
    mutate(sig_label = case_when(
      p.FDR < 0.05 ~ "**",
      pval  < 0.05 ~ "*",
      TRUE         ~ ""
    ))
  
  if (nrow(df_cat) == 0) {
    cat("Sin datos para categoría:", cat, "\n")
    return(NULL)
  }
  
  cat("Categoría:", cat, "| n asociaciones:", nrow(df_cat), "\n")
  
  taxa_order_cat <- df_cat %>%
    group_by(Taxa_label) %>%
    summarise(mean_abs = mean(abs(rho)), .groups = "drop") %>%
    arrange(desc(mean_abs)) %>%
    pull(Taxa_label)
  
  df_cat$Taxa_label <- factor(df_cat$Taxa_label,
                              levels = rev(taxa_order_cat))
  
  ggplot(df_cat, aes(x = Variable, y = Taxa_label, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sig_label),
              size = 4.5, color = "black", vjust = 0.8) +
    facet_grid(. ~ Strain) +
    scale_fill_gradient2(
      low    = gradient_low,
      mid    = gradient_mid,
      high   = gradient_high,
      name   = "Spearman ρ",
      limits = c(-1, 1)
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y      = element_text(size = 11),
      strip.text       = element_text(size = 13, face = "bold"),
      strip.background = element_rect(fill = "grey95", color = NA),
      plot.title       = element_text(size = 14, face = "bold", hjust = 0),
      plot.subtitle    = element_text(size = 10, color = "grey40"),
      legend.position  = "right",
      plot.margin      = margin(10, 10, 10, 10)
    ) +
    labs(
      title    = paste0(cat, " — microbiota associations (24 weeks)"),
      subtitle = "* p<0.05  ** FDR<0.05 | Spearman | Prevalence ≥ 30%",
      x        = "Clinical variable",
      y        = NULL
    )
})

# Quit NULLs
plots_por_categoria <- Filter(Negate(is.null), plots_por_categoria)

cat("\nCategorías con datos:", length(plots_por_categoria), "\n")

if (length(plots_por_categoria) == 0) {
  cat("⚠️ Sin correlaciones significativas en ninguna categoría\n")
} else {
  
  # ──  panel ────────────────────────────────────────
  n_taxa_por_panel <- sapply(categorias, function(cat) {
    df_cat <- cor_df %>% filter(pval < 0.05, Category == cat)
    length(unique(df_cat$Taxa_label))
  })
  n_taxa_por_panel <- n_taxa_por_panel[n_taxa_por_panel > 0]
  
  alto_total <- sum(pmax(n_taxa_por_panel * 0.5, 4)) + 2
  
  Figure_Cor_byCategory <- wrap_plots(plots_por_categoria, ncol = 1) +
    plot_annotation(
      title    = "Microbiota–clinical correlations by biological category (24 weeks)",
      subtitle = "* p<0.05  ** FDR<0.05 | Spearman rank correlation | Prevalence ≥ 30% | |ρ| < 0.99",
      theme    = theme(
        plot.title    = element_text(size = 15, face = "bold",    hjust = 0.5),
        plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0.5)
      )
    )
  
  ggsave(paste0(Figures, "FigureS_Correlations_24w_byCategory.png"),
         Figure_Cor_byCategory,
         width  = 16,
         height = alto_total,
         dpi    = 300)
  
  # separate panel save
  for (i in seq_along(plots_por_categoria)) {
    cat_name <- gsub(" ", "_", categorias[categorias %in%
                                            sapply(plots_por_categoria,
                                                   function(p) p$labels$title)][i])
    ggsave(
      paste0(Figures, "FigureS_Cor_", i, "_",
             gsub(" ", "_", categorias[i]), ".png"),
      plots_por_categoria[[i]],
      width  = max(8, length(unique(
        cor_df %>%
          filter(pval < 0.05, Category == categorias[i]) %>%
          pull(Variable)
      )) * 0.8),
      height = max(4, length(unique(
        cor_df %>%
          filter(pval < 0.05, Category == categorias[i]) %>%
          pull(Taxa_label)
      )) * 0.5),
      dpi = 300
    )
  }
  
  cat("✔ Figures by category saved\n")
  cat("  Combined: FigureS_Correlations_24w_byCategory.png\n")
  cat("  Separate: FigureS_Cor_1_Metabolic.png, etc.\n")
}
