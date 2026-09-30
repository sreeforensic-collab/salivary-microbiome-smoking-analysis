# ============================================================
# Salivary Microbiome Smoking-Status Analysis
# Kraken2 + Bracken -> R downstream analysis
# 
#
#
# Upstream:
# FASTQ -> FastQC -> Trimmomatic -> Bowtie2 -> Kraken2 -> Bracken
#
# This script starts from Bracken genus reports and generates:
# 1. Top-10 dominant genera
# 2. Shannon and Simpson diversity
# 3. Bray-Curtis PCoA + PERMANOVA
# 4. DESeq2 pairwise differential abundance
# 5. Differential-abundance heatmap
# 6. Selected-genus relative-abundance plots
# 7. Prevalence plot on Selected genera
# ============================================================

# ---------- 0. Packages ----------
required <- c("tidyverse", "vegan", "DESeq2", "pheatmap")

for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "DESeq2") {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("DESeq2", ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}

library(tidyverse)
library(vegan)
library(DESeq2)
library(pheatmap)

# ---------- 1. INPUT PATHS ----------
bracken_dir <- "analysis_27samples/bracken_reports"
metadata_file <- "sample_metadata.csv"

out_dir <- "R_poster_analysis"
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 2. METADATA ----------
# Required columns:
# sample_id
# smoking_status
#
# smoking_status:
# Current_smoker / Former_smoker / Never_smoker

metadata <- read.csv("metadata.csv", stringsAsFactors = FALSE,
                     check.names = FALSE)

metadata$sample_id <- as.character(metadata$sample_id)

metadata$smoking_status <- factor(
  metadata$smoking_status,
  levels = c("Current_smoker", "Former_smoker", "Never_smoker")
)

print(table(metadata$smoking_status))

# ---------- 3. READ BRACKEN GENUS REPORTS ----------
# Expected columns include:
# name, new_est_reads, fraction_total_reads

files <- list.files(
  "bracken_reports",
  pattern = "_genus\\.txt$",
  full.names = TRUE
)

if (length(files) == 0)
  stop("No Bracken genus files found. Check bracken_dir.")

# ============================================================
# 3. READ BRACKEN GENUS REPORTS
# ============================================================

read_bracken <- function(file) {
  
  x <- read.delim(
    file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Check required columns
  required_cols <- c("name", "new_est_reads")
  
  if (!all(required_cols %in% colnames(x))) {
    stop(
      paste(
        "Required Bracken columns missing in:",
        basename(file)
      )
    )
  }
  
  sample_id <- sub(
    "_genus\\.txt$",
    "",
    basename(file)
  )
  
  x %>%
    dplyr::select(
      name,
      new_est_reads
    ) %>%
    dplyr::rename(
      Genus = name,
      Count = new_est_reads
    ) %>%
    dplyr::mutate(
      Count = as.numeric(Count),
      sample_id = sample_id
    )
}

bracken_long <- purrr::map_dfr(files, read_bracken)

# ---------- 4. GENUS COUNT MATRIX ----------
genus_matrix_df <- bracken_long %>%
  group_by(Genus, sample_id) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sample_id, values_from = Count, values_fill = 0)

genus_matrix <- as.data.frame(genus_matrix_df)
rownames(genus_matrix) <- genus_matrix$Genus
genus_matrix$Genus <- NULL
genus_matrix <- as.matrix(genus_matrix)
mode(genus_matrix) <- "numeric"

# Keep only samples with metadata
common_samples <- intersect(colnames(genus_matrix), metadata$sample_id)
genus_matrix <- genus_matrix[, common_samples, drop = FALSE]

metadata <- metadata[
  match(common_samples, metadata$sample_id),
  ,
  drop = FALSE
]
rownames(metadata) <- metadata$sample_id

stopifnot(identical(colnames(genus_matrix), rownames(metadata)))

write.csv(genus_matrix,
          file.path(out_dir, "genus_count_matrix.csv"))

# ---------- 5. RELATIVE ABUNDANCE ----------
relative_abundance <- sweep(
  genus_matrix, 2, colSums(genus_matrix), "/"
) * 100

write.csv(relative_abundance,
          file.path(out_dir, "relative_abundance_genus.csv"))

# ============================================================
# 6. TOP 10 DOMINANT GENERA
# ============================================================

mean_abundance <- sort(
  rowMeans(relative_abundance, na.rm = TRUE),
  decreasing = TRUE
)

top10_genera <- names(mean_abundance)[1:10]

top10_long <- as.data.frame(
  relative_abundance[top10_genera, , drop = FALSE]
) %>%
  mutate(Genus = rownames(.)) %>%
  pivot_longer(
    cols = -Genus,
    names_to = "SampleID",
    values_to = "Relative_abundance"
  ) %>%
  left_join(
    metadata %>% select(sample_id, smoking_status),
    by = c("SampleID" = "sample_id")
  )

top10_group <- top10_long %>%
  group_by(smoking_status, Genus) %>%
  summarise(
    Mean_abundance = mean(Relative_abundance, na.rm = TRUE),
    .groups = "drop"
  )

top10_group$smoking_status <- factor(
  top10_group$smoking_status,
  levels = c("Current_smoker", "Former_smoker", "Never_smoker")
)

top10_plot <- ggplot(
  top10_group,
  aes(smoking_status, Mean_abundance, fill = Genus)
) +
  geom_col(width = 0.65) +
  labs(
    title = "Top 10 Dominant Bacterial Genera",
    subtitle = "Oral microbiome composition across smoking-status groups",
    x = NULL,
    y = "Mean relative abundance (%)",
    fill = "Genus"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

ggsave(
  file.path(fig_dir, "Top10_Genera_Smoking_Status.png"),
  top10_plot, width = 9, height = 6, dpi = 300
)

write.csv(
  data.frame(
    Genus = top10_genera,
    Overall_mean = mean_abundance[top10_genera]
  ),
  file.path(out_dir, "Top10_Genera_Summary.csv"),
  row.names = FALSE
)

# ============================================================
# 7. ALPHA DIVERSITY
# ============================================================

shannon <- diversity(t(genus_matrix), index = "shannon")
simpson <- diversity(t(genus_matrix), index = "simpson")

alpha_data <- metadata %>%
  select(sample_id, smoking_status) %>%
  mutate(
    Shannon = shannon[sample_id],
    Simpson = simpson[sample_id]
  )

write.csv(alpha_data,
          file.path(out_dir, "alpha_diversity.csv"),
          row.names = FALSE)

p_shannon <- ggplot(
  alpha_data,
  aes(smoking_status, Shannon, fill = smoking_status)
) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.10, size = 2) +
  labs(title = "Shannon Diversity", x = NULL, y = "Shannon index") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

p_simpson <- ggplot(
  alpha_data,
  aes(smoking_status, Simpson, fill = smoking_status)
) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.10, size = 2) +
  labs(title = "Simpson Diversity", x = NULL, y = "Simpson index") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "Shannon_Diversity.png"),
       p_shannon, width = 5, height = 4, dpi = 300)

ggsave(file.path(fig_dir, "Simpson_Diversity.png"),
       p_simpson, width = 5, height = 4, dpi = 300)

shannon_kw <- kruskal.test(Shannon ~ smoking_status, data = alpha_data)
simpson_kw <- kruskal.test(Simpson ~ smoking_status, data = alpha_data)

capture.output(
  shannon_kw, simpson_kw,
  file = file.path(out_dir, "alpha_diversity_statistics.txt")
)

# ============================================================
# 8. BRAY-CURTIS PCoA + PERMANOVA
# ============================================================

bray_dist <- vegdist(t(relative_abundance), method = "bray")

pcoa <- cmdscale(bray_dist, k = 2, eig = TRUE)

pcoa_df <- data.frame(
  SampleID = rownames(pcoa$points),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
) %>%
  left_join(
    metadata %>% select(sample_id, smoking_status),
    by = c("SampleID" = "sample_id")
  )

eig <- pcoa$eig
PC1_percent <- round(100 * eig[1] / sum(eig[eig > 0]), 1)
PC2_percent <- round(100 * eig[2] / sum(eig[eig > 0]), 1)

pcoa_plot <- ggplot(
  pcoa_df,
  aes(PC1, PC2, fill = smoking_status)
) +
  geom_point(shape = 21, size = 4, colour = "black") +
  labs(
    title = "Bray-Curtis PCoA",
    x = paste0("PCoA1 (", PC1_percent, "%)"),
    y = paste0("PCoA2 (", PC2_percent, "%)")
  ) +
  theme_classic(base_size = 14)

ggsave(file.path(fig_dir, "Bray_Curtis_PCoA.png"),
       pcoa_plot, width = 6, height = 5, dpi = 300)

permanova <- adonis2(
  bray_dist ~ smoking_status,
  data = metadata,
  permutations = 999
)

write.table(
  as.data.frame(permanova),
  file.path(out_dir, "PERMANOVA_Bray_Curtis.txt"),
  sep = "\t", quote = FALSE
)

# ============================================================
# 9. DESEQ2 DIFFERENTIAL ABUNDANCE
# ============================================================

dds <- DESeqDataSetFromMatrix(
  countData = round(genus_matrix),
  colData = metadata,
  design = ~ smoking_status
)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds$smoking_status <- relevel(
  dds$smoking_status,
  ref = "Never_smoker"
)

dds <- DESeq(dds)

res_current_vs_never <- results(
  dds,
  contrast = c("smoking_status", "Current_smoker", "Never_smoker")
)

res_former_vs_never <- results(
  dds,
  contrast = c("smoking_status", "Former_smoker", "Never_smoker")
)

res_current_vs_former <- results(
  dds,
  contrast = c("smoking_status", "Current_smoker", "Former_smoker")
)

write.csv(
  as.data.frame(res_current_vs_never),
  file.path(out_dir, "DESeq2_Current_vs_Never.csv")
)

write.csv(
  as.data.frame(res_former_vs_never),
  file.path(out_dir, "DESeq2_Former_vs_Never.csv")
)

write.csv(
  as.data.frame(res_current_vs_former),
  file.path(out_dir, "DESeq2_Current_vs_Former.csv")
)

# ============================================================
# 10. DIFFERENTIAL-ABUNDANCE HEATMAP
# ============================================================

get_sig <- function(res) {
  as.data.frame(res) %>%
    rownames_to_column("Genus") %>%
    filter(!is.na(padj), padj < 0.05) %>%
    arrange(padj) %>%
    pull(Genus)
}

sig_taxa <- unique(c(
  get_sig(res_current_vs_never),
  get_sig(res_former_vs_never),
  get_sig(res_current_vs_former)
))

sig_taxa <- sig_taxa[1:min(20, length(sig_taxa))]

if (length(sig_taxa) > 0) {

  hm <- relative_abundance[
    intersect(sig_taxa, rownames(relative_abundance)),
    , drop = FALSE
  ]

  hm_z <- t(scale(t(hm)))

  annotation_col <- data.frame(
    Smoking_status = metadata$smoking_status
  )
  rownames(annotation_col) <- metadata$sample_id

  pheatmap(
    hm_z,
    annotation_col = annotation_col,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize = 9,
    main = "Differentially Abundant Genera",
    filename = file.path(fig_dir, "Differential_Abundance_Heatmap.png"),
    width = 8, height = 6
  )
}

# ============================================================
# 11. SELECTED GENERA
# ============================================================

selected_genera <- c(
  "Blautia",
  "Filifactor",
  "Gemella",
  "Neisseria"
)

selected_genera <- intersect(
  selected_genera,
  rownames(relative_abundance)
)

selected_long <- as.data.frame(
  relative_abundance[selected_genera, , drop = FALSE]
) %>%
  mutate(Genus = rownames(.)) %>%
  pivot_longer(
    cols = -Genus,
    names_to = "SampleID",
    values_to = "Relative_abundance"
  ) %>%
  left_join(
    metadata %>% select(sample_id, smoking_status),
    by = c("SampleID" = "sample_id")
  )

selected_plot <- ggplot(
  selected_long,
  aes(smoking_status, Relative_abundance, fill = smoking_status)
) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.10, size = 1.8) +
  facet_wrap(~Genus, scales = "free_y") +
  labs(
    title = "Relative Abundance of Selected Genera",
    x = NULL,
    y = "Relative abundance (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "Selected_Genera_Relative_Abundance.png"),
  selected_plot, width = 8, height = 6, dpi = 300
)
# ============================================================
# 8. PREVALENCE ANALYSIS
# ============================================================
# Prevalence = percentage of samples in each smoking-status
# group in which a genus is detected.
#
# Detection criterion:
# Genus is considered present when Count > 0.
# ============================================================

prevalence_genus <- bracken_long %>%
  left_join(
    metadata %>%
      select(sample_id, smoking_status),
    by = "sample_id"
  ) %>%
  group_by(Genus, smoking_status) %>%
  summarise(
    prevalence_percent = mean(Count >= 20) * 100,
    .groups = "drop"
  )


# ============================================================
# Selected genera for prevalence plot
# ============================================================

selected_genera <- c(
  "Blautia",
  "Filifactor",
  "Gemella",
  "Neisseria"
)

prevalence_selected <- prevalence_genus %>%
  filter(Genus %in% selected_genera)


# ============================================================
# Prevalence plot
# ============================================================

prevalence_plot <- ggplot(
  prevalence_selected,
  aes(
    x = Genus,
    y = prevalence_percent,
    fill = smoking_status
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  labs(
    title = "Prevalence of Selected Oral Bacterial Genera",
    x = "Genus",
    y = "Prevalence (%)",
    fill = "Smoking status"
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold"
    )
  )


# Display plot
print(prevalence_plot)


# ============================================================
# Save prevalence plot
# ============================================================

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "results/figures/prevalence_selected_genera.png",
  prevalence_plot,
  width = 8,
  height = 5,
  dpi = 300
)

# ============================================================
# 12. SUMMARY
# ============================================================

summary_lines <- c(
  "Salivary microbiome — Kraken2/Bracken R analysis",
  "",
  "Smoking-status groups:",
  capture.output(table(metadata$smoking_status)),
  "",
  "Top 10 genera:",
  paste(seq_along(top10_genera), top10_genera, sep = ". "),
  "",
  paste("Shannon Kruskal-Wallis p =", signif(shannon_kw$p.value, 4)),
  paste("Simpson Kruskal-Wallis p =", signif(simpson_kw$p.value, 4)),
  "",
  "PERMANOVA:",
  capture.output(permanova)
)

writeLines(
  summary_lines,
  file.path(out_dir, "poster_analysis_summary.txt")
)

cat("\nAnalysis complete. Results are in:", out_dir, "\n")
