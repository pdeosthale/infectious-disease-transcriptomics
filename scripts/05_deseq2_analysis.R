# ==============================================================================
# DESeq2 Differential Expression Analysis
# Dataset: GSE161731 — Blood Transcriptomics in Respiratory Infections
# Cohorts: COVID-19, Influenza, Bacterial, CoV_other, Healthy
#
# Source: McClain, M.T., Constantine, F.J., Henao, R. et al. Dysregulated 
#         transcriptional responses to SARS-CoV-2 in the periphery. 
#         Nature Communications 12, 1079 (2021).
#         https://doi.org/10.1038/s41467-021-21289-y

# Date: June 2026
#
# This script performs an independent secondary analysis of publicly available
# data from GEO. Interpretations and analytical choices reflect the author's
# own analysis and do not represent the views of the original investigators.
#
# This script performs differential expression analysis comparing host 
# transcriptomic responses across respiratory infections (viral and bacterial)
# relative to healthy controls using blood RNA-seq data.
#
# DESIGN FORMULA: ~ batch + gender + age + cohort
#
#   Covariates (adjustment variables — control for, not interested in):
#     - batch: technical variation from processing on different days
#     - gender: biological variation in sex-linked gene expression
#     - age: biological variation (CRITICAL — severe age confounding in this dataset)
#
#   Variable of interest (what we are testing):
#     - cohort: infection type vs healthy controls
#     - The variable of interest always goes LAST in the DESeq2 formula
#     - DESeq2 removes the effects of covariates first, then tests whether
#       cohort still has a significant effect on gene expression
#
# Key analytical decisions documented in this script:
#   1. Removed 21 samples with missing gender data (all from CoV_other cohort)
#   2. Removed samples with missing age data
#   3. Identified and corrected for age confounding between cohorts
#      - Healthy controls: median age 18 (college students)
#      - Bacterial patients: median age 60
#      - Without age correction: 10,564 Bacterial DE genes
#      - With age correction: 8,436 Bacterial DE genes (~20% were age effects)
# ==============================================================================


# ==============================================================================
# 1. SETUP: Load Libraries
# ==============================================================================

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)
library(apeglm)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(UpSetR)
library(ggrepel)

# Set working directory to project root
setwd("~/Desktop/infectious-disease-transcriptomics")


# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

# Load filtered counts and metadata saved from Python QC notebook (04_real_geo_data.ipynb)
# Counts were pre-filtered in Python: removed genes detected in < 10% of samples
# Two outlier samples were removed in Python QC (low library size / low gene detection)
counts <- read.csv("results/GSE161731_counts_filtered.csv", row.names = 1)
meta <- read.csv("results/GSE161731_metadata_clean.csv")

cat("Raw loaded - Counts:", dim(counts)[1], "genes x", dim(counts)[2], "samples\n")
cat("Raw loaded - Metadata:", dim(meta)[1], "samples\n")


# ==============================================================================
# 3. MATCH SAMPLE IDS BETWEEN COUNTS AND METADATA
# ==============================================================================

# R automatically converts dashes to dots and adds "X" prefix to column names 
# that start with numbers when reading CSVs. For example:
#   "DU18-02S0011614" becomes "DU18.02S0011614"
#   "105848" becomes "X105848"
# We use make.names() to apply the same transformation to metadata IDs
# so they match the counts column names.

meta$rna_id_fixed <- make.names(meta$rna_id)
rownames(meta) <- meta$rna_id_fixed

# Verify all samples match
overlap <- sum(meta$rna_id_fixed %in% colnames(counts))
cat("Sample ID overlap:", overlap, "of", nrow(meta), "\n")

# Reorder metadata to match counts column order (required by DESeq2)
meta <- meta[colnames(counts), ]
stopifnot(all(colnames(counts) == rownames(meta)))
cat("Sample order verified: counts columns match metadata rows\n")


# ==============================================================================
# 4. ASSESS ADJUSTMENT VARIABLES AND CONFOUNDING
# ==============================================================================
# Before building the DESeq2 model, we check the variables that will be used
# as covariates (adjustment variables). Covariates are variables that affect
# gene expression but are NOT what we're studying. We adjust for them so they
# don't confuse our results about the variable of interest (cohort/infection type).
#
# Important distinction:
#   COVARIATES = batch, gender, age (adjust for these, not interested in them)
#   VARIABLE OF INTEREST = cohort (this is what we're actually testing)

# --- 4a. Check gender completeness ---
# 21 samples have missing gender data, all from CoV_other cohort.
# Decision: Remove these samples to keep gender as an adjustment variable.
# Rationale: Gender affects immune gene expression (X-linked genes, hormonal effects).
# If we don't adjust for gender, sex differences could be mistaken for infection effects.
# CoV_other drops from 60 to 39 samples — still sufficient for analysis.

cat("\nGender completeness:\n")
cat("Missing gender:", sum(is.na(meta$gender) | meta$gender == ""), "\n")
meta %>% filter(is.na(gender) | gender == "") %>% count(cohort) %>% print()

# Remove samples with missing gender
meta <- meta %>% filter(gender != "" & !is.na(gender))
counts <- counts[, rownames(meta)]
cat("After gender filter:", nrow(meta), "samples\n")

# --- 4b. Check age confounding ---
# CRITICAL FINDING: Severe age confounding between cohorts.
# Healthy controls are ~18 years old (college students), while Bacterial
# patients are ~60 years old. Without adjusting for age, thousands of
# age-related gene expression changes would be falsely attributed to infection.
# Age MUST be included as an adjustment variable in the model.

cat("\nAge distribution by cohort:\n")
meta %>%
  group_by(cohort) %>%
  summarize(
    n = n(),
    mean_age = mean(age, na.rm = TRUE),
    median_age = median(age, na.rm = TRUE),
    min_age = min(age, na.rm = TRUE),
    max_age = max(age, na.rm = TRUE)
  ) %>%
  print()

# --- 4c. Check batch confounding ---
# Batch effects are moderate. Most cohorts are predominantly Batch 1.
# COVID-19 and CoV_other have the most Batch 2 samples.
# Batch is included as an adjustment variable to remove technical variation.

cat("\nBatch distribution by cohort:\n")
print(table(meta$cohort, meta$batch))

# --- 4d. Remove samples with missing age ---
cat("\nMissing age:", sum(is.na(meta$age)), "\n")
meta <- meta %>% filter(!is.na(age))
counts <- counts[, rownames(meta)]
cat("After age filter:", nrow(meta), "samples\n")


# ==============================================================================
# 5. FINAL COHORT SUMMARY
# ==============================================================================

cat("\n=== Final Cohort Summary ===\n")
cat("Total samples:", nrow(meta), "\n")
cat("Total genes:", nrow(counts), "\n\n")
print(table(meta$cohort))

meta %>%
  group_by(cohort) %>%
  summarize(
    n = n(),
    median_age = median(age),
    pct_female = round(sum(gender == "Female") / n() * 100, 1),
    pct_batch1 = round(sum(batch == 1) / n() * 100, 1)
  ) %>%
  print()


# ==============================================================================
# 6. CREATE DESeqDataSet
# ==============================================================================

# Set up factor levels
# - cohort: healthy is reference (all comparisons are infection vs healthy)
# - Clean cohort names to remove special characters (R-safe names)
meta$cohort <- recode(meta$cohort,
                      "COVID-19" = "COVID19",
                      "CoV other" = "CoV_other")
meta$cohort <- factor(meta$cohort, 
                      levels = c("healthy", "COVID19", "Influenza", "Bacterial", "CoV_other"))
meta$batch <- factor(meta$batch)
meta$gender <- factor(meta$gender)
meta$age <- as.numeric(meta$age)

cat("Reference level:", levels(meta$cohort)[1], "\n")
cat("All comparisons will be relative to healthy controls\n")

# Design formula: ~ batch + gender + age + cohort
#
# How to read this formula:
#   "Model gene expression as a function of batch, gender, age, and cohort."
#
# COVARIATES (adjustment variables — listed first):
#   - batch: technical variation from samples processed on different days.
#            We don't care about batch differences — we want to remove them.
#   - gender: biological variation in sex-linked gene expression.
#             We don't care about male vs female differences here — we want 
#             to remove them so they don't confuse the infection comparison.
#   - age: biological variation in age-related gene expression.
#          CRITICAL in this dataset because healthy controls (median 18) are 
#          much younger than bacterial patients (median 60). Without adjusting 
#          for age, ~20% of Bacterial DE genes were actually age effects.
#
# VARIABLE OF INTEREST (listed last):
#   - cohort: the infection type. This is what we're actually testing.
#     DESeq2 first removes the effects of batch, gender, and age, then asks:
#     "After accounting for all covariates, is there still a significant
#     difference in gene expression between infection groups and healthy?"
#
# The variable of interest MUST go last in the DESeq2 formula.
# Covariates go before it. Order among covariates does not matter.

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ batch + gender + age + cohort
)

cat("DESeqDataSet created:", nrow(dds), "genes x", ncol(dds), "samples\n")
cat("Design: ~ batch + gender + age + cohort\n")


# ==============================================================================
# 7. RUN DESeq2
# ==============================================================================

# This performs three steps internally:
#   1. estimateSizeFactors — normalizes for library size differences
#   2. estimateDispersions — models gene-level variance across samples
#   3. nbinomWaldTest — tests for differential expression using negative binomial model

dds <- DESeq(dds)
cat("DESeq2 complete\n")
cat("Available comparisons:\n")
print(resultsNames(dds))


# ==============================================================================
# 8. GENE ID MAPPING (Ensembl → Gene Symbol)
# ==============================================================================

# The count matrix uses Ensembl IDs (ENSG...). We need gene symbols (ISG15, IFI27, etc.)
# for biological interpretation. org.Hs.eg.db provides the mapping.

ensembl_ids <- rownames(dds)
gene_map <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = ensembl_ids,
                                  keytype = "ENSEMBL",
                                  columns = "SYMBOL")

# Some Ensembl IDs map to multiple symbols — keep first match only
gene_map <- gene_map %>% distinct(ENSEMBL, .keep_all = TRUE)
cat("Mapped", nrow(gene_map), "of", length(ensembl_ids), "Ensembl IDs to gene symbols\n")


# ==============================================================================
# 9. EXTRACT AND ANNOTATE RESULTS — ALL COMPARISONS
# ==============================================================================

# Extract results for all four comparisons with log2FC shrinkage using apeglm.
#
# WHAT IS APEGLM?
# apeglm (Approximate Posterior Estimation for Generalized Linear Models) is a 
# method for shrinking log2 fold change estimates. Here's why we use it:
#
# The problem: Genes with very low counts produce unreliable fold change estimates.
# A gene with 2 counts in healthy and 0 counts in COVID has a huge fold change, 
# but that's based on almost no data — it's noise, not signal.
#
# What apeglm does: It uses a Bayesian approach to "shrink" fold changes. Genes 
# with low information (low counts, high variance) get pulled toward zero, while 
# genes with strong evidence (high counts, consistent across replicates) keep their 
# original fold change values largely intact.
#
# Why it matters: Without shrinkage, your top fold-change genes would be dominated 
# by low-count noise. With apeglm shrinkage, ranking genes by fold change gives 
# biologically meaningful results. This is critical for:
#   - Ranking genes for GSEA (which uses the full ranked list)
#   - Selecting top genes for heatmaps and figures
#   - Reporting fold changes in publications
#
# Note: apeglm does NOT change the p-values or significance calls. It only changes 
# the fold change estimates. The number of significant genes (padj < 0.05) stays 
# the same whether or not you use shrinkage.
#
# Alternatives: "ashr" and "normal" shrinkage are also available in DESeq2.
# apeglm is recommended by the DESeq2 authors as the default choice.

comparisons <- c("cohort_COVID19_vs_healthy", 
                 "cohort_Influenza_vs_healthy",
                 "cohort_Bacterial_vs_healthy", 
                 "cohort_CoV_other_vs_healthy")

all_results <- list()

for (comp in comparisons) {
  # Shrink log2FC estimates
  res <- lfcShrink(dds, coef = comp, type = "apeglm")
  
  # Annotate with gene symbols
  res_annotated <- as.data.frame(res) %>%
    rownames_to_column("ENSEMBL") %>%
    left_join(gene_map, by = "ENSEMBL") %>%
    filter(!is.na(padj))
  
  # Store results
  comp_name <- gsub("cohort_", "", gsub("_vs_healthy", "", comp))
  all_results[[comp_name]] <- res_annotated
  
  # Print summary
  n_sig <- sum(res_annotated$padj < 0.05)
  n_up <- sum(res_annotated$padj < 0.05 & res_annotated$log2FoldChange > 0)
  n_down <- sum(res_annotated$padj < 0.05 & res_annotated$log2FoldChange < 0)
  cat(comp_name, ": ", n_sig, " significant (", n_up, " up, ", n_down, " down)\n", sep = "")
  
  # Save to CSV
  write.csv(res_annotated, 
            paste0("results/DE_", comp_name, "_vs_healthy_age_corrected.csv"), 
            row.names = FALSE)
}


# ==============================================================================
# 10. CHECK KNOWN INTERFERON-STIMULATED GENES (ISGs)
# ==============================================================================

# ISGs are the hallmark of antiviral immune response. Checking their behavior
# validates that the analysis is capturing real biology.
# Note: This is blood transcriptomics — ISG signal is weaker in blood than at
# the infection site (e.g., nasal swabs) because the signal is diluted across
# all circulating cell types.

isgs <- c("IFI27", "ISG15", "MX1", "OAS1", "IFIT1", "RSAD2",
          "IFI44L", "SIGLEC1", "CXCL10", "IFITM3", "IFI44", "OAS2")

cat("\n=== ISG Check: COVID19 vs Healthy ===\n")
all_results$COVID19 %>%
  filter(SYMBOL %in% isgs) %>%
  arrange(padj) %>%
  dplyr::select(SYMBOL, log2FoldChange, padj) %>%
  print()


# ==============================================================================
# 11. VISUALIZATION — VOLCANO PLOTS
# ==============================================================================

# Custom volcano plot function for all comparisons.
# Design choices:
#   - Only FDR < 0.05 genes are colored as significant (no fold change cutoff 
#     for significance — FC cutoff lines are shown for visual reference only)
#   - Orange = upregulated in infection, Skyblue = upregulated in healthy
#   - Only top 5 genes by p-value and top 5 by fold change are labeled per direction
#   - This avoids overplotting while highlighting the most important genes

library(ggrepel)

make_volcano <- function(res_df, comparison_name, filename) {
  # Classify genes by significance and direction
  res_df$color <- "NS"
  res_df$color[res_df$padj < 0.05 & res_df$log2FoldChange > 0] <- "Up_Infection"
  res_df$color[res_df$padj < 0.05 & res_df$log2FoldChange < 0] <- "Up_Healthy"
  
  # Select top 5 genes to label for each direction:
  # Top 5 by lowest p-value (most statistically significant)
  top_pval_up <- res_df %>%
    filter(color == "Up_Infection") %>%
    arrange(padj) %>%
    head(5)
  
  top_pval_down <- res_df %>%
    filter(color == "Up_Healthy") %>%
    arrange(padj) %>%
    head(5)
  
  # Top 5 by highest absolute fold change (largest effect size)
  top_fc_up <- res_df %>%
    filter(color == "Up_Infection") %>%
    arrange(desc(log2FoldChange)) %>%
    head(5)
  
  top_fc_down <- res_df %>%
    filter(color == "Up_Healthy") %>%
    arrange(log2FoldChange) %>%
    head(5)
  
  # Combine and deduplicate
  genes_to_label <- bind_rows(top_pval_up, top_pval_down, top_fc_up, top_fc_down) %>%
    distinct(ENSEMBL, .keep_all = TRUE)
  
  # Add label column
  res_df$label <- ""
  res_df$label[res_df$ENSEMBL %in% genes_to_label$ENSEMBL] <- 
    res_df$SYMBOL[res_df$ENSEMBL %in% genes_to_label$ENSEMBL]
  
  # Count significant genes for subtitle
  n_up <- sum(res_df$color == "Up_Infection")
  n_down <- sum(res_df$color == "Up_Healthy")
  
  # Build plot
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = color), size = 1.2, alpha = 0.7) +
    scale_color_manual(
      values = c("Up_Infection" = "#E67E22", "Up_Healthy" = "#5DADE2", "NS" = "#BFBFBF"),
      labels = c("Up_Infection" = paste0("Up in ", comparison_name, " (", n_up, ")"),
                 "Up_Healthy" = paste0("Up in Healthy (", n_down, ")"),
                 "NS" = "Not Significant"),
      breaks = c("Up_Infection", "Up_Healthy", "NS"),
      name = ""
    ) +
    geom_text_repel(
      data = res_df %>% filter(label != ""),
      aes(label = label),
      size = 4,
      max.overlaps = 20,
      segment.color = "grey30",
      segment.size = 0.4,
      fontface = "bold.italic",
      box.padding = 0.6
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.6) +
    labs(
      title = paste0(comparison_name, " vs Healthy Controls"),
      subtitle = "GSE161731 Blood Transcriptomics | Model: ~ batch + gender + age + cohort",
      x = expression(bold(log[2]~"Fold Change")),
      y = expression(bold(-log[10]~"(FDR)"))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 11),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  ggsave(filename, p, width = 10, height = 8, dpi = 300)
  cat("Saved:", filename, "\n")
}

# Generate volcano plots for all four comparisons
make_volcano(all_results$COVID19, "COVID-19", "figures/GSE161731/05a_volcano_covid_vs_healthy.png")
make_volcano(all_results$Influenza, "Influenza", "figures/GSE161731/05b_volcano_influenza_vs_healthy.png")
make_volcano(all_results$Bacterial, "Bacterial", "figures/GSE161731/05c_volcano_bacterial_vs_healthy.png")
make_volcano(all_results$CoV_other, "CoV Other", "figures/GSE161731/05d_volcano_cov_other_vs_healthy.png")


# ==============================================================================
# 12. VISUALIZATION — DE GENE COUNTS COMPARISON
# ==============================================================================

# Bar plot comparing number of DE genes across all infections
summary_list <- lapply(names(all_results), function(name) {
  res <- all_results[[name]]
  data.frame(
    comparison = name,
    direction = c("Upregulated", "Downregulated"),
    count = c(
      sum(res$padj < 0.05 & res$log2FoldChange > 0),
      -sum(res$padj < 0.05 & res$log2FoldChange < 0)
    )
  )
})
summary_df <- bind_rows(summary_list)

ggplot(summary_df, aes(x = comparison, y = count, fill = direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Downregulated" = "#3498db", "Upregulated" = "#e74c3c")) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  labs(title = "Differentially Expressed Genes vs Healthy Controls (Age-Corrected)",
       subtitle = "GSE161731 Blood Transcriptomics | padj < 0.05 | Model: ~ batch + gender + age + cohort",
       x = "", y = "Number of DE Genes", fill = "") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/GSE161731/06_de_gene_counts.png", width = 8, height = 6, dpi = 300)


# ==============================================================================
# 13. OVERLAP ANALYSIS — SHARED VS UNIQUE DE GENES
# ==============================================================================

# Identify which DE genes are shared across infections (common host response)
# versus unique to each infection (pathogen-specific response)

sig_genes <- lapply(all_results, function(res) {
  genes <- res %>% filter(padj < 0.05) %>% pull(SYMBOL) %>% na.omit()
  as.character(unique(genes))
})

# Verify gene lists are populated
cat("\nDE genes per comparison:\n")
for (name in names(sig_genes)) {
  cat("  ", name, ":", length(sig_genes[[name]]), "genes\n")
}

# Summary of overlaps
shared_all <- Reduce(intersect, sig_genes)
cat("\n=== Gene Overlap Summary ===\n")
cat("Shared across ALL four infections:", length(shared_all), "genes\n")
cat("Shared between COVID19 and Influenza:", 
    length(intersect(sig_genes$COVID19, sig_genes$Influenza)), "genes\n")

cat("\nUnique to each infection:\n")
for (name in names(sig_genes)) {
  others <- sig_genes[names(sig_genes) != name]
  all_others <- Reduce(union, others)
  unique_genes <- setdiff(sig_genes[[name]], all_others)
  cat("  ", name, ":", length(unique_genes), "unique genes\n")
}

# --- UpSet plot ---
# UpSetR requires a data frame where each column is a set (1 = gene is in set, 0 = not)
# fromList() converts our named list of gene vectors into this format automatically

upset_input <- fromList(sig_genes)
rownames(upset_input) <- unique(unlist(sig_genes))

png("figures/GSE161731/07_upset_de_genes.png", width = 10, height = 6, units = "in", res = 300)
upset(upset_input,
      sets = c("COVID19", "Influenza", "Bacterial", "CoV_other"),
      keep.order = TRUE,
      order.by = "freq",
      sets.bar.color = c("#e74c3c", "#3498db", "#e67e22", "#9b59b6"),
      main.bar.color = "gray30",
      text.scale = 1.3,
      mainbar.y.label = "Number of Shared DE Genes",
      sets.x.label = "Total DE Genes",
      nsets = 4,
      nintersects = 15)
dev.off()
cat("UpSet plot saved\n")


# ==============================================================================
# 14. PATHWAY ANALYSIS — GO ENRICHMENT
# ==============================================================================

# --- 14a. COVID-19 unique genes ---
# Genes DE only in COVID-19 and no other infection — COVID-specific host response

covid_unique <- setdiff(sig_genes$COVID19, Reduce(union, sig_genes[names(sig_genes) != "COVID19"]))
cat("\nCOVID-19 unique DE genes:", length(covid_unique), "\n")

covid_entrez <- bitr(covid_unique, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_covid <- enrichGO(
  gene = covid_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

cat("COVID-19 unique: significant GO terms:", nrow(go_covid), "\n")

dotplot(go_covid, showCategory = 15, title = "COVID-19 Unique Genes: GO Biological Process") +
  theme(axis.text.y = element_text(size = 9))
ggsave("figures/GSE161731/08_go_covid_unique.png", width = 10, height = 8, dpi = 300)

# --- 14b. Shared infection response genes ---
# Genes DE in ALL four infections — universal host response to respiratory infection

shared_entrez <- bitr(shared_all, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_shared <- enrichGO(
  gene = shared_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

cat("Shared infection response: significant GO terms:", nrow(go_shared), "\n")

dotplot(go_shared, showCategory = 15, title = "Shared Infection Response (All 4 Groups): GO Biological Process") +
  theme(axis.text.y = element_text(size = 9))
ggsave("figures/GSE161731/09_go_shared_infection.png", width = 10, height = 8, dpi = 300)


# ==============================================================================
# 15. DOCUMENTATION — AGE CONFOUNDING ANALYSIS
# ==============================================================================

# This section documents the impact of including age as an adjustment variable.
# This analysis was prompted by observing that Bacterial infection showed an
# unexpectedly large number of DE genes (10,564 without age adjustment).
# Investigation revealed severe age confounding:
#   - Healthy controls: median age 18
#   - Bacterial patients: median age 60
#
# Without age in the model, DESeq2 cannot distinguish between genes that change
# because of infection vs genes that change because of aging. By adding age as 
# a covariate (adjustment variable), DESeq2 first accounts for age-related 
# expression changes, then tests what remains for infection effects.
#
# After adding age to the model, Bacterial DE genes dropped by ~20%,
# confirming that a substantial fraction of the original signal was age-related,
# not infection-related.

cat("\n=== Age Correction Impact ===\n\n")
cat("Model without age: ~ batch + gender + cohort\n")
cat("Model with age:    ~ batch + gender + age + cohort\n\n")
cat("Comparison           | Without Age | With Age | Removed | % Change\n")
cat("---------------------+-------------+----------+---------+---------\n")
cat("Bacterial vs Healthy |      10,564 |    8,436 |   2,128 |   20.1%\n")
cat("COVID19 vs Healthy   |       3,816 |    2,820 |     996 |   26.1%\n")
cat("\nConclusion: Age correction removed 20-26% of DE genes.\n")
cat("The age-corrected model is used for all final results.\n")


# ==============================================================================
# 16. LIMMA-VOOM PARALLEL ANALYSIS
# ==============================================================================

# limma-voom is an alternative to DESeq2 for differential expression.
#
# KEY DIFFERENCES:
#   DESeq2: models raw counts directly using a negative binomial distribution.
#   limma-voom: transforms counts using a precision-weighted approach, then
#   applies limma's linear modeling framework (originally built for microarrays).
#
# WHEN TO USE WHICH:
#   - DESeq2: great for small sample sizes (3-5 per group), more conservative
#   - limma-voom: excels with complex designs (interactions, blocking, many factors),
#     faster with large sample sizes, more flexible model specification
#   - Both give similar results for straightforward comparisons
#   - Running both and comparing builds confidence in your findings
#
# WHY COMPARE: Genes significant in BOTH methods are high-confidence hits.
# Genes significant in only one method may be borderline or method-sensitive.

library(limma)
library(edgeR)

# --- 16a. Create DGEList object ---
# edgeR's DGEList is the input format for limma-voom
# We use the same filtered counts and metadata as DESeq2

dge <- DGEList(counts = counts, group = meta$cohort)

# Calculate normalization factors (TMM method)
# TMM = Trimmed Mean of M-values, similar purpose to DESeq2's size factors
dge <- calcNormFactors(dge, method = "TMM")

cat("Normalization factors (first 5 samples):\n")
print(head(dge$samples))

# --- 16b. Build design matrix ---
# In limma, the design matrix is explicit (unlike DESeq2's formula interface)
# model.matrix() converts the formula into a numeric matrix of 0s and 1s
#
# This is the SAME model as DESeq2: ~ batch + gender + age + cohort
# Covariates (adjustment variables): batch, gender, age
# Variable of interest: cohort

design <- model.matrix(~ batch + gender + age + cohort, data = meta)

cat("Design matrix columns:\n")
print(colnames(design))

# --- 16c. voom transformation ---
# voom estimates the mean-variance relationship of the log-counts,
# then assigns a precision weight to each observation.
#
# Genes with high variance get lower weights (less influence on the model).
# This accounts for the fact that RNA-seq count data is heteroscedastic
# (variance depends on the mean), which violates the assumptions of
# standard linear models.
#
# The voom plot shows the mean-variance trend — it should show a
# decreasing trend (high-count genes have lower relative variance).

png("figures/GSE161731/10a_voom_mean_variance.png", width = 8, height = 6, units = "in", res = 300)
v <- voom(dge, design, plot = TRUE)
dev.off()

v <- voom(dge, design, plot = TRUE)
cat("voom transformation complete\n")

# --- 16d. Fit linear model ---
# lmFit fits a linear model for each gene using the precision weights from voom.
# eBayes applies empirical Bayes moderation — it borrows information across genes
# to produce more stable variance estimates (similar concept to DESeq2's dispersion
# shrinkage, but implemented differently).

fit <- lmFit(v, design)
fit <- eBayes(fit)

cat("Model coefficients available:\n")
print(colnames(fit$coefficients))

# --- 16e. Extract results for all comparisons ---
limma_all_results <- list()

for (coef_name in c("cohortCOVID19", "cohortInfluenza", "cohortBacterial", "cohortCoV_other")) {
  res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  res$ENSEMBL <- rownames(res)
  res <- res %>% left_join(gene_map, by = "ENSEMBL") %>% filter(!is.na(adj.P.Val))
  
  comp_name <- gsub("cohort", "", coef_name)
  limma_all_results[[comp_name]] <- res
  
  n_sig <- sum(res$adj.P.Val < 0.05)
  n_up <- sum(res$adj.P.Val < 0.05 & res$logFC > 0)
  n_down <- sum(res$adj.P.Val < 0.05 & res$logFC < 0)
  cat(comp_name, ": ", n_sig, " significant (", n_up, " up, ", n_down, " down)\n", sep = "")
  
  # Save to CSV
  write.csv(res, paste0("results/DE_limma_", comp_name, "_vs_healthy.csv"), row.names = FALSE)
}


# ==============================================================================
# 17. DESeq2 vs LIMMA-VOOM COMPARISON
# ==============================================================================

# Compare significance calls and fold changes between the two methods.
# High concordance = robust findings. Low concordance = borderline genes near
# the significance threshold that different statistical models handle differently.

# --- 17a. Concordance table across all comparisons ---
cat("\n=== Method Comparison Across All Comparisons ===\n\n")
cat("Comparison        | DESeq2  | limma   | Both    | Concordance\n")
cat("------------------+---------+---------+---------+------------\n")

deseq_names <- c("COVID19", "Influenza", "Bacterial", "CoV_other")

for (name in deseq_names) {
  d_sig <- all_results[[name]] %>% filter(padj < 0.05) %>% pull(ENSEMBL)
  l_sig <- limma_all_results[[name]] %>% filter(adj.P.Val < 0.05) %>% pull(ENSEMBL)
  both_sig <- intersect(d_sig, l_sig)
  concordance <- round(length(both_sig) / length(union(d_sig, l_sig)) * 100, 1)
  
  cat(sprintf("%-18s| %7d | %7d | %7d | %5.1f%%\n",
              name, length(d_sig), length(l_sig), length(both_sig), concordance))
}

# --- 17b. Fold change correlation plot (COVID19) ---
# Compare the actual fold change values between methods.
# High correlation = methods agree on direction and magnitude of change,
# even if they disagree on the significance cutoff.

comparison <- all_results$COVID19 %>%
  dplyr::select(ENSEMBL, SYMBOL, deseq2_log2FC = log2FoldChange, deseq2_padj = padj) %>%
  inner_join(
    limma_all_results$COVID19 %>% dplyr::select(ENSEMBL, limma_logFC = logFC, limma_fdr = adj.P.Val),
    by = "ENSEMBL"
  )

cor_val <- cor(comparison$deseq2_log2FC, comparison$limma_logFC, use = "complete.obs")

comparison$category <- case_when(
  comparison$deseq2_padj < 0.05 & comparison$limma_fdr < 0.05 ~ "Both significant",
  comparison$deseq2_padj < 0.05 ~ "DESeq2 only",
  comparison$limma_fdr < 0.05 ~ "limma only",
  TRUE ~ "Neither"
)

p <- ggplot(comparison, aes(x = deseq2_log2FC, y = limma_logFC, color = category)) +
  geom_point(size = 0.8, alpha = 0.5) +
  scale_color_manual(
    values = c("Both significant" = "#8E44AD", "DESeq2 only" = "#E67E22",
               "limma only" = "#5DADE2", "Neither" = "#BFBFBF"),
    name = ""
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(
    title = "DESeq2 vs limma-voom: Fold Change Comparison",
    subtitle = paste0("COVID-19 vs Healthy | Pearson r = ", round(cor_val, 3)),
    x = expression(bold(DESeq2 ~ log[2] ~ "Fold Change")),
    y = expression(bold(limma-voom ~ log[2] ~ "Fold Change"))
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p)
ggsave("figures/GSE161731/10b_deseq2_vs_limma_fc.png", p, width = 8, height = 7, dpi = 300)

cat("\nFold change correlation (Pearson r):", round(cor_val, 3), "\n")
cat("Interpretation: r > 0.8 indicates strong agreement between methods.\n")
cat("The methods agree on the biology; differences are at the significance threshold.\n")


# ==============================================================================
# 18. SESSION INFO
# ==============================================================================

cat("\n=== Session Info ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("R version:", R.version.string, "\n")
cat("DESeq2:", as.character(packageVersion("DESeq2")), "\n")
cat("apeglm:", as.character(packageVersion("apeglm")), "\n")
cat("limma:", as.character(packageVersion("limma")), "\n")
cat("edgeR:", as.character(packageVersion("edgeR")), "\n")
cat("clusterProfiler:", as.character(packageVersion("clusterProfiler")), "\n")
sessionInfo()