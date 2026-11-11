# ==============================================================================
# Differential Expression Analysis: DESeq2 and limma-voom
# Dataset: GSE161731 — Blood Transcriptomics in Respiratory Infections
# Cohorts: COVID-19, Influenza, Bacterial, CoV_other, Healthy
#
# Source: McClain, M.T., Constantine, F.J., Henao, R. et al. Dysregulated 
#         transcriptional responses to SARS-CoV-2 in the periphery. 
#         Nature Communications 12, 1079 (2021).
#         https://doi.org/10.1038/s41467-021-21289-y
#
# Date: June 2026
#
# This script performs an independent secondary analysis of publicly available
# data from GEO. Interpretations and analytical choices reflect the author's
# own analysis and do not represent the views of the original investigators.
#
# DESIGN FORMULA: ~ batch + gender + age + cohort
#
#   Covariates (adjustment variables to control for, not interested in):
#     - batch: technical variation from processing on different days
#     - gender: biological variation in sex-linked gene expression
#     - age: biological variation
#
#   Variable of interest (what we are testing):
#     - cohort: infection type vs healthy controls
#     - The variable of interest always goes LAST in the DESeq2 formula
#     - DESeq2 removes the effects of covariates first, then tests whether
#       cohort still has a significant effect on gene expression
#
# Key analytical decisions documented in this script:
#   1. Removed 2 cross-cohort subjects appearing as both healthy and CoV_other
#   2. Removed repeated subjects : kept earliest timepoint per subject
#   3. Removed samples with missing gender (all from CoV_other cohort)
#   4. Removed samples with missing age
#   5. Corrected for age confounding between cohorts
#   6. Ran both DESeq2 and limma-voom on identical independent samples
#
# PSEUDOREPLICATION NOTE:
#   The original dataset contains 14 COVID-19 patients sampled at multiple
#   timepoints (early/middle/late disease). Without correction, these repeated 
#   samples are treated as independent observations, inflating statistical 
#   significance. We remove repeated subjects and keep only the earliest 
#   timepoint to ensure all samples are independent.
#
#   ALTERNATIVE APPROACH: The original authors (McClain et al. 2021) used
#   limma::duplicateCorrelation() to account for within-subject correlation
#   while retaining all timepoints. This is statistically valid but:
#     1. Cannot be applied in DESeq2, preventing fair method comparison
#     2. Assumes a single correlation value for all genes across subjects
#     3. Retains non-independent observations rather than removing them
#   We chose the more conservative approach: one sample per subject,
#   keeping the earliest timepoint for longitudinal COVID-19 patients
#   to match the acute presentation timepoint of the other infection cohorts,
#   who were enrolled at emergency department visits.
# ==============================================================================
# ==============================================================================
# 1. SETUP: Load Libraries
# ==============================================================================

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(apeglm)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(ReactomePA)
library(UpSetR)
library(ggrepel)
library(limma)
library(edgeR)
library(RColorBrewer)
library(patchwork)

# Set working directory to project root
setwd("~/Desktop/github_repos/infectious-disease-transcriptomics")
dir.create("figures/GSE161731/deseq2_comparison", showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

counts <- read.csv("results/GSE161731_counts_filtered.csv", row.names = 1)
meta   <- read.csv("results/GSE161731_metadata_clean.csv")

cat("Raw loaded - Counts:", dim(counts)[1], "genes x", dim(counts)[2], "samples\n")
cat("Raw loaded - Metadata:", dim(meta)[1], "samples\n")


# ==============================================================================
# QC SUMMARY: ALL FILTERING STEPS APPLIED TO THIS DATASET
# ==============================================================================
#
# STEP 1 — GENE FILTERING (Python)
#   Removed genes detected in fewer than 10% of samples. 60,675 -> 28,996 genes.
#
# STEP 2 — OUTLIER SAMPLE REMOVAL (Python)
#   2 extreme outliers removed via PCA (very low library size and gene detection).
#   198 -> 196 samples.
#
# STEP 3 — CROSS-COHORT SUBJECT REMOVAL (R, Section 4a)
#   Two subjects (A34BBC, E77BE3) appeared in both healthy and CoV_other.
#   Removed both subjects entirely.
#
# STEP 4 — PSEUDOREPLICATION REMOVAL (R, Section 4b)
#   14 COVID-19 patients sampled at multiple timepoints.
#   Kept only the earliest timepoint per subject.
#   Removing pseudoreplication reduced COVID-19 DE genes by ~68%.
#
# STEP 5 — MISSING GENDER REMOVAL (R, Section 4c)
#   Removed samples with missing gender (all from CoV_other cohort).
#
# STEP 6 — MISSING AGE REMOVAL (R, Section 4d)
#   Removed samples with NA age values.
#
# STEP 7 — COHORT NAME CLEANING (R, Section 5)
#   Renamed "COVID-19" -> "COVID19" and "CoV other" -> "CoV_other".
#
# STEP 8 — DESeq2 INDEPENDENT FILTERING (automatic)
#   Genes with very low mean counts excluded automatically by results().
#
# WHAT WE DID NOT FILTER:
#   - No adapter trimming: GEO data is pre-processed counts, not raw FASTQ.
#   - No library size cutoff: outliers identified by PCA, not arbitrary threshold.
#   - No fold change cutoff: documented in Section 8 with rationale.
#   - No mitochondrial percentage filter: bulk RNA-seq only. Mitochondrial
#     filtering applies to single-cell RNA-seq where high mitochondrial content
#     identifies dead/damaged individual cells. In bulk RNA-seq, each sample
#     represents millions of cells and mitochondrial content reflects normal
#     tissue metabolic activity, not sample quality.


# ==============================================================================
# 3. MATCH SAMPLE IDS BETWEEN COUNTS AND METADATA
# ==============================================================================

meta$rna_id_fixed <- make.names(meta$rna_id)
rownames(meta) <- meta$rna_id_fixed

overlap <- sum(meta$rna_id_fixed %in% colnames(counts))
cat("Sample ID overlap:", overlap, "of", nrow(meta), "\n")

meta <- meta[colnames(counts), ]
stopifnot(all(colnames(counts) == rownames(meta)))
cat("Sample order verified\n")


# ==============================================================================
# 4. SAMPLE FILTERING AND DEDUPLICATION
# ==============================================================================

# --- 4a. Remove cross-cohort subjects ---
cross_cohort_subjects <- c("A34BBC", "E77BE3")
meta   <- meta %>% filter(!subject_id %in% cross_cohort_subjects)
counts <- counts[, rownames(meta)]
cat("After removing cross-cohort subjects:", nrow(meta), "samples\n")

# --- 4b. Remove repeated subjects (pseudoreplication) ---
subject_counts_dup <- meta %>% count(subject_id) %>% filter(n > 1)
cat("Subjects with multiple samples:", nrow(subject_counts_dup), "\n")

time_rank    <- c("early" = 1, "middle" = 2, "late" = 3)
meta$time_rank <- ifelse(
  is.na(meta$time_since_onset) | meta$time_since_onset == "",
  0, time_rank[meta$time_since_onset]
)

meta <- meta %>%
  group_by(subject_id) %>%
  arrange(time_rank) %>%
  slice(1) %>%
  ungroup() %>%
  as.data.frame()

rownames(meta) <- meta$rna_id_fixed
counts <- counts[, rownames(meta)]
cat("After deduplication:", nrow(meta), "independent samples\n")

# --- 4c. Remove samples with missing gender ---
meta   <- meta %>% filter(gender != "" & !is.na(gender))
counts <- counts[, rownames(meta)]
cat("After gender filter:", nrow(meta), "samples\n")

# --- 4d. Remove samples with missing age ---
meta   <- meta %>% filter(!is.na(age))
counts <- counts[, rownames(meta)]
cat("After age filter:", nrow(meta), "samples\n")


# ==============================================================================
# 5. PREPARE FACTORS AND FINAL COHORT SUMMARY
# ==============================================================================

meta$cohort <- recode(meta$cohort,
                      "COVID-19"  = "COVID19",
                      "CoV other" = "CoV_other")
meta$cohort  <- factor(meta$cohort,
                       levels = c("healthy","COVID19","Influenza","Bacterial","CoV_other"))
meta$batch   <- factor(meta$batch)
meta$gender  <- factor(meta$gender)
meta$age     <- as.numeric(meta$age)

cat("\n=== Final Independent Cohort Summary ===\n")
print(table(meta$cohort))
cat("All samples independent:", nrow(meta) == length(unique(meta$subject_id)), "\n")
cat("Match:", all(colnames(counts) == rownames(meta)), "\n")


# ==============================================================================
# 6. CREATE DESeqDataSet AND RUN DESeq2
# ==============================================================================

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ batch + gender + age + cohort
)
dds <- DESeq(dds)
cat("DESeq2 complete\n")
print(resultsNames(dds))


# ==============================================================================
# 7. GENE ID MAPPING (Ensembl -> Gene Symbol)
# ==============================================================================

ensembl_ids <- rownames(dds)
gene_map <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys     = ensembl_ids,
                                  keytype  = "ENSEMBL",
                                  columns  = "SYMBOL")
gene_map <- gene_map %>% distinct(ENSEMBL, .keep_all = TRUE)
cat("Mapped", nrow(gene_map), "of", length(ensembl_ids), "Ensembl IDs\n")


# ==============================================================================
# 8. EXTRACT AND ANNOTATE RESULTS — DESeq2
# ==============================================================================
#
# apeglm shrinks unreliable log2FC estimates toward zero without changing
# p-values or significance calls. Essential for ranking genes reliably.
#
# We use FDR < 0.05 only (no |log2FC| cutoff) because:
#   1. DESeq2 authors recommend against hard FC cutoffs
#   2. Small FC genes (TFs, signaling) can be biologically important
#   3. Volcano dashed lines at |log2FC|=1 provide visual context without
#      discarding data

comparisons <- c("cohort_COVID19_vs_healthy",
                 "cohort_Influenza_vs_healthy",
                 "cohort_Bacterial_vs_healthy",
                 "cohort_CoV_other_vs_healthy")

deseq2_results <- list()

cat("\n=== DESeq2 Results (Independent Samples) ===\n")
for (comp in comparisons) {
  res <- lfcShrink(dds, coef = comp, type = "apeglm")
  res_ann <- as.data.frame(res) %>%
    rownames_to_column("ENSEMBL") %>%
    left_join(gene_map, by = "ENSEMBL") %>%
    filter(!is.na(padj))
  comp_name <- gsub("cohort_", "", gsub("_vs_healthy", "", comp))
  deseq2_results[[comp_name]] <- res_ann
  n_sig  <- sum(res_ann$padj < 0.05)
  n_up   <- sum(res_ann$padj < 0.05 & res_ann$log2FoldChange > 0)
  n_down <- sum(res_ann$padj < 0.05 & res_ann$log2FoldChange < 0)
  cat(comp_name, ": ", n_sig, " (", n_up, " up, ", n_down, " down)\n", sep="")
  write.csv(res_ann,
            paste0("results/DE_deseq2_", comp_name, "_vs_healthy.csv"),
            row.names = FALSE)
}


# ==============================================================================
# 9. LIMMA-VOOM ANALYSIS
# ==============================================================================

dge    <- DGEList(counts = counts, group = meta$cohort)
dge    <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~ batch + gender + age + cohort, data = meta)

png("figures/GSE161731/deseq2_comparison/10a_voom_mean_variance.png",
    width=8, height=6, units="in", res=300)
v <- voom(dge, design, plot=TRUE)
dev.off()
v <- voom(dge, design, plot=TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)

limma_results <- list()
cat("\n=== limma-voom Results (Independent Samples) ===\n")
for (coef_name in c("cohortCOVID19","cohortInfluenza",
                    "cohortBacterial","cohortCoV_other")) {
  res <- topTable(fit, coef=coef_name, number=Inf, sort.by="none")
  res$ENSEMBL <- rownames(res)
  res <- res %>% left_join(gene_map, by="ENSEMBL") %>% filter(!is.na(adj.P.Val))
  comp_name <- gsub("cohort","",coef_name)
  limma_results[[comp_name]] <- res
  n_sig  <- sum(res$adj.P.Val < 0.05)
  n_up   <- sum(res$adj.P.Val < 0.05 & res$logFC > 0)
  n_down <- sum(res$adj.P.Val < 0.05 & res$logFC < 0)
  cat(comp_name, ": ", n_sig, " (", n_up, " up, ", n_down, " down)\n", sep="")
  write.csv(res, paste0("results/DE_limma_", comp_name, "_vs_healthy.csv"),
            row.names=FALSE)
}


# ==============================================================================
# 10. DESeq2 vs LIMMA-VOOM COMPARISON
# ==============================================================================

cat("\n=== Method Comparison  ===\n\n")
cat("Comparison        | DESeq2  | limma   | Both    | Concordance\n")
cat("------------------+---------+---------+---------+------------\n")
for (name in c("COVID19","Influenza","Bacterial","CoV_other")) {
  d_sig <- deseq2_results[[name]] %>% filter(padj < 0.05) %>% pull(ENSEMBL)
  l_sig <- limma_results[[name]]  %>% filter(adj.P.Val < 0.05) %>% pull(ENSEMBL)
  both  <- intersect(d_sig, l_sig)
  conc  <- round(length(both)/length(union(d_sig,l_sig))*100, 1)
  cat(sprintf("%-18s| %7d | %7d | %7d | %5.1f%%\n",
              name, length(d_sig), length(l_sig), length(both), conc))
}

comp_df <- deseq2_results$COVID19 %>%
  dplyr::select(ENSEMBL, SYMBOL, deseq2_log2FC=log2FoldChange, deseq2_padj=padj) %>%
  inner_join(limma_results$COVID19 %>%
               dplyr::select(ENSEMBL, limma_logFC=logFC, limma_fdr=adj.P.Val),
             by="ENSEMBL")
cor_val <- cor(comp_df$deseq2_log2FC, comp_df$limma_logFC, use="complete.obs")

comp_df$category <- case_when(
  comp_df$deseq2_padj<0.05 & comp_df$limma_fdr<0.05 ~ "Both significant",
  comp_df$deseq2_padj<0.05 ~ "DESeq2 only",
  comp_df$limma_fdr<0.05   ~ "limma only",
  TRUE ~ "Neither")

p_fc <- ggplot(comp_df, aes(x=deseq2_log2FC, y=limma_logFC, color=category)) +
  geom_point(size=0.8, alpha=0.5) +
  scale_color_manual(
    values=c("Both significant"="#8E44AD","DESeq2 only"="#E67E22",
             "limma only"="#5DADE2","Neither"="#BFBFBF"), name="") +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="black", linewidth=0.5) +
  labs(title="DESeq2 vs limma-voom: Fold Change Comparison",
       subtitle=paste0("COVID-19 vs Healthy | Pearson r = ", round(cor_val,3)),
       x=expression(bold(DESeq2~log[2]~"Fold Change")),
       y=expression(bold(limma-voom~log[2]~"Fold Change"))) +
  theme_minimal(base_size=13) +
  theme(legend.position="top", plot.title=element_text(face="bold",size=14),
        axis.title=element_text(face="bold"), panel.grid.minor=element_blank())
print(p_fc)
ggsave("figures/GSE161731/deseq2_comparison/10b_deseq2_vs_limma_fc.png", p_fc,
       width=8, height=7, dpi=300)
cat("\nFold change correlation (Pearson r):", round(cor_val,3), "\n")


# ==============================================================================
# 11. CHECK KNOWN INTERFERON-STIMULATED GENES (ISGs)
# ==============================================================================

isgs <- c("IFI27","ISG15","MX1","OAS1","IFIT1","RSAD2",
          "IFI44L","SIGLEC1","CXCL10","IFITM3","IFI44","OAS2")
cat("\n=== ISG Check: COVID19 vs Healthy ===\n")
deseq2_results$COVID19 %>%
  filter(SYMBOL %in% isgs) %>%
  arrange(padj) %>%
  dplyr::select(SYMBOL, log2FoldChange, padj) %>%
  print()


# ==============================================================================
# 12. VISUALIZATION — VOLCANO PLOTS
# ==============================================================================

make_volcano <- function(res_df, comparison_name, filename) {
  res_df$color <- "NS"
  res_df$color[res_df$padj<0.05 & res_df$log2FoldChange>0] <- "Up_Infection"
  res_df$color[res_df$padj<0.05 & res_df$log2FoldChange<0] <- "Up_Healthy"
  
  top_pu <- res_df %>% filter(color=="Up_Infection") %>% arrange(padj) %>% head(5)
  top_pd <- res_df %>% filter(color=="Up_Healthy")   %>% arrange(padj) %>% head(5)
  top_fu <- res_df %>% filter(color=="Up_Infection") %>% arrange(desc(log2FoldChange)) %>% head(5)
  top_fd <- res_df %>% filter(color=="Up_Healthy")   %>% arrange(log2FoldChange) %>% head(5)
  
  lab <- bind_rows(top_pu,top_pd,top_fu,top_fd) %>% distinct(ENSEMBL,.keep_all=TRUE)
  res_df$label <- ""
  res_df$label[res_df$ENSEMBL %in% lab$ENSEMBL] <- res_df$SYMBOL[res_df$ENSEMBL %in% lab$ENSEMBL]
  
  n_up   <- sum(res_df$color=="Up_Infection")
  n_down <- sum(res_df$color=="Up_Healthy")
  
  p <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(aes(color=color), size=1.2, alpha=0.7) +
    scale_color_manual(
      values=c("Up_Infection"="#E67E22","Up_Healthy"="#5DADE2","NS"="#BFBFBF"),
      labels=c("Up_Infection"=paste0("Upregulated in ",comparison_name," (",n_up,")"),
               "Up_Healthy"  =paste0("Upregulated in Healthy (",n_down,")"),
               "NS"="Not Significant"),
      breaks=c("Up_Infection","Up_Healthy","NS"), name="") +
    geom_text_repel(data=res_df %>% filter(label!=""), aes(label=label),
                    size=4, max.overlaps=20, segment.color="grey30",
                    segment.size=0.4, fontface="bold.italic", box.padding=0.6) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black", linewidth=0.6) +
    geom_vline(xintercept=c(-1,1),       linetype="dashed", color="black", linewidth=0.6) +
    labs(title=paste0(comparison_name," vs Healthy Controls"),
         subtitle="GSE161731 | Model: ~ batch + gender + age + cohort",
         x=expression(bold(log[2]~"Fold Change")),
         y=expression(bold(-log[10]~"(FDR)"))) +
    theme_minimal(base_size=14) +
    theme(legend.position="top",
          legend.text=element_text(size=12,face="bold"),
          plot.title=element_text(face="bold",size=16),
          plot.subtitle=element_text(size=11),
          axis.title=element_text(face="bold",size=13),
          panel.grid.minor=element_blank())
  print(p)
  ggsave(filename, p, width=10, height=8, dpi=300)
  cat("Saved:", filename, "\n")
}

make_volcano(deseq2_results$COVID19,   "COVID-19",  "figures/GSE161731/deseq2_comparison/05a_volcano_covid_vs_healthy.png")
make_volcano(deseq2_results$Influenza, "Influenza",  "figures/GSE161731/deseq2_comparison/05b_volcano_influenza_vs_healthy.png")
make_volcano(deseq2_results$Bacterial, "Bacterial",  "figures/GSE161731/deseq2_comparison/05c_volcano_bacterial_vs_healthy.png")
make_volcano(deseq2_results$CoV_other, "CoV Other",  "figures/GSE161731/deseq2_comparison/05d_volcano_cov_other_vs_healthy.png")


# ==============================================================================
# 13. VISUALIZATION — DE GENE COUNTS
# ==============================================================================

summary_list <- lapply(names(deseq2_results), function(name) {
  res <- deseq2_results[[name]]
  data.frame(
    comparison = name,
    direction  = c("Upregulated","Downregulated"),
    count      = c(sum(res$padj<0.05 & res$log2FoldChange>0),
                   -sum(res$padj<0.05 & res$log2FoldChange<0)))
})
summary_df <- bind_rows(summary_list)

ggplot(summary_df, aes(x=comparison, y=count, fill=direction)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("Downregulated"="#3498db","Upregulated"="#e74c3c")) +
  geom_hline(yintercept=0, linewidth=0.5) +
  labs(title="DE Genes vs Healthy Controls",
       subtitle="GSE161731 | padj < 0.05 | Model: ~ batch + gender + age + cohort",
       x="", y="Number of DE Genes", fill="") +
  theme_minimal(base_size=13) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("figures/GSE161731/deseq2_comparison/06_de_gene_counts.png", width=8, height=6, dpi=300)


# ==============================================================================
# 14. OVERLAP ANALYSIS — SHARED VS UNIQUE DE GENES
# ==============================================================================

sig_genes <- lapply(deseq2_results, function(res) {
  as.character(unique(na.omit(res$SYMBOL[res$padj < 0.05])))
})

cat("\nDE genes per comparison:\n")
for (name in names(sig_genes)) cat("  ",name,":",length(sig_genes[[name]]),"genes\n")

shared_all <- Reduce(intersect, sig_genes)
cat("\nShared across ALL four infections:", length(shared_all), "genes\n")

cat("\nUnique to each infection:\n")
for (name in names(sig_genes)) {
  others <- sig_genes[names(sig_genes)!=name]
  cat("  ",name,":",length(setdiff(sig_genes[[name]],Reduce(union,others))),"unique genes\n")
}

# UpSet plot
# Build binary matrix manually — more reliable than fromList() for PNG output
all_genes   <- unique(unlist(sig_genes))
upset_input <- data.frame(
  COVID19   = as.integer(all_genes %in% sig_genes$COVID19),
  Influenza = as.integer(all_genes %in% sig_genes$Influenza),
  Bacterial = as.integer(all_genes %in% sig_genes$Bacterial),
  CoV_other = as.integer(all_genes %in% sig_genes$CoV_other)
)

dir.create("figures/GSE161731/deseq2_comparison", showWarnings=FALSE, recursive=TRUE)
svg("figures/GSE161731/deseq2_comparison/07_upset_de_genes.svg", width=14, height=8)
upset(upset_input,
      sets           = c("COVID19","Influenza","Bacterial","CoV_other"),
      keep.order     = TRUE,
      order.by       = "freq",
      sets.bar.color = c("#e74c3c","#3498db","#e67e22","#9b59b6"),
      main.bar.color = "gray30",
      text.scale     = 1.5,
      mainbar.y.label = "Number of Shared DE Genes",
      sets.x.label   = "Total DE Genes",
      nsets          = 4,
      nintersects    = 15)
print(upset_input)
dev.off()
cat("UpSet plot saved\n")

# ==============================================================================
# 15: edgeR DIFFERENTIAL EXPRESSION
# ==============================================================================
#
# edgeR uses the negative binomial distribution like DESeq2, but differs in
# how it estimates dispersion:
#
# DESeq2:  shrinks gene-wise dispersions toward a fitted trend across all genes
# edgeR:   estimates three levels of dispersion:
#            - Common dispersion: one value shared across all genes
#            - Trended dispersion: dispersion as a function of mean expression
#            - Tagwise dispersion: gene-specific, shrunk toward the trended value
#
# Both use empirical Bayes shrinkage but implement it differently.
# edgeR uses quasi-likelihood F-tests (glmQLFTest) which are more conservative
# and better control type I error than likelihood ratio tests.
#
# INPUT:  counts_qc, meta_qc (same as DESeq2 and limma-voom)
# OUTPUT: DE results per comparison, saved as CSV
# ==============================================================================

cat("\n========== edgeR Analysis ==========\n")

# ------------------------------------------------------------------------------
# Build DGEList and normalize
# ------------------------------------------------------------------------------
# Re-use dge_qc from limma-voom section (already TMM normalized)
# Build DGEList — same input as limma-voom
dge_qc <- DGEList(counts = counts, group = meta$cohort)
dge_qc <- calcNormFactors(dge_qc, method = "TMM")
design <- model.matrix(~ batch + gender + age + cohort, data = meta)
meta_edger <- meta
meta_edger$age <- scale(meta_edger$age)[, 1]
design <- model.matrix(~ batch + gender + age + cohort, data = meta_edger)

# Estimate dispersions
# Three levels: common → trended → tagwise (each refines the previous)
dge_edger <- estimateDisp(dge_qc, design, robust = TRUE)

cat("Common dispersion (BCV):",
    round(sqrt(dge_edger$common.dispersion), 3), "\n")
# BCV (biological coefficient of variation) = sqrt(dispersion)
# Typical values: 0.1-0.4 for well-controlled experiments
# Higher BCV = more biological variability between samples

# ------------------------------------------------------------------------------
# Fit quasi-likelihood model
# ------------------------------------------------------------------------------
# glmQLFit: fits a negative binomial GLM with quasi-likelihood dispersion
# More conservative than glmLRT (likelihood ratio test)
# Recommended for most RNA-seq experiments
fit_edger <- glmQLFit(dge_edger, design, robust = TRUE)

# ------------------------------------------------------------------------------
# DE testing — one comparison per infection group vs healthy
# ------------------------------------------------------------------------------
comparisons <- c("cohortCOVID19", "cohortInfluenza",
                 "cohortBacterial", "cohortCoV_other")

comp_labels <- c(
  "cohortCOVID19"   = "COVID19_vs_healthy",
  "cohortInfluenza" = "Influenza_vs_healthy",
  "cohortBacterial" = "Bacterial_vs_healthy",
  "cohortCoV_other" = "CoV_other_vs_healthy"
)

edger_results <- list()

for (coef in comparisons) {
  label <- comp_labels[coef]
  
  # glmQLFTest: quasi-likelihood F-test for the specified coefficient
  qlf <- glmQLFTest(fit_edger, coef = coef)
  
  # Extract results — topTags returns sorted by p-value
  res <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "PValue")$table
  res$ENSEMBL <- rownames(res)
  
  # Count significant genes
  n_sig <- sum(res$FDR < 0.05, na.rm = TRUE)
  cat(label, "— significant genes (FDR < 0.05):", n_sig, "\n")
  
  edger_results[[label]] <- res
  
  # Save results
  write.csv(res,
            file = paste0("results/DE_edger_", label, ".csv"),
            row.names = FALSE)
}

# ==============================================================================
# 15. GSEA REACTOME PATHWAY ANALYSIS — DESeq2 vs limma-voom
# ==============================================================================
#
# WHY GSEA INSTEAD OF ORA?
#   ORA takes a binary gene list that discards fold change information and is sensitive
#   to arbitrary cutoffs. GSEA uses ALL genes ranked by fold change, testing
#   whether pathway genes cluster at the top or bottom of the ranked list.
#   More sensitive and not dependent on a significance cutoff.
#
# WHY REACTOME?
#   More comprehensive and better curated than GO for human immune biology.
#
# WHY RUN GSEA FROM BOTH METHODS?
#   The core question of this repo is whether DESeq2 and limma-voom lead to
#   the same biological conclusions. Gene level concordance (Section 10) shows
#   agreement on individual genes. Pathway level concordance asks a more
#   important question: do both methods identify the same biological processes?
#   Pathway level concordance is expected to be higher than gene level because
#   pathways are more robust : two methods may disagree on borderline genes
#   but still rank the same biological pathways as most enriched.
#
# RANKING METRIC: sign(log2FC) * -log10(pvalue)
#   Weights both direction and statistical confidence.
#   DESeq2: uses log2FoldChange and pvalue columns
#   limma-voom: uses logFC and P.Value columns
#
# PLOT DESIGN:
#   Side by side: top 5 activated + top 5 suppressed from each method.
#   Shared pathways highlighted to show pathway level concordance.

# --- Helper: build ranked list from DESeq2 results ---
make_ranked_list_deseq2 <- function(res_df) {
  ranked <- res_df %>%
    filter(!is.na(log2FoldChange) & !is.na(pvalue) & !is.na(SYMBOL)) %>%
    mutate(rank_metric = sign(log2FoldChange) * -log10(pvalue)) %>%
    arrange(desc(rank_metric))
  entrez_map <- bitr(ranked$SYMBOL, fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ranked %>%
    left_join(entrez_map, by = "SYMBOL") %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    { setNames(.$rank_metric, .$ENTREZID) }
}

# --- Helper: build ranked list from limma-voom results ---
# limma uses logFC and P.Value (note: capital V in P.Value)
make_ranked_list_limma <- function(res_df) {
  ranked <- res_df %>%
    filter(!is.na(logFC) & !is.na(P.Value) & !is.na(SYMBOL)) %>%
    mutate(rank_metric = sign(logFC) * -log10(P.Value)) %>%
    arrange(desc(rank_metric))
  entrez_map <- bitr(ranked$SYMBOL, fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ranked %>%
    left_join(entrez_map, by = "SYMBOL") %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    { setNames(.$rank_metric, .$ENTREZID) }
}

# --- Helper: build ranked list from edgeR results ---
# edgeR uses logFC and PValue (capital P, capital V — different from limma)
make_ranked_list_edger <- function(res_df) {
  ranked <- res_df %>%
    filter(!is.na(logFC) & !is.na(PValue)) %>%
    mutate(rank_metric = sign(logFC) * -log10(PValue)) %>%
    arrange(desc(rank_metric))
  
  entrez_map <- bitr(res_df$ENSEMBL, fromType = "ENSEMBL",
                     toType   = "ENTREZID",
                     OrgDb    = org.Hs.eg.db)
  
  ranked %>%
    left_join(entrez_map, by = "ENSEMBL") %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    { setNames(.$rank_metric, .$ENTREZID) }
}

# --- Helper: run GSEA safely ---
run_gsea <- function(gene_list, label) {
  set.seed(42)
  tryCatch(
    gsePathway(geneList = gene_list, organism = "human",
               pvalueCutoff = 0.05, pAdjustMethod = "BH",
               verbose = FALSE, eps = 0),
    error = function(e) {
      cat("  Error for", label, ":", e$message, "\n")
      return(NULL)
    }
  )
}

# --- Updated helper: three-panel GSEA comparison plot ---
make_gsea_comparison_plot <- function(gsea_d2, gsea_lm, gsea_er,
                                      comparison_name, filename) {
  
  if (is.null(gsea_d2) && is.null(gsea_lm) && is.null(gsea_er)) {
    cat("No significant pathways for any method:", comparison_name, "\n")
    return(NULL)
  }
  
  # Extract top 5 activated + 5 suppressed from each method
  get_top <- function(gsea_res, method_label) {
    if (is.null(gsea_res) || nrow(gsea_res@result) == 0) return(NULL)
    top_act <- gsea_res@result %>% filter(NES > 0) %>%
      arrange(p.adjust) %>% head(5) %>% mutate(direction = "Activated")
    top_sup <- gsea_res@result %>% filter(NES < 0) %>%
      arrange(p.adjust) %>% head(5) %>% mutate(direction = "Suppressed")
    bind_rows(top_act, top_sup) %>% mutate(method = method_label)
  }
  
  d2_data <- get_top(gsea_d2, "DESeq2")
  lm_data <- get_top(gsea_lm, "limma-voom")
  er_data <- get_top(gsea_er, "edgeR")
  
  all_data_list <- list(d2_data, lm_data, er_data)
  all_data_list <- all_data_list[!sapply(all_data_list, is.null)]
  
  if (length(all_data_list) == 0) return(NULL)
  
  # Identify pathways shared across ALL available methods
  path_lists <- lapply(all_data_list, function(x) x$Description)
  shared_all <- Reduce(intersect, path_lists)
  
  # Pathways shared between at least two methods
  shared_any <- unique(unlist(lapply(
    combn(path_lists, 2, simplify = FALSE),
    function(x) intersect(x[[1]], x[[2]])
  )))
  
  cat("  Pathways shared across all methods:", length(shared_all), "\n")
  cat("  Pathways shared across any two methods:", length(shared_any), "\n")
  if (length(shared_all) > 0) {
    cat("  Shared all:", paste(str_trunc(shared_all, 40),
                               collapse = "\n              "), "\n")
  }
  
  all_data <- bind_rows(all_data_list) %>%
    mutate(
      Description_short = str_trunc(Description, 45),
      shared_all = Description %in% shared_all,
      shared_any = Description %in% shared_any,
      # Shape: diamond if shared across all, triangle if shared by two, circle if unique
      shape_cat = case_when(
        shared_all ~ "Shared (all 3)",
        shared_any ~ "Shared (2 methods)",
        TRUE       ~ "Method-specific"
      ),
      # Reorder by NES within each method panel
      Description_short = factor(
        Description_short,
        levels = unique(Description_short[order(NES)])
      ),
      # Set method factor order
      method = factor(method, levels = c("DESeq2", "limma-voom", "edgeR"))
    )
  
  p <- ggplot(all_data,
              aes(x = NES, y = Description_short,
                  size = setSize, color = p.adjust)) +
    geom_point(aes(shape = shape_cat), alpha = 0.85) +
    scale_shape_manual(
      values = c("Shared (all 3)"    = 18,   # diamond
                 "Shared (2 methods)"= 17,   # triangle
                 "Method-specific"   = 16),  # circle
      name = ""
    ) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    scale_color_gradient(low = "#C0392B", high = "#2980B9", name = "FDR") +
    scale_size_continuous(range = c(2, 8), name = "Gene set\nsize") +
    facet_wrap(~ method, scales = "free_y", ncol = 3) +
    labs(
      title    = paste0("GSEA Reactome: ", comparison_name,
                        " vs Healthy — DESeq2 | limma-voom | edgeR"),
      subtitle = paste0("GSE161731 | FDR < 0.05 | Top 5 activated + top 5 suppressed | ",
                        "Diamond = shared all 3 | Triangle = shared 2 methods"),
      x = "Normalized Enrichment Score (NES)",
      y = ""
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 9, color = "grey40"),
      strip.text       = element_text(face = "bold", size = 11),
      axis.text.y      = element_text(size = 8),
      axis.title.x     = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )
  
  print(p)
  ggsave(filename, p, width = 18, height = 7, dpi = 300)
  cat("Saved:", filename, "\n")
}


# ==============================================================================
#  GSEA LOOP — runs DESeq2, limma-voom, AND edgeR
# ==============================================================================

gsea_deseq2_results <- list()
gsea_limma_results  <- list()
gsea_edger_results  <- list()

comp_labels <- c(COVID19   = "COVID-19",
                 Influenza = "Influenza",
                 Bacterial = "Bacterial",
                 CoV_other = "CoV Other")

comp_files <- c(
  COVID19   = "figures/GSE161731/deseq2_comparison/08a_gsea_reactome_covid.png",
  Influenza = "figures/GSE161731/deseq2_comparison/08b_gsea_reactome_influenza.png",
  Bacterial = "figures/GSE161731/deseq2_comparison/08c_gsea_reactome_bacterial.png",
  CoV_other = "figures/GSE161731/deseq2_comparison/08d_gsea_reactome_cov_other.png"
)

# Map comparison names to edgeR results list keys
edger_comp_keys <- c(
  COVID19   = "COVID19_vs_healthy",
  Influenza = "Influenza_vs_healthy",
  Bacterial = "Bacterial_vs_healthy",
  CoV_other = "CoV_other_vs_healthy"
)

for (name in names(comp_labels)) {
  cat("\n=== GSEA Reactome:", name, "===\n")
  
  # --- DESeq2 ranked list ---
  gl_d2  <- make_ranked_list_deseq2(deseq2_results[[name]])
  cat("DESeq2 ranked genes:", length(gl_d2), "\n")
  gsea_d2 <- run_gsea(gl_d2, paste(name, "DESeq2"))
  gsea_deseq2_results[[name]] <- gsea_d2
  if (!is.null(gsea_d2)) {
    cat("DESeq2:", nrow(gsea_d2), "significant pathways\n")
    write.csv(gsea_d2@result,
              paste0("results/GSEA_deseq2_", name, "_vs_healthy.csv"),
              row.names = FALSE)
  }
  
  # --- limma-voom ranked list ---
  gl_lm  <- make_ranked_list_limma(limma_results[[name]])
  cat("limma-voom ranked genes:", length(gl_lm), "\n")
  gsea_lm <- run_gsea(gl_lm, paste(name, "limma"))
  gsea_limma_results[[name]] <- gsea_lm
  if (!is.null(gsea_lm)) {
    cat("limma-voom:", nrow(gsea_lm), "significant pathways\n")
    write.csv(gsea_lm@result,
              paste0("results/GSEA_limma_", name, "_vs_healthy.csv"),
              row.names = FALSE)
  }
  
  # --- edgeR ranked list ---
  er_key <- edger_comp_keys[[name]]
  er_res <- edger_results[[er_key]]
  gl_er  <- make_ranked_list_edger(er_res)
  cat("edgeR ranked genes:", length(gl_er), "\n")
  gsea_er <- run_gsea(gl_er, paste(name, "edgeR"))
  gsea_edger_results[[name]] <- gsea_er
  if (!is.null(gsea_er)) {
    cat("edgeR:", nrow(gsea_er), "significant pathways\n")
    write.csv(gsea_er@result,
              paste0("results/GSEA_edger_", name, "_vs_healthy.csv"),
              row.names = FALSE)
  }
  
  # --- Three-panel comparison plot ---
  make_gsea_comparison_plot(gsea_d2, gsea_lm, gsea_er,
                            comp_labels[[name]], comp_files[[name]])
}

# --- Summary: pathway-level concordance across all three methods ---
cat("\n=== GSEA Pathway Concordance (three-way) ===\n\n")
cat(sprintf("%-18s | %7s | %7s | %7s | %8s | %8s | %8s | %6s\n",
            "Comparison", "DESeq2", "limma", "edgeR",
            "D2∩limma", "D2∩edgeR", "lm∩edgeR", "All 3"))
cat(strrep("-", 90), "\n")

for (name in names(comp_labels)) {
  d2_paths <- if (!is.null(gsea_deseq2_results[[name]]))
    gsea_deseq2_results[[name]]@result$Description else character(0)
  lm_paths <- if (!is.null(gsea_limma_results[[name]]))
    gsea_limma_results[[name]]@result$Description  else character(0)
  er_paths <- if (!is.null(gsea_edger_results[[name]]))
    gsea_edger_results[[name]]@result$Description  else character(0)
  
  shared_dl  <- length(intersect(d2_paths, lm_paths))
  shared_de  <- length(intersect(d2_paths, er_paths))
  shared_le  <- length(intersect(lm_paths, er_paths))
  shared_all <- length(Reduce(intersect, list(d2_paths, lm_paths, er_paths)))
  
  cat(sprintf("%-18s | %7d | %7d | %7d | %8d | %8d | %8d | %6d\n",
              name, length(d2_paths), length(lm_paths), length(er_paths),
              shared_dl, shared_de, shared_le, shared_all))
}


# ==============================================================================
# 16. HEATMAP — DESeq2 VST vs limma-voom log2-CPM
# ==============================================================================
#
# Top 50 DE genes from COVID19 vs Healthy shown across ALL samples.
# DESeq2 VST: sophisticated variance stabilization that is recommended for visualization.
# limma-voom log2-CPM: the E matrix from voom() : same values used in the model.
# Both produce log-scale stabilized values. Z-scoring per gene centers each
# gene so colors reflect relative expression, not absolute levels.
#
# NOTE: dev.off() called before each PNG to ensure clean graphics device.

top50_covid <- deseq2_results$COVID19 %>%
  filter(!is.na(padj) & !is.na(SYMBOL)) %>%
  arrange(padj) %>% head(50)
cat("\nTop 50 genes for heatmap, padj range:", range(top50_covid$padj), "\n")

col_annotation <- data.frame(
  Cohort = meta$cohort,
  Batch  = meta$batch,
  row.names = rownames(meta))
annotation_colors <- list(
  Cohort = c(healthy="#2ecc71",COVID19="#e74c3c",Influenza="#3498db",
             Bacterial="#e67e22",CoV_other="#9b59b6"),
  Batch  = c("1"="#95a5a6","2"="#2c3e50"))
heatmap_colors <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
heatmap_breaks <- seq(-2.5, 2.5, length.out=101)

# DESeq2 VST heatmap
vst_counts  <- assay(vst(dds, blind=FALSE))
vst_top50   <- vst_counts[top50_covid$ENSEMBL, ]
rownames(vst_top50) <- top50_covid$SYMBOL
vst_top50_z <- t(scale(t(vst_top50)))

if (dev.cur() > 1) dev.off()  # close any open device safely
png("figures/GSE161731/deseq2_comparison/09a_heatmap_deseq2_vst.png",
    width=14, height=12, units="in", res=300)
pheatmap(vst_top50_z,
         annotation_col=col_annotation, annotation_colors=annotation_colors,
         show_colnames=FALSE, cluster_rows=TRUE, cluster_cols=TRUE,
         clustering_method="ward.D2", color=heatmap_colors, breaks=heatmap_breaks,
         fontsize_row=8, border_color=NA,
         main="Top 50 DE Genes: COVID-19 vs Healthy\nDESeq2 VST normalization (Z-scored, ward.D2 clustering)")
dev.off()
cat("DESeq2 VST heatmap saved\n")

# limma-voom log2-CPM heatmap
voom_top50   <- v$E[top50_covid$ENSEMBL, ]
rownames(voom_top50) <- top50_covid$SYMBOL
voom_top50_z <- t(scale(t(voom_top50)))

if (dev.cur() > 1) dev.off()  # close any open device safely
png("figures/GSE161731/deseq2_comparison/09b_heatmap_limma_voom.png",
    width=14, height=12, units="in", res=300)
pheatmap(voom_top50_z,
         annotation_col=col_annotation, annotation_colors=annotation_colors,
         show_colnames=FALSE, cluster_rows=TRUE, cluster_cols=TRUE,
         clustering_method="ward.D2", color=heatmap_colors, breaks=heatmap_breaks,
         fontsize_row=8, border_color=NA,
         main="Top 50 DE Genes: COVID-19 vs Healthy\nlimma-voom log2-CPM normalization (Z-scored, ward.D2 clustering)")
dev.off()
cat("limma-voom heatmap saved\n")

cat("\nHeatmap verification:\n")
cat("VST z-score range:", round(range(vst_top50_z,  na.rm=TRUE), 2), "\n")
cat("Voom z-score range:", round(range(voom_top50_z, na.rm=TRUE), 2), "\n")
cat("Correlation:", round(cor(as.vector(vst_top50_z),
                              as.vector(voom_top50_z), use="complete.obs"), 3), "\n")

# ==============================================================================
# edgeR HEATMAP
# ==============================================================================
#
# edgeR uses TMM-normalized log2-CPM values for visualization.
# Same top 50 genes as DESeq2 and limma-voom heatmaps for direct comparison.
# log2-CPM computed from the TMM-normalized DGEList object (dge_edger).

# Compute log2-CPM from TMM-normalized edgeR object
log_cpm_edger <- cpm(dge_edger, log = TRUE, prior.count = 1)

# Subset to top 50 COVID-19 DE genes
edger_top50   <- log_cpm_edger[top50_covid$ENSEMBL, ]
rownames(edger_top50) <- top50_covid$SYMBOL

# Z-score per gene (same as DESeq2 and limma-voom heatmaps)
edger_top50_z <- t(scale(t(edger_top50)))

cat("\nedgeR heatmap z-score range:",
    round(range(edger_top50_z, na.rm = TRUE), 2), "\n")
cat("Correlation with DESeq2 VST:",
    round(cor(as.vector(vst_top50_z),
              as.vector(edger_top50_z), use = "complete.obs"), 3), "\n")
cat("Correlation with limma-voom:",
    round(cor(as.vector(voom_top50_z),
              as.vector(edger_top50_z), use = "complete.obs"), 3), "\n")

if (dev.cur() > 1) dev.off()
png("figures/GSE161731/deseq2_comparison/09c_heatmap_edger.png",
    width = 14, height = 12, units = "in", res = 300)
pheatmap(edger_top50_z,
         annotation_col    = col_annotation,
         annotation_colors = annotation_colors,
         show_colnames     = FALSE,
         cluster_rows      = TRUE,
         cluster_cols      = TRUE,
         clustering_method = "ward.D2",
         color             = heatmap_colors,
         breaks            = heatmap_breaks,
         fontsize_row      = 8,
         border_color      = NA,
         main = "Top 50 DE Genes: COVID-19 vs Healthy\nedgeR TMM log2-CPM normalization (Z-scored, ward.D2 clustering)")
dev.off()
cat("edgeR heatmap saved\n")

# ==============================================================================
# THREE-WAY COMPARISON: DESeq2 vs limma-voom vs edgeR
# ==============================================================================
#
# Gene-level concordance: proportion of significant genes shared between methods
# Fold change correlation: Pearson r across all genes with results in all three
#
# edgeR output columns:
#   logFC    — log2 fold change (not shrunk — same as limma-voom)
#   logCPM   — average log2 CPM
#   F        — quasi-likelihood F statistic
#   PValue   — raw p-value
#   FDR      — Benjamini-Hochberg adjusted p-value

cat("\n========== Three-Way Method Comparison ==========\n")

# Comparisons to run
comp_names <- c("COVID19_vs_healthy", "Influenza_vs_healthy",
                "Bacterial_vs_healthy", "CoV_other_vs_healthy")

# Store concordance results
concordance_3way <- data.frame()

for (comp in comp_names) {
  
  # Load DESeq2 results
  deseq2_res <- read.csv(paste0("results/DE_deseq2_", comp, ".csv"))
  deseq2_sig <- deseq2_res$ENSEMBL[deseq2_res$padj < 0.05 &
                                     !is.na(deseq2_res$padj)]
  
  # Load limma results
  limma_res <- read.csv(paste0("results/DE_limma_", comp, ".csv"))
  limma_sig <- limma_res$ENSEMBL[limma_res$adj.P.Val < 0.05 &
                                   !is.na(limma_res$adj.P.Val)]
  
  # edgeR results
  edger_res <- edger_results[[comp]]
  edger_sig <- edger_res$ENSEMBL[edger_res$FDR < 0.05 &
                                   !is.na(edger_res$FDR)]
  
  # Gene concordance — pairwise and all three
  d_vs_l <- length(intersect(deseq2_sig, limma_sig))
  d_vs_e <- length(intersect(deseq2_sig, edger_sig))
  l_vs_e <- length(intersect(limma_sig, edger_sig))
  all_3  <- length(Reduce(intersect, list(deseq2_sig, limma_sig, edger_sig)))
  
  cat("\n", comp, "\n")
  cat("  DESeq2:", length(deseq2_sig), "| limma:", length(limma_sig),
      "| edgeR:", length(edger_sig), "\n")
  cat("  DESeq2 ∩ limma:", d_vs_l,
      "(", round(100*d_vs_l/max(length(deseq2_sig),1)), "%)\n")
  cat("  DESeq2 ∩ edgeR:", d_vs_e,
      "(", round(100*d_vs_e/max(length(deseq2_sig),1)), "%)\n")
  cat("  limma ∩ edgeR:", l_vs_e,
      "(", round(100*l_vs_e/max(length(limma_sig),1)), "%)\n")
  cat("  All three:", all_3, "\n")
  
  # Fold change correlation — all genes present in all three methods
  common_genes <- Reduce(intersect, list(deseq2_res$ENSEMBL,
                                         limma_res$ENSEMBL,
                                         edger_res$ENSEMBL))
  
  d_fc <- deseq2_res$log2FoldChange[match(common_genes, deseq2_res$ENSEMBL)]
  l_fc <- limma_res$logFC[match(common_genes, limma_res$ENSEMBL)]
  e_fc <- edger_res$logFC[match(common_genes, edger_res$ENSEMBL)]
  
  r_dl <- round(cor(d_fc, l_fc, use = "complete.obs"), 3)
  r_de <- round(cor(d_fc, e_fc, use = "complete.obs"), 3)
  r_le <- round(cor(l_fc, e_fc, use = "complete.obs"), 3)
  
  cat("  FC correlation DESeq2 vs limma:", r_dl, "\n")
  cat("  FC correlation DESeq2 vs edgeR:", r_de, "\n")
  cat("  FC correlation limma vs edgeR:", r_le, "\n")
  
  concordance_3way <- bind_rows(concordance_3way, data.frame(
    comparison   = comp,
    n_deseq2     = length(deseq2_sig),
    n_limma      = length(limma_sig),
    n_edger      = length(edger_sig),
    n_all3       = all_3,
    pct_d_vs_l   = round(100*d_vs_l/max(length(deseq2_sig),1)),
    pct_d_vs_e   = round(100*d_vs_e/max(length(deseq2_sig),1)),
    pct_l_vs_e   = round(100*l_vs_e/max(length(limma_sig),1)),
    r_deseq2_limma = r_dl,
    r_deseq2_edger = r_de,
    r_limma_edger  = r_le
  ))
}

write.csv(concordance_3way,
          "results/method_concordance_3way.csv",
          row.names = FALSE)
cat("\nThree-way concordance saved to results/method_concordance_3way.csv\n")

# ==============================================================================
# THREE-WAY FC SCATTER PLOT
# ==============================================================================

# Use COVID19 vs Healthy for visualization
comp <- "COVID19_vs_healthy"

deseq2_res <- read.csv(paste0("results/DE_deseq2_", comp, ".csv"))
limma_res  <- read.csv(paste0("results/DE_limma_", comp, ".csv"))
edger_res  <- edger_results[[comp]]

common_genes <- Reduce(intersect, list(deseq2_res$ENSEMBL,
                                       limma_res$ENSEMBL,
                                       edger_res$ENSEMBL))

fc_df <- data.frame(
  ENSEMBL  = common_genes,
  deseq2   = deseq2_res$log2FoldChange[match(common_genes, deseq2_res$ENSEMBL)],
  limma    = limma_res$logFC[match(common_genes, limma_res$ENSEMBL)],
  edger    = edger_res$logFC[match(common_genes, edger_res$ENSEMBL)]
)

# Add significance category to fc_df
fc_df <- fc_df %>%
  left_join(deseq2_res %>% dplyr::select(ENSEMBL, padj), by = "ENSEMBL") %>%
  left_join(limma_res  %>% dplyr::select(ENSEMBL, adj.P.Val), by = "ENSEMBL") %>%
  left_join(edger_res  %>% dplyr::select(ENSEMBL, FDR), by = "ENSEMBL") %>%
  mutate(
    cat_dl = case_when(
      padj < 0.05 & adj.P.Val < 0.05 ~ "Both significant",
      padj < 0.05                     ~ "DESeq2 only",
      adj.P.Val < 0.05                ~ "limma only",
      TRUE                            ~ "Neither"),
    cat_de = case_when(
      padj < 0.05 & FDR < 0.05 ~ "Both significant",
      padj < 0.05               ~ "DESeq2 only",
      FDR < 0.05                ~ "edgeR only",
      TRUE                      ~ "Neither"),
    cat_le = case_when(
      adj.P.Val < 0.05 & FDR < 0.05 ~ "Both significant",
      adj.P.Val < 0.05               ~ "limma only",
      FDR < 0.05                     ~ "edgeR only",
      TRUE                           ~ "Neither")
  )

col_both <- c("Both significant" = "#8E44AD",
              "DESeq2 only"      = "#E67E22",
              "limma only"       = "#5DADE2",
              "edgeR only"       = "#2ecc71",
              "Neither"          = "#BFBFBF")
r_dl <- round(cor(fc_df$deseq2, fc_df$limma,  use = "complete.obs"), 2)
r_de <- round(cor(fc_df$deseq2, fc_df$edger,  use = "complete.obs"), 2)
r_le <- round(cor(fc_df$limma,  fc_df$edger,  use = "complete.obs"), 2)
p_dl <- ggplot(fc_df, aes(x = deseq2, y = limma, color = cat_dl)) +
  geom_point(size = 0.8, alpha = 0.5) +
  scale_color_manual(values = col_both, name = "") +
  geom_abline(slope=1, intercept=0, linetype="dashed",
              color="black", linewidth=0.5) +
  annotate("text", x=-Inf, y=Inf,
           label=paste0("r = ", r_dl),
           hjust=-0.3, vjust=1.5, size=4, fontface="bold") +
  labs(title="DESeq2 vs limma-voom",
       x="DESeq2 log2FC (apeglm-shrunk)",
       y="limma-voom logFC") +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold", size=11),
        legend.position="bottom",
        legend.text=element_text(size=8))

p_de <- ggplot(fc_df, aes(x = deseq2, y = edger, color = cat_de)) +
  geom_point(size = 0.8, alpha = 0.5) +
  scale_color_manual(values = col_both, name = "") +
  geom_abline(slope=1, intercept=0, linetype="dashed",
              color="black", linewidth=0.5) +
  annotate("text", x=-Inf, y=Inf,
           label=paste0("r = ", r_de),
           hjust=-0.3, vjust=1.5, size=4, fontface="bold") +
  labs(title="DESeq2 vs edgeR",
       x="DESeq2 log2FC (apeglm-shrunk)",
       y="edgeR logFC") +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold", size=11),
        legend.position="bottom",
        legend.text=element_text(size=8))

p_le <- ggplot(fc_df, aes(x = limma, y = edger, color = cat_le)) +
  geom_point(size = 0.8, alpha = 0.5) +
  scale_color_manual(values = col_both, name = "") +
  geom_abline(slope=1, intercept=0, linetype="dashed",
              color="black", linewidth=0.5) +
  annotate("text", x=-Inf, y=Inf,
           label=paste0("r = ", r_le),
           hjust=-0.3, vjust=1.5, size=4, fontface="bold") +
  labs(title="limma-voom vs edgeR",
       x="limma-voom logFC",
       y="edgeR logFC") +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold", size=11),
        legend.position="bottom",
        legend.text=element_text(size=8))

p_3way <- p_dl + p_de + p_le +
  plot_annotation(
    title    = "Three-Way Fold Change Correlation — COVID-19 vs Healthy",
    subtitle = paste0("DESeq2 log2FC is apeglm-shrunk | limma and edgeR logFC are unshrunk | ",
                      nrow(fc_df), " genes"),
    theme = theme(plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 10))
  )

ggsave("figures/GSE161731/deseq2_comparison/11_threeway_fc_scatter.png",
       p_3way, width = 14, height = 5, dpi = 300)
cat("Three-way FC scatter plot saved\n")

# ==============================================================================
# THREE-WAY UPSET PLOT — one file per comparison
# ==============================================================================
install.packages("gridGraphics")
library(gridGraphics)
for (comp in comp_names) {
  deseq2_res <- read.csv(paste0("results/DE_deseq2_", comp, ".csv"))
  limma_res  <- read.csv(paste0("results/DE_limma_", comp, ".csv"))
  edger_res  <- edger_results[[comp]]
  
  deseq2_sig <- deseq2_res$ENSEMBL[deseq2_res$padj < 0.05 &
                                     !is.na(deseq2_res$padj)]
  limma_sig  <- limma_res$ENSEMBL[limma_res$adj.P.Val < 0.05 &
                                    !is.na(limma_res$adj.P.Val)]
  edger_sig  <- edger_res$ENSEMBL[edger_res$FDR < 0.05 &
                                    !is.na(edger_res$FDR)]
  
  all_genes <- unique(c(deseq2_sig, limma_sig, edger_sig))
  
  upset_mat <- data.frame(
    DESeq2 = as.integer(all_genes %in% deseq2_sig),
    limma  = as.integer(all_genes %in% limma_sig),
    edgeR  = as.integer(all_genes %in% edger_sig),
    row.names = all_genes
  )
  
  comp_label <- gsub("_vs_healthy", " vs Healthy", gsub("_", " ", comp))
  filename <- paste0("figures/GSE161731/deseq2_comparison/12_upset_threeway_", comp, ".svg")
  library(gridGraphics)
  
  svg(filename, width = 14, height = 8)
  grid.echo(function() {
    upset(upset_mat,
          sets           = c("DESeq2", "limma", "edgeR"),
          order.by       = "freq",
          main.bar.color = "#2c3e50",
          sets.bar.color = c("#e74c3c", "#3498db", "#2ecc71"),
          text.scale     = 1.2,
          mainbar.y.label = "Genes in Intersection",
          sets.x.label   = "DE Genes per Method",
          mb.ratio       = c(0.6, 0.4))
  })
  dev.off()
  cat("UpSet saved:", filename, "\n")
}

# ==============================================================================
# 17. SESSION INFO
# ==============================================================================

cat("\n=== Session Info ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("R version:", R.version.string, "\n")
cat("DESeq2:", as.character(packageVersion("DESeq2")), "\n")
cat("apeglm:", as.character(packageVersion("apeglm")), "\n")
cat("limma:", as.character(packageVersion("limma")), "\n")
cat("edgeR:", as.character(packageVersion("edgeR")), "\n")
cat("ReactomePA:", as.character(packageVersion("ReactomePA")), "\n")
cat("clusterProfiler:", as.character(packageVersion("clusterProfiler")), "\n")
sessionInfo()
