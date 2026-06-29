# ==============================================================================
# Infectious Disease Transcriptomics: Host Blood Transcriptome Analysis
# Dataset: GSE161731 — Blood RNA-seq in Respiratory Infections
# Cohorts: COVID-19, Influenza, CoV_other, Bacterial, Healthy
#
# Source: McClain, M.T., Constantine, F.J., Henao, R. et al.
#         Dysregulated transcriptional responses to SARS-CoV-2 in the periphery.
#         Nature Communications 12, 1079 (2021).
#         https://doi.org/10.1038/s41467-021-21289-y
#
# Date: June 2026
#
# This script performs an independent secondary analysis of publicly available
# data from GEO. Interpretations and analytical choices reflect the author's
# own analysis and do not represent the views of the original investigators.
#
# METHOD: limma-voom with TMM normalization
#   Consistent with the original McClain et al. paper which used limma GLMs.
#   Key difference from original: one sample per subject (earliest timepoint)
#   rather than duplicateCorrelation() to keep the analysis simple and
#   ensure all observations are independent.
#
# DESIGN FORMULA: ~ batch + gender + age + cohort
#   Covariates: batch, gender, age
#   Variable of interest: cohort (infection type vs healthy)
#
# BIOLOGICAL QUESTION:
#   What transcriptomic signatures are shared vs unique across four respiratory
#   infections? We identify pan-infection, viral-specific, and pathogen-specific
#   gene sets and characterize their pathway enrichment.
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

library(limma)
library(edgeR)
library(tidyverse)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(ReactomePA)
library(UpSetR)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(VennDiagram)
library(grid)

setwd("~/Desktop/infectious-disease-transcriptomics")
dir.create("figures/GSE161731", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE)

palette_cohort <- c(
  healthy   = "#2ecc71",
  COVID19   = "#e74c3c",
  Influenza = "#3498db",
  Bacterial = "#e67e22",
  CoV_other = "#9b59b6"
)


# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

counts <- read.csv("results/GSE161731_counts_filtered.csv", row.names = 1)
meta   <- read.csv("results/GSE161731_metadata_clean.csv")

cat("Loaded:", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("Metadata:", nrow(meta), "samples\n")


# ==============================================================================
# 3. MATCH SAMPLE IDs
# ==============================================================================

meta$rna_id_fixed <- make.names(meta$rna_id)
rownames(meta)    <- meta$rna_id_fixed

meta   <- meta[colnames(counts), ]
stopifnot(all(colnames(counts) == rownames(meta)))
cat("Sample order verified\n")


# ==============================================================================
# 4. QC AND SAMPLE FILTERING
# ==============================================================================
#
# STEP 1: Remove cross-cohort subjects (appeared as both healthy and CoV_other)
# STEP 2: Keep earliest timepoint per subject (14 COVID-19 patients had
#         multiple timepoints). Original authors used duplicateCorrelation();
#         we use one sample per subject to ensure independence.
# STEP 3: Remove samples with missing gender (required as covariate)
# STEP 4: Remove samples with missing age (required as covariate)

# Step 1
cross_cohort <- c("A34BBC", "E77BE3")
meta   <- meta %>% filter(!subject_id %in% cross_cohort)
counts <- counts[, rownames(meta)]
cat("After cross-cohort removal:", nrow(meta), "samples\n")

# Step 2
time_rank      <- c("early" = 1, "middle" = 2, "late" = 3)
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

# Step 3
meta   <- meta %>% filter(gender != "" & !is.na(gender))
counts <- counts[, rownames(meta)]
cat("After gender filter:", nrow(meta), "samples\n")

# Step 4
meta   <- meta %>% filter(!is.na(age))
counts <- counts[, rownames(meta)]
cat("After age filter:", nrow(meta), "samples\n")


# ==============================================================================
# 5. PREPARE FACTORS
# ==============================================================================

meta$cohort <- recode(meta$cohort,
                      "COVID-19"  = "COVID19",
                      "CoV other" = "CoV_other")
meta$cohort <- factor(meta$cohort,
                      levels = c("healthy","COVID19","Influenza",
                                 "Bacterial","CoV_other"))
meta$batch  <- factor(meta$batch)
meta$gender <- factor(meta$gender)
meta$age    <- as.numeric(meta$age)

cat("\n=== Final Sample Summary ===\n")
print(table(meta$cohort))
cat("Independent samples:", nrow(meta) == length(unique(meta$subject_id)), "\n")


# ==============================================================================
# 6. LIMMA-VOOM DIFFERENTIAL EXPRESSION
# ==============================================================================

# TMM normalization
dge <- DGEList(counts = counts, group = meta$cohort)
dge <- calcNormFactors(dge, method = "TMM")

# Design matrix
design <- model.matrix(~ batch + gender + age + cohort, data = meta)

# Voom transformation — estimates mean-variance relationship and assigns
# precision weights to each observation, stabilizing variance across
# expression levels before linear modeling
png("figures/GSE161731/00_voom_mean_variance.png",
    width = 8, height = 6, units = "in", res = 300)
v <- voom(dge, design, plot = TRUE)
title("")
dev.off()
v <- voom(dge, design, plot = FALSE)

# Fit linear model and apply empirical Bayes moderation
# eBayes borrows information across genes to stabilize variance estimates —
# particularly important for small sample sizes
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Extract results for each infection vs healthy
coef_map <- c(
  COVID19   = "cohortCOVID19",
  Influenza = "cohortInfluenza",
  Bacterial = "cohortBacterial",
  CoV_other = "cohortCoV_other"
)

limma_results <- list()

cat("\n=== limma-voom DE Results ===\n")
for (name in names(coef_map)) {
  res <- topTable(fit, coef = coef_map[[name]], number = Inf, sort.by = "none")
  res$ENSEMBL <- rownames(res)
  n_sig  <- sum(res$adj.P.Val < 0.05)
  n_up   <- sum(res$adj.P.Val < 0.05 & res$logFC > 0)
  n_down <- sum(res$adj.P.Val < 0.05 & res$logFC < 0)
  cat(name, ":", n_sig, "(", n_up, "up,", n_down, "down)\n")
  limma_results[[name]] <- res
  write.csv(res,
            paste0("results/DE_limma_", name, "_vs_healthy.csv"),
            row.names = FALSE)
}


# ==============================================================================
# 7. GENE ID MAPPING
# ==============================================================================

all_ensembl <- rownames(v$E)
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = all_ensembl,
  keytype = "ENSEMBL",
  columns = "SYMBOL"
) %>% distinct(ENSEMBL, .keep_all = TRUE)

cat("Mapped", sum(!is.na(gene_map$SYMBOL)), "of",
    nrow(gene_map), "Ensembl IDs to gene symbols\n")

# Add symbols to results
for (name in names(limma_results)) {
  limma_results[[name]] <- limma_results[[name]] %>%
    left_join(gene_map, by = "ENSEMBL")
}


# ==============================================================================
# 8. VOLCANO PLOTS
# ==============================================================================

make_volcano <- function(res_df, comparison_name, filename) {
  res_df <- res_df %>% filter(!is.na(adj.P.Val) & !is.na(logFC))
  res_df$color <- "NS"
  res_df$color[res_df$adj.P.Val < 0.05 & res_df$logFC > 0] <- "Up"
  res_df$color[res_df$adj.P.Val < 0.05 & res_df$logFC < 0] <- "Down"
  
  top_up <- res_df %>% filter(color == "Up") %>%
    arrange(adj.P.Val) %>% head(5)
  top_dn <- res_df %>% filter(color == "Down") %>%
    arrange(adj.P.Val) %>% head(5)
  top_fc_up <- res_df %>% filter(color == "Up") %>%
    arrange(desc(logFC)) %>% head(5)
  top_fc_dn <- res_df %>% filter(color == "Down") %>%
    arrange(logFC) %>% head(5)
  
  lab_genes <- bind_rows(top_up, top_dn, top_fc_up, top_fc_dn) %>%
    distinct(ENSEMBL, .keep_all = TRUE)
  res_df$label <- ifelse(res_df$ENSEMBL %in% lab_genes$ENSEMBL,
                         res_df$SYMBOL, "")
  
  n_up   <- sum(res_df$color == "Up")
  n_down <- sum(res_df$color == "Down")
  
  p <- ggplot(res_df, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
    geom_point(size = 1.2, alpha = 0.7) +
    scale_color_manual(
      values = c("Up" = "#E67E22", "Down" = "#5DADE2", "NS" = "#BFBFBF"),
      labels = c("Up"   = paste0("Up in ", comparison_name, " (", n_up, ")"),
                 "Down" = paste0("Up in Healthy (", n_down, ")"),
                 "NS"   = "Not Significant"),
      breaks = c("Up", "Down", "NS"),
      name   = ""
    ) +
    geom_text_repel(data = res_df %>% filter(label != ""),
                    aes(label = label), size = 3.5, max.overlaps = 20,
                    segment.color = "grey40", fontface = "bold.italic",
                    box.padding = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "black", linewidth = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed",
               color = "black", linewidth = 0.5) +
    labs(
      title    = paste0(comparison_name, " vs Healthy Controls"),
      subtitle = "GSE161731 | limma-voom | FDR < 0.05 | Dashed lines at |log2FC| = 1",
      x        = expression(bold(log[2]~"Fold Change")),
      y        = expression(bold(-log[10]~"(FDR)"))
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position  = "top",
          legend.text      = element_text(size = 11, face = "bold"),
          plot.title       = element_text(face = "bold", size = 14),
          panel.grid.minor = element_blank())
  ggsave(filename, p, width = 10, height = 8, dpi = 300)
  cat("Saved:", filename, "\n")
}

make_volcano(limma_results$COVID19,
             "COVID-19",
             "figures/GSE161731/05a_volcano_covid_vs_healthy.png")
make_volcano(limma_results$Influenza,
             "Influenza",
             "figures/GSE161731/05b_volcano_influenza_vs_healthy.png")
make_volcano(limma_results$Bacterial,
             "Bacterial",
             "figures/GSE161731/05c_volcano_bacterial_vs_healthy.png")
make_volcano(limma_results$CoV_other,
             "CoV Other",
             "figures/GSE161731/05d_volcano_cov_other_vs_healthy.png")


# ==============================================================================
# 9. DE GENE COUNTS
# ==============================================================================

count_df <- map_dfr(names(limma_results), function(name) {
  res <- limma_results[[name]]
  data.frame(
    comparison  = name,
    direction   = c("Upregulated", "Downregulated"),
    count       = c(
      sum(res$adj.P.Val < 0.05 & res$logFC > 0, na.rm = TRUE),
      -sum(res$adj.P.Val < 0.05 & res$logFC < 0, na.rm = TRUE)
    )
  )
})

p_counts <- ggplot(count_df, aes(x = comparison, y = count, fill = direction)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("Upregulated" = "#e74c3c",
                               "Downregulated" = "#3498db"),
                    name = "") +
  labs(title    = "DE Genes vs Healthy Controls",
       subtitle  = "GSE161731 | limma-voom | FDR < 0.05",
       x = "", y = "Number of DE Genes") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        plot.title   = element_text(face = "bold"))

ggsave("figures/GSE161731/06_de_gene_counts.png",
       p_counts, width = 8, height = 6, dpi = 300)
cat("Plot 06 saved: DE gene counts\n")


# ==============================================================================
# 10. INFECTION-SPECIFIC GENE SIGNATURES
# ==============================================================================
#
# Extract significant genes per infection and identify:
#   (a) Infection-specific genes — unique to one infection
#   (b) Shared viral genes — shared among COVID19, Influenza, CoV_other
#   (c) Pan-infection genes — shared across all four infections
#
# These gene sets are used for GSEA to characterize the biological processes
# driving each type of response.

sig_genes <- lapply(limma_results, function(res) {
  res %>%
    filter(adj.P.Val < 0.05 & !is.na(SYMBOL)) %>%
    pull(ENSEMBL) %>%
    unique()
})

cat("\n=== Significant Genes per Infection ===\n")
for (name in names(sig_genes)) {
  cat(" ", name, ":", length(sig_genes[[name]]), "genes\n")
}

# Define biologically meaningful gene sets
covid_only     <- setdiff(sig_genes$COVID19,
                          Reduce(union, sig_genes[c("Influenza","Bacterial","CoV_other")]))
influenza_only <- setdiff(sig_genes$Influenza,
                          Reduce(union, sig_genes[c("COVID19","Bacterial","CoV_other")]))
bacterial_only <- setdiff(sig_genes$Bacterial,
                          Reduce(union, sig_genes[c("COVID19","Influenza","CoV_other")]))
cov_other_only <- setdiff(sig_genes$CoV_other,
                          Reduce(union, sig_genes[c("COVID19","Influenza","Bacterial")]))

# Shared viral — in all three viral groups but not bacterial
shared_viral   <- Reduce(intersect,
                         list(sig_genes$COVID19, sig_genes$Influenza,
                              sig_genes$CoV_other)) %>%
  setdiff(sig_genes$Bacterial)

# Pan-infection — shared across all four
shared_all     <- Reduce(intersect,
                         list(sig_genes$COVID19, sig_genes$Influenza,
                              sig_genes$Bacterial, sig_genes$CoV_other))

# COVID19 vs CoV_other — what distinguishes SARS-CoV-2 from seasonal coronaviruses
covid_vs_cov   <- intersect(sig_genes$COVID19, sig_genes$CoV_other) %>%
  setdiff(union(sig_genes$Influenza, sig_genes$Bacterial))

cat("\n=== Infection-Specific Gene Sets ===\n")
cat("  COVID-19 specific:        ", length(covid_only), "\n")
cat("  Influenza specific:       ", length(influenza_only), "\n")
cat("  Bacterial specific:       ", length(bacterial_only), "\n")
cat("  CoV_other specific:       ", length(cov_other_only), "\n")
cat("  Shared viral (no bact):   ", length(shared_viral), "\n")
cat("  Pan-infection (all four): ", length(shared_all), "\n")
cat("  COVID19 + CoV_other only: ", length(covid_vs_cov), "\n")


# ==============================================================================
# 11. UPSET PLOT — INFECTION OVERLAP
# ==============================================================================

all_sig <- unique(unlist(sig_genes))
upset_mat <- data.frame(
  COVID19   = as.integer(all_sig %in% sig_genes$COVID19),
  Influenza = as.integer(all_sig %in% sig_genes$Influenza),
  Bacterial = as.integer(all_sig %in% sig_genes$Bacterial),
  CoV_other = as.integer(all_sig %in% sig_genes$CoV_other),
  row.names = all_sig
)

svg("figures/GSE161731/07_upset_infection_overlap.svg", width = 14, height = 8)
p_upset <- upset(upset_mat,
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
print(p_upset)
dev.off()
cat("Plot 07 saved: UpSet infection overlap\n")


# ==============================================================================
# 12. VENN DIAGRAM — FOUR-WAY
# ==============================================================================
#
# Four-way Venn showing overlap between COVID19, Influenza, Bacterial, CoV_other
# significant genes. Highlights infection-specific and shared signatures.

venn_colors <- c("#e74c3c80","#3498db80","#e67e2280","#9b59b680")

png("figures/GSE161731/08_venn_four_way.png",
    width = 10, height = 10, units = "in", res = 300)
grid.newpage()
venn_plot <- draw.quad.venn(
  area1         = length(sig_genes$COVID19),
  area2         = length(sig_genes$Influenza),
  area3         = length(sig_genes$Bacterial),
  area4         = length(sig_genes$CoV_other),
  n12           = length(intersect(sig_genes$COVID19, sig_genes$Influenza)),
  n13           = length(intersect(sig_genes$COVID19, sig_genes$Bacterial)),
  n14           = length(intersect(sig_genes$COVID19, sig_genes$CoV_other)),
  n23           = length(intersect(sig_genes$Influenza, sig_genes$Bacterial)),
  n24           = length(intersect(sig_genes$Influenza, sig_genes$CoV_other)),
  n34           = length(intersect(sig_genes$Bacterial, sig_genes$CoV_other)),
  n123          = length(Reduce(intersect, sig_genes[c("COVID19","Influenza","Bacterial")])),
  n124          = length(Reduce(intersect, sig_genes[c("COVID19","Influenza","CoV_other")])),
  n134          = length(Reduce(intersect, sig_genes[c("COVID19","Bacterial","CoV_other")])),
  n234          = length(Reduce(intersect, sig_genes[c("Influenza","Bacterial","CoV_other")])),
  n1234         = length(shared_all),
  category      = c("COVID-19","Influenza","Bacterial","CoV Other"),
  fill          = venn_colors,
  col           = c("#e74c3c","#3498db","#e67e22","#9b59b6"),
  fontfamily    = "sans",
  cat.fontface  = "bold",
  cat.fontfamily = "sans",
  cex           = 1.2,
  cat.cex       = 1.3
)
dev.off()
cat("Plot 08 saved: four-way Venn diagram\n")


# ==============================================================================
# 13. GSEA REACTOME — INFECTION-SPECIFIC SIGNATURES
# ==============================================================================
#
# Run GSEA on ranked gene lists for each infection to characterize:
#   - What pathways are uniquely activated in each infection
#   - What pathways are shared across viral infections
#   - What pathways define the pan-infection host response
#
# RANKING METRIC: sign(logFC) * -log10(P.Value)
#   Uses raw P.Value (not adjusted) for ranking — GSEA handles multiple testing
#   internally through permutation. P.Value preserves more signal than padj.

make_ranked_list <- function(res_df) {
  res_df %>%
    filter(!is.na(logFC) & !is.na(P.Value) & !is.na(SYMBOL)) %>%
    mutate(rank_metric = sign(logFC) * -log10(P.Value)) %>%
    arrange(desc(rank_metric)) %>%
    { left_join(., bitr(.$SYMBOL, fromType = "SYMBOL",
                        toType = "ENTREZID", OrgDb = org.Hs.eg.db),
                by = "SYMBOL",
                relationship = "many-to-many") } %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    { setNames(.$rank_metric, .$ENTREZID) }
}

run_gsea <- function(gene_list, label) {
  set.seed(42)
  tryCatch(
    gsePathway(geneList     = gene_list,
               organism     = "human",
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               verbose      = FALSE,
               eps          = 0),
    error = function(e) {
      cat("  GSEA error for", label, ":", e$message, "\n")
      return(NULL)
    }
  )
}

make_gsea_plot <- function(gsea_res, title, filename, n_top = 10) {
  if (is.null(gsea_res) || nrow(gsea_res@result) == 0) {
    cat("No significant pathways for:", title, "\n")
    return(NULL)
  }
  
  plot_df <- gsea_res@result %>%
    mutate(direction = ifelse(NES > 0, "Activated", "Suppressed")) %>%
    group_by(direction) %>%
    arrange(p.adjust) %>%
    slice_head(n = n_top) %>%
    ungroup() %>%
    mutate(Description_short = str_trunc(Description, 50))
  
  p <- ggplot(plot_df,
              aes(x = NES, y = reorder(Description_short, NES),
                  size = setSize, color = p.adjust)) +
    geom_point(alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    scale_color_gradient(low = "#C0392B", high = "#2980B9", name = "FDR") +
    scale_size_continuous(range = c(2, 8), name = "Gene set size") +
    labs(title    = title,
         subtitle = paste0("GSE161731 | limma-voom | Reactome GSEA | FDR < 0.05 | ",
                           "Top ", n_top, " activated + suppressed"),
         x = "Normalized Enrichment Score (NES)",
         y = "") +
    theme_minimal(base_size = 11) +
    theme(plot.title       = element_text(face = "bold", size = 12),
          axis.text.y      = element_text(size = 9),
          panel.grid.minor = element_blank(),
          legend.position  = "right")
  
  ggsave(filename, p, width = 12, height = 7, dpi = 300)
  cat("Saved:", filename, "\n")
  return(gsea_res)
}

# Run GSEA for each infection vs healthy
gsea_results <- list()
for (name in names(limma_results)) {
  cat("\n=== GSEA:", name, "===\n")
  gl <- make_ranked_list(limma_results[[name]])
  cat("Ranked genes:", length(gl), "\n")
  gsea_res <- run_gsea(gl, name)
  gsea_results[[name]] <- gsea_res
  if (!is.null(gsea_res)) {
    cat("Significant pathways:", nrow(gsea_res), "\n")
    write.csv(gsea_res@result,
              paste0("results/GSEA_", name, "_vs_healthy.csv"),
              row.names = FALSE)
    make_gsea_plot(
      gsea_res,
      title    = paste0(name, " vs Healthy — Reactome Pathways"),
      filename = paste0("figures/GSE161731/09a_gsea_", tolower(name), ".png")
    )
  }
}

# ==============================================================================
# 14. GSEA — INFECTION-SPECIFIC GENE SETS
# ==============================================================================
#
# Instead of running GSEA on the full ranked list, we now characterize the
# biological processes specifically enriched in:
#   (a) Infection-specific genes (unique to each infection)
#   (b) Shared viral signature (COVID19 + Influenza + CoV_other, not Bacterial)
#   (c) Pan-infection signature (all four infections)
#   (d) Coronavirus-specific (COVID19 + CoV_other only)
#
# For gene-set-level ORA (overrepresentation analysis) rather than GSEA,
# we use enrichPathway() on each gene set.

run_ora <- function(ensembl_ids, label, universe_ids) {
  entrez_map <- bitr(ensembl_ids, fromType = "ENSEMBL",
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  universe_entrez <- bitr(universe_ids, fromType = "ENSEMBL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (nrow(entrez_map) == 0) {
    cat("  No Entrez IDs for", label, "\n")
    return(NULL)
  }
  tryCatch(
    enrichPathway(gene         = entrez_map$ENTREZID,
                  universe     = universe_entrez$ENTREZID,
                  organism     = "human",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  readable     = TRUE),
    error = function(e) {
      cat("  ORA error for", label, ":", e$message, "\n")
      return(NULL)
    }
  )
}

make_ora_plot <- function(ora_res, title, filename, n_top = 15) {
  if (is.null(ora_res) || nrow(ora_res@result) == 0) {
    cat("No significant pathways for:", title, "\n")
    return(NULL)
  }
  
  plot_df <- ora_res@result %>%
    arrange(p.adjust) %>%
    head(n_top) %>%
    mutate(Description_short = str_trunc(Description, 55),
           GeneRatio_num = sapply(GeneRatio, function(x) {
             parts <- strsplit(x, "/")[[1]]
             as.numeric(parts[1]) / as.numeric(parts[2])
           }))
  
  p <- ggplot(plot_df,
              aes(x = GeneRatio_num,
                  y = reorder(Description_short, GeneRatio_num),
                  size = Count, color = p.adjust)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(low = "#C0392B", high = "#2980B9", name = "FDR") +
    scale_size_continuous(range = c(3, 10), name = "Gene count") +
    labs(title    = title,
         subtitle = paste0("GSE161731 | Reactome ORA | FDR < 0.05 | Top ", n_top),
         x = "Gene Ratio",
         y = "") +
    theme_minimal(base_size = 11) +
    theme(plot.title       = element_text(face = "bold", size = 12),
          axis.text.y      = element_text(size = 9),
          panel.grid.minor = element_blank(),
          legend.position  = "right")
  
  ggsave(filename, p, width = 12, height = 7, dpi = 300)
  cat("Saved:", filename, "\n")
}

universe <- unique(unlist(sig_genes))

gene_sets <- list(
  covid_specific    = list(genes = covid_only,
                           label = "COVID-19 Specific Genes",
                           file  = "figures/GSE161731/09b_ora_covid_specific.png"),
  influenza_specific = list(genes = influenza_only,
                            label = "Influenza Specific Genes",
                            file  = "figures/GSE161731/09c_ora_influenza_specific.png"),
  bacterial_specific = list(genes = bacterial_only,
                            label = "Bacterial Specific Genes",
                            file  = "figures/GSE161731/09d_ora_bacterial_specific.png"),
  cov_other_specific = list(genes = cov_other_only,
                            label = "CoV Other Specific Genes",
                            file  = "figures/GSE161731/09e_ora_cov_other_specific.png"),
  shared_viral       = list(genes = shared_viral,
                            label = "Shared Viral Signature (COVID19 + Influenza + CoV_other)",
                            file  = "figures/GSE161731/09f_ora_shared_viral.png"),
  shared_all         = list(genes = shared_all,
                            label = "Pan-Infection Signature (All Four Infections)",
                            file  = "figures/GSE161731/09g_ora_shared_all.png"),
  coronavirus        = list(genes = covid_vs_cov,
                            label = "Coronavirus Signature (COVID-19 + CoV_other only)",
                            file  = "figures/GSE161731/09h_ora_coronavirus.png")
)

ora_results <- list()
for (set_name in names(gene_sets)) {
  gs <- gene_sets[[set_name]]
  cat("\n=== ORA:", gs$label, "| n =", length(gs$genes), "===\n")
  if (length(gs$genes) < 5) {
    cat("  Too few genes, skipping\n")
    next
  }
  ora_res <- run_ora(gs$genes, gs$label, universe)
  ora_results[[set_name]] <- ora_res
  if (!is.null(ora_res) && nrow(ora_res@result) > 0) {
    cat("  Significant pathways:", nrow(ora_res@result), "\n")
    write.csv(ora_res@result,
              paste0("results/ORA_", set_name, ".csv"),
              row.names = FALSE)
    make_ora_plot(ora_res, gs$label, gs$file)
  }
}


# ==============================================================================
# 15. HEATMAP — TOP 50 GENES
# ==============================================================================
#
# Show top 50 DE genes from COVID-19 vs Healthy across all samples.
# VST equivalent for limma-voom is the log2-CPM E matrix from voom().
# Z-score per gene so colors reflect relative expression across samples.

top50 <- limma_results$COVID19 %>%
  filter(!is.na(adj.P.Val) & !is.na(SYMBOL)) %>%
  arrange(adj.P.Val) %>%
  head(50)

cat("\nTop 50 genes padj range:", range(top50$adj.P.Val), "\n")

ann_col <- data.frame(
  Cohort = meta$cohort,
  Batch  = meta$batch,
  row.names = rownames(meta)
)
ann_colors <- list(
  Cohort = palette_cohort,
  Batch  = c("1" = "#95a5a6", "2" = "#2c3e50")
)

expr_top50   <- v$E[top50$ENSEMBL, ]
rownames(expr_top50) <- top50$SYMBOL
expr_top50_z <- t(scale(t(expr_top50)))

heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
heatmap_breaks <- seq(-2.5, 2.5, length.out = 101)

cohort_order          <- order(meta$cohort)
expr_top50_ordered    <- expr_top50_z[, cohort_order]
ann_col_ordered       <- ann_col[cohort_order, , drop = FALSE]
cohort_sizes          <- table(meta$cohort[cohort_order])
gaps                  <- cumsum(cohort_sizes)[1:(length(cohort_sizes) - 1)]

if (dev.cur() > 1) dev.off()
png("figures/GSE161731/10_heatmap_top50.png",
    width = 14, height = 12, units = "in", res = 300)
pheatmap(expr_top50_ordered,
         annotation_col    = ann_col_ordered,
         annotation_colors = ann_colors,
         show_colnames     = FALSE,
         cluster_rows      = TRUE,
         cluster_cols      = FALSE,
         gaps_col          = as.integer(gaps),
         color             = heatmap_colors,
         breaks            = heatmap_breaks,
         fontsize_row      = 8,
         border_color      = NA,
         main = paste0("Top 50 DE Genes: COVID-19 vs Healthy\n",
                       "limma-voom log2-CPM (Z-scored | Samples ordered by cohort)"))
dev.off()
cat("Plot 10 saved: heatmap top 50\n")

# ==============================================================================
# 16. CORONAVIRUS-SPECIFIC mRNA SPLICING SIGNATURE
# ==============================================================================

cat("\n=== Coronavirus mRNA Splicing Signature ===\n")

ora_corona <- ora_results$coronavirus

if (!is.null(ora_corona) && nrow(ora_corona@result) > 0) {
  
  # ----------------------------------------------------------------------------
  # Step 1: Extract splicing genes from ORA result
  # ----------------------------------------------------------------------------
  splicing_pathways <- c("mRNA Splicing - Major Pathway",
                         "mRNA Splicing",
                         "Processing of Capped Intron-Containing Pre-mRNA")
  
  splicing_genes <- ora_corona@result %>%
    filter(Description %in% splicing_pathways) %>%
    pull(geneID) %>%
    strsplit("/") %>%
    unlist() %>%
    unique()
  
  cat("Genes driving coronavirus mRNA splicing signature:\n")
  print(splicing_genes)
  
  # ----------------------------------------------------------------------------
  # Step 2: Check direction in each cohort
  # ----------------------------------------------------------------------------
  cat("\nCOVID-19 expression:\n")
  limma_results$COVID19 %>%
    filter(SYMBOL %in% splicing_genes) %>%
    arrange(adj.P.Val) %>%
    dplyr::select(SYMBOL, logFC, adj.P.Val) %>%
    print()
  
  cat("\nCoV_other expression:\n")
  limma_results$CoV_other %>%
    filter(SYMBOL %in% splicing_genes) %>%
    arrange(adj.P.Val) %>%
    dplyr::select(SYMBOL, logFC, adj.P.Val) %>%
    print()
  
  cat("\nAbsence in non-coronavirus infections:\n")
  for (name in c("Influenza", "Bacterial")) {
    sig <- limma_results[[name]] %>%
      filter(SYMBOL %in% splicing_genes & adj.P.Val < 0.05)
    cat(" ", name, "— splicing genes significant:", nrow(sig), "\n")
  }
  
  # ----------------------------------------------------------------------------
  # Step 3: Build expression data frame for plotting
  # ----------------------------------------------------------------------------
  splicing_plot_df <- map_dfr(splicing_genes, function(gene) {
    gm  <- gene_map
    ens <- gm$ENSEMBL[gm$SYMBOL == gene]
    ens <- ens[!is.na(ens)][1]
    cat("Gene:", gene, "Ensembl:", ens, "\n")
    if (is.na(ens) || !ens %in% rownames(v$E)) return(NULL)
    tibble(
      gene   = gene,
      expr   = as.numeric(v$E[ens, ]),
      cohort = as.character(meta[colnames(v$E), "cohort"]),
      sample = colnames(v$E)
    )
  })
  
  splicing_plot_df$cohort <- factor(
    splicing_plot_df$cohort,
    levels = c("healthy","COVID19","Influenza","Bacterial","CoV_other")
  )
  
  cat("splicing_plot_df rows:", nrow(splicing_plot_df), "\n")
  
  # ----------------------------------------------------------------------------
  # Step 4: Build significance labels (FDR from limma vs healthy)
  # ----------------------------------------------------------------------------
  make_padj_label <- function(p) {
    if (is.na(p))   return("ns")
    if (p < 0.001)  return("***")
    if (p < 0.01)   return("**")
    if (p < 0.05)   return("*")
    return("ns")
  }
  
  padj_labels <- map_dfr(splicing_genes, function(gene) {
    map_dfr(names(limma_results), function(name) {
      res <- limma_results[[name]] %>% filter(SYMBOL == gene)
      if (nrow(res) == 0) return(NULL)
      tibble(
        gene      = gene,
        cohort    = name,
        adj.P.Val = res$adj.P.Val[1],
        label     = make_padj_label(res$adj.P.Val[1])
      )
    })
  })
  
  padj_labels$cohort <- factor(
    padj_labels$cohort,
    levels = c("healthy","COVID19","Influenza","Bacterial","CoV_other")
  )
  
  # ----------------------------------------------------------------------------
  # Step 5: Build bracket positions
  # ----------------------------------------------------------------------------
  # x positions: healthy=1, COVID19=2, Influenza=3, Bacterial=4, CoV_other=5
  x_map <- c("healthy"=1,"COVID19"=2,"Influenza"=3,"Bacterial"=4,"CoV_other"=5)
  
  y_range <- splicing_plot_df %>%
    group_by(gene) %>%
    summarize(y_max   = max(expr, na.rm = TRUE),
              y_min   = min(expr, na.rm = TRUE),
              y_range = y_max - y_min,
              .groups = "drop")
  
  bracket_labels <- padj_labels %>%
    left_join(y_range, by = "gene") %>%
    mutate(
      x1           = 1,
      x2           = x_map[as.character(cohort)],
      bracket_rank = x_map[as.character(cohort)] - 1,
      y_pos        = y_max + y_range * 0.08 * bracket_rank,
      y_tick       = y_pos - y_range * 0.02,
      y_text       = y_pos + y_range * 0.015,
      x_mid        = (1 + x_map[as.character(cohort)]) / 2
    ) %>%
    filter(cohort != "healthy")
  
  # ----------------------------------------------------------------------------
  # Step 6: Boxplot — healthy first, brackets vs healthy
  # ----------------------------------------------------------------------------
  p_splicing <- ggplot(splicing_plot_df,
                       aes(x = cohort, y = expr, fill = cohort)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.4, color = "grey30") +
    # Horizontal bracket line
    geom_segment(data = bracket_labels,
                 aes(x = x1, xend = x2,
                     y = y_pos, yend = y_pos,
                     fill = NULL),
                 color = "black", linewidth = 0.4,
                 inherit.aes = FALSE) +
    # Left tick
    geom_segment(data = bracket_labels,
                 aes(x = x1, xend = x1,
                     y = y_tick, yend = y_pos,
                     fill = NULL),
                 color = "black", linewidth = 0.4,
                 inherit.aes = FALSE) +
    # Right tick
    geom_segment(data = bracket_labels,
                 aes(x = x2, xend = x2,
                     y = y_tick, yend = y_pos,
                     fill = NULL),
                 color = "black", linewidth = 0.4,
                 inherit.aes = FALSE) +
    # Significance label
    geom_text(data = bracket_labels,
              aes(x = x_mid, y = y_text, label = label, fill = NULL),
              size = 3.5, fontface = "bold",
              inherit.aes = FALSE) +
    scale_fill_manual(values = palette_cohort, name = "") +
    facet_wrap(~ gene, scales = "free_y", nrow = 1) +
    labs(
      title    = "Coronavirus-Specific mRNA Splicing Genes",
      subtitle = paste0(
        "GSE161731 | limma-voom log2-CPM | Brackets show FDR vs healthy controls\n",
        "* FDR<0.05  ** FDR<0.01  *** FDR<0.001  ns = not significant"
      ),
      x = "",
      y = "log2-CPM expression"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(size = 9, color = "grey40"),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
      strip.text       = element_text(face = "bold.italic", size = 12),
      legend.position  = "none",
      panel.grid.minor = element_blank()
    )
  
  ggsave("figures/GSE161731/11_boxplot_coronavirus_splicing.png",
         p_splicing, width = 14, height = 7, dpi = 300)
  cat("Plot 11 saved: coronavirus splicing boxplots\n")
  
} 
# ==============================================================================
# 16. SESSION INFO
# ==============================================================================

cat("\n=== Analysis Complete ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("R version:", R.version.string, "\n")
cat("limma:", as.character(packageVersion("limma")), "\n")
cat("edgeR:", as.character(packageVersion("edgeR")), "\n")
cat("ReactomePA:", as.character(packageVersion("ReactomePA")), "\n")
cat("clusterProfiler:", as.character(packageVersion("clusterProfiler")), "\n")
sessionInfo()