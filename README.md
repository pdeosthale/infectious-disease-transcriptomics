# Infectious Disease Transcriptomics

Computational analysis of host immune responses to respiratory infections using public blood transcriptomic data. This project applies differential expression analysis to characterize pathogen-specific and shared transcriptomic signatures across viral and bacterial respiratory infections.

## Biological Question

How does the blood transcriptome differ across respiratory infections (COVID-19, Influenza, other coronaviruses, bacterial pneumonia) compared to healthy controls? Are there shared host response genes across all infections, or does each pathogen trigger a distinct transcriptomic program?

## Disclaimer

This is an independent secondary analysis of publicly available data from GEO (GSE161731) for educational and portfolio purposes. The original study was designed and conducted by McClain et al. (2021). Any interpretations, limitations identified, or analytical choices reflect the author's own learning process and do not represent the views or conclusions of the original investigators.

## Dataset

**GSE161731** — Blood RNA-seq from patients with respiratory infections.

Source: McClain, M.T., Constantine, F.J., Henao, R. et al. Dysregulated transcriptional responses to SARS-CoV-2 in the periphery. *Nature Communications* 12, 1079 (2021). https://doi.org/10.1038/s41467-021-21289-y

| Cohort | Samples (after QC) | Median Age |
|--------|-------------------|------------|
| COVID-19 | 45 | 33 |
| CoV other | 36 | 18 |
| Influenza | 17 | 18 |
| Bacterial | 23 | 60 |
| Healthy | 17 | 18 |

**Final sample count: 136 independent samples** after all QC and deduplication steps.

**QC steps applied:**

- Gene filtering: removed genes detected in fewer than 10% of samples (60,675 to 28,996 genes)
- Outlier removal: 2 samples removed via PCA (very low library size and gene detection)
- Cross-cohort subjects: 2 subjects (A34BBC, E77BE3) removed — appeared as both healthy and CoV_other
- Pseudoreplication: 14 COVID-19 patients had multiple timepoints; kept earliest timepoint only
- Missing gender: 21 samples removed (all CoV_other) to include gender as adjustment variable
- Missing age: samples with NA age removed to include age as adjustment variable

## Key Analytical Choices

### 1. One sample per subject for independent comparisons

The dataset includes 14 COVID-19 patients sampled at multiple timepoints (early/middle/late disease). The original authors accounted for this using `limma::duplicateCorrelation()`, which is a statistically valid approach that retains all timepoints. We took a more conservative approach — restricting to one sample per subject (earliest timepoint) — to enable a fair comparison between DESeq2 and limma-voom on identical data. DESeq2 has no equivalent to `duplicateCorrelation()`, so including repeated samples would have made the two-method comparison unequal. This reduced COVID-19 DE genes from 2,820 to 874, reflecting the reduced sample size and removal of within-subject correlation, not a correction of an error in the original study.

### 2. Age adjustment in a multi-cohort design

The study enrolled participants from different sources — healthy controls from a university undergraduate cohort (median age 18) and infected patients from emergency departments and clinical studies (median age 33–60 depending on cohort). This age variation across cohorts is an inherent feature of the multi-cohort design. In our reanalysis, we included age as an adjustment variable to isolate infection-specific transcriptomic effects from age-related gene expression differences.

| Comparison | Without Age Adjustment | With Age Adjustment | Difference | % Change |
|---|---|---|---|---|
| Bacterial vs Healthy | 10,564 | 7,607 | 2,957 | 28% |
| COVID-19 vs Healthy | 3,816 | 874* | — | — |

*After both one-sample-per-subject restriction and age adjustment.

## Key Findings

### 1. Bacterial infection triggers the largest transcriptomic shift

In our age-adjusted, one-sample-per-subject analysis, bacterial infection produces the most DE genes (approximately 7,600) — roughly 3x more than any viral group. This reflects fundamental biological differences: bacterial infections activate neutrophil-dominated innate immunity through TLR2/4 pattern recognition and acute phase responses, while viral infections primarily trigger more targeted interferon-mediated pathways.

### 2. IFI27 is a robust COVID-19 blood biomarker

IFI27 (interferon-stimulated gene 27) was the most significant upregulated gene in COVID-19 vs Healthy (log2FC = 3.95, padj = 4.2e-06) and remained significant in our age-adjusted, one-sample-per-subject analysis. This is consistent with published literature identifying IFI27 as one of the most reliable COVID-19 blood transcriptomic biomarkers.

### 3. Cell cycle and B cell activation are COVID-19 unique pathways

GSEA Reactome analysis identified two distinct biological themes in genes upregulated in COVID-19:

- **Adaptive immune activation:** B cell mediated immunity, immunoglobulin-mediated immune response, lymphocyte mediated immunity — reflecting strong antibody and B cell responses
- **Cell cycle/mitosis:** Mitotic prometaphase, cell cycle checkpoints, resolution of sister chromatid cohesion — reflecting lymphocyte proliferation and clonal expansion

### 4. Translation machinery is suppressed across all infections

GSEA Reactome identified suppression of Eukaryotic Translation Elongation, Peptide Chain Elongation, rRNA processing, and oxidative phosphorylation across all infection groups. This is consistent with host shutdown of translation machinery during acute infection — redirecting cellular resources toward immune defense.

### 5. DESeq2 and limma-voom show strong agreement, especially at the pathway level

Both methods were run on the same independent samples with the same design formula. The comparison operates at three levels:

**Gene-level fold change correlation:** Pearson r = 0.81, indicating strong agreement on direction and magnitude of individual gene expression changes.

**Gene-level significance concordance:** Ranges from 48–76% depending on comparison. Lower concordance reflects borderline genes near the significance threshold, not biological disagreement.

**Pathway-level GSEA Reactome concordance:** Higher than gene-level concordance for Influenza, Bacterial, and CoV_other, confirming that both methods converge on the same biological processes even when disagreeing on borderline genes. COVID-19 is an exception — see note below.

| Comparison | Gene concordance | Pathway concordance |
|---|---|---|
| COVID-19 | approx. 55% | 36.4% (see note) |
| Influenza | approx. 58% | 74.3% |
| Bacterial | approx. 75% | 73.2% |
| CoV other | approx. 58% | 73.4% |

Note on COVID-19: Pathway concordance is lower than gene concordance due to asymmetric pathway detection (DESeq2: 143 pathways, limma-voom: 322 pathways). Reduced statistical power after one-sample-per-subject restriction (n=45, 874 DE genes) likely contributes — with fewer DE genes, small differences in gene ranking push different pathways across the significance threshold.

**Key conclusion:** For three of four comparisons (Influenza, Bacterial, CoV_other), both methods identify the same biological pathways as significantly enriched, suggesting pathway-level conclusions are more robust to method choice than individual gene calls.

**Voom mean-variance trend:** The voom plot shows the expected pattern for well-behaved RNA-seq data — low-count genes (left) have high relative variance, while high-count genes (right) show decreasing variance. The smooth red trend line is used by voom to assign precision weights: high-variance genes get lower weights in the linear model, preventing noisy low-count genes from dominating the results. This is analogous to DESeq2's dispersion shrinkage but implemented through a different statistical framework.

**Practical recommendation:** For downstream analyses (biomarker selection, pathway analysis), using the intersection of DESeq2 and limma-voom results provides the most robust gene set. For exploratory analysis, either method alone is appropriate.

## Methods

### Statistical Model

```
Design formula: ~ batch + gender + age + cohort
```

- **Covariates (adjustment variables):** batch, gender, age
- **Variable of interest:** cohort (infection type vs healthy controls)
- Reference level: healthy controls
- DESeq2 with apeglm log2 fold change shrinkage
- limma-voom with TMM normalization and empirical Bayes moderation
- Both methods applied to identical independent samples for direct comparison
- Significance threshold: FDR < 0.05

**Fold change filtering:** FDR < 0.05 only — no |log2FC| cutoff applied. The DESeq2 authors recommend against hard FC cutoffs; apeglm shrinkage handles noisy low-count genes. Dashed lines at |log2FC| = 1 on volcano plots provide visual context without discarding data.

### Analysis Pipeline

1. Data download from GEO (Python: GEOparse)
2. Exploratory QC: library size, gene detection, PCA, outlier identification (Python: pandas, scikit-learn)
3. Sample deduplication: cross-cohort subject removal, one sample per subject (R)
4. Differential expression: DESeq2 and limma-voom with age-adjusted model (R)
5. Gene-level comparison: concordance table and fold change correlation (Pearson r)
6. GSEA Reactome: run independently from DESeq2 and limma-voom ranked lists, pathway concordance compared
7. Visualization: volcano plots, UpSet overlap plot, dual heatmaps (R: ggplot2, ggrepel, pheatmap, UpSetR)
8. Pathway comparison figures: side-by-side DESeq2 vs limma-voom top pathways per comparison

## Project Structure

```
scripts/
├── 04_real_geo_data.ipynb
└── 05_deseq2_analysis.R

figures/
└── GSE161731/
    ├── 01_cohort_overview.png
    ├── 02_qc_library_size.png
    ├── 03_pca_before_qc.png
    ├── 04_pca_after_qc.png
    ├── 05a_volcano_covid_vs_healthy.png
    ├── 05b_volcano_influenza_vs_healthy.png
    ├── 05c_volcano_bacterial_vs_healthy.png
    ├── 05d_volcano_cov_other_vs_healthy.png
    ├── 06_de_gene_counts.png
    ├── 07_upset_de_genes.png
    ├── 08a_gsea_reactome_covid.png
    ├── 08b_gsea_reactome_influenza.png
    ├── 08c_gsea_reactome_bacterial.png
    ├── 08d_gsea_reactome_cov_other.png
    ├── 09a_heatmap_deseq2_vst.png
    ├── 09b_heatmap_limma_voom.png
    ├── 10a_voom_mean_variance.png
    └── 10b_deseq2_vs_limma_fc.png

results/
├── DE_deseq2_COVID19_vs_healthy.csv
├── DE_deseq2_Influenza_vs_healthy.csv
├── DE_deseq2_Bacterial_vs_healthy.csv
├── DE_deseq2_CoV_other_vs_healthy.csv
├── DE_limma_COVID19_vs_healthy.csv
├── DE_limma_Influenza_vs_healthy.csv
├── DE_limma_Bacterial_vs_healthy.csv
├── DE_limma_CoV_other_vs_healthy.csv
├── GSEA_deseq2_COVID19_vs_healthy.csv
├── GSEA_deseq2_Influenza_vs_healthy.csv
├── GSEA_deseq2_Bacterial_vs_healthy.csv
├── GSEA_deseq2_CoV_other_vs_healthy.csv
├── GSEA_limma_COVID19_vs_healthy.csv
├── GSEA_limma_Influenza_vs_healthy.csv
├── GSEA_limma_Bacterial_vs_healthy.csv
└── GSEA_limma_CoV_other_vs_healthy.csv

data/                    # Raw downloads (not tracked by git)
```

## Tools

**Python:** NumPy, pandas, matplotlib, seaborn, scikit-learn, GEOparse, SciPy

**R:** DESeq2, apeglm, limma, edgeR, tidyverse, ggplot2, ggrepel, ReactomePA, clusterProfiler, org.Hs.eg.db, pheatmap, RColorBrewer, UpSetR, patchwork

## How to Reproduce

1. Clone this repository
2. Download raw data: run `scripts/04_real_geo_data.ipynb` to download GSE161731 from GEO
3. Run `scripts/05_deseq2_analysis.R` for all differential expression, GSEA, and visualization

## Limitations

- **Multi-cohort age variation:** Age varies systematically across cohorts due to different recruitment sources (university students vs emergency department patients). We include age as an adjustment variable, but residual confounding cannot be fully excluded in a cross-cohort observational design.
- **Blood transcriptomics:** Captures systemic immune response, not local infection site biology. Nasal or respiratory samples would show stronger pathogen-specific signatures, particularly for interferon-stimulated genes.
- **Unequal group sizes:** COVID-19 n=45 vs Influenza n=17 after deduplication — affects statistical power for Influenza comparisons.
- **One sample per subject:** Original authors used `duplicateCorrelation()` to retain all timepoints while accounting for within-subject correlation. We used one sample per subject (earliest timepoint) to enable fair DESeq2 vs limma-voom comparison. Results for COVID-19 may differ from the original paper due to reduced sample size.

## References

McClain, M.T., Constantine, F.J., Henao, R. et al. Dysregulated transcriptional responses to SARS-CoV-2 in the periphery. *Nature Communications* 12, 1079 (2021). https://doi.org/10.1038/s41467-021-21289-y

Stupnikov, A. et al. Robustness of differential gene expression analysis of RNA-seq. *Computational and Structural Biotechnology Journal* 19, 3470-3481 (2021). https://doi.org/10.1016/j.csbj.2021.05.040

Costa-Silva, J., Domingues, D. & Lopes, F.M. RNA-Seq differential expression analysis: An extended review and a software tool. *PLOS ONE* 12, e0190152 (2017). https://doi.org/10.1371/journal.pone.0190152