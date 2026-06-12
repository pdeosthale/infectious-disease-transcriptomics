# Infectious Disease Transcriptomics

Computational analysis of host immune responses to respiratory infections using public blood transcriptomic data. This project applies differential expression analysis to characterize pathogen-specific and shared transcriptomic signatures across viral and bacterial respiratory infections.

## Biological Question

How does the blood transcriptome differ across respiratory infections (COVID-19, Influenza, other coronaviruses, bacterial pneumonia) compared to healthy controls? Are there shared host response genes across all infections, or does each pathogen trigger a distinct transcriptomic program?

## Dataset

**GSE161731** Blood RNA-seq from patients with respiratory infections.

Source: McClain, M.T., Constantine, F.J., Henao, R. et al. Dysregulated transcriptional responses to SARS-CoV-2 in the periphery. *Nature Communications* 12, 1079 (2021). https://doi.org/10.1038/s41467-021-21289-y

| Cohort | Samples (after QC) | Median Age |
|--------|-------------------|------------|
| COVID-19 | 77 | 33 |
| CoV other | 39 | 18 |
| Influenza | 17 | 18 |
| Bacterial | 23 | 60 |
| Healthy | 19 | 18 |

- 28,996 genes after filtering (detected in ≥10% of samples)
- 2 outlier samples removed during QC (low library size and gene detection)
- 21 samples removed due to missing gender data (all from CoV other cohort)

## Disclaimer

This is an independent secondary analysis of publicly available data from GEO (GSE161731) for educational and portfolio purposes. The original study was designed and conducted by McClain et al. (2021). Any interpretations, limitations identified, or analytical choices in this project reflect the author's own learning process and do not represent the views or conclusions of the original investigators. Observations about cohort characteristics (such as age distribution) describe features of the dataset as used in this reanalysis and are not intended as criticism of the original study design.

## Key Findings

### 1. Bacterial infection triggers the largest transcriptomic shift

Even after correcting for age confounding, bacterial infection produced ~8,400 DE genes; roughly 3x more than any viral infection. This reflects fundamental biological differences: bacterial infections activate neutrophil dominated innate immunity through TLR2/4 pattern recognition and acute phase responses, while viral infections primarily trigger more targeted interferon mediated pathways.

### 2. Age confounding substantially inflated the bacterial signal

Healthy controls in this dataset are young (median 18), while bacterial patients are elderly (median 60). Without adjusting for age, bacterial infection appeared to have 10,564 DE genes. Adding age as an adjustment variable removed ~2,128 genes (20%), and COVID-19 dropped by ~996 genes (26%). This demonstrates that clinical cohort studies with unbalanced demographics can produce misleading results without proper covariate modeling.

| Comparison | Without Age Adjustment | With Age Adjustment | Removed | % Change |
|---|---|---|---|---|
| Bacterial vs Healthy | 10,564 | 8,436 | 2,128 | 20.1% |
| COVID-19 vs Healthy | 3,816 | 2,820 | 996 | 26.1% |

### 3. Viral infections share more with each other than with bacterial

COVID-19, Influenza, and CoV other show substantial overlap in DE genes, consistent with shared type I/III interferon responses and antiviral ISG programs. The shared core across all four infections represents a universal host response to respiratory infection — likely general inflammation, stress response, and innate immune activation.

### 4. COVID-19 has a distinct transcriptomic signature

IFI27 was the strongest upregulated gene in COVID-19 (log2FC = 4.15), consistent with published literature identifying it as a COVID-19 blood biomarker. Immunoglobulin genes (IGLV3-19, IGHV1-69-2) reflect strong B cell activation. GDF15 (growth differentiation factor 15), a stress cytokine associated with COVID severity, was also among the top upregulated genes.

### 5. The interferon response is muted in blood

Only 3 of 12 classic interferon-stimulated genes (ISGs) were significantly upregulated in blood (IFI27, IFITM3, OAS1). This is expected since the interferon response is strongest at the mucosal infection site (nasal epithelium, lungs). In blood, the signal is diluted across all circulating cell types. Nasal swab RNAseq would show much stronger ISG activation.

### 6. Downregulation dominates across all infections

More genes are downregulated than upregulated in every infection group. This is consistent with the known biology of acute infection: the immune system redirects transcriptional resources toward defense, suppressing baseline functions like ribosomal biogenesis and translation. Ribosomal proteins (RPS24, RPL6, RPL11) were among the most significantly downregulated genes.

### 7. DESeq2 and Limma-voom show agreement

Both methods were run on the same data with the same model formula (`~ batch + gender + age + cohort`). This dual method approach builds confidence in the findings that genes significant in both are high confidence hits, while method specific genes are borderline.

**Fold change correlation:** Pearson r = 0.81, indicating strong agreement on the direction and magnitude of expression changes. The methods agree on the biology; differences occur at the significance threshold where borderline genes (padj ≈ 0.05) get pushed to different sides of the cutoff by different statistical models.

**Concordance in significance calls:**

| Comparison | DESeq2 | Limma-voom | Both | Concordance |
|---|---|---|---|---|
| COVID-19 | 2,820 | 2,821 | 1,843 | 48.5% |
| Influenza | 3,299 | 2,661 | 2,192 | 58.2% |
| Bacterial | 8,436 | 8,253 | 7,187 | 75.6% |
| CoV other | 2,946 | 2,746 | 2,088 | 57.9% |

**Key observations:**

- Bacterial has the highest concordance (75.6%) : when the biological signal is strong, both methods agree.
- COVID-19 has the lowest concordance (48.5%) : many COVID DE genes are borderline, hovering near the significance threshold. The 1,843 genes significant in both methods represent the highest-confidence COVID-19 signature.
- Total DE gene counts are remarkably similar between methods (e.g., COVID-19: 2,820 vs 2,821), confirming that overall sensitivity is comparable.
- Published benchmarking studies show that concordance between DE methods varies substantially depending on the dataset, effect size, and sample characteristics (Costa-Silva et al., 2017; Stupnikov et al., 2021). Clean experimental data with strong effects typically shows 80-90%+ overlap, while complex clinical datasets with weaker signals and confounders tend to show lower concordance, consistent with our observations.

**Voom mean-variance trend:** The voom plot shows the expected pattern for good RNAseq data where low-count genes (left) cluster with variable variance, while high-count genes (right) show decreasing variance. The smooth red trend line is used by voom to assign precision weights: high-variance genes get lower weights in the linear model, preventing noisy low-count genes from dominating the results. This is analogous to DESeq2's dispersion shrinkage but implemented through a different statistical framework.

**Practical recommendation:** For downstream analyses (biomarker selection, pathway analysis), using the intersection of DESeq2 and Limma-voom results provides the most robust gene set. For exploratory analysis, either method alone is appropriate.

## Methods

### Statistical Model

```
Design formula: ~ batch + gender + age + cohort
```

- **Covariates (adjustment variables):** batch (technical), gender (biological), age (biological)
- **Variable of interest:** cohort (infection type vs healthy)
- Reference level: healthy controls
- **DESeq2** with apeglm log2 fold change shrinkage
- **Limma-voom** with TMM normalization and empirical Bayes moderation
- Both methods applied with the same design formula for direct comparison
- Significance threshold: FDR < 0.05

### Analysis Pipeline

1. Data download from GEO (Python: GEOparse)
2. Exploratory QC: library size, gene detection, PCA, outlier identification (Python: pandas, scikit-learn)
3. Differential expression with DESeq2: age-corrected model, apeglm shrinkage (R)
4. Differential expression with Limma-voom: TMM normalization, empirical Bayes (R)
5. Methods comparison: concordance analysis and fold change correlation between DESeq2 and Limma-voom
6. Visualization: custom volcano plots, UpSet overlap plots (R: ggplot2, ggrepel, UpSetR)
7. Pathway analysis: GO enrichment of COVID-19 unique and shared infection response genes (R: clusterProfiler)

## Project Structure

```
scripts/                          # Analysis notebooks and scripts
├── 01_python_foundations.ipynb    # Python basics with simulated expression data
├── 02_python_programming.ipynb   # Functions, loops, file I/O, pairwise analysis
├── 03_numpy_and_data_wrangling.ipynb  # NumPy, z-scores, correlation, PCA
├── 04_real_geo_data.ipynb        # GEO data download, QC, outlier detection
└── 05_deseq2_analysis.R          # Full DESeq2 pipeline with covariates

figures/
├── practice/                     # Learning exercises with simulated data
│   ├── 01_cohort_overview.png
│   ├── 02_expression_heatmap.png
│   ├── 03_isg_analysis.png
│   ├── 04_pairwise_log2fc.png
│   ├── 05_gene_correlation.png
│   ├── 06_pca_analysis.png
│   └── 07_immune_panel.png
└── GSE161731/                    # Publication-quality figures from real data
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
    ├── 08_go_covid_unique.png
    ├── 09_go_shared_infection.png
    ├── 10_voom_mean_variance.png
    └── 11_deseq2_vs_limma_fc.png

results/                          # Processed data and DE results
├── GSE161731_counts_filtered.csv
├── GSE161731_metadata_clean.csv
├── DE_COVID19_vs_healthy_age_corrected.csv
├── DE_Influenza_vs_healthy_age_corrected.csv
├── DE_Bacterial_vs_healthy_age_corrected.csv
├── DE_CoV_other_vs_healthy_age_corrected.csv
├── DE_limma_COVID19_vs_healthy.csv
├── DE_limma_Influenza_vs_healthy.csv
├── DE_limma_Bacterial_vs_healthy.csv
└── DE_limma_CoV_other_vs_healthy.csv

data/                             # Raw downloads (not tracked by git)
```

## Tools

**Python:** NumPy, pandas, matplotlib, seaborn, scikit-learn, GEOparse, SciPy

**R:** DESeq2, apeglm, Limma, edgeR, tidyverse, ggplot2, ggrepel, clusterProfiler, org.Hs.eg.db, pheatmap, UpSetR

## How to Reproduce

1. Clone this repository
2. Download raw data: run `scripts/04_real_geo_data.ipynb` to download from GEO
3. Run QC and filtering in Python notebooks (01–04)
4. Run `scripts/05_deseq2_analysis.R` for differential expression and Limma-voom analysis (runs both methods and generates comparison)

## Limitations

- Severe age imbalance between cohorts (healthy controls are younger). Adjusted for statistically but cannot fully eliminate confounding.
- Blood transcriptomics captures systemic immune response, not local infection site biology. Nasal or respiratory samples would show stronger pathogen-specific signatures.
- Unequal group sizes (COVID-19 n=77 vs Influenza n=17) affect statistical power for some comparisons.
- 21 CoV other samples excluded due to missing gender data, reducing that group from 60 to 39.
- DESeq2 vs Limma-voom concordance for COVID-19 is 48.5%, meaning roughly half of the significant genes are method specific. Downstream analyses should consider using the intersection (high-confidence) or union (exploratory) depending on the goal.

## References

- McClain, M.T., Constantine, F.J., Henao, R. et al. Dysregulated transcriptional responses to SARS-CoV-2 in the periphery. *Nature Communications* 12, 1079 (2021). https://doi.org/10.1038/s41467-021-21289-y
- Stupnikov, A. et al. Robustness of differential gene expression analysis of RNA-seq. Computational and Structural Biotechnology Journal 19, 3470-3481 (2021). https://doi.org/10.1016/j.csbj.2021.05.040
- Costa-Silva, J., Domingues, D. & Lopes, F.M. RNA-Seq differential expression analysis: An extended review and a software tool. PLOS ONE 12, e0190152 (2017). https://doi.org/10.1371/journal.pone.0190152


