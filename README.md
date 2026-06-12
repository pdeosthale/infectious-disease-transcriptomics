# infectious-disease-transcriptomics
Computational analysis of host-pathogen interactions using public transcriptomic datasets
# Infectious Disease Transcriptomics

Computational analysis of host immune responses to respiratory infections using public blood transcriptomic data. This project applies differential expression analysis to characterize pathogen-specific and shared transcriptomic signatures across viral and bacterial respiratory infections.

## Biological Question

How does the blood transcriptome differ across respiratory infections (COVID-19, Influenza, other coronaviruses, bacterial pneumonia) compared to healthy controls? Are there shared host response genes across all infections, or does each pathogen trigger a distinct transcriptomic program?

## Dataset

**GSE161731** — Blood RNA-seq from patients with respiratory infections, published in *Nature Communications*.

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

## Key Findings

### 1. Bacterial infection triggers the largest transcriptomic shift

Even after correcting for age confounding, bacterial infection produced ~8,400 DE genes — roughly 3x more than any viral infection. This reflects fundamental biological differences: bacterial infections activate neutrophil-dominated innate immunity through TLR2/4 pattern recognition and acute phase responses, while viral infections primarily trigger more targeted interferon-mediated pathways.

### 2. Age confounding substantially inflated the bacterial signal

Healthy controls in this dataset are college-age students (median 18), while bacterial patients are elderly (median 60). Without adjusting for age, bacterial infection appeared to have 10,564 DE genes. Adding age as an adjustment variable removed ~2,128 genes (20%), and COVID-19 dropped by ~996 genes (26%). This demonstrates that clinical cohort studies with unbalanced demographics can produce misleading results without proper covariate modeling.

| Comparison | Without Age Adjustment | With Age Adjustment | Removed | % Change |
|---|---|---|---|---|
| Bacterial vs Healthy | 10,564 | 8,436 | 2,128 | 20.1% |
| COVID-19 vs Healthy | 3,816 | 2,820 | 996 | 26.1% |

### 3. Viral infections share more with each other than with bacterial

COVID-19, Influenza, and CoV other show substantial overlap in DE genes, consistent with shared type I/III interferon responses and antiviral ISG programs. The shared core across all four infections represents a universal host response to respiratory infection — likely general inflammation, stress response, and innate immune activation.

### 4. COVID-19 has a distinct transcriptomic signature

IFI27 was the strongest upregulated gene in COVID-19 (log2FC = 4.15), consistent with published literature identifying it as a COVID-19 blood biomarker. Immunoglobulin genes (IGLV3-19, IGHV1-69-2) reflect strong B cell activation. GDF15 (growth differentiation factor 15), a stress cytokine associated with COVID severity, was also among the top upregulated genes.

### 5. The interferon response is muted in blood

Only 3 of 12 classic interferon-stimulated genes (ISGs) were significantly upregulated in blood (IFI27, IFITM3, OAS1). This is expected — the interferon response is strongest at the mucosal infection site (nasal epithelium, lungs). In blood, the signal is diluted across all circulating cell types. Nasal swab RNA-seq would show much stronger ISG activation.

### 6. Downregulation dominates across all infections

More genes are downregulated than upregulated in every infection group. This is consistent with the known biology of acute infection: the immune system redirects transcriptional resources toward defense, suppressing baseline functions like ribosomal biogenesis and translation. Ribosomal proteins (RPS24, RPL6, RPL11) were among the most significantly downregulated genes.

## Methods

### Statistical Model

```
Design formula: ~ batch + gender + age + cohort
```

- **Covariates (adjustment variables):** batch (technical), gender (biological), age (biological)
- **Variable of interest:** cohort (infection type vs healthy)
- Reference level: healthy controls
- DESeq2 with apeglm log2 fold change shrinkage
- Significance threshold: FDR < 0.05

### Analysis Pipeline

1. Data download from GEO (Python: GEOparse)
2. Exploratory QC: library size, gene detection, PCA, outlier identification (Python: pandas, scikit-learn)
3. Differential expression: DESeq2 with age-corrected model (R)
4. Visualization: custom volcano plots, UpSet overlap plots (R: ggplot2, ggrepel, UpSetR)
5. Pathway analysis: GO enrichment of COVID-19 unique and shared infection response genes (R: clusterProfiler)

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
    └── 09_go_shared_infection.png

results/                          # Processed data and DE results
├── GSE161731_counts_filtered.csv
├── GSE161731_metadata_clean.csv
├── DE_COVID19_vs_healthy_age_corrected.csv
├── DE_Influenza_vs_healthy_age_corrected.csv
├── DE_Bacterial_vs_healthy_age_corrected.csv
└── DE_CoV_other_vs_healthy_age_corrected.csv

data/                             # Raw downloads (not tracked by git)
```

## Tools

**Python:** NumPy, pandas, matplotlib, seaborn, scikit-learn, GEOparse, SciPy

**R:** DESeq2, apeglm, limma, tidyverse, ggplot2, ggrepel, clusterProfiler, org.Hs.eg.db, pheatmap, UpSetR

## How to Reproduce

1. Clone this repository
2. Download raw data: run `scripts/04_real_geo_data.ipynb` to download from GEO
3. Run QC and filtering in Python notebooks (01–04)
4. Run `scripts/05_deseq2_analysis.R` for differential expression analysis

## Limitations

- Severe age imbalance between cohorts (healthy controls are college students). Adjusted for statistically but cannot fully eliminate confounding.
- Blood transcriptomics captures systemic immune response, not local infection site biology. Nasal or respiratory samples would show stronger pathogen-specific signatures.
- Unequal group sizes (COVID-19 n=77 vs Influenza n=17) affect statistical power for some comparisons.
- 21 CoV other samples excluded due to missing gender data, reducing that group from 60 to 39.

## References

Mick et al. (2020). Upper airway gene expression reveals suppressed immune responses to SARS-CoV-2 compared with other respiratory viruses. *Nature Communications*, 11, 5854.

## Status

Module 1 (Python Foundations) — Complete
Module 4 (Bulk RNA-seq / DESeq2) — Complete