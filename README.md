# Infectious Disease Transcriptomics

Host blood transcriptome analysis across four respiratory infections using publicly available RNA-seq data. This project characterizes pathogen-specific and shared transcriptomic signatures, identifies infection-specific gene sets, and uses pathway analysis to interpret the biological processes driving each host response.

## Biological Question

What transcriptomic signatures are shared vs unique across COVID-19, Influenza, bacterial pneumonia, and seasonal coronavirus infections? Can we identify a pan-infection host response, a shared viral signature, and pathogen-specific gene sets — and what biological pathways do each represent?

## Disclaimer

This is an independent secondary analysis of publicly available data from GEO (GSE161731) for educational and portfolio purposes. The original study was designed and conducted by McClain et al. (2021). Analytical choices and interpretations reflect the author's own learning process and do not represent the views of the original investigators.

## Dataset

**GSE161731** — Blood RNA-seq from patients with acute respiratory infections and healthy controls.

Source: McClain, M.T., Constantine, F.J., Henao, R. et al. Dysregulated transcriptional responses to SARS-CoV-2 in the periphery. *Nature Communications* 12, 1079 (2021). https://doi.org/10.1038/s41467-021-21289-y

| Cohort | Samples (after QC) | Median Age |
|---|---|---|
| Healthy | 17 | 18 |
| COVID-19 | 45 | 33 |
| Influenza | 17 | 18 |
| Bacterial | 23 | 60 |
| CoV other | 36 | 18 |

**Final sample count: 136 independent samples** after all QC and deduplication steps.

**QC steps applied:**

- Gene filtering: removed genes detected in fewer than 10% of samples (60,675 to 28,996 genes)
- Outlier removal: 2 samples removed via PCA (very low library size and gene detection)
- Cross-cohort subjects: 2 subjects (A34BBC, E77BE3) removed — appeared as both healthy and CoV_other
- Pseudoreplication: 14 COVID-19 patients had multiple timepoints; kept earliest timepoint only
- Missing gender: 21 samples removed (all CoV_other) to include gender as adjustment variable
- Missing age: samples with NA age removed to include age as adjustment variable

## Method

**limma-voom** with TMM normalization — consistent with the original McClain et al. paper which used limma GLMs.

limma-voom is not a generalized linear model (GLM). Rather than modeling count variance through a distributional assumption (as DESeq2 does with the negative binomial), voom estimates the mean-variance relationship non-parametrically and converts it into per-observation precision weights. These weights are fed into a standard weighted least squares linear model, enabling the full flexibility of limma's empirical Bayes moderation.

**Design formula:** `~ batch + gender + age + cohort`

- Covariates: batch (technical), gender (biological), age (biological — varies systematically across cohorts due to different recruitment sources)
- Variable of interest: cohort (infection type vs healthy controls)
- Reference level: healthy controls
- Significance threshold: FDR < 0.05

**Key difference from original paper:** The original authors used `duplicateCorrelation()` to retain all timepoints while accounting for within-subject correlation. We used one sample per subject (earliest timepoint) to ensure all observations are independent and simplify the analysis.

## Key Analytical Choices

### 1. One sample per subject

14 COVID-19 patients were sampled at multiple timepoints (early/middle/late disease). We kept only the earliest timepoint to ensure independence. This is more conservative than `duplicateCorrelation()` and means COVID-19 results reflect early disease only — relevant for interpreting effect sizes relative to other cohorts sampled at peak illness.

### 2. Age adjustment

Age varies systematically across cohorts (healthy median 18, bacterial median 60) due to different recruitment sources. We included age as an adjustment variable to isolate infection-specific transcriptomic effects from age-related gene expression differences.

## Key Findings

### 1. Bacterial infection triggers the largest transcriptomic shift

Bacterial pneumonia produced the most DE genes, with neutrophil degranulation as the top activated pathway. This reflects fundamental biological differences: bacterial infections activate neutrophil-dominated innate immunity through TLR2/4 pattern recognition, while viral infections primarily trigger interferon-mediated pathways.

### 2. Translation machinery is universally suppressed across all four infections

The pan-infection gene set (shared across COVID-19, Influenza, Bacterial, and CoV_other) is strongly enriched for translation machinery suppression — ribosomal subunits, peptide chain elongation, eukaryotic translation initiation, and nonsense-mediated decay pathways (FDR approximately 10^-50). This universal host response likely reflects redirection of cellular resources toward immune defense during acute infection.

### 3. COVID-19-specific signature is enriched for cell cycle and mitotic pathways

Genes uniquely differentially expressed in COVID-19 (not in Influenza, Bacterial, or CoV_other) are enriched for cell cycle, M phase, mitotic checkpoints, and cell cycle checkpoints. This is consistent with lymphocyte proliferation and B cell clonal expansion — a hallmark of SARS-CoV-2 adaptive immune activation documented in the original paper.

### 4. Influenza drives a stronger interferon response than COVID-19

Influenza shows the clearest interferon signature of all four infections — Interferon Signaling, Interferon gamma signaling, and Interferon alpha/beta signaling are all among the top activated pathways. COVID-19 notably lacks a prominent interferon signature in the GSEA results, consistent with published evidence that SARS-CoV-2 actively suppresses type I interferon responses through multiple mechanisms (NSP1, NSP16).

### 5. Coronavirus-specific suppression of mRNA splicing machinery

ORA on the coronavirus-specific gene set (COVID-19 and CoV_other only, absent in Influenza and Bacterial) identified enrichment of mRNA splicing pathways. The driving genes are LSM6, LSM8, SLU7, and SF3A1 — components of the spliceosome LSm2-8 ring and U2 snRNP complex.

Expression analysis confirms:

| Gene | Role | COVID-19 logFC | CoV_other logFC | Influenza | Bacterial |
|---|---|---|---|---|---|
| LSM8 | LSm ring, U6 snRNA stabilization | -0.90 * | -0.93 *** | ns | ns |
| LSM6 | LSm ring, U6 snRNA stabilization | -0.47 * | -0.57 *** | ns | ns |
| SLU7 | 3' splice site recognition | -0.54 * | -0.62 ** | ns | ns |
| SF3A1 | U2 snRNP, branch point recognition | +0.51 * | +0.46 ** | ns | ns |

LSM6, LSM8, and SLU7 are downregulated in both coronavirus infections while SF3A1 is upregulated — a structured bidirectional pattern suggesting specific disruption of spliceosome assembly rather than nonspecific transcriptional noise.

This is consistent with published evidence that SARS-CoV-2 NSP16 binds U1/U2 snRNA to suppress global host mRNA splicing (Bhatt et al. 2021, Cell). Crucially, this pattern is present in CoV_other (seasonal coronaviruses) as well as COVID-19, confirming it is a coronavirus class effect rather than SARS-CoV-2-specific.

**Important caveat:** Effect sizes are modest (|log2FC| 0.47-0.93) with high inter-individual variability. The stronger effect in CoV_other relative to COVID-19 likely reflects the analytical choice to use only the earliest COVID-19 timepoint — peak splicing dysregulation may occur later in disease. This is an exploratory finding requiring validation in larger cohorts. A longitudinal analysis is planned as a follow-up in the covid19-host-response-dynamics repository.

### 6. CoV_other-specific signature shows unexpected neuronal pathway enrichment

Genes uniquely differentially expressed in CoV_other show enrichment for nervous system development, axon guidance, and neuronal system pathways (FDR 0.2-0.6). These findings have high FDR and should be considered hypothesis-generating only. The neuronal enrichment may reflect cohort-specific biology, the younger age of CoV_other subjects, or neurotropic effects of seasonal coronaviruses.

## Analysis Pipeline

1. Data download: Python — GEOparse, count matrix and key file from GEO FTP
2. Exploratory QC: Python — library size, gene detection, PCA, outlier identification
3. Sample QC and deduplication: R — cross-cohort subject removal, one sample per subject
4. Differential expression: R — limma-voom, TMM normalization, empirical Bayes moderation
5. Gene set definition: R — infection-specific, shared viral, pan-infection, coronavirus gene sets
6. Pathway analysis: R — GSEA (Reactome, per infection) and ORA (Reactome, per gene set)
7. Visualization: R — volcano plots, UpSet, Venn, GSEA dot plots, ORA dot plots, heatmap, boxplots

## Project Structure

```
scripts/
├── 01_data_download_qc.ipynb        # Python: data download and QC
└── 02_limma_voom_analysis.R         # R: DE analysis, pathway analysis, visualization

figures/
└── GSE161731/
    ├── 00_voom_mean_variance.png
    ├── 01_cohort_overview.png
    ├── 02_qc_library_size.png
    ├── 03_pca_before_qc.png
    ├── 04_pca_after_qc.png
    ├── 05a_volcano_covid_vs_healthy.png
    ├── 05b_volcano_influenza_vs_healthy.png
    ├── 05c_volcano_bacterial_vs_healthy.png
    ├── 05d_volcano_cov_other_vs_healthy.png
    ├── 06_de_gene_counts.png
    ├── 07_upset_infection_overlap.svg
    ├── 08_venn_four_way.png
    ├── 09a_gsea_covid19.png
    ├── 09a_gsea_influenza.png
    ├── 09a_gsea_bacterial.png
    ├── 09a_gsea_cov_other.png
    ├── 09b_ora_covid_specific.png
    ├── 09c_ora_influenza_specific.png
    ├── 09d_ora_bacterial_specific.png
    ├── 09e_ora_cov_other_specific.png
    ├── 09f_ora_shared_viral.png
    ├── 09g_ora_shared_all.png
    ├── 09h_ora_coronavirus.png
    ├── 10_heatmap_top50.png
    └── 11_boxplot_coronavirus_splicing.png

results/
├── DE_limma_COVID19_vs_healthy.csv
├── DE_limma_Influenza_vs_healthy.csv
├── DE_limma_Bacterial_vs_healthy.csv
├── DE_limma_CoV_other_vs_healthy.csv
├── GSEA_COVID19_vs_healthy.csv
├── GSEA_Influenza_vs_healthy.csv
├── GSEA_Bacterial_vs_healthy.csv
├── GSEA_CoV_other_vs_healthy.csv
├── ORA_covid_specific.csv
├── ORA_influenza_specific.csv
├── ORA_bacterial_specific.csv
├── ORA_cov_other_specific.csv
├── ORA_shared_viral.csv
├── ORA_shared_all.csv
└── ORA_coronavirus.csv

data/                                # Downloaded from GEO on first run (not tracked)
```

## Tools

**Python:** GEOparse, NumPy, pandas, matplotlib, seaborn, scikit-learn

**R:** limma, edgeR, tidyverse, ggplot2, ggrepel, pheatmap, RColorBrewer, ReactomePA, clusterProfiler, org.Hs.eg.db, AnnotationDbi, UpSetR, VennDiagram, patchwork

## How to Reproduce

1. Clone this repository
2. Run `scripts/01_data_download_qc.ipynb` — downloads GSE161731 from GEO and performs initial QC
3. Run `scripts/02_limma_voom_analysis.R` — all DE analysis, pathway analysis, and figures

## Limitations

- **Blood transcriptomics:** Captures systemic immune response, not local infection site biology. Pathogen-specific signatures would be stronger in nasal or respiratory samples.
- **Early timepoint only for COVID-19:** 14 COVID-19 patients had longitudinal samples; only the earliest timepoint was retained. Results for COVID-19 reflect early disease and may underestimate effects that develop at peak viral load.
- **Age imbalance across cohorts:** Age varies systematically by cohort due to different recruitment sources. Age is included as a covariate but residual confounding cannot be fully excluded.
- **Small group sizes:** Influenza n=17 and healthy n=17 limit statistical power for those comparisons.
- **Exploratory ORA findings:** Several infection-specific ORA results have FDR above 0.1 (notably Influenza-specific and CoV_other-specific gene sets). These should be interpreted as hypothesis-generating.
- **Coronavirus splicing finding:** Effect sizes are modest with high inter-individual variability. Requires validation in larger cohorts.

## References

McClain, M.T., Constantine, F.J., Henao, R. et al. Dysregulated transcriptional responses to SARS-CoV-2 in the periphery. *Nature Communications* 12, 1079 (2021). https://doi.org/10.1038/s41467-021-21289-y

Law, C.W., Chen, Y., Shi, W. & Smyth, G.K. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biology* 15, R29 (2014). https://doi.org/10.1186/gb-2014-15-2-r29

Bhatt, D.L. et al. SARS-CoV-2 Disrupts Splicing, Translation, and Protein Trafficking to Suppress Host Defenses. *Cell* 183, 1325-1339 (2020). https://doi.org/10.1016/j.cell.2020.10.004

Liao, Y., Smyth, G.K. & Shi, W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics* 30, 923-930 (2014). https://doi.org/10.1093/bioinformatics/btt656