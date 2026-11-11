# Infectious Disease Blood Transcriptomics

Differential gene expression and pathway analysis of blood RNA-seq from patients with COVID-19, influenza, bacterial infection, and other coronavirus infections versus healthy controls (GSE161731, McClain et al. 2021, *Nat Commun*).

---

## Biological Question

Do respiratory infections leave distinct transcriptional signatures in blood, and can shared vs infection-specific programs be identified? Specifically: which pathways are universally dysregulated across all four infection types, which are viral-specific, and which are unique to coronavirus infections?

---

## Dataset

| Accession | Description | Cohorts | n |
|---|---|---|---|
| [GSE161731](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161731) | McClain et al. 2021, *Nat Commun* | COVID-19, Influenza, Bacterial, CoV other, Healthy | 136 samples (analysis set; see QC below) |

**Design:** 5 cohorts, 2 batches. Age is a significant covariate — bacterial cohort median age ~60 vs ~20 for other groups. Model: `~ batch + gender + age + cohort`.

---

## Repository Structure

```
infectious-disease-transcriptomics/
│
├── scripts/
│   ├── 01_data_download_qc.ipynb        # GEO download, cohort overview, QC, PCA (Python)
│   ├── 02_limma_analysis.R              # Main analysis: limma-voom DE, DEG overlaps (UpSet, Venn),
│   │                                    #   Reactome GSEA/ORA, top-50 heatmap, splicing boxplots
│   └── 02_deseq2_limma_analysis.R       # Supplementary: DESeq2 vs limma-voom vs edgeR comparison
│
├── data/                                # Downloaded files (not tracked)
├── results/                             # Filtered counts, metadata, DE results (not tracked)
│
└── figures/GSE161731/
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
```

---

## Methods

### Quality Control
- Library size and gene detection assessed per cohort
- PCA-based outlier removal: 2 extreme outliers identified on PC1 and removed (198 → 196 samples)
- Sample filtering before DE: 2 subjects appearing in more than one cohort removed (→ 192); earliest timepoint kept per subject (→ 151); samples missing gender (→ 137) or age (→ 136) removed
- Final analysis set: 136 independent samples — healthy 17, COVID-19 45, Influenza 17, Bacterial 23, CoV other 34

### Differential Expression
limma-voom with model `~ batch + gender + age + cohort`. One sample per subject (earliest timepoint). Healthy controls as reference. Four contrasts: COVID-19, Influenza, Bacterial, CoV other vs Healthy. FDR < 0.05, |log2FC| > 1.

### Overlap Analysis
DEG union across contrasts — UpSet plot and four-way Venn diagram to identify infection-specific and shared gene sets.

### Pathway Analysis
- **GSEA**: Reactome gene sets, ranked by sign(log2FC) × −log10(adj.p), per infection type
- **ORA**: Reactome, on infection-specific and shared DEG subsets defined by overlap analysis

---

## Key Results

### 1. Bacterial infection produces the largest transcriptional response

Bacterial infection shows the most extensive DEG signature (3,382 upregulated / 4,975 downregulated vs healthy), followed by Influenza (1,360 / 1,460), CoV other (978 / 1,473), and COVID-19 (250 / 636). COVID-19 produces a notably modest response relative to the other infections at this timepoint.

### 2. Ribosomal gene suppression is a universal infection signature

373 genes are differentially expressed across all four infections. Pathway analysis of this pan-infection signature identifies ribosomal biogenesis and translation as the dominant theme — eukaryotic translation initiation/elongation, NMD, peptide chain elongation all strongly suppressed. This pattern is visible in all four volcano plots and confirmed in the top 50 DE gene heatmap, which is dominated by RPL/RPS family genes.

### 3. Interferon signaling distinguishes viral from bacterial infections

Viral infections (Influenza, CoV other) show strong upregulation of interferon signaling (IFN-α/β, IFN-γ) by GSEA — absent in bacterial infection. Influenza shows the clearest IFN signature (Interferon Signaling NES ~2.3, Interferon gamma signaling NES ~2.2, Interferon alpha/beta signaling NES ~2.1). Bacterial infection instead activates neutrophil degranulation and Rho GTPase/cytoskeletal pathways.

### 4. Coronavirus infections share an mRNA splicing signature absent in other infections

Genes shared between COVID-19 and CoV other but not Influenza or Bacterial infections are enriched for mRNA splicing pathways — "Processing of Capped Intron-Containing Pre-mRNA", "mRNA Splicing — Major Pathway", "mRNA Splicing" (all FDR < 0.05, gene ratio ~0.22). Individual splicing factor genes (LSM6, LSM8, SF3A1, SLU7) are significantly downregulated in COVID-19 and CoV other vs healthy, but not in Influenza or Bacterial infections.

**Important caveat:** The coronavirus group (CoV other) in this dataset consists predominantly of young adults (~19 years median age) while COVID-19 patients are older and more severely ill. Age and disease severity confounding cannot be fully excluded despite model adjustment. This finding warrants validation in age-matched cohorts.

### 5. COVID-19-specific genes are enriched for cell cycle and mitotic pathways

145 genes are specific to COVID-19 (not shared with other infections). Reactome ORA identifies Cell Cycle, Cell Cycle Mitotic, M Phase, and Cell Cycle Checkpoints as top enriched pathways — consistent with coronavirus-driven dysregulation of host cell cycle machinery and B cell proliferative responses during adaptive immune activation.

---

## Limitations

- **Multi-cohort age variation:** Age varies systematically across cohorts due to different recruitment sources (university students vs emergency department patients). We include age as an adjustment variable, but residual confounding cannot be fully excluded in a cross-cohort observational design.
- **Blood transcriptomics:** Captures systemic immune response, not local infection site biology. Nasal or respiratory samples would show stronger pathogen-specific signatures, particularly for interferon-stimulated genes.
- **Unequal group sizes:** COVID-19 n=45 vs Influenza n=17 after deduplication — affects statistical power for Influenza comparisons.
- **One sample per subject:** Original authors used `duplicateCorrelation()` to retain all timepoints while accounting for within-subject correlation. We used one sample per subject (earliest timepoint) to enable fair DESeq2 vs limma-voom comparison. Results for COVID-19 may differ from the original paper due to reduced sample size.

---

## Tools and Languages

**Python:** GEOparse, pandas, matplotlib, seaborn (QC, cohort overview, PCA)

**R:** limma, edgeR, clusterProfiler, ReactomePA, ggplot2, ComplexHeatmap, ggVennDiagram, UpSetR, tidyverse

---

## Reference

McClain MT, Constantine FJ, Henao R, et al. Dysregulated transcriptional responses to SARS-CoV-2 in the periphery. *Nat Commun.* 2021;12:1079. https://doi.org/10.1038/s41467-021-21289-y

## License

MIT License — see [LICENSE](LICENSE) for details.
