# Validation of a Machine Learning-Derived Methylation-Based Diagnostic Classifier for Sinonasal Tumors

## Project Overview

This project externally validates the Jurmeister et al. (2022) 18-class DNA-methylation-based classifier for sinonasal tumors in an independent, retrospective archival cohort at RCSI, as part of the SPELCASTER Work Package 1b decentralised validation framework.

The classifier had not previously been evaluated outside its original development cohort. This study assesses:

- Diagnostic concordance between classifier predictions and gold-standard histopathological diagnosis
- The relationship between sample quality (bisulfite conversion efficiency, probe failure rate) and classifier yield/confidence
- A reproducible EPICv2 methylation array processing pipeline for classifier-ready data

**Note:** This project validated and applied an existing, pre-trained classifier (Jurmeister et al., 2022); model development was not part of the scope.

---

## Dataset

- **Samples:** 38 archival FFPE sinonasal tumour samples, RCSI pathology archive (Beaumont Hospital)
- **Platform:** Illumina MethylationEPIC v2 (EPICv2) array
- **Ethics:** Approved by Beaumont Hospital Ethics Committee, 15 Nov 2024 (Ref: 24/53), TRANSCAN2023-1858-077
- **Reference:** Jurmeister et al. (2022) development/reference cohort (n=453), used for joint dimensionality reduction

---

## Workflow Summary

```
Raw IDATs (EPICv2) → SeSAMe QCDPB preprocessing → EPICv2-to-EPICv1 liftover (mLiftOver)
  → QC stratification (GCT score, fraction NA) → SVM classification + confidence calibration
  → CpG sensitivity analysis → HiTIMED tumour purity deconvolution
  → UMAP/t-SNE + CNV profiling → Differential methylation & pathway enrichment
  → Classifier-histopathology concordance analysis
```

---

## Analysis Pipeline

### 1. Preprocessing and Probe Harmonisation
- SeSAMe (v1.24.0): quality masking, Infinium I channel inference, dye bias correction, pOOBAH masking, noob background correction
- EPICv2 → EPICv1 probe-space liftover via `mLiftOver()`
- Final consensus beta matrix: 866,553 probes × 38 samples

### 2. Quality Control
- Samples stratified as PASS/FLAGGED by GCT bisulfite conversion score (>1.5) and fraction NA (>0.20)
- 20 PASS, 18 FLAGGED (8 by GCT alone, 9 by fraction NA alone, 1 auto-flagged)

### 3. SVM Classification and Confidence Calibration
- Jurmeister et al. (2022) 18-class SVM classifier, 20,000-CpG feature set
- Multinomial ridge regression calibration
- Confidence tiers per Capper et al. (2018): Strong (>0.90), Acceptable (0.84–0.90), Suggestive (0.50–0.84), Non-informative (<0.50)

### 4. CpG Sensitivity Analysis
- Numerical perturbation (±0.02) of top 300 variance-ranked CpGs per sample; finite-difference derivative (dp/dx) to identify influential features

### 5. HiTIMED Tumour Purity Deconvolution
- Hierarchical level 6, HNSC reference (closest available approximation; no sinonasal-specific reference exists)

### 6. Dimensionality Reduction and CNV Profiling
- UMAP and t-SNE, RCSI-only and jointly with Jurmeister reference cohort
- Genome-wide CNV via `cnSegmentation()` (SeSAMe)

### 7. Differential Methylation and Pathway Enrichment
- limma with empirical Bayes moderation, PASS vs. FLAGGED (n=35)
- gometh/missMethyl and clusterProfiler for probe-bias corrected GO/KEGG enrichment

### 8. Classifier-Histopathology Concordance
- Outcome framework (Jurmeister et al., 2022): Diagnosis Confirmed / Diagnosis Reclassified / Non-classifiable

---

## Key Results

- 16 of 38 samples (42.1%) received a confident class assignment; 22 returned Unknown
- **75% concordance** with histopathological diagnosis among classified samples (12 of 16), spanning 7 tumour types
- 3 samples (2 diagnosed ONB, 1 SNUC) were independently reclassified to the **NEC-like SMARCA4/ARID1A** class — replicating a key finding from the original Jurmeister et al. (2022) development cohort
- HiTIMED tumour purity: median 13%, range 7.8–87.9%, no significant association with classifier confidence (Spearman's ρ = −0.069, p = 0.681)
- No CpGs met FDR-corrected significance between FLAGGED and PASS samples; 1,793 nominally significant CpGs (69.4% hypomethylated in FLAGGED) taken forward for exploratory enrichment
- Overall conclusion: classifier yield was primarily limited by sample-level technical degradation and low tumour purity — not classifier validity. The classifier performs reliably when input quality is adequate.

---

## Repo Structure

```
scripts/
    sesame.R                     - EPICv2 preprocessing, QC, EPICv2-to-EPICv1 liftover
    classifier.R                 - SVM classification and confidence calibration
    concordance_analysis.R       - Classifier-histopathology concordance analysis
    enrichment_analysis.R        - Differential methylation and GO/KEGG enrichment
    annotation_visualisation.R   - Summary annotation and visualisation

data/
    (archival FFPE methylation array data - not included; patient data restrictions.
    Raw beta value matrices and per-sample imputed feature values are excluded from
    this repository. See Reproducibility Notes below.)

results/
    classifier/                  - Predictions, confidence scores, CpG contribution scores
    cnv/                         - Genome-wide copy number variation profiles
    dimensionality_reduction/    - UMAP/t-SNE plots and coordinates
    enrichment_analysis/         - Differential methylation results, GO/KEGG enrichment
    pathology_concordance/       - Classifier-histopathology concordance results
    qc/                          - Quality control metrics and per-sample QC plots
    tumour_purity/               - HiTIMED tumour purity estimates
```

---

## Reproducibility Notes

- All scripts used in the validation pipeline are provided in the `scripts/` directory.
- The SVM classifier itself (Jurmeister et al., 2022) is pre-trained and not part of this repository; this work covers external validation and application only.
- Raw archival methylation data are not publicly deposited due to patient data restrictions (Beaumont Hospital Ethics Ref: 24/53).
- Full methodological detail, results, and discussion are provided in the project report (not yet included in this repo).

---

## Author Details

**Dara Ó Murchú**
MSc Bioinformatics and Genomic Data Science, RCSI
SPELCASTER Work Package 1b

---

## Key Reference

Jurmeister P, Gloss S, Roller R, et al. DNA methylation-based classification of sinonasal tumors. *Nat Commun*. 2022;13(1):7148.
