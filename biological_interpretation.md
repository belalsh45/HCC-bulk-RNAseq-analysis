# Biological interpretation — TCGA-LIHC paired bulk RNA-seq

**Belal Abdellah | Exploratory interpretation, 10 matched tumor–normal pairs**

## Summary

The paired DESeq2 comparison retained 19,782 features and reported 3,782 features meeting both BH-adjusted p < 0.05 and absolute log2FC > 1: 1,574 higher and 2,208 lower in tumors. These numbers were independently recomputed from the supplied CSV; the cohort remains a small, nonrandom exploratory sample.

## Tumor-enriched proliferation programs

The Hallmark NES plot and accompanying E2F-targets running-enrichment curve show positive enrichment of **E2F TARGETS**, **G2M CHECKPOINT**, and **MITOTIC SPINDLE** in the tumor direction. The E2F figure itself displays nominal p ≈ 5.706 × 10⁻⁴⁵ and adjusted p ≈ 2.853 × 10⁻⁴³; these numbers should be cross-checked against the local GSEA results table before formal scientific reuse.

In the paired VST-difference heatmap, tumor-increased genes include **FOXM1, E2F2, CCNB2, CENPF, SKA3, ANLN** and other proliferation-related features. Together these patterns are *consistent with* transcriptional proliferation and mitotic programs in sampled tumor tissue. They do **not** establish cell-cycle activity by functional assays or prove causal drivers.

## Normal-enriched tissue, metabolism and immune signatures

The NES bar plot shows normal-direction enrichment for multiple Hallmarks, including **FATTY ACID METABOLISM, XENOBIOTIC METABOLISM, BILE ACID METABOLISM, OXIDATIVE PHOSPHORYLATION, INFLAMMATORY RESPONSE, IL6 JAK STAT3 SIGNALING, TNFA SIGNALING VIA NFKB**, and others shown in the figure. These are *relative gene-expression enrichment directions within the tumor-versus-normal contrast*, **not demonstrations of pathway inhibition in malignant hepatocytes**.

The paired heatmap also shows negative tumor-minus-normal expression differences for **HAMP, HGFAC, FCN3, CLEC4M, STAB2** and other genes connected to hepatic physiology, innate immune recognition, vascular/stromal niches, or tissue composition. Some may reflect the proportion of normal nonmalignant cell types rather than only tumor-cell-autonomous regulation. Treat interpretations of immune signatures particularly cautiously because bulk RNA-seq measures mixed cellular populations.

## Patient heterogeneity and clustering

The patient-level figure depicts **VST(tumor) − VST(normal)** per selected gene. Most patients follow broadly similar directions for the displayed features; two patient columns (TCGA-BC-A10X and TCGA-BC-A110) appear to have smaller-amplitude changes for many selected genes. Hierarchical clustering is based on **1 − Pearson correlation** with average linkage, which considers profile similarity rather than simply amplitude. Ten patients and 40 post hoc selected genes cannot establish stable tumor subtypes.

## What can and cannot be concluded

- **Supported:** Within this ten-pair dataset, paired differential-expression analysis identifies tumor-associated transcript changes; pathway-level enrichment shows strong cell-cycle signatures and opposite-direction tissue-function/immune-related signatures.
- **Not established:** Causality, validated cancer biomarkers, universal TCGA-LIHC generalizability, cell-type-specific changes, functional pathway activation or clinically meaningful prognostic associations.
- **Important caveat:** The paired heatmap displays **VST differences**, not the DESeq2 coefficient or a row-wise Z-score. The color scale is clipped near the extreme percentiles to improve visualization.

## Reasonable follow-up checks

Use all eligible matched pairs (50 in the original metadata query); inspect library QC and VST PCA before selecting markers; verify that Hallmark NES directions and multiple-testing-adjusted significance are stable. Consider gene-level annotations and cell-type composition deconvolution only if the question calls for them and requisite QC is possible. Validate interpretations in independent cohorts before drawing general biological conclusions.

## Note on verification status

The DESeq2 CSV and four figures are present in the initial portfolio package. Full Hallmark GSEA numeric results and gene-set version metadata must still be exported from the RStudio analysis for detailed numerical confirmation; figure appearance alone does not substitute for an analysis table.
