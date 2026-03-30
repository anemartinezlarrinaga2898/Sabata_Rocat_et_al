# Sabata_Rocat_et_al

**Context-dependent response of endothelial cells to PIK3CA mutation**

Welcome to this repository, where we delve into the code used to reproduce the analysis and generate all necessary R objects for the study: *"Context-dependent response of endothelial cells to PIK3CA mutation*"*.

**Abstract**

Cancer-associated mutations in PIK3CA cause congenital disorders characterised by tissue overgrowth. The endothelium is among the most frequently affected tissues, displaying aberrant vascular overgrowth in the form of malformations. PIK3CA-driven vascular phenotypes predominantly affect veins and capillaries but rarely arteries for reasons that remain unclear. Here, using lineage tracing, we show that expression of mutated PIK3CAH1047R in endothelial cells leads to marked clonal expansions in capillary and venous endothelial populations. In contrast, mature arterial endothelial cells are resistant to physiological levels of PIK3CA signalling. Moreover, PIK3CAH1047R expression in arterial precursors interrupts arterial differentiation and induces the acquisition of venous features. Impaired arterial differentiation, therefore, provides an additional layer of protection against arterial damage in response to PIK3CA genetic perturbation. Mechanistically, this process is mediated by two complementary downstream arms of PI3K signalling: upregulation of the venous-specifying transcription factor NR2F2 (COUP-TFII), which suppresses arterial differentiation, and increased proliferation activity, which promotes vascular overgrowth. Our findings reveal that pathogenic responses to PIK3CAH1047R depend on the differentiation stage and fate trajectory of the targeted cell. Arteries are thus shielded against PIK3CA mutations, explaining the rarity of PIK3CA-associated arterial malformations in patients.

## 🧬 Data description

This repository includes single-cell RNA sequencing (scRNA-seq) data from two experimental systems:

### 1. In vivo retinal endothelial cells

- Model: Esm1(BAC)-CreERT2; R26-mTmG ± Pik3caH1047R
- Induction: 4-OHT at P1–P2
- Collection: P6 retinas
- Enrichment: CD31+ magnetic sorting (Miltenyi Biotec)
- Platform: 10x Genomics (Chromium Controller)
- Chemistry: 3’ v4
- Replicates:
  - Exp1: n = 5 control, n = 3 mutant
  - Exp2: n = 4 control, n = 4 mutant

Special processing:
- Ambient RNA correction using CellBender
- GFP sequence added to mm10 reference for lineage tracing

---

### 2. In vitro lung endothelial cells

- Model: Pdgfb-CreERT2; Pik3caH1047R/+
- Induction: 4-OHT (2 μM, overnight)
- Collection: 48h post induction
- Replicates: 3 biological replicates per condition (pooled)
- Platform: 10x Genomics
- Chemistry: 3’ v3

---

### Sequencing

- Platform: Illumina NovaSeq X
- Reads: 2 × 150 bp
- Alignment: Cell Ranger (v7–8)
- Reference genome: mm10

## ⚙️ Computational analysis

All datasets (in vivo and in vitro) were processed using a unified pipeline in R with the Seurat package.

### Preprocessing

- Filtering:
  - nFeature > 250
  - percent.mt < 20%
  - complexity < 0.8
- Removal of genes with zero counts
- Doublet detection using DoubletFinder (v2.0.3)
  - Doublets were annotated but not removed

### Normalization and feature selection

- Log-normalization (scale factor = 10,000)
- Highly variable genes (HVGs):
  - Defined as 3rd quartile + 100 genes

### Dimensionality reduction

- PCA
- Batch correction using Harmony (batch variable: sample ID)
- UMAP embedding (for visualization)

### Clustering

- Graph-based clustering:
  - SNN graph (FindNeighbors)
  - Louvain algorithm (FindClusters)
- Tested resolutions: 0.1, 0.3, 0.5
- Selected resolution: 0.3

### Endothelial cell selection

- Subset based on:
  - Pecam1
  - Cdh5

- Re-analysis performed on EC subset:
  - HVG selection
  - PCA
  - Harmony
  - Clustering

### Annotation

- Marker-based annotation using FindAllMarkers
- Statistical test:
  - Wilcoxon rank-sum
  - Bonferroni correction

## 🧪 Differential expression analysis

Differential expression was performed comparing:

- GFP⁺ MUT vs GFP⁺ WT (in vivo)
- MUT vs WT (in vitro)

### Methods:

1. Single-cell approach:
   - Seurat FindMarkers (Wilcoxon test)

2. Pseudobulk approach:
   - Aggregation by sample (AggregateExpression)
   - Analysis using DESeq2

### Validation:

- Correlation of log2 fold changes
- Overlap of significant genes (adj. p < 0.05)

Wilcoxon results were used for visualization due to higher sensitivity.
Pseudobulk results were used for validation.

## 🧬 GFP-based classification (in vivo only)

- GFP sequence added to mm10 reference genome
- Expression quantified during alignment (Cell Ranger)

Cell classification:
- GFP⁺: log-normalized expression > 1.5
- GFP⁻: below threshold

This enables lineage tracing of recombined endothelial cells.

## 🔁 Downsampling analysis

To control for unequal cell numbers:

- 1,000 bootstrap iterations
- Random subsampling to equalize group sizes
- DE recomputed at each iteration

Metrics:
- Frequency of significance (adj. p < 0.05)
- Consistency of log2FC direction

- ## 📥 Data availability

Data are available at GEO:

- Accession IN VIVO: XXXXX
- Accession IN VITRO: XXXXX
