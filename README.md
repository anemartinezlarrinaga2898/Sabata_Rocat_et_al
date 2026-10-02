# Sabata_Rocat_et_al

## Context-dependent response of endothelial cells to PIK3CA mutation

This repository contains the R and Python code used to reproduce the single-cell transcriptomic analyses associated with the study **“Context-dependent response of endothelial cells to PIK3CA mutation”**.

The code covers the main analysis branches used for the retinal endothelial-cell dataset, including quality control, doublet estimation, Seurat preprocessing, Harmony integration, endothelial-cell selection and re-clustering, WT and MUT integration, eGFP-based lineage classification, endothelial population annotation, differential-expression analyses, downsampling, pseudobulk analysis, and RNA-velocity analysis using scVelo and CellRank.

Script numbering reflects the historical execution order within each analysis branch. Some numbering gaps are retained intentionally to preserve correspondence with the original analysis workflow.

---

## Abstract

Cancer-associated mutations in PIK3CA cause congenital disorders characterised by tissue overgrowth. The endothelium is among the most frequently affected tissues, displaying aberrant vascular overgrowth in the form of malformations. PIK3CA-driven vascular phenotypes predominantly affect veins and capillaries but rarely arteries for reasons that remain unclear. Here, using lineage tracing, we show that expression of mutated PIK3CAH1047R in endothelial cells leads to marked clonal expansions in capillary and venous endothelial populations. In contrast, mature arterial endothelial cells are resistant to physiological levels of PIK3CA signalling. Moreover, PIK3CAH1047R expression in arterial precursors interrupts arterial differentiation and induces the acquisition of venous features. Impaired arterial differentiation, therefore, provides an additional layer of protection against arterial damage in response to PIK3CA genetic perturbation. Mechanistically, this process is mediated by two complementary downstream arms of PI3K signalling: upregulation of the venous-specifying transcription factor NR2F2 (COUP-TFII), which suppresses arterial differentiation, and increased proliferation activity, which promotes vascular overgrowth. Our findings reveal that pathogenic responses to PIK3CAH1047R depend on the differentiation stage and fate trajectory of the targeted cell. Arteries are thus shielded against PIK3CA mutations, explaining the rarity of PIK3CA-associated arterial malformations in patients.

---

## 🧬 Data description

### 1. In vivo retinal endothelial cells

- **Model:** Esm1(BAC)-CreERT2; R26-mTmG ± Pik3caH1047R
- **Induction:** 4-OHT at P1–P2
- **Collection:** P6 retinas
- **Enrichment:** CD31+ magnetic sorting (Miltenyi Biotec)
- **Platform:** 10x Genomics Chromium Controller
- **Chemistry:** 3′ v4
- **Replicates:**
  - Exp1: n = 5 control, n = 3 mutant
  - Exp2: n = 4 control, n = 4 mutant

Additional processing:

- Ambient RNA correction using CellBender.
- An eGFP sequence was included in the mm10 reference used for alignment to enable lineage-tracing analyses.

### 2. In vitro lung endothelial cells

- **Model:** Pdgfb-CreERT2; Pik3caH1047R/+
- **Induction:** 4-OHT (2 μM, overnight)
- **Collection:** 48 h post induction
- **Replicates:** 3 biological replicates per condition (pooled)
- **Platform:** 10x Genomics
- **Chemistry:** 3′ v3

### Sequencing

- **Platform:** Illumina NovaSeq X
- **Reads:** 2 × 150 bp
- **Alignment:** Cell Ranger (v7–8)
- **Reference genome:** mm10, with eGFP added where required for lineage tracing

---

## ⚙️ Computational analysis

Core scRNA-seq processing and clustering were performed in **R using Seurat**. Batch correction was performed using **Harmony**. The RNA-velocity branch was performed in **Python** using **Scanpy, scVelo and CellRank**.

The in vivo analysis was organised into four connected branches:

```text
WT preprocessing ──────┐
                       ├──> WT + MUT merged EC dataset ──> downstream analyses
MUT preprocessing ─────┘                │
                                        └──> eGFP+ WT + MUT subset
                                                  │
                                                  └──> RNA velocity / CellRank
```

---

## 1. Initial quality control

Seurat objects were generated from corrected count matrices and sample metadata.

The following cell-level metrics were calculated:

- `nFeature_RNA`
- `nCount_RNA`
- mitochondrial percentage (`percent.mt`)
- ribosomal percentage (`percent.rb`)
- library complexity (`log10GenesPerUMI`)

Cells were retained using the following thresholds:

```text
nFeature_RNA >= 250
percent.mt < 20
log10GenesPerUMI > 0.8
```

WT and MUT samples were subsequently processed through separate preprocessing branches before being merged for the combined analysis.

---

## 2. Doublet estimation

Doublets were estimated independently for each sample using **DoubletFinder**.

The workflow was:

1. Split the Seurat object by sample ID.
2. Run the preprocessing workflow independently for each sample.
3. Estimate doublet classifications using DoubletFinder.
4. Add DoubletFinder results back to the metadata of the corresponding Seurat objects.
5. Merge the sample-specific objects.

Doublet classifications were retained as metadata and were **not removed during the main preprocessing workflow**.

Individual downstream analyses may explicitly restrict the dataset to singlets. In particular, the pseudobulk workflow performs this filtering before sample-level aggregation.

---

## 3. Seurat preprocessing and dimensionality reduction

The preprocessing workflow uses shared Seurat utility functions to perform the core processing steps required for analysis of the complete dataset and endothelial-cell subsets.

These steps include:

- log-normalization,
- highly variable gene selection,
- scaling,
- PCA,
- nearest-neighbor graph construction,
- UMAP generation,
- graph-based clustering.

Highly variable genes were dynamically selected based on the distribution of genes detected per cell.

Principal components used for downstream analysis were selected using a shared PCA-dimension selection function.

Initial graph-based clustering was evaluated over multiple resolutions, including:

```text
0.1, 0.3, 0.5, 1, 1.5, 2
```

Following Harmony integration, additional resolutions were commonly evaluated:

```text
0.1, 0.3, 0.5, 0.7, 0.9, 1
```

Higher resolutions were also explored during intermediate stages of the WT and MUT analyses.

Therefore, **no single clustering resolution was applied to every stage of the workflow**. Different resolutions were examined during endothelial-cell selection, refinement and final annotation.

---

## 4. Harmony integration

Batch correction was performed using **Harmony**, with sample ID (`ID`) as the integration variable.

The Harmony embedding was used to:

- construct the nearest-neighbor graph,
- generate UMAP embeddings,
- evaluate alternative clustering resolutions,
- inspect sample mixing,
- inspect phenotype distribution.

WT and MUT datasets were first analysed independently.

Their curated endothelial-cell objects were subsequently merged and reprocessed to generate the combined WT + MUT endothelial dataset.

---

## 5. Endothelial-cell selection and refinement

Endothelial populations were identified by inspection of cluster-level transcriptional profiles and canonical vascular markers, including:

- `Pecam1`
- `Cdh5`
- `Vwf`
- arterial markers
- venous markers
- capillary markers
- angiogenic and tip-cell markers
- proliferative markers

Candidate endothelial clusters were subsetted and reprocessed.

Additional non-endothelial, low-quality or ambiguous populations were removed iteratively, followed by repeated Seurat preprocessing and Harmony integration.

Because WT and MUT datasets were initially analysed independently, cluster IDs used during these filtering steps are dataset-specific and should not be interpreted as shared biological labels across branches.

---

## 6. WT + MUT combined endothelial dataset

The curated WT and MUT endothelial objects were merged and processed again as a single dataset.

The combined workflow includes:

1. Merge of curated WT and MUT endothelial objects.
2. Definition of eGFP-positive and eGFP-negative cells.
3. Reprocessing with the Seurat subset pipeline.
4. Harmony integration using sample ID.
5. Evaluation of multiple clustering resolutions.
6. Marker analysis.
7. Endothelial-population annotation.

The repository retains historical and revised annotation scripts to preserve the analytical provenance of the study.

The final annotation strategy uses information from more than one clustering resolution to refine endothelial identities rather than relying on a single clustering resolution.

Annotated endothelial populations include:

- Angiogenic EC
- Proliferative EC states
- Venous-capillary EC
- Capillary EC
- Artery-capillary EC
- Arterial EC
- Tip EC, where applicable

---

## 7. eGFP-based lineage classification

For the in vivo retinal dataset, eGFP expression was used to identify recombined lineage-traced cells.

Cells were classified using log-normalized eGFP expression:

```text
eGFP+ : eGFP expression > 1.5
eGFP− : eGFP expression <= 1.5
```

This classification was stored in the Seurat metadata and used for:

- visualization of eGFP distribution,
- comparison of endothelial population composition,
- selection of eGFP-positive cells,
- eGFP+ MUT versus eGFP+ WT differential-expression analysis,
- downsampling analyses,
- pseudobulk analyses,
- RNA-velocity analysis of the eGFP-positive subset.

---

## 8. eGFP-positive WT + MUT subset

The combined endothelial dataset was subsetted to retain eGFP-positive cells and subsequently reprocessed independently.

The eGFP-positive workflow includes:

1. Selection of eGFP-positive endothelial cells.
2. Seurat reprocessing.
3. Harmony integration by sample ID.
4. Marker estimation.
5. Endothelial population annotation.

The final eGFP-positive annotation workflow uses Harmony clustering at resolution `0.7` and groups cells into endothelial states including:

- Angiogenic EC
- Proliferative 1
- Proliferative 2
- Venous-capillary EC
- Capillary
- Artery-capillary EC
- Arterial EC

---

## 🧪 Differential-expression analysis

### Single-cell differential expression

For the in vivo dataset, differential expression was performed independently within annotated endothelial populations using Seurat `FindMarkers`.

The principal lineage-traced comparison implemented in the combined analysis is:

```text
eGFP+ MUT vs eGFP+ WT
```

Comparisons were performed independently within each `Layer_1` endothelial population.

Populations lacking sufficient cells in both groups were skipped.

Differential-expression tables contain the statistics returned by Seurat, including:

- average log2 fold change,
- P value,
- adjusted P value,
- percentage of expressing cells.

---

## Pseudobulk differential expression

A complementary sample-level pseudobulk analysis was performed using **DESeq2**.

The pseudobulk workflow:

1. Restricts the analysis to cells classified as singlets.
2. Defines eGFP-positive WT and MUT groups.
3. Aggregates raw RNA counts by endothelial population and biological sample using `AggregateExpression`.
4. Constructs a sample-level count matrix for each endothelial population.
5. Uses WT as the DESeq2 reference condition.
6. Tests MUT versus WT using a `~ condition` design.
7. Filters genes with fewer than 10 total counts before model fitting.

The pseudobulk branch also contains downstream:

- Gene Ontology enrichment analysis
- Gene Set Enrichment Analysis (GSEA)

based on the DESeq2 results.

---

## 🔁 Downsampling analysis

Repeated random downsampling was used to evaluate whether unequal numbers of WT and MUT cells could influence differential-expression results.

For each endothelial population:

1. Only eGFP-positive cells are retained.
2. WT and MUT cells are identified.
3. Both groups are randomly sampled to the size of the smaller group.
4. Differential expression is recomputed using Seurat `FindMarkers` with:

```text
logfc.threshold = 0
min.pct = 0
```

5. The procedure is repeated **1,000 times**.

Per iteration, the workflow records:

- number of significant genes,
- number of non-significant genes,
- number of upregulated genes,
- number of downregulated genes,
- number of significant upregulated genes,
- number of significant downregulated genes.

This analysis represents **repeated random downsampling without replacement**, rather than bootstrap resampling.

---

## 🧭 RNA velocity and trajectory analysis

RNA velocity was performed on the **eGFP-positive WT + MUT endothelial subset** using Python.

### Seurat-to-AnnData conversion

Count matrices, cell metadata, gene names and PCA coordinates exported from the Seurat workflow were converted into sample-specific AnnData objects.

UMAP coordinates were transferred when available.

### Addition of splicing layers

For each sample, the cleaned and annotated AnnData object was matched to a corresponding AnnData object containing RNA-velocity layers.

Cells and genes common to both objects were retained.

The following layers were transferred when available:

- `spliced`
- `unspliced`
- `ambiguous`

The sample-specific AnnData objects were subsequently combined into a single AnnData object for velocity analysis.

### scVelo preprocessing

The combined AnnData object was processed using Scanpy and scVelo.

Initial preprocessing included:

```text
Minimum genes per cell: 200
Minimum cells per gene: 10
Library normalization target: 10,000
PCA components calculated: 50
Neighbors: 15
PCs used for neighbor graph: 30
```

WT and MUT cells were then separated and RNA-velocity dynamics were estimated independently for each condition.

For each condition, the scVelo workflow included:

```text
filter_and_normalize(
    min_shared_counts = 10,
    n_top_genes = 2000
)

recover_dynamics()
velocity(mode = "dynamical")
velocity_graph()
latent_time()
velocity_pseudotime()
```

The processed WT and MUT velocity objects were saved independently.

### CellRank

CellRank was applied separately to WT and MUT using the **VelocityKernel**.

Velocity-derived transition matrices were calculated and stored together with the scVelo-derived velocity information.

The resulting objects were used to visualize:

- velocity streamlines,
- latent time,
- velocity pseudotime.

---

## 📊 Marker analysis and visualization

Marker-based characterization across the repository includes:

- `FindMarkers`
- `FindAllMarkers`
- FeaturePlots
- violin plots
- DotPlots
- average-expression heatmaps
- log2FC heatmaps
- differential-expression dot plots
- endothelial population proportion plots

Marker sets were used to support biological interpretation of:

- arterial states,
- venous states,
- capillary states,
- artery-capillary states,
- angiogenic and tip-cell states,
- proliferative endothelial states.

---

## 📁 Repository structure

The analysis is organised according to the biological and computational progression of the study:

```text
Sabata_Rocat_et_al/
│
├── InVivo_Retina/
│   │
│   ├── 01_WT/
│   │   └── scripts/
│   │
│   ├── 02_MUT/
│   │   └── scripts/
│   │
│   ├── 03_WT_MUT_Merged/
│   │   └── scripts/
│   │
│   └── 04_eGFPpos_WT_MUT/
│       └── scripts/
│
├── InVitro_Lung/
│   └── scripts/
│
├── utils/
│   └── shared helper functions
│
└── README.md
```

Within each analysis branch, numbered scripts reflect the historical order of execution.

---

## 🔬 Analysis workflow

The main in vivo workflow can be summarized as:

```text
Raw / corrected count matrices
          │
          ▼
   Seurat object creation
          │
          ▼
      Quality control
          │
          ▼
    Doublet estimation
          │
          ├───────────────┐
          ▼               ▼
     WT analysis       MUT analysis
          │               │
          ▼               ▼
   EC identification   EC identification
          │               │
          ▼               ▼
   EC refinement       EC refinement
          │               │
          └───────┬───────┘
                  ▼
        WT + MUT EC merge
                  │
                  ▼
           Harmony integration
                  │
                  ▼
         EC population annotation
                  │
          ┌───────┴─────────┐
          ▼                 ▼
  Differential          eGFP+ subset
   expression                │
          │                  ▼
          │             Re-clustering
          │                  │
          │                  ▼
          │              Annotation
          │                  │
          │                  ▼
          │             RNA velocity
          │                  │
          │                  ▼
          │               CellRank
          │
          ├── Downsampling
          │
          └── Pseudobulk / DESeq2
```

---

## Reproducibility notes

The repository preserves the analytical history of the study. Consequently, some scripts represent intermediate or alternative annotation and filtering strategies that were evaluated during analysis.

Important points include:

- WT and MUT were processed independently before construction of the merged WT + MUT dataset.
- Multiple clustering resolutions were evaluated at different stages of the workflow.
- There is no universal clustering resolution of `0.3` for the complete analysis.
- Earlier and revised combined-annotation scripts are retained for provenance.
- Some downstream scripts originate from different historical annotation branches (`Annotations` and `New_Annotations`).
- Doublet predictions were retained during the main preprocessing workflow.
- The pseudobulk branch explicitly restricts the analysis to cells classified as singlets.
- RNA-velocity analysis requires AnnData objects containing `spliced`, `unspliced` and `ambiguous` layers.
- Shared utility functions sourced by the R scripts are required for full reproduction of the Seurat preprocessing workflow and PCA-dimension selection.

---

## 💻 Main software

### R

The analysis uses packages including:

- Seurat
- Harmony
- DoubletFinder
- UCell
- DESeq2
- clusterProfiler
- ComplexHeatmap
- tidyverse

### Python

The RNA-velocity workflow uses packages including:

- Scanpy
- AnnData
- scVelo
- CellRank
- pandas
- NumPy
- SciPy
- matplotlib

Exact software versions are reported in the accompanying software documentation and Reporting Summary associated with the manuscript.

---

## 📥 Data availability

The scRNA-seq datasets generated for this study were deposited in the Gene Expression Omnibus (GEO):

- **In vivo retinal scRNA-seq:** `GSE326941`
- **In vitro lung endothelial scRNA-seq:** `GSE287780`

The repository contains the analysis code and does not replace the primary sequencing-data records deposited in GEO.

---

## 💾 Code availability

All original code used for the analyses described in the manuscript is deposited in this repository.

No new algorithms were developed specifically for this study.

Repository:

`https://github.com/anemartinezlarrinaga2898/Sabata_Rocat_et_al.git`

---

## 📖 Citation

If you use this code or the associated datasets, please cite the corresponding article:

**Sabata-Rocat et al. _Context-dependent response of endothelial cells to PIK3CA mutation_.**
