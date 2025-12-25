# Bacterial Transcriptomic Analysis (Dual RNA-seq)

## Project overview
This directory contains the **bacterial-side analysis** of a dual RNA-seq study investigating
how a **bacterial gene knockout (KO)** affects bacterial transcriptional dynamics during infection,
in comparison to the **wild-type (WT)** strain, across **two time points (12h and 24h)**.

The goal is to assess:
- the **impact of the bacterial genotype (WT vs KO)** on bacterial gene expression,
- the **temporal transcriptional dynamics** within each bacterial strain,
- and whether the knockout alters the **magnitude or structure of the bacterial response over time**.

---

## Experimental design (bacteria)

- **Organism**: Organism: *Helicobacter pylori* (Gram-negative gastric pathogen)
- **Strains**:
  - WT: wild-type bacterial strain
  - KO: bacterial strain with a targeted gene knockout
- **Time points**:
  - 12 hours post-infection
  - 24 hours post-infection
- **Replicates**: 3 biological replicates per condition
- **Data type**: raw read counts (used for DESeq2)

---

## Analysis strategy

All analyses were performed using **raw integer counts** and **DESeq2**, following a
minimal and justified modeling strategy to avoid overfitting.
A global PCA was performed as a quality control step prior to differential expression analyses.

### Comparisons performed

#### Genotype effect (within time point)
- WT vs KO at 12h
- WT vs KO at 24h

#### Time effect (within strain)
- WT: 24h vs 12h
- KO: 24h vs 12h

This design allows disentangling:
- the **direct effect of the bacterial knockout**, and
- its influence on **temporal transcriptional adaptation**.

---

## Directory structure
Bacteria/
├── data/
│   ├── raw/          # Original expression tables
│   ├── processed/    # Count matrices used for DESeq2
│   └── metadata/     # Sample metadata
│
├── scripts/
│   └── *.R           # Reproducible analysis scripts
│
├── results/
│   ├── tables/
│   │   ├── DESeq2/       # Full DESeq2 results
│   │   └── key_genes/    # Robust, condition-dependent bacterial genes
│   │
│   └── figures/
│       ├── pca_*.png     # Global PCA (QC)
│       └── volcano_*.png # Differential expression volcano plots
│
└── README.md

---

## Identification of condition-dependent bacterial genes

To focus on biologically meaningful signals, a subset of **key genes** was extracted using
the following criteria:

- adjusted p-value (`padj`) < 0.05  
- |log2 fold change| ≥ 1  
- `baseMean` ≥ 50  

These genes are:
- sufficiently expressed (robust signal),
- and differentially regulated across conditions.

They are stored in:
results/tables/key_genes/

---

## Key findings (summary)

- The **bacterial knockout strain** shows a **reduced transcriptional response** compared to WT.
- In the WT strain, bacterial gene expression undergoes a **strong time-dependent remodeling**
  between 12h and 24h.
- This temporal dynamic is **attenuated in the KO strain**, suggesting that the deleted gene
  contributes to bacterial transcriptional adaptation during infection.
- A global PCA performed after count matrix preparation shows a strong separation driven by time post-infection (PC1),
with an additional genotype effect (WT vs KO) captured on PC2, confirming both data quality and biological relevance.

---

## Notes on reproducibility

- All analyses rely on **explicit scripts** (no manual steps).
- No additional packages beyond those already used in the project were introduced.
- The same analytical framework was applied across all bacterial comparisons to ensure
  methodological consistency.

---

## Relation to the host analysis

This bacterial analysis complements the **host-side transcriptomic analysis** performed in
the `Host/` directory. Together, both components form a complete **dual RNA-seq workflow**
allowing independent yet integrated interpretation of host and bacterial responses.
