# Host Transcriptomic Analysis (Dual RNA-seq)

## Project overview
This directory contains the **host-side transcriptomic analysis** of a dual RNA-seq study
investigating the host cellular response to infection by *Helicobacter pylori*.

The analysis focuses on how host gene expression is modulated:
- by **infection status**,
- by **bacterial genotype** (WT vs KO),
- and by **time post-infection (12h vs 24h)**.

The host and bacterial analyses were conducted independently but using a
**consistent analytical framework**, enabling integrated interpretation.

---

## Experimental design (host)

- **Organism**: *Homo sapiens* (human host cells)
- **Conditions**:
  - Control (non-infected)
  - Infected
- **Bacterial strains involved**:
  - WT: wild-type *H. pylori*
  - KO: *H. pylori* knockout strain
- **Time points**:
  - 12 hours post-infection
  - 24 hours post-infection
- **Replicates**: biological replicates per condition
- **Data type**: raw read counts (used for DESeq2)

---

## Analysis strategy

Differential expression analyses were performed using **DESeq2** on **raw integer counts**.
Simple and biologically justified models were chosen to avoid overfitting.

### Comparisons performed

- Infected vs Control (global host response)
- WT vs KO infection at 12h
- WT vs KO infection at 24h
- Temporal effects:
  - WT: 24h vs 12h
  - KO: 24h vs 12h

These comparisons allow disentangling:
- the **host response to infection**,
- the **impact of bacterial genotype on the host**,
- and the **temporal evolution of host transcriptional responses**.

---

## Directory structure
```text
Host/
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
│   │   └── key_genes/    # Robust, condition-dependent host genes
│   │
│   └── figures/
│       ├── pca_*.png     # Global PCA (QC)
│       └── volcano_*.png # Differential expression volcano plots
│
└── README.md

---

## Identification of condition-dependent host genes

To focus on biologically meaningful host responses, **key genes** were extracted using
the following criteria:

- adjusted p-value (`padj`) < 0.05  
- |log2 fold change| ≥ 1  
- `baseMean` ≥ 50  

This approach prioritizes genes that are:
- robustly expressed,
- and specifically regulated according to infection conditions.

Key gene tables are stored in:
results/tables/key_genes/

---

## Key findings (summary)

- Host cells exhibit a **clear transcriptional response to infection**.
- The magnitude and structure of this response vary depending on the
  **bacterial genotype (WT vs KO)**.
- Temporal analyses indicate that **host transcriptional responses evolve between 12h and 24h**,
  with differences in intensity depending on the infecting bacterial strain.

These results highlight a **dynamic host response shaped by both time and bacterial genetic factors**.

---

## Relation to the bacterial analysis

This host analysis complements the **bacterial-side transcriptomic analysis** located in
the `Bacteria/` directory. Together, both components form a complete and coherent
**dual RNA-seq workflow**, enabling independent yet integrated interpretation of host–pathogen interactions.

---

## Notes on reproducibility

- All analyses are script-based and fully reproducible.
- No unnecessary model complexity was introduced.
- The same methodological standards were applied across host and bacterial analyses.

