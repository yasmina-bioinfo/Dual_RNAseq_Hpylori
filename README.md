# Dual RNA-seq Analysis of *Helicobacter pylori* Infection in Humans

## Project overview
This repository presents a **dual RNA-seq analysis** investigating host–pathogen interactions
between **human host cells (*Homo sapiens*)** and the bacterial pathogen
***Helicobacter pylori***.

The project independently analyzes:
- the **host transcriptional response** to infection,
- the **bacterial transcriptional dynamics** during infection,

while using a **consistent, reproducible analytical framework** on both sides.
This approach enables integrated interpretation without conflating host and bacterial signals.

---

## Biological context
*Helicobacter pylori* is a Gram-negative gastric pathogen associated with chronic inflammation
and gastric disease in humans. Understanding how both the host and the bacterium
adjust their gene expression during infection is essential to decipher host–pathogen interactions.

This project focuses on:
- the effect of **bacterial genotype** (wild-type vs knockout strain),
- the role of **time post-infection** (12h vs 24h),
- and how these factors shape transcriptional responses on both sides.

---

## Experimental design (summary)

### Host
- **Organism**: *Homo sapiens*
- **Conditions**:
  - Control (non-infected)
  - Infected
- **Bacterial strains**:
  - WT: wild-type *H. pylori*
  - KO: *H. pylori* knockout strain
- **Time points**: 12h and 24h post-infection
- **Data type**: raw RNA-seq read counts

### Bacterium
- **Organism**: *Helicobacter pylori*
- **Strains**:
  - WT (wild-type)
  - KO (gene knockout)
- **Time points**: 12h and 24h post-infection
- **Data type**: raw RNA-seq read counts

---

## Analysis strategy

All analyses were conducted using **DESeq2** on **raw integer counts**.
To ensure robustness and interpretability:

- simple and biologically justified models were used,
- no unnecessary interactions were introduced,
- the same methodological principles were applied across host and bacterial analyses.

Comparisons include:
- genotype effects (WT vs KO),
- temporal effects (24h vs 12h),
- and condition-dependent transcriptional responses.

---

## Repository structure

Dual_RNAseq_Hpylori/
├── Host/
│ ├── data/
│ ├── scripts/
│ ├── results/
│ └── README.md
│
├── Bacteria/
│ ├── data/
│ ├── scripts/
│ ├── results/
│ └── README.md
│
└── README.md

---

Each subdirectory contains:
- raw and processed data,
- fully reproducible analysis scripts,
- result tables and figures,
- and a dedicated README describing the analytical choices.

---

## Identification of condition-dependent genes

For both host and bacterial analyses, **key genes** were defined as genes that are:
- significantly differentially expressed (`padj` < 0.05),
- show substantial expression changes (|log2FC| ≥ 1),
- and are robustly expressed (`baseMean` ≥ 50).

This strategy prioritizes **biologically meaningful transcriptional changes**
over low-expression statistical artefacts.

---

## Key conclusions (high-level)

- Human host cells display a **dynamic transcriptional response** to *H. pylori* infection.
- The **bacterial genotype (WT vs KO)** influences both host and bacterial gene expression.
- The wild-type bacterium shows a **strong time-dependent transcriptional remodeling**,
  which is **attenuated in the knockout strain**.
- Together, these results highlight a **functional role of the bacterial gene knockout**
  in shaping transcriptional dynamics during infection.

---

## Reproducibility and scope

- All results are generated via explicit scripts (no manual steps).
- The pipeline is designed for **clarity, reproducibility, and methodological consistency**.
- This repository represents a **methodological and analytical mini-project**, not a full
  mechanistic study, and avoids over-interpretation.

---

## Author’s note
This project was developed as part of a **PhD preparation portfolio**, with an emphasis on:
- sound experimental reasoning,
- clean computational workflows,
- and defensible biological interpretation.
