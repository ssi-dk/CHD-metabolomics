# Metabolomics Preprocessing Pipeline

Preprocessing scripts for untargeted metabolomics data from the matched case-control CHD study (Nastou, et al, 2026). The pipeline takes raw feature tables from MZmine and produces a normalized, batch-corrected feature table ready for downstream statistical analysis.

## Scripts

Run in order:

| Step | Script | Language | Description |
|------|--------|----------|-------------|
| 1 | `1_CreateDummyFeatureTable.R` | R | Creates an anonymized feature table (zeroed intensities, generic sample names) for public deposition e.g. on GNPS |
| 2 | `2_filter_IIMN_adducts.py` | Python | Filters redundant adduct ions identified by Ion Identity Molecular Networking (IIMN), retaining the most abundant adduct per metabolite |
| 3 | `3_BlankRemoval_Imputation_Wave.R` | R | Blank removal, missing value filtering (75%/25% thresholds), RSD-based QC filtering, LOD imputation, and WaveICA2.0 batch correction |
| 4 | `4_StatisticalAnalysis.R` | R | Conditional logistic regression (per metabolite) for CHD case/control association, CHD subtype interaction analysis, multiple testing correction (FDR), and export of results tables |

## Dependencies

**R packages:** `tidyverse`, `ggsci`, `ggpubr`, `cowplot`, `vegan`, `WaveICA2.0`, `survival`, `mice`, `gt`, `writexl`

**Python packages:** `pandas`, `networkx`

## Usage

Each script has a short **Paths** section at the top — set your input and output paths there before running. No command-line arguments are needed.

## Data

Raw data is not included in this repository. The scripts expect:
- A feature table CSV exported from MZmine
- A metadata CSV with sample information (sample type, batch, run order, QC flags)
- An ion identity network edge table from IIMN (for step 2)