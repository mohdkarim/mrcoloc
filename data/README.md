# Data Requirements

This directory contains small reference files. Large data files must be downloaded separately.

## Files included in this repository

### Gene family lists (`gene_lists/`)
- `ab_tractable.tsv` - Antibody-tractable targets
- `sm_tractable.tsv` - Small molecule-tractable targets
- `kinases.tsv` - Kinase gene list
- `enzymes.tsv` - Enzyme gene list
- `ion_channels.tsv` - Ion channel gene list
- `nuclear_receptors.tsv` - Nuclear receptor gene list
- `catalytic_receptors.tsv` - Catalytic receptor gene list
- `rhodop_gpcr.tsv` - Rhodopsin-class GPCR gene list
- `transporters.tsv` - Transporter gene list

### Reference mappings
- `areas.tsv` - Therapeutic area definitions
- `indic.tsv` - Indication annotations with genetic insight classifications
- `indic_topl_match.tsv` - Indication to top-level therapeutic area mappings

## Files to download

### 1. Minikel et al. therapeutic index data

**Source**: Supplementary materials from [Minikel et al. Nature 2024](https://doi.org/10.1038/s41586-024-07316-0)

**Required files**:
- `merge2.tsv.gz` - Main therapeutic index with genetic support annotations

**Download instructions**:
1. Go to the Nature paper supplementary information
2. Download the supplementary data files
3. Place `merge2.tsv.gz` in `genetic_support-main/data/`

### 2. pQTL MR-coloc results

**Source**: Available from [link pending - will be deposited upon publication]

**Required files** (place in `data_raw/`):
- `pqtl_mrcoloc_2025.rds` - Filtered MR-coloc results
- `ukb_ppp_mr_coloc_results.rds` - UKB-PPP MR-coloc results
- `mr_prot_filtered_dataset_v1_v2.rds` - Legacy MR results (filtered)
- `mr_prot_unfiltered_dataset_v1_v2_without_egger_with_transcoloc.rds` - Legacy MR results (unfiltered)

### 3. Supporting data files

**Required files** (place in `data_raw/`):
- `chembl.rds` - ChEMBL drug-target annotations
- `panukb.rds` - Pan-UK Biobank trait mappings
- `trans_genes.rds` - Trans-coloc gene annotations
- `ttpairs_tested.rds` - Background target-trait pairs tested

## Directory structure after setup

```
data/
├── README.md (this file)
├── gene_lists/
│   ├── ab_tractable.tsv
│   ├── sm_tractable.tsv
│   ├── kinases.tsv
│   ├── enzymes.tsv
│   ├── ion_channels.tsv
│   ├── nuclear_receptors.tsv
│   ├── catalytic_receptors.tsv
│   ├── rhodop_gpcr.tsv
│   └── transporters.tsv
├── areas.tsv
├── indic.tsv
└── indic_topl_match.tsv

data_raw/  (create this directory, not tracked by git)
├── pqtl_mrcoloc_2025.rds
├── ukb_ppp_mr_coloc_results.rds
├── mr_prot_filtered_dataset_v1_v2.rds
├── mr_prot_unfiltered_dataset_v1_v2_without_egger_with_transcoloc.rds
├── chembl.rds
├── panukb.rds
├── trans_genes.rds
└── ttpairs_tested.rds

genetic_support-main/data/  (Minikel et al. data)
└── merge2.tsv.gz
```

## Verification

After downloading all files, you can verify your setup by running:

```r
# Check required files exist
required_files <- c(
  "data_raw/pqtl_mrcoloc_2025.rds",
  "genetic_support-main/data/merge2.tsv.gz",
  "data/areas.tsv",
  "data/indic.tsv"
)

missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  cat("Missing files:\n")
  cat(paste("-", missing, collapse = "\n"))
} else {
  cat("All required files present!")
}
```