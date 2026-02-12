# Impact of proteogenomic evidence on clinical success

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18451758.svg)](https://doi.org/10.5281/zenodo.18451758)

Code and data to reproduce analyses from:

> **Impact of proteogenomic evidence on clinical success**  
> Karim MA*, Hukku A, Ariano B, Holzinger E, Tsepilov Y, Hayhurst J, Buniello A, McDonagh EM, Castel SE, Nelson MR, Maranville J, Yerges-Armstrong L, Ghoussaini M*  
> *Manuscript in preparation*

## Key Findings

- pQTL-supported target-indication pairs show **4.7× higher** probability of clinical success (Phase I → Launch)
- This exceeds the 2.6× enrichment from human genetic evidence lacking pQTL support
- pQTL-based enrichment is prominent in druggable protein families (enzymes, kinases) that show limited enrichment from genetic evidence alone

## Interactive Browser

Browse all MR results (FDR < 0.05): **https://mk31.shinyapps.io/pqtl_mr_fdr05/**

## Quick Start
```bash
# Clone the repository
git clone https://github.com/mohdkarim/mrcoloc.git
cd mrcoloc

# Download all data (~1 GB)
Rscript scripts/download_data.R

# Create derived datasets
Rscript scripts/create_derived_data.R

# Generate main figures (Figure 1a-c)
Rscript scripts/mrcoloc_paper_2025_main_figures.R

# Generate supplementary figures (Figures S2-S4)
Rscript scripts/mrcoloc_paper_2025_supp_figures.R

# Generate supplementary tables (ST1-ST17)
Rscript scripts/generate_mrcoloc_supplement.R
```

Or run the full setup in one command:
```bash
Rscript scripts/setup.R
```

## Data Availability

All data is automatically downloaded by `scripts/download_data.R`:

| Source | Files | Size | Description |
|--------|-------|------|-------------|
| [Zenodo](https://doi.org/10.5281/zenodo.18451758) | 5 | ~922 MB | MR-coloc results, ChEMBL annotations |
| [Minikel et al. GitHub](https://github.com/ericminikel/genetic_support) | 30 | ~120 MB | Therapeutic index, drug phase data |
| [Gene lists](https://github.com/ericminikel/genetic_support) | 9 | ~100 KB | Protein family annotations |

### Data Sources
- **Zenodo DOI**: [10.5281/zenodo.18451758](https://doi.org/10.5281/zenodo.18451758)
- **Minikel et al. Nature 2024**: [10.1038/s41586-024-07316-0](https://doi.org/10.1038/s41586-024-07316-0)

## Repository Structure
```
mrcoloc/
├── scripts/                    # Analysis scripts
│   ├── download_data.R                 # Downloads all required data
│   ├── create_derived_data.R           # Creates derived datasets
│   ├── setup.R                         # Full setup (download + derive)
│   ├── mrcoloc_paper_2025_main_figures.R    # Figure 1a-c
│   ├── mrcoloc_paper_2025_supp_figures.R    # Figures S2-S4
│   ├── generate_mrcoloc_supplement.R   # Supplementary Tables ST1-ST17
│   └── flowchart_numbers.R             # Flowchart statistics
├── R/                          # Helper functions
│   ├── advancement_rr.R        # Relative success calculations
│   ├── genes_tbl.R             # Gene annotation utilities
│   ├── pipeline_best.R         # T-I pair processing
│   └── triangulate_ot.R        # Open Targets triangulation
├── data/                       # Reference data (downloaded)
│   ├── minikel/                # Core Minikel files (merge2, assoc)
│   ├── gene_lists/             # Protein family annotations
│   └── *.tsv                   # Therapeutic area mappings, etc.
├── data_raw/                   # Analysis data (downloaded from Zenodo)
├── figures/                    # Generated figures
└── output/                     # Generated supplementary tables
```

## Prerequisites
```r
# Install required packages
install.packages(c(
  "tidyverse", "data.table", "openxlsx", "DescTools",
  "ggplot2", "UpSetR", "googlesheets4", "ckbplotr"
))

# Bioconductor packages
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db"))
```

## Outputs

### Main Figures
- **Figure 1a**: Forest plot showing relative success by genetic evidence source
- **Figure 1b**: UpSet plot of MR-coloc overlap for launched T-I pairs
- **Figure 1c**: Gene family enrichment comparing L2G vs L2G+pQTL

### Supplementary Figures
- **Figure S2**: Relative success by therapeutic area
- **Figure S3**: Relative success by phase transition
- **Figure S4**: pQTL enrichment by background universe (sensitivity analysis)

### Supplementary Tables
| Table | Description |
|-------|-------------|
| ST1 | Outcome GWAS datasets |
| ST2 | Proteomic GWAS datasets |
| ST3 | Excluded traits |
| ST4 | Figure 1a enrichment data |
| ST5 | Figure 1b UpSet plot data |
| ST6 | Figure 1c gene family enrichment |
| ST7 | Chi-square distribution of pQTL across TAs |
| ST8-9 | Therapeutic area enrichment (Universe 1 & 2) |
| ST10 | Spearman correlation summary |
| ST11-12 | Per-TA relative success (Universe 1 & 2) |
| ST13 | Breslow-Day heterogeneity test |
| ST14-15 | Leave-one-out sensitivity (Universe 1 & 2) |
| ST16 | All Bonferroni-significant MR target-trait pairs |
| ST17 | Successful pQTL-supported T-I pairs |

## Citation
```bibtex
@article{karim2025proteogenomic,
  title={Impact of proteogenomic evidence on clinical success},
  author={Karim, Mohd Anisul and Hukku, Abhay and Ariano, Bruno and others},
  journal={Manuscript in preparation},
  year={2025},
  doi={}
}
```

## Related Resources

- **Open Targets Genetics**: https://genetics.opentargets.org/
- **Minikel et al. (2024)**: https://doi.org/10.1038/s41586-024-07316-0

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) for details.

## Contact

- Mohd Anisul Karim (mohd@variantbio.com)
- Maya Ghoussaini (mayaghoussainy@hotmail.com)
