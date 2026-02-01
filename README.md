# Impact of proteogenomic evidence on clinical success

Code and data to reproduce analyses from:

> **Impact of proteogenomic evidence on clinical success**  
> Karim MA*, Hukku A, Ariano B, Holzinger E, Tsepilov Y, Hayhurst J, Buniello A, McDonagh EM, Castel SE, Nelson MR, Maranville J, Yerges-Armstrong L, Ghoussaini M*  
> *Nature Genetics* (2025)

## Key Findings

- pQTL-supported target-indication pairs show **4.7× higher** probability of clinical success (Phase I → Launch)
- This exceeds the 2.6× enrichment from human genetic evidence lacking pQTL support
- pQTL-based enrichment is prominent in druggable protein families (enzymes, kinases) that show limited enrichment from genetic evidence alone

## Interactive Browser

Browse all MR results (FDR < 0.05): **https://mk31.shinyapps.io/pqtl_mr_fdr05/**

## Repository Structure

```
mrcoloc/
├── R/                          # Helper functions
│   ├── advancement_rr.R        # Relative success calculations
│   ├── genes_tbl.R             # Gene annotation utilities
│   ├── pipeline_best.R         # T-I pair processing
│   └── triangulate_ot.R        # Open Targets triangulation
├── scripts/                    # Analysis scripts
│   ├── generate_mrcoloc_supplement.R   # Supplementary Tables ST1-ST17
│   ├── flowchart_numbers.R             # Flowchart statistics
│   ├── mrcoloc_paper_2025_main_figures.R    # Figure 1a-c
│   └── mrcoloc_paper_2025_supp_figures.R    # Figures S2-S3
├── data/                       # Reference data files
│   ├── gene_lists/             # Protein family annotations
│   ├── areas.tsv               # Therapeutic area mappings
│   └── indic.tsv               # Indication annotations
├── figures/                    # Publication figures
└── output/                     # Generated supplementary tables
```

## Data Requirements

### Files included in this repository
- Gene family lists (kinases, enzymes, etc.)
- Therapeutic area mappings
- Indication annotations

### Files to download (see data/README.md)
1. **Minikel et al. therapeutic index** - From [Nature 2024 supplement](https://doi.org/10.1038/s41586-024-07316-0)
2. **pQTL summary statistics** - From [GWAS Catalog](https://www.ebi.ac.uk/gwas/downloads/summary-statistics)
3. **MR-coloc results** - Available from [link pending]

## Reproducing the Analysis

### Prerequisites

```r
# Install required packages
install.packages(c(
  "tidyverse", "data.table", "openxlsx", "DescTools",
  "ggplot2", "UpSetR", "DiagrammeR", "googlesheets4"
))

# Bioconductor packages
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db", "biomaRt"))
```

### Running the analysis

```r
# Set project root
Sys.setenv(PQTL_ENRICH_ROOT = "/path/to/mrcoloc")

# Generate supplementary tables (ST1-ST17)
source("scripts/generate_mrcoloc_supplement.R")

# Generate main figures
source("scripts/mrcoloc_paper_2025_main_figures.R")

# Generate supplementary figures
source("scripts/mrcoloc_paper_2025_supp_figures.R")
```

## Outputs

### Main Figures
- **Figure 1a**: Forest plot showing relative success by genetic evidence source
- **Figure 1b**: UpSet plot of MR-coloc overlap for launched T-I pairs
- **Figure 1c**: Gene family enrichment comparing L2G vs L2G+pQTL

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
  journal={Nature Genetics},
  year={2025},
  doi={pending}
}
```

## Related Resources

- **Open Targets Genetics**: https://genetics.opentargets.org/
- **Minikel et al. (2024)**: https://doi.org/10.1038/s41586-024-07316-0
- **Previous preprint**: https://doi.org/10.1101/2023.06.01.23290252

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) for details.

## Contact

- Mohd Anisul Karim (mohd@variantbio.com)
- Maya Ghoussaini (mayaghoussainy@hotmail.com)