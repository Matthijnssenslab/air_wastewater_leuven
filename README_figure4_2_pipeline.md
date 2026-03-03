# Figure 4 Part 2: Influenza A Resistance Analysis

`scripts/figure4_2_pipeline.R` performs consensus sequence generation from filtered BAM files, resistance mutation calling, and generates coverage plots with mutation annotations for four key influenza A genes (HA, NA, M2, PA).

The R script is designed to skip BAM file processing if intermediate files already exist — these are included in this repository, so you can reproduce the figures without needing the original BAM files.

## Dependencies

**R packages:**
```r
install.packages(c("dplyr", "ggplot2", "cowplot", "readr", "jsonlite"))
```

**Python packages:**
```bash
pip install pandas biopython
```

**System tools:** `samtools`, `bcftools` (must be available on PATH or configured in the script).

## Key parameters

- Mapping quality filter: > 30 (primary alignments only)
- Minimum depth for consensus: 100×
- Minimum allele frequency for consensus: 90%
- Positions not meeting these thresholds are masked as 'N'

## Usage

Set your working directory to the repository root in R, then:

```r
source("scripts/figure4_2_pipeline.R")
```

The script skips consensus generation if intermediate files already exist in `resistance_analysis/`, so subsequent runs are faster.