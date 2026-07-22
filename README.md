
# sammyR

<!-- badges: start -->
  <!-- badges: end -->
  
  R package for SAMMY-seq chromatin compartment and differential solubility analysis.

## Installation

You can install the development version of sammyR from [GitHub](https://github.com/daisymut/sammyR) with:
  
  ``` r
# install.packages("remotes")
remotes::install_github("daisymut/sammyR")
```


## Overview

sammyR provides a workflow for importing SAMMY-seq signal, binning it, and calling A/B compartments and subcompartments along the genome and perform differential solubility analysis between groups of samples.

Main functions:
  
  - import_and_rebin__bw() — import bigWig tracks and rebin to a chosen resolution
- call_subcompartments_sammy() — main call A/B (sub)compartments function
- Bins_selector() — main select and filter genomic bins for differential solubility analysis
- add_metadata() — attach sample/patient metadata
- generate_files() — export compartments results
- setup_gene_annotation() — prepare gene annotation 
- save_bins_data() — save binned data objects

## Example

This is a basic example :
  
  ``` r
library(sammyR)

patient    <- "PZ01"
chromosome <- "chr2"
binsize    <- 100000L
input_file <- "samples_fractions_bws.csv"
chrom_bed  <- "chr2_binned.bed"
gene_gtf   <- "genes.gtf"

### Load data
comp_df <- fread(input_file, data.table = FALSE)
comp_df_replica <- comp_df[comp_df$Patient_name == patient, ]

bins_gr <- import(chrom_bed, format = "BED")
genes_gr <- import(gene_gtf, format = "GTF")

sub2_colors <- c("B" = "#4575b4", "A" = "#d73027")
subs_file <- paste0(patient, "_compartment___", chromosome, '_', binsize, ".Rdata")

### Run SAMMY compartments call per chromosome
sub_objs <- call_subcompartments_sammy(
  patients = patient,
  tracks_db = comp_df_replica,
  bins_gr = bins_gr,
  subs_file = subs_file,
  binsize = binsize,
  chr = chromosome,
  genes_gr = genes_gr,
  keeping_bins1 = "all",
  sublevel = "sub.2",
  sub_colors = sub2_colors
)

### create output files
generate_files(sub_objs, chromosome)

```

## Authors

- [Margherita Mutarelli](https://github.com/daisymut) (maintainer) — margherita.mutarelli@cnr.it
- <Coautore 2> —
- <Coautore 3> —

## License

MIT © 2026
