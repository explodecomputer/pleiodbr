# pleiodbr

R reader for [pleiodb](https://github.com/explodecomputer/pleiodb) — a binary
database format for storing GWAS z-scores across thousands of traits at
genome-wide scale.

`pleiodbr` reads `.pleiodb` databases directly from R with no Python
dependency. It supports PheWAS, GWAS, top-hits, arbitrary variant×trait block
queries, and phenotypic correlation (rho) queries. All queries return tidy
tibbles.


## Installation

```r
# install.packages("remotes")
remotes::install_github("explodecomputer/pleiodbr")
```

The `zstdlite` dependency is on R-universe (not yet on CRAN). If you hit a
`zstdlite` installation error, install it first:

```r
install.packages("zstdlite", repos = c("https://coolbutuseless.r-universe.dev",
                                        "https://cloud.r-project.org"))
```

**Dependencies**: `zstdlite`, `tibble`, `dplyr`, `jsonlite`.

## Quick start

```r
library(pleiodbr)

# Open a database
db <- open_pleiodb("/path/to/main.pleiodb")
print(db)
#> pleiodb database
#>   path:           /path/to/main.pleiodb
#>   format version: 3
#>   variants (V):   95,378
#>   traits (T):     4,159
#>   chunk shape:    512×512

# PheWAS — one variant across all traits
phewas(db, "16:53800954_C_T")      # exact ALID
phewas(db, "16:53e6-54e6")         # genomic region

# GWAS — all variants for one trait
gwas(db, "ukb-b-19953")            # Body mass index

# Top hits
tophits(db, traits = c("ukb-b-19953", "ebi-a-GCST006867"), pval = 5e-8)

# Arbitrary variant × trait block
associations(
  db,
  variants = c("16:53800954_C_T", "19:45412079_C_T"),
  traits   = c("ukb-b-19953", "ebi-a-GCST006867")
)

# Phenotypic correlation
rho(db, "ukb-b-19953", "ebi-a-GCST006867")
```

## Output columns

All query functions except `rho()` return a tibble with:

| Column | Type | Description |
|--------|------|-------------|
| `variant_id` | chr | ALID (`chrom:pos_REF_ALT`) |
| `trait_id` | chr | OpenGWAS trait ID |
| `z` | dbl | Z-score |
| `beta` | dbl | Effect estimate (reconstructed from z and Neff) |
| `se` | dbl | Standard error |
| `pval` | dbl | Two-sided p-value |
| `eaf` | dbl | Effect allele frequency |
| `n` | dbl | Effective sample size |
| `imputed` | lgl | Whether z-score was imputed by reference completion |

`rho()` returns a long-format tibble with columns `trait_id_1`, `trait_id_2`,
`rho`.

## What is a pleiodb database?

A `.pleiodb` directory stores a V×T matrix of GWAS z-scores (variants ×
traits) as zstd-compressed 512×512 int16 chunks. It also stores effective
sample sizes (Neff), pairwise phenotypic correlations (rho), and optional
imputed z-scores produced by LD reference completion. The Python `pleiodb`
package builds these databases from OpenGWAS VCF files.

## Vignettes

- `vignette("pleiodbr")` — full walkthrough with Manhattan plots, PheWAS
  plots, and a Mendelian randomisation example.
- `vignette("benchmarks")` — wall-clock timings for every query type and
  a guide to what drives performance.

### Building vignettes locally

Vignettes require a live `.pleiodb` database.  The easiest way to set up
the environment is via [mamba](https://mamba.readthedocs.io/):

```bash
mamba env create -f environment.yml
mamba run -n pleiodbr Rscript setup.R
mamba run -n pleiodbr Rscript -e "devtools::build_vignettes()"
```
