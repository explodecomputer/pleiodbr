# Top hits query

Returns significant variant-trait associations below a p-value threshold
for a specified set of traits.

## Usage

``` r
tophits(db, traits, pval = 5e-08)
```

## Arguments

- db:

  A `pleiodb` object from
  [`open_pleiodb()`](https://explodecomputer.github.io/pleiodbr/reference/open_pleiodb.md).

- traits:

  Character vector of OpenGWAS trait IDs. Required.

- pval:

  P-value threshold (default `5e-8`). The database stores pre-built
  masks for `5e-8` and `1e-5`; other thresholds trigger a full-column
  scan.

## Value

A tibble with columns `variant_id`, `trait_id`, `z`, `beta`, `se`,
`pval`, `eaf`, `n`, `imputed`.

## Examples

``` r
if (FALSE) { # \dontrun{
db <- open_pleiodb("/path/to/main.pleiodb")
tophits(db, traits = "ukb-b-19953")
tophits(db, traits = c("ukb-b-19953", "ebi-a-GCST006867"), pval = 1e-5)
} # }
```
