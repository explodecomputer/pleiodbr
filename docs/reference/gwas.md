# GWAS query

Returns association statistics for all variants for a single trait.

## Usage

``` r
gwas(db, trait_id)
```

## Arguments

- db:

  A `pleiodb` object from
  [`open_pleiodb()`](https://explodecomputer.github.io/pleiodbr/reference/open_pleiodb.md).

- trait_id:

  OpenGWAS trait ID string.

## Value

A tibble with columns `variant_id`, `trait_id`, `z`, `beta`, `se`,
`pval`, `eaf`, `n`, `imputed`. Rows with `z = NA` are dropped.

## Examples

``` r
if (FALSE) { # \dontrun{
db  <- open_pleiodb("/path/to/main.pleiodb")
gwas(db, "ukb-b-19953")
} # }
```
