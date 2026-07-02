# PheWAS query

Returns association statistics for a single variant (or all variants in
a genomic region) across every trait in the database.

## Usage

``` r
phewas(db, variant)
```

## Arguments

- db:

  A `pleiodb` object from
  [`open_pleiodb()`](https://explodecomputer.github.io/pleiodbr/reference/open_pleiodb.md).

- variant:

  Either an ALID string (`"1:103574777_C_G"`) or a region string
  (`"1:103e6-104e6"`, positions in base-pairs).

## Value

A tibble with columns `variant_id`, `trait_id`, `z`, `beta`, `se`,
`pval`, `eaf`, `n`, `imputed`. Rows with `z = NA` are dropped.

## Examples

``` r
if (FALSE) { # \dontrun{
db <- open_pleiodb("/path/to/main.pleiodb")
phewas(db, "16:53800954_C_T")
phewas(db, "16:53e6-54e6")
} # }
```
