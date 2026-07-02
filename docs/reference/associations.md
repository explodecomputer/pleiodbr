# Associations block query

Returns association statistics for an arbitrary set of variants ×
traits.

## Usage

``` r
associations(db, variants, traits)
```

## Arguments

- db:

  A `pleiodb` object from
  [`open_pleiodb()`](https://explodecomputer.github.io/pleiodbr/reference/open_pleiodb.md).

- variants:

  Character vector of ALID strings.

- traits:

  Character vector of OpenGWAS trait IDs.

## Value

A tibble with columns `variant_id`, `trait_id`, `z`, `beta`, `se`,
`pval`, `eaf`, `n`, `imputed`. Rows with `z = NA` are dropped.

## Examples

``` r
if (FALSE) { # \dontrun{
db <- open_pleiodb("/path/to/main.pleiodb")
associations(db,
  variants = c("16:53800954_C_T", "19:45412079_C_T"),
  traits   = c("ukb-b-19953", "ebi-a-GCST006867"))
} # }
```
