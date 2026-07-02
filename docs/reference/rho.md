# Rho (phenotypic correlation) query

Returns pairwise phenotypic correlations between sets of traits from the
sample-overlap-weighted rho matrix.

## Usage

``` r
rho(db, traits_1, traits_2)
```

## Arguments

- db:

  A `pleiodb` object from
  [`open_pleiodb()`](https://explodecomputer.github.io/pleiodbr/reference/open_pleiodb.md).

- traits_1:

  Character vector of OpenGWAS trait IDs.

- traits_2:

  Character vector of OpenGWAS trait IDs.

## Value

A tibble with columns `trait_id_1`, `trait_id_2`, `rho` covering all
cross-product pairs. Use
[`tidyr::pivot_wider()`](https://tidyr.tidyverse.org/reference/pivot_wider.html)
for matrix form.

## Examples

``` r
if (FALSE) { # \dontrun{
db <- open_pleiodb("/path/to/main.pleiodb")
rho(db, "ukb-b-19953", "ebi-a-GCST006867")
rho(db, c("ukb-b-19953", "ebi-a-GCST90018961"), c("ebi-a-GCST006867"))
} # }
```
