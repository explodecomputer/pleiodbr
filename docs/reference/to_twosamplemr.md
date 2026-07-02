# Convert pleiodbr results to TwoSampleMR format

Reformats a tibble from
[`phewas()`](https://explodecomputer.github.io/pleiodbr/reference/phewas.md),
[`gwas()`](https://explodecomputer.github.io/pleiodbr/reference/gwas.md),
[`tophits()`](https://explodecomputer.github.io/pleiodbr/reference/tophits.md),
or
[`associations()`](https://explodecomputer.github.io/pleiodbr/reference/associations.md)
into the format expected by
[`TwoSampleMR::mr()`](https://mrcieu.github.io/TwoSampleMR/reference/mr.html).

## Usage

``` r
to_twosamplemr(dat, type = c("exposure", "outcome"), trait_name = NULL)
```

## Arguments

- dat:

  A tibble returned by a pleiodbr query.

- type:

  `"exposure"` or `"outcome"`. Controls which column suffix set is
  produced (`beta.exposure`/`se.exposure`/... vs `beta.outcome`/...).

- trait_name:

  Optional character scalar. If supplied, used as the `exposure` or
  `outcome` label; otherwise `trait_id` is used.

## Value

A data frame in TwoSampleMR format, ready to pass to
[`TwoSampleMR::harmonise_data()`](https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.html).

## Examples

``` r
if (FALSE) { # \dontrun{
db  <- open_pleiodb("/path/to/main.pleiodb")
exp <- gwas(db, "ukb-b-19953") |> to_twosamplemr("exposure")
out <- gwas(db, "ebi-a-GCST006867") |> to_twosamplemr("outcome")
dat <- TwoSampleMR::harmonise_data(exp, out)
TwoSampleMR::mr(dat)
} # }
```
