# Manhattan plot

Produces a genome-wide Manhattan plot from the output of
[`gwas()`](https://explodecomputer.github.io/pleiodbr/reference/gwas.md)
or a filtered
[`phewas()`](https://explodecomputer.github.io/pleiodbr/reference/phewas.md)
call.

## Usage

``` r
manhattan_plot(
  dat,
  threshold = 5e-08,
  suggestive = 1e-05,
  highlight_imputed = TRUE,
  title = NULL
)
```

## Arguments

- dat:

  A tibble with at least `variant_id`, `pval`, and optionally `imputed`
  (logical) columns.

- threshold:

  Genome-wide significance threshold (default `5e-8`). Drawn as a dashed
  horizontal line.

- suggestive:

  Suggestive significance threshold (default `1e-5`). Drawn as a dotted
  line. Set to `NULL` to omit.

- highlight_imputed:

  Logical. Colour imputed variants in a distinct shade (default `TRUE`).
  Has no effect if `dat` has no `imputed` column.

- title:

  Plot title string (default `NULL` = no title).

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
db  <- open_pleiodb("/path/to/main.pleiodb")
res <- gwas(db, "ukb-b-19953")
manhattan_plot(res, title = "Body mass index")
} # }
```
