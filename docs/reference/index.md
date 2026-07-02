# Package index

## Open a database

- [`open_pleiodb()`](https://explodecomputer.github.io/pleiodbr/reference/open_pleiodb.md)
  : Open a pleiodb database

## Queries

All return a tibble with columns: variant_id, trait_id, z, beta, se,
pval, eaf, n, imputed.

- [`phewas()`](https://explodecomputer.github.io/pleiodbr/reference/phewas.md)
  : PheWAS query
- [`gwas()`](https://explodecomputer.github.io/pleiodbr/reference/gwas.md)
  : GWAS query
- [`tophits()`](https://explodecomputer.github.io/pleiodbr/reference/tophits.md)
  : Top hits query
- [`associations()`](https://explodecomputer.github.io/pleiodbr/reference/associations.md)
  : Associations block query

## Phenotypic correlation

- [`rho()`](https://explodecomputer.github.io/pleiodbr/reference/rho.md)
  : Rho (phenotypic correlation) query

## Helpers

Utilities for formatting outputs and building plots.

- [`parse_alid()`](https://explodecomputer.github.io/pleiodbr/reference/parse_alid.md)
  : Parse ALID strings into components
- [`to_twosamplemr()`](https://explodecomputer.github.io/pleiodbr/reference/to_twosamplemr.md)
  : Convert pleiodbr results to TwoSampleMR format
- [`manhattan_plot()`](https://explodecomputer.github.io/pleiodbr/reference/manhattan_plot.md)
  : Manhattan plot
