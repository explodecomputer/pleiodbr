# Parse ALID strings into components

Splits one or more ALID strings (`"chrom:pos_EA_OA"`) into their four
components. Useful for building LD clumping inputs or annotation tables.

## Usage

``` r
parse_alid(alid)
```

## Arguments

- alid:

  Character vector of ALID strings.

## Value

A tibble with columns `alid`, `chrom`, `pos` (integer), `ea` (effect
allele), `oa` (other allele).

## Examples

``` r
parse_alid(c("16:53800954_C_T", "19:45412079_C_T"))
#> # A tibble: 2 × 5
#>   alid            chrom      pos ea    oa   
#>   <chr>           <chr>    <int> <chr> <chr>
#> 1 16:53800954_C_T 16    53800954 C     T    
#> 2 19:45412079_C_T 19    45412079 C     T    
```
