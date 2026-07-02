# Open a pleiodb database

Reads metadata from a `.pleiodb` directory and returns a connection
object. All query functions (`phewas`, `gwas`, `tophits`,
`associations`, `rho`) accept this object as their first argument.

## Usage

``` r
open_pleiodb(path)
```

## Arguments

- path:

  Path to the `.pleiodb` directory.

## Value

An S3 object of class `"pleiodb"`.

## Examples

``` r
if (FALSE) { # \dontrun{
db <- open_pleiodb("/path/to/main.pleiodb")
print(db)
} # }
```
