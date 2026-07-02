# Read a rectangular block from a chunked matrix file.

Read a rectangular block from a chunked matrix file.

## Usage

``` r
.get_block(db, name, v_start, v_end, t_start, t_end)
```

## Arguments

- db:

  pleiodb S3 object

- name:

  matrix name without extension ("zscore", "neff", "rho")

- v_start, v_end:

  row range \[v_start, v_end)

- t_start, t_end:

  column range \[t_start, t_end)

## Value

numeric matrix (v_end-v_start) × (t_end-t_start), decoded
