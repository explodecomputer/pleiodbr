#' PheWAS query
#'
#' Returns association statistics for a single variant (or all variants in a
#' genomic region) across every trait in the database.
#'
#' @param db A `pleiodb` object from [open_pleiodb()].
#' @param variant Either an ALID string (`"1:103574777_C_G"`) or a region
#'   string (`"1:103e6-104e6"`, positions in base-pairs).
#' @return A tibble with columns `variant_id`, `trait_id`, `z`, `beta`, `se`,
#'   `pval`, `eaf`, `n`, `imputed`. Rows with `z = NA` are dropped.
#' @export
#' @examples
#' \dontrun{
#' db <- open_pleiodb("/path/to/main.pleiodb")
#' phewas(db, "16:53800954_C_T")
#' phewas(db, "16:53e6-54e6")
#' }
phewas <- function(db, variant) {
  stopifnot(inherits(db, "pleiodb"))
  stopifnot(is.character(variant), length(variant) == 1L)

  v_idx <- .resolve_variant(db, variant)
  if (length(v_idx) == 0L)
    stop("No variants found for: ", variant)

  rows <- vector("list", length(v_idx))

  for (i in seq_along(v_idx)) {
    vi    <- v_idx[i]
    z_row <- .get_block(db, "zscore", vi, vi + 1L, 0L, db$T)[1L, ]
    keep  <- which(!is.na(z_row))
    if (length(keep) == 0L) next

    rows[[i]] <- .build_tibble(
      v_idx  = rep(vi, length(keep)),
      t_idx  = keep - 1L,
      z_vals = z_row[keep],
      db     = db
    )
  }

  dplyr::bind_rows(rows)
}

# ---- helpers ----------------------------------------------------------------

.resolve_variant <- function(db, variant) {
  if (grepl("-", variant, fixed = TRUE)) {
    .parse_region(db, variant)
  } else {
    idx <- match(variant, db$variants$alid)
    if (is.na(idx))
      stop("Variant not found in database: ", variant)
    idx - 1L   # 0-based
  }
}

.parse_region <- function(db, region_str) {
  # Accepts "chrom:start-end" where start/end are bp (scientific notation ok)
  parts <- strsplit(region_str, ":", fixed = TRUE)[[1L]]
  if (length(parts) != 2L)
    stop("Invalid region string (expected chrom:start-end): ", region_str)
  chrom <- parts[1L]
  se    <- strsplit(parts[2L], "-", fixed = TRUE)[[1L]]
  if (length(se) != 2L)
    stop("Invalid region string (expected chrom:start-end): ", region_str)
  start <- as.numeric(se[1L])
  end   <- as.numeric(se[2L])
  if (any(is.na(c(start, end))))
    stop("Could not parse positions in region string: ", region_str)

  hits <- which(
    db$variants$chrom == chrom &
    db$variants$pos   >= start &
    db$variants$pos   <= end
  )
  if (length(hits) == 0L)
    warning("No variants found in region: ", region_str)
  hits - 1L   # 0-based
}
