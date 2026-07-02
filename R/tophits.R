#' Top hits query
#'
#' Returns significant variant-trait associations below a p-value threshold
#' for a specified set of traits.
#'
#' @param db A `pleiodb` object from [open_pleiodb()].
#' @param traits Character vector of OpenGWAS trait IDs. Required.
#' @param pval P-value threshold (default `5e-8`). The database stores
#'   pre-built masks for `5e-8` and `1e-5`. Any threshold that is more
#'   stringent than an available mask (e.g. `1e-6` uses the `1e-5` mask with
#'   a subsequent filter). Thresholds less stringent than all masks (e.g.
#'   `1e-4`) trigger a full-column scan.
#' @return A tibble with columns `variant_id`, `trait_id`, `z`, `beta`, `se`,
#'   `pval`, `eaf`, `n`, `imputed`.
#' @export
#' @examples
#' \dontrun{
#' db <- open_pleiodb("/path/to/main.pleiodb")
#' tophits(db, traits = "ukb-b-19953")
#' tophits(db, traits = c("ukb-b-19953", "ebi-a-GCST006867"), pval = 1e-5)
#' }
tophits <- function(db, traits, pval = 5e-8) {
  stopifnot(inherits(db, "pleiodb"))
  if (missing(traits))
    stop("`traits` is required. Pass a character vector of trait IDs.")
  stopifnot(is.character(traits))

  t_idx <- match(traits, db$traits$trait_id)
  missing_t <- traits[is.na(t_idx)]
  if (length(missing_t) > 0L)
    stop("Traits not found in database: ", paste(missing_t, collapse = ", "))
  t_idx <- t_idx - 1L   # 0-based

  coo <- .load_sig_coo(db, pval)

  if (!is.null(coo)) {
    keep <- coo[, 2L] %in% t_idx
    coo  <- coo[keep, , drop = FALSE]
    if (nrow(coo) == 0L)
      return(.empty_tibble())

    z_vals <- .fetch_for_pairs(db, "zscore", coo[, 1L], coo[, 2L])

    keep2 <- !is.na(z_vals) & abs(z_vals) >= qnorm(pval / 2, lower.tail = FALSE)
    if (!any(keep2)) return(.empty_tibble())

    .build_tibble(coo[keep2, 1L], coo[keep2, 2L], z_vals[keep2], db)
  } else {
    z_thresh <- qnorm(pval / 2, lower.tail = FALSE)
    rows <- vector("list", length(t_idx))
    for (i in seq_along(t_idx)) {
      ti    <- t_idx[i]
      z_col <- .get_block(db, "zscore", 0L, db$V, ti, ti + 1L)[, 1L]
      keep  <- which(!is.na(z_col) & abs(z_col) >= z_thresh)
      if (length(keep) == 0L) next
      rows[[i]] <- .build_tibble(keep - 1L, rep(ti, length(keep)),
                                 z_col[keep], db)
    }
    dplyr::bind_rows(rows)
  }
}

.empty_tibble <- function() {
  tibble::tibble(
    variant_id = character(), trait_id = character(),
    z = numeric(), beta = numeric(), se = numeric(),
    pval = numeric(), eaf = numeric(), n = numeric(), imputed = logical()
  )
}
