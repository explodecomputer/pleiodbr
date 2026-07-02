#' Associations block query
#'
#' Returns association statistics for an arbitrary set of variants × traits.
#'
#' @param db A `pleiodb` object from [open_pleiodb()].
#' @param variants Character vector of ALID strings.
#' @param traits Character vector of OpenGWAS trait IDs.
#' @return A tibble with columns `variant_id`, `trait_id`, `z`, `beta`, `se`,
#'   `pval`, `eaf`, `n`, `imputed`. Rows with `z = NA` are dropped.
#' @export
#' @examples
#' \dontrun{
#' db <- open_pleiodb("/path/to/main.pleiodb")
#' associations(db,
#'   variants = c("16:53800954_C_T", "19:45412079_C_T"),
#'   traits   = c("ukb-b-19953", "ebi-a-GCST006867"))
#' }
associations <- function(db, variants, traits) {
  stopifnot(inherits(db, "pleiodb"))
  stopifnot(is.character(variants), is.character(traits))

  v_idx <- match(variants, db$variants$alid)
  missing_v <- variants[is.na(v_idx)]
  if (length(missing_v) > 0L)
    stop("Variants not found in database: ", paste(missing_v, collapse = ", "))

  t_idx <- match(traits, db$traits$trait_id)
  missing_t <- traits[is.na(t_idx)]
  if (length(missing_t) > 0L)
    stop("Traits not found in database: ", paste(missing_t, collapse = ", "))

  v_idx <- v_idx - 1L   # 0-based
  t_idx <- t_idx - 1L

  pairs  <- expand.grid(vi = v_idx, ti = t_idx, KEEP.OUT.ATTRS = FALSE)
  z_vals <- .fetch_for_pairs(db, "zscore", pairs$vi, pairs$ti)

  keep <- !is.na(z_vals)
  if (!any(keep)) return(.empty_tibble())

  .build_tibble(pairs$vi[keep], pairs$ti[keep], z_vals[keep], db)
}
