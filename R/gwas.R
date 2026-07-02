#' GWAS query
#'
#' Returns association statistics for all variants for a single trait.
#'
#' @param db A `pleiodb` object from [open_pleiodb()].
#' @param trait_id OpenGWAS trait ID string.
#' @return A tibble with columns `variant_id`, `trait_id`, `z`, `beta`, `se`,
#'   `pval`, `eaf`, `n`, `imputed`. Rows with `z = NA` are dropped.
#' @export
#' @examples
#' \dontrun{
#' db  <- open_pleiodb("/path/to/main.pleiodb")
#' gwas(db, "ukb-b-19953")
#' }
gwas <- function(db, trait_id) {
  stopifnot(inherits(db, "pleiodb"))
  stopifnot(is.character(trait_id), length(trait_id) == 1L)

  ti <- match(trait_id, db$traits$trait_id)
  if (is.na(ti))
    stop("Trait not found in database: ", trait_id)
  ti <- ti - 1L   # 0-based

  z_col <- .get_block(db, "zscore", 0L, db$V, ti, ti + 1L)[, 1L]
  keep  <- which(!is.na(z_col))
  if (length(keep) == 0L)
    return(tibble::tibble(
      variant_id = character(), trait_id = character(),
      z = numeric(), beta = numeric(), se = numeric(),
      pval = numeric(), eaf = numeric(), n = numeric(), imputed = logical()
    ))

  imp_coo <- .load_imputed_coo(db)
  .build_tibble(
    v_idx       = keep - 1L,
    t_idx       = rep(ti, length(keep)),
    z_vals      = z_col[keep],
    db          = db,
    imputed_coo = imp_coo
  )
}
