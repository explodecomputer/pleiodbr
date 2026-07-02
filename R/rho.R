#' Rho (phenotypic correlation) query
#'
#' Returns pairwise phenotypic correlations between sets of traits from the
#' sample-overlap-weighted rho matrix.
#'
#' @param db A `pleiodb` object from [open_pleiodb()].
#' @param traits_1 Character vector of OpenGWAS trait IDs.
#' @param traits_2 Character vector of OpenGWAS trait IDs.
#' @return A tibble with columns `trait_id_1`, `trait_id_2`, `rho` covering
#'   all cross-product pairs. Use `tidyr::pivot_wider()` for matrix form.
#' @export
#' @examples
#' \dontrun{
#' db <- open_pleiodb("/path/to/main.pleiodb")
#' rho(db, "ukb-b-19953", "ebi-a-GCST006867")
#' rho(db, c("ukb-b-19953", "ebi-a-GCST90018961"), c("ebi-a-GCST006867"))
#' }
rho <- function(db, traits_1, traits_2) {
  stopifnot(inherits(db, "pleiodb"))
  stopifnot(is.character(traits_1), is.character(traits_2))

  rho_path <- file.path(db$path, "rho.bin")
  if (!file.exists(rho_path))
    stop("No rho matrix found in this database (rho.bin missing).")

  t1_idx <- match(traits_1, db$traits$trait_id)
  missing1 <- traits_1[is.na(t1_idx)]
  if (length(missing1) > 0L)
    stop("Traits not found in database: ", paste(missing1, collapse = ", "))

  t2_idx <- match(traits_2, db$traits$trait_id)
  missing2 <- traits_2[is.na(t2_idx)]
  if (length(missing2) > 0L)
    stop("Traits not found in database: ", paste(missing2, collapse = ", "))

  t1_idx <- t1_idx - 1L
  t2_idx <- t2_idx - 1L

  pairs <- expand.grid(t1 = t1_idx, t2 = t2_idx, KEEP.OUT.ATTRS = FALSE)

  # rho matrix is T×T using same chunk layout but v-dim = T, t-dim = T
  # Temporarily override T for the rho matrix read
  rho_db <- db
  rho_db$V <- db$T
  rho_db$n_v_chunks <- db$n_t_chunks

  rho_vals <- vapply(seq_len(nrow(pairs)), function(i) {
    .get_block(rho_db, "rho",
               pairs$t1[i], pairs$t1[i] + 1L,
               pairs$t2[i], pairs$t2[i] + 1L)[1L, 1L]
  }, numeric(1L))

  tibble::tibble(
    trait_id_1 = db$traits$trait_id[pairs$t1 + 1L],
    trait_id_2 = db$traits$trait_id[pairs$t2 + 1L],
    rho        = rho_vals
  )
}
