#' Open a pleiodb database
#'
#' Reads metadata from a `.pleiodb` directory and returns a connection object.
#' All query functions (`phewas`, `gwas`, `tophits`, `associations`, `rho`)
#' accept this object as their first argument.
#'
#' @param path Path to the `.pleiodb` directory.
#' @return An S3 object of class `"pleiodb"`.
#' @export
#' @examples
#' \dontrun{
#' db <- open_pleiodb("/path/to/main.pleiodb")
#' print(db)
#' }
open_pleiodb <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)

  meta_path <- file.path(path, "meta.json")
  if (!file.exists(meta_path))
    stop("Not a pleiodb directory (meta.json not found): ", path)

  meta <- jsonlite::fromJSON(meta_path)

  supported_version <- 3L
  fv <- as.integer(meta$format_version %||% 1L)
  if (fv > supported_version)
    stop(
      "This database uses pleiodb format version ", fv,
      " but pleiodbr only supports up to version ", supported_version,
      ". Please update the pleiodbr package."
    )

  # Load variants
  variants <- utils::read.table(
    file.path(path, "variants.tsv"), header = TRUE, sep = "\t",
    stringsAsFactors = FALSE, quote = "", comment.char = ""
  )

  # Load traits
  traits <- utils::read.table(
    file.path(path, "traits.tsv"), header = TRUE, sep = "\t",
    stringsAsFactors = FALSE, quote = "", comment.char = ""
  )

  V  <- as.integer(meta$V)
  T_ <- as.integer(meta$T)
  CS <- as.integer(meta$chunk_shape)   # [CV, CT]
  CV <- CS[1L]; CT <- CS[2L]

  # Parse chrom and pos from ALID string ("chrom:pos_REF_ALT")
  alid_split     <- strsplit(variants$alid, ":", fixed = TRUE)
  variants$chrom <- vapply(alid_split, `[[`, "", 1L)
  pos_allele     <- vapply(alid_split, `[[`, "", 2L)
  variants$pos   <- as.integer(sub("_.*", "", pos_allele))

  # Pre-load imputed COO and encode as integer keys for fast per-query lookup.
  # Encoding: key = v_idx * T + t_idx (fits in double exactly since V*T < 2^53).
  # This replaces the per-query paste()+%in% approach over ~2M pairs.
  imp_path <- file.path(path, "imputed.coo.zst")
  if (file.exists(imp_path)) {
    raw <- readBin(imp_path, "raw", n = file.size(imp_path))
    dec <- zstdlite::zstd_decompress(raw)
    imp_coo <- matrix(
      readBin(dec, "integer", n = length(dec) %/% 4L,
              size = 4L, endian = "little"),
      ncol = 2L, byrow = TRUE
    )
    imp_keys <- as.double(imp_coo[, 1L]) * T_ + as.double(imp_coo[, 2L])
  } else {
    imp_keys <- double(0L)
  }

  db <- structure(
    list(
      path            = path,
      V               = V,
      T               = T_,
      CV              = CV,
      CT              = CT,
      n_v_chunks      = ceiling(V / CV),
      n_t_chunks      = ceiling(T_ / CT),
      format_version  = fv,
      pval_thresholds = as.numeric(meta$pval_thresholds %||% c(5e-8, 1e-5)),
      variants        = variants,
      traits          = traits,
      eaf             = as.numeric(variants$eaf),
      imp_keys        = imp_keys,
      cidx_cache      = new.env(hash = TRUE, parent = emptyenv())
    ),
    class = "pleiodb"
  )

  db
}

#' @export
print.pleiodb <- function(x, ...) {
  cat(
    "pleiodb database\n",
    "  path:           ", x$path, "\n",
    "  format version: ", x$format_version, "\n",
    "  variants (V):   ", format(x$V, big.mark = ","), "\n",
    "  traits (T):     ", format(x$T, big.mark = ","), "\n",
    "  chunk shape:    ", x$CV, "×", x$CT, "\n",
    sep = ""
  )
  invisible(x)
}

# Null-coalescing operator (avoids importing rlang for one helper)
`%||%` <- function(a, b) if (!is.null(a)) a else b
