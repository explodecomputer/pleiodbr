# Internal format-reading helpers — not exported.
# .pleiodb format: 512×512 zstd-compressed chunks, ti-major ordering.
# cidx is a uint64 array (length n_chunks+1); cidx[i]/cidx[i+1] = byte range.

Z_NA      <- -32768L
NEFF_NA   <- 0xFFFFL
Z_SCALE   <- 100.0
NEFF_FRAC <- 2048.0   # 2^11

# ---- chunk addressing -------------------------------------------------------

.chunk_id <- function(vi, ti, n_v_chunks) {
  ti * n_v_chunks + vi
}

.chunk_bounds <- function(vi, ti, CV, CT, V, T) {
  list(
    v0 = vi * CV,
    v1 = min((vi + 1L) * CV, V),
    t0 = ti * CT,
    t1 = min((ti + 1L) * CT, T)
  )
}

# ---- low-level I/O ----------------------------------------------------------

.read_cidx <- function(path) {
  # uint64 stored as 8-byte little-endian; read pairs as two int32 and
  # reconstruct. R does not support unsigned size-4 reads, but all practical
  # .bin files are < 2 GB so the signed values equal the unsigned values.
  raw_bytes <- readBin(path, what = "raw", n = file.size(path))
  n <- length(raw_bytes) %/% 8L
  lo <- readBin(raw_bytes, what = "integer", n = n * 2L, size = 4L,
                endian = "little")[seq(1L, n * 2L, 2L)]
  hi <- readBin(raw_bytes, what = "integer", n = n * 2L, size = 4L,
                endian = "little")[seq(2L, n * 2L, 2L)]
  # Combine as double (safe to 2^53 bytes)
  as.double(lo) + as.double(hi) * 4294967296.0
}

# .read_raw_chunk is kept for any external callers; internal code now uses the
# single-connection variants inside .get_block / .fetch_for_pairs.
.read_raw_chunk <- function(bin_path, cidx, chunk_id) {
  start <- cidx[chunk_id + 1L]
  end   <- cidx[chunk_id + 2L]
  n     <- as.integer(end - start)
  con   <- file(bin_path, open = "rb")
  on.exit(close(con))
  seek(con, where = start, origin = "start")
  readBin(con, what = "raw", n = n)
}

# Read (and cache) the chunk-index for a named matrix.
# db$cidx_cache is an environment (reference type) so the cache survives even
# when db is passed by value.
.get_cidx <- function(db, name) {
  cidx <- db$cidx_cache[[name]]
  if (is.null(cidx)) {
    cidx <- .read_cidx(file.path(db$path, paste0(name, ".cidx")))
    db$cidx_cache[[name]] <- cidx
  }
  cidx
}

.decompress <- function(raw_bytes) {
  zstdlite::zstd_decompress(raw_bytes)
}

# ---- decode -----------------------------------------------------------------

.decode_z <- function(int16_vec) {
  out <- as.double(int16_vec) / Z_SCALE
  out[int16_vec == Z_NA] <- NA_real_
  out
}

.decode_neff <- function(uint16_vec) {
  out <- 2.0 ^ (as.double(uint16_vec) / NEFF_FRAC)
  out[uint16_vec == NEFF_NA] <- NA_real_
  out
}

.decode_f16 <- function(uint16_vec) {
  # IEEE 754 half-precision → double.
  # sign=bit15, exp=bits14-10 (bias 15), mantissa=bits9-0
  u <- as.integer(uint16_vec)
  sign  <- bitwShiftR(bitwAnd(u, 0x8000L), 15L)
  exp   <- bitwShiftR(bitwAnd(u, 0x7C00L), 10L)
  mant  <- bitwAnd(u, 0x03FFL)
  val <- ifelse(
    exp == 0L,
    # subnormal
    (2.0 ^ -14) * mant / 1024.0,
    ifelse(
      exp == 31L,
      ifelse(mant == 0L, Inf, NaN),
      # normal
      2.0 ^ (exp - 15L) * (1.0 + mant / 1024.0)
    )
  )
  val * ifelse(sign == 1L, -1.0, 1.0)
}

# ---- block reader -----------------------------------------------------------

#' Read a rectangular block from a chunked matrix file.
#'
#' @param db  pleiodb S3 object
#' @param name  matrix name without extension ("zscore", "neff", "rho")
#' @param v_start,v_end  row range [v_start, v_end)
#' @param t_start,t_end  column range [t_start, t_end)
#' @return numeric matrix (v_end-v_start) × (t_end-t_start), decoded
#' @keywords internal
.get_block <- function(db, name, v_start, v_end, t_start, t_end) {
  v_end <- min(v_end, db$V)
  t_end <- min(t_end, db$T)
  CV <- db$CV; CT <- db$CT
  n_v_chunks <- db$n_v_chunks

  cidx     <- .get_cidx(db, name)
  bin_path <- file.path(db$path, paste0(name, ".bin"))

  vi_lo <- v_start %/% CV
  vi_hi <- (v_end - 1L) %/% CV
  ti_lo <- t_start %/% CT
  ti_hi <- (t_end - 1L) %/% CT

  nrows <- v_end - v_start
  ncols <- t_end - t_start
  out   <- matrix(NA_real_, nrow = nrows, ncol = ncols)

  # Open the .bin file once for all chunk reads in this block.
  fh <- file(bin_path, open = "rb")
  on.exit(close(fh))

  for (ti in ti_lo:ti_hi) {
    for (vi in vi_lo:vi_hi) {
      cid    <- .chunk_id(vi, ti, n_v_chunks)
      bounds <- .chunk_bounds(vi, ti, CV, CT, db$V, db$T)

      start <- cidx[cid + 1L]
      seek(fh, where = start, origin = "start")
      dec <- .decompress(readBin(fh, what = "raw", n = as.integer(cidx[cid + 2L] - start)))

      chunk_nrows <- bounds$v1 - bounds$v0
      chunk_ncols <- bounds$t1 - bounds$t0

      # byrow=TRUE: binary format is row-major (trait varies fastest within a
      # chunk); R's matrix() default is column-major, which would transpose.
      if (name == "zscore") {
        vals <- matrix(
          .decode_z(readBin(dec, "integer", n = chunk_nrows * chunk_ncols,
                            size = 2L, signed = TRUE, endian = "little")),
          nrow = chunk_nrows, ncol = chunk_ncols, byrow = TRUE
        )
      } else if (name == "neff") {
        vals <- matrix(
          .decode_neff(readBin(dec, "integer", n = chunk_nrows * chunk_ncols,
                               size = 2L, signed = FALSE, endian = "little")),
          nrow = chunk_nrows, ncol = chunk_ncols, byrow = TRUE
        )
      } else if (name == "rho") {
        vals <- matrix(
          .decode_f16(readBin(dec, "integer", n = chunk_nrows * chunk_ncols,
                              size = 2L, signed = FALSE, endian = "little")),
          nrow = chunk_nrows, ncol = chunk_ncols, byrow = TRUE
        )
      } else {
        stop("Unknown matrix name: ", name)
      }

      r0  <- max(v_start, bounds$v0) - v_start + 1L
      r1  <- min(v_end,   bounds$v1) - v_start
      c0  <- max(t_start, bounds$t0) - t_start + 1L
      c1  <- min(t_end,   bounds$t1) - t_start
      sr0 <- max(v_start, bounds$v0) - bounds$v0 + 1L
      sr1 <- min(v_end,   bounds$v1) - bounds$v0
      sc0 <- max(t_start, bounds$t0) - bounds$t0 + 1L
      sc1 <- min(t_end,   bounds$t1) - bounds$t0

      out[r0:r1, c0:c1] <- vals[sr0:sr1, sc0:sc1]
    }
  }
  out
}

# ---- beta/SE reconstruction -------------------------------------------------

.reconstruct_beta_se <- function(z, neff, eaf) {
  # se = 1 / sqrt(neff * 2 * eaf * (1 - eaf))
  # beta = z * se
  denom <- neff * 2.0 * eaf * (1.0 - eaf)
  se    <- ifelse(denom > 0, 1.0 / sqrt(denom), NA_real_)
  beta  <- z * se
  list(beta = beta, se = se)
}

# ---- imputed COO loader -----------------------------------------------------

.load_imputed_coo <- function(db) {
  path <- file.path(db$path, "imputed.coo.zst")
  if (!file.exists(path)) return(matrix(integer(0), ncol = 2L))
  raw  <- readBin(path, "raw", n = file.size(path))
  dec  <- zstdlite::zstd_decompress(raw)
  # Data is interleaved pairs: v0,t0,v1,t1,... so byrow=TRUE gives correct columns
  matrix(
    readBin(dec, "integer", n = length(dec) %/% 4L,
            size = 4L, endian = "little"),
    ncol = 2L, byrow = TRUE
  )  # col1=v_idx, col2=t_idx, 0-based
}

# ---- significance COO loader ------------------------------------------------

.load_sig_coo <- function(db, pval) {
  thresholds <- db$pval_thresholds
  # Use the most stringent mask whose threshold is >= pval (i.e. a superset of
  # all pval-significant pairs).  This lets any custom threshold below 1e-5 use
  # the 1e-5 mask with a subsequent filter, avoiding a full column scan.
  # If pval is less stringent than every available mask, fall back to scanning.
  candidates <- thresholds[thresholds >= pval]
  if (length(candidates) == 0L) return(NULL)
  match_thr <- min(candidates)

  fname <- sprintf("masks/%.0e.coo.zst", match_thr)
  fname <- gsub("e-0*(\\d+)", "e-0\\1", fname)  # ensure two-digit exponent
  path  <- file.path(db$path, fname)
  if (!file.exists(path)) return(NULL)
  raw <- readBin(path, "raw", n = file.size(path))
  dec <- zstdlite::zstd_decompress(raw)
  matrix(
    readBin(dec, "integer", n = length(dec) %/% 4L,
            size = 4L, endian = "little"),
    ncol = 2L, byrow = TRUE
  )  # col1=v_idx, col2=t_idx, 0-based
}

# ---- chunk-batched pair reader ----------------------------------------------
# Reads a named matrix (zscore / neff / rho) for an arbitrary set of
# (v_idx, t_idx) pairs by grouping them into their natural chunk boundaries.
# One chunk is decompressed per unique (v_chunk, t_chunk) combination, and all
# pairs in that chunk are extracted in a single matrix subscript.  This
# eliminates the O(n_pairs) file seeks that arise from calling .get_block once
# per pair.

.fetch_for_pairs <- function(db, name, v_idx, t_idx) {
  if (length(v_idx) == 0L) return(numeric(0L))
  CV <- db$CV; CT <- db$CT
  n_v_chunks <- db$n_v_chunks

  cidx     <- .get_cidx(db, name)
  bin_path <- file.path(db$path, paste0(name, ".bin"))

  vc  <- v_idx %/% CV
  tc  <- t_idx %/% CT
  cid <- tc * n_v_chunks + vc   # ti-major ordering

  vals <- numeric(length(v_idx))

  # Open the .bin file once for all chunk reads.
  fh <- file(bin_path, open = "rb")
  on.exit(close(fh))

  for (cid_val in unique(cid)) {
    sel      <- which(cid == cid_val)
    vi_chunk <- cid_val %%  n_v_chunks
    ti_chunk <- cid_val %/% n_v_chunks

    v0 <- vi_chunk * CV;  v1 <- min((vi_chunk + 1L) * CV, db$V)
    t0 <- ti_chunk * CT;  t1 <- min((ti_chunk + 1L) * CT, db$T)
    nr <- v1 - v0;        nc <- t1 - t0

    start <- cidx[cid_val + 1L]
    seek(fh, where = start, origin = "start")
    dec <- .decompress(readBin(fh, what = "raw", n = as.integer(cidx[cid_val + 2L] - start)))

    chunk_vals <- if (name == "zscore") {
      matrix(.decode_z(readBin(dec, "integer", n = nr * nc,
                               size = 2L, signed = TRUE,  endian = "little")),
             nrow = nr, ncol = nc, byrow = TRUE)
    } else if (name == "neff") {
      matrix(.decode_neff(readBin(dec, "integer", n = nr * nc,
                                  size = 2L, signed = FALSE, endian = "little")),
             nrow = nr, ncol = nc, byrow = TRUE)
    } else if (name == "rho") {
      matrix(.decode_f16(readBin(dec, "integer", n = nr * nc,
                                 size = 2L, signed = FALSE, endian = "little")),
             nrow = nr, ncol = nc, byrow = TRUE)
    } else {
      stop("Unknown matrix name: ", name)
    }

    vals[sel] <- chunk_vals[cbind(v_idx[sel] - v0 + 1L, t_idx[sel] - t0 + 1L)]
  }
  vals
}

# ---- shared tibble builder --------------------------------------------------

.build_tibble <- function(v_idx, t_idx, z_vals, db) {
  eaf  <- db$eaf[v_idx + 1L]
  neff <- .fetch_for_pairs(db, "neff", v_idx, t_idx)

  bs   <- .reconstruct_beta_se(z_vals, neff, eaf)
  pval <- 2 * pnorm(-abs(z_vals))

  # Imputed flag: db$imp_keys holds pre-encoded (v*T + t) doubles for all
  # imputed pairs, so membership is a single vectorised %in% with no string
  # allocation (replacing the old paste()+%in% over the full 1.95M-row COO).
  if (length(db$imp_keys) > 0L) {
    pair_keys <- as.double(v_idx) * db$T + as.double(t_idx)
    is_imp    <- pair_keys %in% db$imp_keys
  } else {
    is_imp <- rep(FALSE, length(v_idx))
  }

  tibble::tibble(
    variant_id = db$variants$alid[v_idx + 1L],
    trait_id   = db$traits$trait_id[t_idx + 1L],
    z          = z_vals,
    beta       = bs$beta,
    se         = bs$se,
    pval       = pval,
    eaf        = eaf,
    n          = neff,
    imputed    = is_imp
  )
}
