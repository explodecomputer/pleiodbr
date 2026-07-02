#' Parse ALID strings into components
#'
#' Splits one or more ALID strings (`"chrom:pos_EA_OA"`) into their four
#' components. Useful for building LD clumping inputs or annotation tables.
#'
#' @param alid Character vector of ALID strings.
#' @return A tibble with columns `alid`, `chrom`, `pos` (integer), `ea`
#'   (effect allele), `oa` (other allele).
#' @export
#' @examples
#' parse_alid(c("16:53800954_C_T", "19:45412079_C_T"))
parse_alid <- function(alid) {
  stopifnot(is.character(alid))
  parts <- strsplit(alid, "[:_]")
  tibble::tibble(
    alid  = alid,
    chrom = vapply(parts, `[[`, "", 1L),
    pos   = as.integer(vapply(parts, `[[`, "", 2L)),
    ea    = vapply(parts, `[[`, "", 3L),
    oa    = vapply(parts, `[[`, "", 4L)
  )
}

#' Convert pleiodbr results to TwoSampleMR format
#'
#' Reformats a tibble from [phewas()], [gwas()], [tophits()], or
#' [associations()] into the format expected by [TwoSampleMR::mr()].
#'
#' @param dat A tibble returned by a pleiodbr query.
#' @param type `"exposure"` or `"outcome"`. Controls which column suffix set is
#'   produced (`beta.exposure`/`se.exposure`/... vs `beta.outcome`/...).
#' @param trait_name Optional character scalar. If supplied, used as the
#'   `exposure` or `outcome` label; otherwise `trait_id` is used.
#' @return A data frame in TwoSampleMR format, ready to pass to
#'   [TwoSampleMR::harmonise_data()].
#' @export
#' @examples
#' \dontrun{
#' db  <- open_pleiodb("/path/to/main.pleiodb")
#' exp <- gwas(db, "ukb-b-19953") |> to_twosamplemr("exposure")
#' out <- gwas(db, "ebi-a-GCST006867") |> to_twosamplemr("outcome")
#' dat <- TwoSampleMR::harmonise_data(exp, out)
#' TwoSampleMR::mr(dat)
#' }
to_twosamplemr <- function(dat, type = c("exposure", "outcome"),
                           trait_name = NULL) {
  type <- match.arg(type)
  stopifnot(
    is.data.frame(dat),
    all(c("variant_id", "trait_id", "beta", "se", "pval", "eaf") %in% names(dat))
  )

  parsed <- parse_alid(dat$variant_id)

  label <- if (!is.null(trait_name)) trait_name else dat$trait_id

  if (type == "exposure") {
    out <- data.frame(
      SNP                         = dat$variant_id,
      beta.exposure               = dat$beta,
      se.exposure                 = dat$se,
      pval.exposure               = dat$pval,
      eaf.exposure                = dat$eaf,
      effect_allele.exposure      = parsed$ea,
      other_allele.exposure       = parsed$oa,
      exposure                    = label,
      id.exposure                 = dat$trait_id,
      mr_keep.exposure            = !is.na(dat$beta) & !is.na(dat$se),
      stringsAsFactors            = FALSE
    )
  } else {
    out <- data.frame(
      SNP                        = dat$variant_id,
      beta.outcome               = dat$beta,
      se.outcome                 = dat$se,
      pval.outcome               = dat$pval,
      eaf.outcome                = dat$eaf,
      effect_allele.outcome      = parsed$ea,
      other_allele.outcome       = parsed$oa,
      outcome                    = label,
      id.outcome                 = dat$trait_id,
      mr_keep.outcome            = !is.na(dat$beta) & !is.na(dat$se),
      stringsAsFactors           = FALSE
    )
  }
  out
}

#' Manhattan plot
#'
#' Produces a genome-wide Manhattan plot from the output of [gwas()] or a
#' filtered [phewas()] call.
#'
#' @param dat A tibble with at least `variant_id`, `pval`, and optionally
#'   `imputed` (logical) columns.
#' @param threshold Genome-wide significance threshold (default `5e-8`).
#'   Drawn as a dashed horizontal line.
#' @param suggestive Suggestive significance threshold (default `1e-5`).
#'   Drawn as a dotted line. Set to `NULL` to omit.
#' @param highlight_imputed Logical. Colour imputed variants in a distinct
#'   shade (default `TRUE`). Has no effect if `dat` has no `imputed` column.
#' @param title Plot title string (default `NULL` = no title).
#' @return A `ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' db  <- open_pleiodb("/path/to/main.pleiodb")
#' res <- gwas(db, "ukb-b-19953")
#' manhattan_plot(res, title = "Body mass index")
#' }
manhattan_plot <- function(dat, threshold = 5e-8, suggestive = 1e-5,
                           highlight_imputed = TRUE, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for manhattan_plot(). ",
         "Install it with: install.packages('ggplot2')")

  stopifnot(
    is.data.frame(dat),
    all(c("variant_id", "pval") %in% names(dat))
  )

  parsed <- parse_alid(dat$variant_id)
  dat$chrom_int <- suppressWarnings(as.integer(parsed$chrom))
  dat$pos       <- parsed$pos
  dat$log10p    <- -log10(dat$pval)

  # Drop unparseable chroms (X/Y/MT)
  dat <- dat[!is.na(dat$chrom_int), ]
  dat <- dat[order(dat$chrom_int, dat$pos), ]

  # Build cumulative genome position
  max_per_chrom <- tapply(dat$pos, dat$chrom_int, max)
  chrom_offset  <- c(0L, cumsum(as.numeric(max_per_chrom[-length(max_per_chrom)])))
  names(chrom_offset) <- names(max_per_chrom)
  dat$gpos <- dat$pos + chrom_offset[as.character(dat$chrom_int)]

  # Chromosome label positions
  chrom_mids <- tapply(dat$gpos, dat$chrom_int, mean)

  # Colour: alternating blues, with imputed in orange
  dat$colour_group <- factor(dat$chrom_int %% 2)
  if (highlight_imputed && "imputed" %in% names(dat)) {
    dat$colour_group <- ifelse(dat$imputed, "imputed", as.character(dat$colour_group))
    colours <- c("0" = "#3A7DC9", "1" = "#7EB6E8", "imputed" = "#E07B39")
  } else {
    colours <- c("0" = "#3A7DC9", "1" = "#7EB6E8")
  }

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = gpos, y = log10p,
                                         colour = colour_group)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.7, shape = 16) +
    ggplot2::geom_hline(yintercept = -log10(threshold),
                        linetype = "dashed", colour = "grey30", linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = colours, guide = "none") +
    ggplot2::scale_x_continuous(
      breaks = chrom_mids,
      labels = names(chrom_mids)
    ) +
    ggplot2::labs(
      x     = "Chromosome",
      y     = expression(-log[10](p)),
      title = title
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor   = ggplot2::element_blank())

  if (!is.null(suggestive)) {
    p <- p + ggplot2::geom_hline(yintercept = -log10(suggestive),
                                 linetype = "dotted", colour = "grey50",
                                 linewidth = 0.4)
  }
  p
}
