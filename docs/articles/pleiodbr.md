# Getting started with pleiodbr

`pleiodbr` reads `.pleiodb` databases — compact binary archives of GWAS
z-scores across thousands of traits — directly from R.

## Opening a database

``` r

library(pleiodbr)
library(ggplot2)
library(dplyr)

db <- open_pleiodb("/path/to/main.pleiodb")
db
```

    pleiodb database
      path:           /path/to/main.pleiodb
      format version: 3
      variants (V):   95,378
      traits (T):     4,159
      chunk shape:    512×512

The `db` object is a lightweight connection: it holds the variant and
trait tables in memory but reads z-scores on demand from the binary
chunks.

------------------------------------------------------------------------

## PheWAS — one variant across all traits

[`phewas()`](https://explodecomputer.github.io/pleiodbr/reference/phewas.md)
fetches z-scores for a single variant across every trait. You can query
by exact ALID (`chrom:pos_REF_ALT`) or by genomic region.

``` r

# rs429358 region (APOE ε4 proxy) — affects LDL, Alzheimer's risk, and more
pw <- phewas(db, "19:45412079_C_T")
pw
```

### PheWAS plot

``` r

# Annotate with trait names from db$traits
pw <- pw |>
  left_join(db$traits[, c("trait_id", "trait_name")], by = "trait_id") |>
  mutate(
    log10p   = -log10(pval),
    highlight = pval < 5e-8
  ) |>
  arrange(log10p)

# Assign a numeric x position per trait (sorted by p-value for clarity)
pw$x <- seq_len(nrow(pw))

ggplot(pw, aes(x = x, y = log10p, colour = highlight)) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"),
                      guide = "none") +
  ggrepel::geom_text_repel(
    data = filter(pw, pval < 1e-20),
    aes(label = trait_name),
    size = 2.5, max.overlaps = 12
  ) +
  labs(
    x = "Trait (ranked by p-value)",
    y = expression(-log[10](p)),
    title = "PheWAS: 19:45412079_C_T"
  ) +
  theme_bw(base_size = 11)
```

------------------------------------------------------------------------

## GWAS — all variants for one trait

[`gwas()`](https://explodecomputer.github.io/pleiodbr/reference/gwas.md)
scans an entire trait column, returning all non-missing associations.

``` r

bmi <- gwas(db, "ukb-b-19953")   # Body mass index (BMI), N ≈ 461,460
bmi
```

### Manhattan plot

``` r

# Parse chromosome and position from ALID
bmi <- bmi |>
  mutate(
    chrom = as.integer(sub(":.*", "", variant_id)),
    pos   = as.numeric(sub(".*:(\\d+)_.*", "\\1", variant_id)),
    log10p = -log10(pval)
  ) |>
  filter(!is.na(chrom)) |>
  arrange(chrom, pos)

# Build x-axis offset per chromosome
chrom_offsets <- bmi |>
  group_by(chrom) |>
  summarise(max_pos = max(pos), .groups = "drop") |>
  arrange(chrom) |>
  mutate(offset = lag(cumsum(as.numeric(max_pos)), default = 0))

bmi <- bmi |>
  left_join(chrom_offsets[, c("chrom", "offset")], by = "chrom") |>
  mutate(gpos = pos + offset)

# Chromosome midpoints for x-axis labels
chrom_labels <- bmi |>
  group_by(chrom) |>
  summarise(mid = mean(gpos), .groups = "drop")

ggplot(bmi, aes(x = gpos, y = log10p,
                colour = factor(chrom %% 2),
                shape  = ifelse(imputed, "imputed", "observed"))) +
  geom_point(size = 0.6, alpha = 0.7) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = c("0" = "#3A7DC9", "1" = "#B0C4DE"), guide = "none") +
  scale_shape_manual(values = c("observed" = 16, "imputed" = 17),
                     name = NULL) +
  scale_x_continuous(breaks = chrom_labels$mid, labels = chrom_labels$chrom) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](p)),
    title = "GWAS Manhattan — Body mass index (BMI)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
```

------------------------------------------------------------------------

## Top hits

[`tophits()`](https://explodecomputer.github.io/pleiodbr/reference/tophits.md)
returns genome-wide-significant associations for a set of traits, using
a pre-built index for speed.

``` r

hits <- tophits(
  db,
  traits = c("ukb-b-19953",       # BMI
             "ebi-a-GCST006867"), # Type 2 diabetes
  pval   = 5e-8
)

hits |>
  left_join(db$traits[, c("trait_id","trait_name")], by = "trait_id") |>
  count(trait_name, sort = TRUE)
```

------------------------------------------------------------------------

## Associations — arbitrary variant × trait block

[`associations()`](https://explodecomputer.github.io/pleiodbr/reference/associations.md)
fetches a specific set of variant–trait pairs. This is the workhorse for
Mendelian randomisation.

``` r

# A handful of LDL instruments (near APOE, PCSK9, LDLR)
ldl_instruments <- c(
  "19:45412079_C_T",   # APOE region
  "1:55496039_G_T",    # PCSK9
  "19:11200038_C_T",   # LDLR region
  "2:21231524_T_C"     # APOB region
)

assoc <- associations(
  db,
  variants = ldl_instruments,
  traits   = c("ebi-a-GCST90018961",  # LDL cholesterol (exposure)
               "ebi-a-GCST003116")    # Coronary artery disease (outcome)
)
assoc
```

------------------------------------------------------------------------

## Mendelian randomisation: LDL → CAD

We use the associations above to run IVW and MR-Egger with
[TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/).

``` r

library(TwoSampleMR)

# Pivot into exposure/outcome data frames
exposure <- assoc |>
  filter(trait_id == "ebi-a-GCST90018961") |>
  transmute(
    SNP    = variant_id,
    beta.exposure = beta,
    se.exposure   = se,
    pval.exposure = pval,
    eaf.exposure  = eaf,
    effect_allele.exposure  = sub(".*_([A-Z])_[A-Z]$", "\\1", variant_id),
    other_allele.exposure   = sub(".*_[A-Z]_([A-Z])$", "\\1", variant_id),
    exposure = "LDL cholesterol",
    id.exposure = "ebi-a-GCST90018961",
    mr_keep.exposure = TRUE
  )

outcome <- assoc |>
  filter(trait_id == "ebi-a-GCST003116") |>
  transmute(
    SNP    = variant_id,
    beta.outcome = beta,
    se.outcome   = se,
    pval.outcome = pval,
    eaf.outcome  = eaf,
    effect_allele.outcome = sub(".*_([A-Z])_[A-Z]$", "\\1", variant_id),
    other_allele.outcome  = sub(".*_[A-Z]_([A-Z])$", "\\1", variant_id),
    outcome = "Coronary artery disease",
    id.outcome = "ebi-a-GCST003116",
    mr_keep.outcome = TRUE
  )

dat <- harmonise_data(exposure, outcome, action = 2)
```

``` r

res <- mr(dat)
res[, c("method", "nsnp", "b", "se", "pval")]
```

### MR scatter plot

``` r

ggplot(dat, aes(x = beta.exposure, y = beta.outcome)) +
  geom_point(
    aes(colour = ifelse(mr_keep, "keep", "excluded")),
    size = 3
  ) +
  geom_errorbar(aes(ymin = beta.outcome - se.outcome,
                    ymax = beta.outcome + se.outcome),
                width = 0, alpha = 0.5) +
  geom_errorbarh(aes(xmin = beta.exposure - se.exposure,
                     xmax = beta.exposure + se.exposure),
                 height = 0, alpha = 0.5) +
  # IVW slope
  geom_abline(
    slope     = res$b[res$method == "Inverse variance weighted"],
    intercept = 0,
    colour    = "steelblue", linewidth = 1
  ) +
  # MR-Egger slope and intercept
  geom_abline(
    slope     = res$b[res$method == "MR Egger"],
    intercept = mr_pleiotropy_test(dat)$egger_intercept,
    colour    = "firebrick", linetype = "dashed", linewidth = 1
  ) +
  scale_colour_manual(
    values = c("keep" = "black", "excluded" = "grey70"),
    guide  = "none"
  ) +
  ggrepel::geom_text_repel(
    aes(label = SNP), size = 2.5, max.overlaps = 10
  ) +
  labs(
    x     = "β LDL cholesterol (per SD)",
    y     = "β Coronary artery disease (log OR)",
    title = "MR: LDL cholesterol → Coronary artery disease",
    caption = "Blue = IVW; Red dashed = MR-Egger"
  ) +
  theme_bw(base_size = 11)
```

------------------------------------------------------------------------

## Phenotypic correlation (rho)

[`rho()`](https://explodecomputer.github.io/pleiodbr/reference/rho.md)
retrieves pairwise phenotypic correlations estimated from sample overlap
in the rho matrix.

``` r

# Correlation between BMI and Type 2 diabetes
rho(db, "ukb-b-19953", "ebi-a-GCST006867")
```

### Correlation heatmap for a trait cluster

``` r

library(tidyr)

metabolic_traits <- c(
  "ukb-b-19953",       # BMI
  "ebi-a-GCST006867",  # Type 2 diabetes
  "ebi-a-GCST90018961",# LDL cholesterol
  "ebi-a-GCST90002412",# LDL (larger study)
  "ebi-a-GCST003116"   # Coronary artery disease
)

rho_mat <- rho(db, metabolic_traits, metabolic_traits)

# Label by trait name
trait_labels <- db$traits |>
  filter(trait_id %in% metabolic_traits) |>
  select(trait_id, trait_name)

rho_mat <- rho_mat |>
  left_join(trait_labels, by = c("trait_id_1" = "trait_id")) |>
  rename(name1 = trait_name) |>
  left_join(trait_labels, by = c("trait_id_2" = "trait_id")) |>
  rename(name2 = trait_name)

ggplot(rho_mat, aes(x = name1, y = name2, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = round(rho, 2)), size = 3) +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                       midpoint = 0, limits = c(-1, 1), name = "ρ") +
  labs(x = NULL, y = NULL, title = "Phenotypic correlations") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
```

------------------------------------------------------------------------

## Working with imputed z-scores

All query functions return an `imputed` column flagging z-scores
produced by LD reference completion rather than direct GWAS. For MR you
may want to check whether instruments are imputed:

``` r

assoc |>
  select(variant_id, trait_id, z, pval, imputed) |>
  arrange(desc(imputed))
```

Imputed z-scores have slightly wider uncertainty because they are
constrained to the observed \|z\| range of the surrounding LD block, so
treat them with the same caution as soft imputed genotypes.

------------------------------------------------------------------------

## Querying a genomic region

To pull all associations in a locus — useful for colocalisation or
regional plots — combine a region-string PheWAS with
[`associations()`](https://explodecomputer.github.io/pleiodbr/reference/associations.md):

``` r

# All variants in the FTO locus (BMI lead region on chr16)
fto_variants <- phewas(db, "16:53e6-54e6") |>
  filter(trait_id == "ukb-b-19953") |>
  arrange(pval)

# Regional Manhattan for BMI at FTO
ggplot(fto_variants, aes(
  x     = as.numeric(sub(".*:(\\d+)_.*", "\\1", variant_id)) / 1e6,
  y     = -log10(pval),
  colour = imputed,
  shape  = imputed
)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(values = c("FALSE" = "steelblue", "TRUE" = "#E07B39"),
                      labels = c("Observed", "Imputed")) +
  scale_shape_manual(values  = c("FALSE" = 16, "TRUE" = 17),
                     labels  = c("Observed", "Imputed")) +
  labs(
    x      = "Position on chr16 (Mb)",
    y      = expression(-log[10](p)),
    colour = NULL, shape = NULL,
    title  = "FTO locus — BMI (chr16:53–54 Mb)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
```
