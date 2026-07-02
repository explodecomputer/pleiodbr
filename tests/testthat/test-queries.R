skip_db <- function() {
  skip_if_not(
    dir.exists("/local-scratch/data/pleiodb/main.pleiodb"),
    "main.pleiodb not available"
  )
}

test_that("phewas returns tibble with correct columns", {
  skip_db()
  db  <- open_pleiodb("/local-scratch/data/pleiodb/main.pleiodb")
  res <- phewas(db, "16:53800954_C_T")
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("variant_id","trait_id","z","beta","se","pval","eaf","n","imputed"))
  expect_true(nrow(res) > 0L)
  expect_true(all(!is.na(res$z)))
  expect_true(is.logical(res$imputed))
})

test_that("phewas region string returns multiple variants", {
  skip_db()
  db  <- open_pleiodb("/local-scratch/data/pleiodb/main.pleiodb")
  res <- phewas(db, "16:53e6-54e6")
  expect_true(length(unique(res$variant_id)) > 1L)
})

test_that("phewas errors on unknown ALID", {
  skip_db()
  db <- open_pleiodb("/local-scratch/data/pleiodb/main.pleiodb")
  expect_error(phewas(db, "99:999999_A_T"), "not found")
})

test_that("gwas returns tibble for known trait", {
  skip_db()
  db  <- open_pleiodb("/local-scratch/data/pleiodb/main.pleiodb")
  res <- gwas(db, "ukb-b-19953")
  expect_s3_class(res, "tbl_df")
  expect_true(nrow(res) > 1000L)
  expect_named(res, c("variant_id","trait_id","z","beta","se","pval","eaf","n","imputed"))
})

test_that("tophits requires traits argument", {
  skip_db()
  db <- open_pleiodb("/local-scratch/data/pleiodb/main.pleiodb")
  expect_error(tophits(db), "required")
})

test_that("tophits returns only significant hits", {
  skip_db()
  db  <- open_pleiodb("/local-scratch/data/pleiodb/main.pleiodb")
  res <- tophits(db, traits = "ukb-b-19953", pval = 5e-8)
  expect_true(all(res$pval <= 5e-8))
  # Regression: chunk data is row-major; byrow=TRUE required to avoid silent
  # transposition.  BMI has ~2000 GWS hits in the mask — a transposed read
  # returns near-zero z-scores for almost all of them, collapsing to ~16.
  expect_gt(nrow(res), 500L)
})

test_that("associations z-scores are at correct (variant, trait) coordinates", {
  skip_db()
  db  <- open_pleiodb("/local-scratch/data/pleiodb/main.pleiodb")
  # APOE/BMI association: known GWS hit.  A transposed read returns a cell
  # from a completely different (variant, trait) position with a near-zero z.
  res <- associations(db,
    variants = "19:45412079_C_T",
    traits   = "ukb-b-19953")
  expect_equal(nrow(res), 1L)
  expect_gt(abs(res$z), 5.0)
})

test_that("associations returns cross-product rows", {
  skip_db()
  db  <- open_pleiodb("/local-scratch/data/pleiodb/main.pleiodb")
  vs  <- c("16:53800954_C_T", "19:45412079_C_T")
  ts  <- c("ukb-b-19953", "ebi-a-GCST006867")
  res <- associations(db, vs, ts)
  expect_s3_class(res, "tbl_df")
  expect_lte(nrow(res), 4L)
  expect_gte(nrow(res), 1L)
})

test_that("rho returns a long-format tibble", {
  skip_db()
  db  <- open_pleiodb("/local-scratch/data/pleiodb/main.pleiodb")
  res <- rho(db, "ukb-b-19953", "ebi-a-GCST006867")
  expect_named(res, c("trait_id_1", "trait_id_2", "rho"))
  expect_equal(nrow(res), 1L)
  expect_true(abs(res$rho) <= 1.0)
})
