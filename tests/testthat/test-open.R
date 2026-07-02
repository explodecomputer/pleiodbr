test_that("open_pleiodb returns correct class and dimensions", {
  skip_if_not(
    dir.exists("/local-scratch/data/pleiodb/main.pleiodb"),
    "main.pleiodb not available"
  )
  db <- open_pleiodb("/local-scratch/data/pleiodb/main.pleiodb")
  expect_s3_class(db, "pleiodb")
  expect_equal(db$V, 95378L)
  expect_equal(db$T, 4159L)
  expect_equal(db$CV, 512L)
  expect_equal(db$CT, 512L)
})

test_that("open_pleiodb errors on future format version", {
  tmp <- withr::local_tempdir()
  writeLines('{"V":10,"T":10,"chunk_shape":[512,512],"format_version":99,
               "pval_thresholds":[5e-8]}', file.path(tmp, "meta.json"))
  expect_error(open_pleiodb(tmp), "format version")
})

test_that("open_pleiodb errors on missing directory", {
  expect_error(open_pleiodb("/nonexistent/path"), class = "error")
})
