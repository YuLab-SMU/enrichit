library(testthat)
library(enrichit)

test_that("calculate_qvalue handles invalid p-values without poisoning all results", {
  skip_if_not_installed("qvalue")

  pvals <- c(0.01, 0.5, NA_real_, -0.1, 1.1)
  qv <- expect_warning(enrichit:::calculate_qvalue(pvals), "Invalid p-values detected")

  expect_length(qv, length(pvals))
  expect_true(is.double(qv))

  expect_false(is.na(qv[1]))
  expect_false(is.na(qv[2]))
  expect_true(qv[1] >= 0 && qv[1] <= 1)
  expect_true(qv[2] >= 0 && qv[2] <= 1)

  expect_true(is.na(qv[3]))
  expect_true(is.na(qv[4]))
  expect_true(is.na(qv[5]))
})
