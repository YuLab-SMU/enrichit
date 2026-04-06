library(testthat)
library(enrichit)

test_that("TERM2NAME maps GSON term ids to names and falls back to ids", {
  x <- gson::gson(
    gsid2gene = data.frame(
      gsid = c("T1", "T1", "T2", "T3"),
      gene = c("g1", "g2", "g3", "g4")
    ),
    gsid2name = data.frame(
      gsid = c("T1", "T2"),
      name = c("Term 1", "Term 2")
    )
  )

  res <- enrichit:::TERM2NAME(c("T1", "T3", "T9"), x)

  expect_type(res, "character")
  expect_named(res, c("T1", "T3", "T9"))
  expect_identical(unname(res), c("Term 1", "T3", "T9"))
})

test_that("TERMID2EXTID returns genes grouped by requested GSON term ids", {
  x <- gson::gson(
    gsid2gene = data.frame(
      gsid = c("T1", "T1", "T2", "T3"),
      gene = c("g1", "g2", "g3", "g4")
    ),
    gsid2name = data.frame(
      gsid = c("T1", "T2"),
      name = c("Term 1", "Term 2")
    )
  )

  res <- enrichit:::TERMID2EXTID(c("T2", "T1", "T9"), x)

  expect_type(res, "list")
  expect_named(res, c("T2", "T1", "T9"))
  expect_true(all(vapply(res, is.character, logical(1))))
  expect_identical(res$T2, "g3")
  expect_identical(res$T1, c("g1", "g2"))
  expect_identical(res$T9, character(0))
})
