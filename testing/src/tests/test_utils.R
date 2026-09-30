library(testthat)
source('../R/utils.R')

intA <- seq(5)
intB <- seq(10)

zeallot::`%<-%`(c(jaccard, overlap, dice, precision, recall), compute_overlap_statistics(intA, intB))

test_that("overlap statistics are being correctly evaluated", {
  expect_equal(jaccard, 5/10)
  expect_equal(overlap, 5/5)
  expect_equal(dice, 2*5 / (5 + 10))
  expect_equal(precision, 5 / 5)
  expect_equal(recall, 5/10)
})
