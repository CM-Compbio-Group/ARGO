library(testthat)
source('../R/run_interval_evaluation.R')

## tests for run_interval_evaluation
test_that("TypeError for N", {
  expect_error(
    run_interval_evaluation(N = 'abc'), 
    regexp = "must be integer"
  )
}) 

test_that("TypeError for N", {
  expect_error(
    run_interval_evaluation(N = 0.123), 
    regexp = "must be integer"
  )
}) 

test_that("TypeError for N", {
  expect_error(
    run_interval_evaluation(N = pi), 
    regexp = "must be integer"
  )
}) 

test_that("TypeError for N", {
  expect_error(
    run_interval_evaluation(N = 1.11e1), 
    regexp = "must be integer"
  )
}) 

test_that("ValueError for N", {
  expect_error(
    run_interval_evaluation(N = 0), 
    regexp = "not suitable for pseudo-bulking"
  )
}) 

test_that("ValueError for N", {
  expect_error(
    run_interval_evaluation(N = -3), 
    regexp = "not suitable for pseudo-bulking"
  )
})


test_that("ValueError for N", {
  expect_error(
    run_interval_evaluation(N = -1e2), 
    regexp = "not suitable for pseudo-bulking"
  )
}) 


test_that("ValueError for N", {
  expect_error(
    run_interval_evaluation(N = 1.25e6), 
    regexp = "not suitable for pseudo-bulking"
  )
}) 

test_that("TypeError for parameter", {
  expect_error(
    run_interval_evaluation(FIT_CONDITIONAL_NOISE = 'FALSE'), 
    regexp = "must be boolean"
  )
}) 

test_that("TypeError for parameter", {
  expect_error(
    run_interval_evaluation(SAVE_FITS = 0), 
    regexp = "must be boolean"
  )
}) 

test_that("TypeError for parameter", {
  expect_error(
    run_interval_evaluation(SAVE_FITS = 'TRUE'), 
    regexp = "must be boolean"
  )
}) 
