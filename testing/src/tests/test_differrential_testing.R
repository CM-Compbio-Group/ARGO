library(testthat)
source('../R/run_differential_testing.R')

test_that("TypeError for expt", {
  expect_error(
    run_differential_testing(expt = 150.32), 
    regexp = "must be integer"
  )
}) 

test_that("TypeError for expt", {
  expect_error(
    run_differential_testing(expt = 'abcx'), 
    regexp = "must be integer"
  )
}) 

test_that("ValueError for expt", {
  expect_error(
    run_differential_testing(expt = 0), 
    regexp = "should be positive"
  )
}) 

test_that("TypeError for N", {
  expect_error(
    run_differential_testing(N = 150.32), 
    regexp = "must be integer"
  )
}) 

test_that("TypeError for N", {
  expect_error(
    run_differential_testing(N = FALSE), 
    regexp = "must be integer"
  )
}) 

test_that("ValueError for N", {
  expect_error(
    run_differential_testing(N = 0), 
    regexp = "not suitable"
  )
}) 

test_that("ValueError for N", {
  expect_error(
    run_differential_testing(N = -2), 
    regexp = "not suitable"
  )
}) 

test_that("ValueError for N", {
  expect_error(
    run_differential_testing(N = 5e3), 
    regexp = "not suitable"
  )
}) 

