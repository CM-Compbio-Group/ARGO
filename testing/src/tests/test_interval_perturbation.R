library(testthat)
source('../R/generate_interval_perturbation.R')

test_that("TypeError for N_GENES", {
  expect_error(
    generate_interval_perturbation(N_GENES = 'abc'), 
    regexp = "must be integer"
  )
}) 

test_that("TypeError for N_GENES", {
  expect_error(
    generate_interval_perturbation(N_GENES = 0.123), 
    regexp = "must be integer"
  )
}) 

test_that("TypeError for N_GENES", {
  expect_error(
    generate_interval_perturbation(N_GENES = FALSE), 
    regexp = "must be integer"
  )
}) 

test_that("TypeError for N_GENES", {
  expect_error(
    generate_interval_perturbation(N_GENES = 1.005e2), 
    regexp = "must be integer"
  )
}) 


test_that("TypeError for N_SAMPLES", {
  expect_error(
    generate_interval_perturbation(N_SAMPLES = 'abc'), 
    regexp = "must be integer"
  )
}) 

test_that("TypeError for N_SAMPLES", {
  expect_error(
    generate_interval_perturbation(N_SAMPLES = -3.23), 
    regexp = "must be integer"
  )
}) 

test_that("TypeError for N_SAMPLES", {
  expect_error(
    generate_interval_perturbation(N_SAMPLES = 1.005e2), 
    regexp = "must be integer"
  )
}) 


test_that("ValueError for N_GENES", {
  expect_error(
    generate_interval_perturbation(N_GENES = 0), 
    regexp = "at least one gene"
  )
}) 

test_that("ValueError for N_GENES", {
  expect_error(
    generate_interval_perturbation(N_GENES = -3), 
    regexp = "at least one gene"
  )
}) 

test_that("ValueError for N_GENES", {
  expect_error(
    generate_interval_perturbation(N_GENES = 1e6), 
    regexp = "seems inflated"
  )
}) 


test_that("ValueError for N_SAMPLES", {
  expect_error(
    generate_interval_perturbation(N_SAMPLES = 0), 
    regexp = "must be positive"
  )
}) 

test_that("ValueError for N_SAMPLES", {
  expect_error(
    generate_interval_perturbation(N_SAMPLES = -10), 
    regexp = "must be positive"
  )
})

test_that("ValueError for N_SAMPLES", {
  expect_error(
    generate_interval_perturbation(N_SAMPLES = 5e2), 
    regexp = "seems inflated"
  )
})


