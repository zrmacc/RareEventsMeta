# Tests for PrepData and MomentEst

test_that("MomentEst returns list with expected names", {
  size_1 <- c(100, 100)
  events_1 <- c(5, 10)
  size_2 <- c(100, 100)
  events_2 <- c(3, 8)
  out <- MomentEst(size_1, events_1, size_2, events_2)
  expect_type(out, "list")
  expect_named(out, c("mu", "mu_se2", "mu_cc", "mu2", "nu"))
  expect_length(out, 5L)
})

test_that("MomentEst mu and mu_se2 are numeric and reasonable", {
  size_1 <- c(50, 50)
  events_1 <- c(2, 4)
  size_2 <- c(50, 50)
  events_2 <- c(1, 3)
  out <- MomentEst(size_1, events_1, size_2, events_2)
  expect_true(is.numeric(out$mu))
  expect_true(out$mu >= 0 && out$mu <= 1)
  expect_true(is.numeric(out$mu_se2) && out$mu_se2 >= 0)
})

test_that("MomentEst with study identifier", {
  size_1 <- c(100)
  events_1 <- c(5)
  size_2 <- c(100)
  events_2 <- c(3)
  out <- MomentEst(size_1, events_1, size_2, events_2, study = 1L)
  expect_equal(out$mu, 5 / (5 + 3), tolerance = 1e-5)
})

test_that("PrepData (internal) produces valid weights when called via MomentEst", {
  # MomentEst calls PrepData internally; we check MomentEst output consistency.
  size_1 <- 10
  events_1 <- 2
  size_2 <- 10
  events_2 <- 1
  out <- MomentEst(size_1, events_1, size_2, events_2)
  expect_true(is.numeric(out$mu) && out$mu >= 0 && out$mu <= 1)
  expect_true(out$mu_se2 >= 0)
})
