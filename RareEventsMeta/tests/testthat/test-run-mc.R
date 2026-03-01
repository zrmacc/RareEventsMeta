# Tests for CalcChi2 and RunMC

test_that("CalcChi2 returns scalar when given scalars", {
  x <- CalcChi2(mu_hat = 0.5, mu = 0.4, se2 = 0.01)
  expect_length(x, 1L)
  expect_true(is.numeric(x) && x >= 0)
})

test_that("CalcChi2 formula", {
  expect_equal(CalcChi2(0.5, 0.5, 0.01), 0)
  expect_equal(CalcChi2(0.6, 0.5, 0.01), (0.1^2) / 0.01)
})

test_that("RunMC returns data.frame when p_only is FALSE", {
  set.seed(42)
  size_1 <- c(50, 50)
  events_1 <- c(3, 5)
  size_2 <- c(50, 50)
  events_2 <- c(2, 4)
  out <- RunMC(
    size_1 = size_1,
    events_1 = events_1,
    size_2 = size_2,
    events_2 = events_2,
    reps = 20L,
    alpha = 5,
    beta = 5
  )
  expect_s3_class(out, "data.frame")
  expect_named(out, c("alpha", "beta", "mu", "nu", "p"))
  expect_true(out$p >= 0 && out$p <= 1)
})

test_that("RunMC with mu instead of alpha, beta", {
  set.seed(42)
  size_1 <- 100
  events_1 <- 10
  size_2 <- 100
  events_2 <- 8
  out <- RunMC(
    size_1 = size_1,
    events_1 = events_1,
    size_2 = size_2,
    events_2 = events_2,
    reps = 15L,
    mu = 0.5
  )
  expect_true("mu" %in% names(out))
  expect_equal(out$mu, 0.5, tolerance = 1e-5)
  expect_true(out$p >= 0 && out$p <= 1)
})

test_that("RunMC p_only returns named vector with p", {
  set.seed(42)
  out <- RunMC(
    size_1 = 50,
    events_1 = 5,
    size_2 = 50,
    events_2 = 3,
    reps = 10L,
    alpha = 5,
    beta = 5,
    p_only = TRUE
  )
  expect_named(out, "p")
  expect_true(out["p"] >= 0 && out["p"] <= 1)
})
