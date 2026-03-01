# Tests for NuSeq

test_that("NuSeq returns data.frame with correct number of rows", {
  out <- NuSeq(alpha = 5, beta = 5, num_nu_vals = 7)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 7L)
  expect_true(all(c("alpha", "beta") %in% names(out)))
})

test_that("NuSeq with explicit mu overrides alpha/beta mean", {
  # When mu is supplied, it is used for the boundary; alpha, beta still required for signature.
  out <- NuSeq(alpha = 5, beta = 5, mu = 0.5, num_nu_vals = 5)
  expect_equal(nrow(out), 5L)
  expect_true("alpha" %in% names(out) && "beta" %in% names(out))
})

test_that("NuSeq alpha and beta are positive", {
  out <- NuSeq(alpha = 3, beta = 7, num_nu_vals = 4)
  expect_true(all(out$alpha > 0))
  expect_true(all(out$beta > 0))
})

test_that("NuSeq with alpha and beta preserves mean", {
  alpha <- 2
  beta <- 8
  mu <- alpha / (alpha + beta)
  out <- NuSeq(alpha = alpha, beta = beta, num_nu_vals = 3)
  for (i in seq_len(nrow(out))) {
    expect_equal(out$alpha[i] / (out$alpha[i] + out$beta[i]), mu, tolerance = 1e-5)
  }
})
