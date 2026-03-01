# Tests for BoundaryAB

test_that("BoundaryAB returns named numeric of length 2", {
  out <- BoundaryAB(mu = 0.5)
  expect_type(out, "double")
  expect_named(out, c("alpha", "beta"))
  expect_length(out, 2L)
})

test_that("BoundaryAB alpha and beta are positive", {
  for (mu in c(0.2, 0.5, 0.8)) {
    out <- BoundaryAB(mu = mu)
    expect_gt(out["alpha"], 0)
    expect_gt(out["beta"], 0)
  }
})

test_that("BoundaryAB with explicit nu", {
  out1 <- BoundaryAB(mu = 0.5, v = NULL)
  out2 <- BoundaryAB(mu = 0.5, v = 0.05)
  expect_named(out2, c("alpha", "beta"))
  expect_gt(out2["alpha"], 0)
  expect_gt(out2["beta"], 0)
  # Different nu gives different (alpha, beta)
  expect_false(isTRUE(all.equal(out1["alpha"], out2["alpha"])))
})

test_that("BoundaryAB mean alpha/(alpha+beta) equals mu", {
  mu <- 0.4
  out <- BoundaryAB(mu = mu)
  alpha <- unname(out["alpha"])
  beta <- unname(out["beta"])
  expect_equal(alpha / (alpha + beta), mu, tolerance = 1e-5)
})
