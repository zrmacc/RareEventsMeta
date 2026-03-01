# Tests for Logit and Expit (internal / used by RunMC)

test_that("Expit and Logit are inverses on (0, 1)", {
  x <- seq(0.01, 0.99, length.out = 20)
  expect_equal(Expit(Logit(x)), x, tolerance = 1e-10)
  expect_equal(Logit(Expit(x)), x, tolerance = 1e-10)
})

test_that("Expit maps R to (0, 1)", {
  expect_equal(Expit(-Inf), 0)
  expect_equal(Expit(Inf), 1)
  expect_equal(Expit(0), 0.5)
  expect_true(all(Expit(c(-10, 0, 10)) > 0 & Expit(c(-10, 0, 10)) < 1))
})

test_that("Logit maps (0, 1) to R", {
  expect_equal(Logit(0.5), 0)
  expect_true(Logit(0.01) < -4)
  expect_true(Logit(0.99) > 4)
})

test_that("Expit and Logit handle vector input", {
  x <- c(0.1, 0.5, 0.9)
  expect_equal(length(Expit(x)), 3L)
  expect_equal(length(Logit(x)), 3L)
  expect_equal(Expit(Logit(x)), x, tolerance = 1e-10)
})
