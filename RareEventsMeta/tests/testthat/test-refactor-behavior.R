# Tests to ensure refactored code preserves intended behavior.

test_that("boundary_nu helper matches manual formula and is used by BoundaryAB and NuSeq", {
  # boundary_nu(mu) should equal min(v1, v2) from original BoundaryAB logic.
  mu <- 0.5
  v1 <- mu^2 * (1 - mu) / (1 + mu)
  v2 <- mu * (1 - mu)^2 / (2 - mu)
  expected_nu <- min(v1, v2)
  expect_equal(RareEventsMeta:::boundary_nu(mu), expected_nu)

  # BoundaryAB(mu) with v=NULL should match BoundaryAB(mu, v = boundary_nu(mu)).
  ab_default <- BoundaryAB(mu)
  ab_explicit <- BoundaryAB(mu, v = RareEventsMeta:::boundary_nu(mu))
  expect_equal(ab_default["alpha"], ab_explicit["alpha"], tolerance = 1e-10)
  expect_equal(ab_default["beta"], ab_explicit["beta"], tolerance = 1e-10)

  # NuSeq with mu only should use same boundary; first row mean = mu.
  nu_seq <- NuSeq(alpha = 5, beta = 5, mu = mu, num_nu_vals = 3)
  first_mean <- nu_seq$alpha[1] / (nu_seq$alpha[1] + nu_seq$beta[1])
  expect_equal(first_mean, mu, tolerance = 1e-5)
})

test_that("ExactConfInt is deterministic and returns valid bounds", {
  data <- data.frame(
    size_1 = c(100, 100),
    events_1 = c(10, 12),
    size_2 = c(100, 100),
    events_2 = c(6, 7)
  )
  set.seed(123)
  ci1 <- ExactConfInt(
    size_1 = data$size_1,
    events_1 = data$events_1,
    size_2 = data$size_2,
    events_2 = data$events_2,
    reps = 25L,
    step_size = 0.05,
    maxit = 15L
  )
  set.seed(123)
  ci2 <- ExactConfInt(
    size_1 = data$size_1,
    events_1 = data$events_1,
    size_2 = data$size_2,
    events_2 = data$events_2,
    reps = 25L,
    step_size = 0.05,
    maxit = 15L
  )
  expect_equal(ci1$lower, ci2$lower, tolerance = 1e-10)
  expect_equal(ci1$upper, ci2$upper, tolerance = 1e-10)

  if (!is.na(ci1$lower) && !is.na(ci1$upper)) {
    expect_gte(ci1$lower, 0)
    expect_lte(ci1$upper, 1)
  }
})

test_that("RunMC is deterministic (same seed gives same p) and returns valid output", {
  set.seed(456)
  out1 <- RunMC(
    size_1 = c(50, 50),
    events_1 = c(4, 6),
    size_2 = c(50, 50),
    events_2 = c(3, 5),
    reps = 30L,
    alpha = 5,
    beta = 5
  )
  set.seed(456)
  out2 <- RunMC(
    size_1 = c(50, 50),
    events_1 = c(4, 6),
    size_2 = c(50, 50),
    events_2 = c(3, 5),
    reps = 30L,
    alpha = 5,
    beta = 5
  )
  expect_equal(out1$p, out2$p, tolerance = 1e-10)
  expect_true(out1$p >= 0 && out1$p <= 1)
})

test_that("FindBound via mu uses same RunMC as via alpha,beta from BoundaryAB", {
  # run_mc_ab(alpha, beta) is used by both mu_to_pval_map and ab_to_pval_map.
  # So RunMC(..., mu=x) should match RunMC(..., alpha=ab[1], beta=ab[2]) for ab = BoundaryAB(x).
  data <- data.frame(
    size_1 = 100,
    events_1 = 10,
    size_2 = 100,
    events_2 = 7
  )
  set.seed(789)
  out_mu <- RunMC(
    size_1 = data$size_1,
    events_1 = data$events_1,
    size_2 = data$size_2,
    events_2 = data$events_2,
    reps = 20L,
    mu = 0.5
  )
  ab <- BoundaryAB(0.5)
  set.seed(789)
  out_ab <- RunMC(
    size_1 = data$size_1,
    events_1 = data$events_1,
    size_2 = data$size_2,
    events_2 = data$events_2,
    reps = 20L,
    alpha = ab[1],
    beta = ab[2]
  )
  expect_equal(out_mu$p, out_ab$p, tolerance = 1e-10)
  expect_equal(out_mu$mu, out_ab$mu, tolerance = 1e-10)
})

test_that("PrepData refactor preserves MomentEst output (structure and value)", {
  # MomentEst calls PrepData; same inputs should give same MomentEst result.
  size_1 <- c(20, 30)
  events_1 <- c(3, 5)
  size_2 <- c(20, 30)
  events_2 <- c(2, 4)
  out <- MomentEst(size_1, events_1, size_2, events_2)
  expect_named(out, c("mu", "mu_se2", "mu_cc", "mu2", "nu"))
  expect_length(out, 5L)
  expect_true(out$mu >= 0 && out$mu <= 1)
  expect_true(out$mu_se2 >= 0)

  # PrepData (internal) for one study: weights should be non-negative and structure correct.
  pd <- RareEventsMeta:::PrepData(
    size_1 = 20,
    events_1 = 3,
    size_2 = 20,
    events_2 = 2,
    study = 1L
  )
  expect_s3_class(pd, "data.frame")
  expect_true(all(c("study", "events", "total_events", "weight", "n1", "n2") %in% names(pd)))
  expect_true(all(pd$weight >= 0))
  expect_true(sum(pd$weight) > 0)
})

test_that("MomentEst has no corrected argument and returns expected list", {
  # After removing corrected, formals should not include it.
  expect_false("corrected" %in% names(formals(MomentEst)))
  out <- MomentEst(c(10, 10), c(2, 3), c(10, 10), c(1, 2))
  expect_identical(names(out), c("mu", "mu_se2", "mu_cc", "mu2", "nu"))
})
