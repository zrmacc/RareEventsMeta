# Smoke tests for FindBound and ExactConfInt (use small reps/coarse step for speed).
# Use fixed data with non-zero events to avoid PrepData/RunMC edge cases.

test_that("FindBound lower returns list with search_results and boundary_region_check", {
  # Fixed data: no zero-event studies to avoid dhyper/PrepData edge cases.
  data <- data.frame(
    size_1 = c(100, 100, 100),
    events_1 = c(8, 12, 10),
    size_2 = c(100, 100, 100),
    events_2 = c(5, 7, 6)
  )
  out <- FindBound(
    size_1 = data$size_1,
    events_1 = data$events_1,
    size_2 = data$size_2,
    events_2 = data$events_2,
    lower_bound = TRUE,
    reps = 25L,
    step_size = 0.05,
    maxit = 20L
  )
  expect_type(out, "list")
  expect_named(out, c("search_results", "boundary_region_check"))
  expect_named(out$search_results, c("bound", "mc", "pval"))
  expect_true(is.numeric(out$search_results$bound))
  expect_true(out$search_results$bound >= 0 && out$search_results$bound <= 1)
})

test_that("FindBound upper returns valid bound", {
  data <- data.frame(
    size_1 = c(100, 100),
    events_1 = c(10, 15),
    size_2 = c(100, 100),
    events_2 = c(6, 8)
  )
  out <- FindBound(
    size_1 = data$size_1,
    events_1 = data$events_1,
    size_2 = data$size_2,
    events_2 = data$events_2,
    lower_bound = FALSE,
    reps = 25L,
    step_size = 0.05,
    maxit = 20L
  )
  expect_true(out$search_results$bound >= 0 && out$search_results$bound <= 1)
})

test_that("ExactConfInt returns data.frame with lower and upper", {
  data <- data.frame(
    size_1 = c(100, 100),
    events_1 = c(10, 12),
    size_2 = c(100, 100),
    events_2 = c(6, 7)
  )
  ci <- ExactConfInt(
    size_1 = data$size_1,
    events_1 = data$events_1,
    size_2 = data$size_2,
    events_2 = data$events_2,
    reps = 20L,
    step_size = 0.05,
    maxit = 15L
  )
  expect_s3_class(ci, "data.frame")
  expect_named(ci, c("lower", "upper"))
  expect_equal(nrow(ci), 1L)
  if (!is.na(ci$lower) && !is.na(ci$upper)) {
    expect_lte(ci$lower, ci$upper)
  }
})
