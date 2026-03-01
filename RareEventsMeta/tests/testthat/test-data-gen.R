# Tests for GenData

test_that("GenData returns data.frame with expected columns", {
  set.seed(123)
  out <- GenData(
    total_studies = 5,
    n1 = 100,
    n2 = 100,
    alpha = 10,
    beta = 10,
    psi = 1
  )
  expect_s3_class(out, "data.frame")
  expect_named(out, c("study", "size_1", "events_1", "size_2", "events_2"))
  expect_equal(nrow(out), 5L)
})

test_that("GenData study sizes match input", {
  set.seed(1)
  out <- GenData(
    total_studies = 3,
    n1 = 50,
    n2 = 75,
    alpha = 5,
    beta = 5,
    psi = 1
  )
  expect_equal(unique(out$size_1), 50)
  expect_equal(unique(out$size_2), 75)
})

test_that("GenData events are non-negative", {
  set.seed(2)
  out <- GenData(
    total_studies = 10,
    n1 = 100,
    n2 = 100,
    alpha = 2,
    beta = 2,
    psi = 1
  )
  expect_gte(min(out$events_1), 0)
  expect_gte(min(out$events_2), 0)
})

test_that("GenData is reproducible with same seed", {
  set.seed(99)
  a <- GenData(5, 50, 50, 2, 2, 1)
  set.seed(99)
  b <- GenData(5, 50, 50, 2, 2, 1)
  expect_identical(a$events_1, b$events_1)
  expect_identical(a$events_2, b$events_2)
})
