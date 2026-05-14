test_that("moc.gapbk returns expected structure on a small problem", {
  set.seed(123)
  x <- matrix(stats::runif(30 * 10, min = -2, max = 2), nrow = 30, ncol = 10)
  d1 <- as.matrix(stats::dist(x, method = "euclidean"))
  d2 <- as.matrix(stats::dist(x, method = "manhattan"))

  res <- moc.gapbk(d1, d2, num_k = 3, generation = 3, pop_size = 6)

  expect_type(res, "list")
  expect_named(res, c("population", "matrix.solutions", "clustering"))
  expect_s3_class(res$population, "data.frame")
  expect_s3_class(res$matrix.solutions, "data.frame")
  expect_type(res$clustering, "list")

  # Each clustering vector should have one entry per object
  expect_true(all(vapply(res$clustering, length, integer(1)) == nrow(d1)))

  # Cluster labels must be within 1..k
  expect_true(all(vapply(res$clustering,
                         function(v) all(v >= 1L & v <= 3L),
                         logical(1))))
})

test_that("input validation rejects invalid arguments", {
  d <- matrix(0, 5, 5)

  # num_k must be > 1
  expect_error(moc.gapbk(d, d, num_k = 1),
               "more than 1")

  # neighborhood must be in [0, 1]
  expect_error(moc.gapbk(d, d, num_k = 2, neighborhood = 1.5),
               "between 0 and 1")

  # mismatched dimensions between d1 and d2
  d_small <- matrix(0, 4, 4)
  expect_error(moc.gapbk(d, d_small, num_k = 2),
               "equal number of rows and columns")

  # non-square matrix
  d_rect <- matrix(0, 5, 4)
  expect_error(moc.gapbk(d_rect, d, num_k = 2),
               "same number of rows and columns")
})
