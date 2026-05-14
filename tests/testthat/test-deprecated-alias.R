test_that("moc.gabk alias warns about deprecation but still works", {
  set.seed(456)
  x <- matrix(stats::runif(20 * 8, min = -1, max = 1), nrow = 20, ncol = 8)
  d1 <- as.matrix(stats::dist(x, method = "euclidean"))
  d2 <- as.matrix(stats::dist(x, method = "manhattan"))

  expect_warning(
    res <- moc.gabk(d1, d2, num_k = 2, generation = 2, pop_size = 4),
    "deprecated"
  )

  expect_named(res, c("population", "matrix.solutions", "clustering"))
})
