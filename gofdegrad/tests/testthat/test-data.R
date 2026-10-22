test_that("gof_data gives the same result for all the data forms", {
  t <- c(0, 1, 2.5, 4)
  X <- rbind(c(0, 0.5, 1.2, 2.0), c(0, 0.3, 1.0, 1.7))
  ref <- gof_data(X, t = t)
  expect_s3_class(ref, "gof_data")
  expect_equal(ref$t, c(1, 2.5, 4))
  expect_equal(ref$X, X[, -1])

  levysim_form <- cbind(t = t, X1 = X[1, ], X2 = X[2, ])
  expect_equal(gof_data(levysim_form), ref)

  long <- data.frame(id = rep(1:2, each = 4), t = rep(t, 2), X = c(X[1, ], X[2, ]))
  expect_equal(gof_data(long), ref)
  expect_equal(gof_data(long[sample(nrow(long)), ]), ref)

  liste <- list(data.frame(t = t, X = X[1, ]), data.frame(t = t, X = X[2, ]))
  expect_equal(gof_data(liste), ref)

  expect_equal(gof_data(X[1, ], t = t)$X, ref$X[1, , drop = FALSE])
})

test_that("gof_data subtracts the level at time 0 and keeps data without time 0", {
  d <- gof_data(data.frame(id = c(1, 1, 1), t = c(0, 1, 2), X = c(0.9, 1.1, 1.6)))
  expect_equal(d$t, c(1, 2))
  expect_equal(as.vector(d$X), c(0.2, 0.7))

  d <- gof_data(matrix(c(0.3, 0.8), nrow = 1), t = c(0.5, 1))
  expect_equal(d$t, c(0.5, 1))
  expect_equal(as.vector(d$X), c(0.3, 0.8))
})

test_that("gof_data accepts the results of levysim and of sim_model", {
  set.seed(1)
  x <- levysim::sim_gamma(0:10, a = 2, b = 3, n_traj = 4)
  d <- gof_data(x)
  expect_equal(dim(d$X), c(4, 10))
  expect_equal(d$t, 1:10)
  y <- sim_model(0:10, "gamma", c(a = 2, b = 3), n_traj = 3)
  expect_equal(dim(gof_data(y)$X), c(3, 10))
})

test_that("gof_data rejects invalid data", {
  expect_error(gof_data(list(data.frame(t = 0:2, X = 0:2), data.frame(t = c(0, 1, 3), X = 0:2))),
               "same times")
  expect_error(gof_data(matrix(1:6, 2), t = c(1, 3, 2)), "strictly increasing")
  expect_error(gof_data(matrix(1:6, 2), t = c(1, 2)), "one inspection time per column")
  expect_error(gof_data(data.frame(a = 1, b = 2)), "id, t and X")
  expect_error(gof_data(matrix(c(1, NA), 1), t = 1:2), "finite")
  expect_error(gof_data(0, t = 0), "after time 0")
})
