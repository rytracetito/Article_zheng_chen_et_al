test_that("power_study returns rejection rates and their standard errors", {
  set.seed(1)
  res <- power_study(function() sim_model(0:10, "gamma", c(a = 2, b = 3), n_traj = 4),
                     "gamma", n_rep = 4, B = 10, n_pool = 50)
  expect_s3_class(res, "gof_power")
  expect_named(res$rate, c("KS_sup", "KS_sum", "CvM_sup", "CvM_sum"))
  expect_true(all(res$rate >= 0 & res$rate <= 1))
  expect_equal(unname(res$se), unname(sqrt(res$rate * (1 - res$rate) / 4)))
  expect_output(print(res), "Rejection rates")
})

test_that("power_study counts or excludes the data sets that cannot be fitted", {
  gen <- function() sim_model(0:10, "wiener", c(mu = 0.1, sigma = 1), n_traj = 4)
  set.seed(2)
  r1 <- power_study(gen, "gamma", n_rep = 3, B = 5, n_pool = 50)
  expect_equal(r1$n_failed, 3)
  expect_equal(unname(r1$rate), rep(1, 4))
  set.seed(2)
  expect_error(power_study(gen, "gamma", n_rep = 3, B = 5, n_pool = 50, failure = "exclude"),
               "could not be fitted to any data set")
})
