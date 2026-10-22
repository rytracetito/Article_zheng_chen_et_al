test_that("the distance matrices follow equations 3 and 4 of the article", {
  set.seed(1)
  X <- matrix(stats::rnorm(5 * 8), 5, 8)
  Y <- matrix(stats::rnorm(13 * 8), 13, 8)
  ks <- cvm <- matrix(0, 5, 13)
  for (i in 1:5) for (k in 1:13) {
    ks[i, k] <- max(abs(X[i, ] - Y[k, ]))
    cvm[i, k] <- sum((X[i, ] - Y[k, ])^2)
  }
  expect_equal(gofdegrad:::distance_matrix(X, Y, "KS"), ks)
  expect_equal(gofdegrad:::distance_matrix(X, Y, "CvM"), cvm)
  s <- similarity_stats(X, Y)
  expect_equal(unname(s), c(max(apply(ks, 1, min)), sum(apply(ks, 1, min)),
                            max(apply(cvm, 1, min)), sum(apply(cvm, 1, min))))
})

test_that("with one trajectory the sup and sum statistics are equal", {
  set.seed(2)
  s <- similarity_stats(matrix(stats::rnorm(8), 1), matrix(stats::rnorm(80), 10))
  expect_equal(s[["KS_sup"]], s[["KS_sum"]])
  expect_equal(s[["CvM_sup"]], s[["CvM_sum"]])
})

test_that("similarity_stats accepts levysim objects", {
  set.seed(3)
  x <- sim_model(0:10, "gamma", c(a = 2, b = 3), n_traj = 3)
  pool <- sim_model(0:10, "gamma", c(a = 2, b = 3), n_traj = 50)
  expect_equal(similarity_stats(x, pool), similarity_stats(gof_data(x)$X, gof_data(pool)$X))
})

test_that("gof_test returns a coherent result", {
  set.seed(4)
  x <- sim_model(0:20, "wiener", c(mu = 0.5, sigma = 0.3), n_traj = 6)
  res <- gof_test(x, "wiener", B = 60, n_pool = 200)
  expect_s3_class(res, "gof_test")
  expect_named(res$estimate, c("mu", "sigma"))
  expect_equal(res$n_pool, 200)
  expect_equal(nrow(res$boot), 60 - res$n_failed)
  expect_true(all(res$p.value >= 0 & res$p.value <= 1))
  expect_equal(res$reject, res$statistic > res$critical)
  expect_equal(unname(res$critical),
               unname(apply(res$boot, 2, stats::quantile, probs = 0.95, names = FALSE)))
  tab <- as.data.frame(res)
  expect_equal(tab$stat, c("KS_sup", "KS_sum", "CvM_sup", "CvM_sum"))
  expect_s3_class(plot(res), "ggplot")
  expect_output(print(res), "H0 rejected by")
})

test_that("gof_test uses n_pool = 100 m by default", {
  set.seed(5)
  x <- sim_model(0:5, "gamma", c(a = 2, b = 3), n_traj = 2)
  expect_equal(gof_test(x, "gamma", B = 5)$n_pool, 200)
})

test_that("gof_test rejects a clearly wrong model", {
  set.seed(6)
  x <- sim_model(1:30, "wiener_power", c(mu = 0.01, sigma = 0.05, r = 2), n_traj = 6)
  res <- gof_test(x, "wiener", B = 60, n_pool = 300)
  expect_true(all(res$reject))
})

test_that("gof_test stops when the model cannot be fitted", {
  set.seed(7)
  x <- sim_model(0:20, "wiener", c(mu = 0.1, sigma = 1), n_traj = 4)
  expect_error(gof_test(x, "gamma", B = 5), "cannot be fitted")
})

test_that("gof_test works when the gamma bootstrap gives zero increments", {
  set.seed(1)
  t <- round(cumsum(runif(50, 0.5, 1.5)), 2)
  set.seed(7)
  x <- sim_model(t, "ig", c(mu = 1.5, eta = 0.1), n_traj = 6)
  res <- suppressWarnings(gof_test(x, "gamma", B = 50, n_pool = 600))
  expect_true(all(is.finite(res$statistic)))
  expect_equal(res$n_failed, 0)
})

test_that("gof_test accepts a custom model", {
  lin <- custom_model("linear path", c("mu", "s_phi", "s_eps"),
    fit = function(X, t) {
      sl <- apply(X, 1, function(x) sum(x * t) / sum(t^2))
      c(mu = mean(sl), s_phi = stats::sd(sl), s_eps = stats::sd(as.vector(X - outer(sl, t))))
    },
    simulate = function(theta, t, n) {
      outer(stats::rnorm(n, theta[["mu"]], theta[["s_phi"]]), t) +
        matrix(stats::rnorm(n * length(t), 0, theta[["s_eps"]]), n)
    })
  set.seed(8)
  x <- sim_path(1:15, "linear", 0.5, 0.1, sigma_eps = 0.2, n_traj = 6)
  res <- gof_test(x, lin, B = 30, n_pool = 200)
  expect_equal(res$model, "linear path")
  expect_named(res$estimate, c("mu", "s_phi", "s_eps"))
})

test_that("compare_models reports the models that cannot be tested", {
  set.seed(9)
  x <- sim_model(0:15, "wiener", c(mu = 0.1, sigma = 1), n_traj = 4)
  res <- compare_models(x, c("wiener", "gamma"), B = 20, n_pool = 100, seed = 1)
  expect_equal(res$model, c("wiener", "gamma"))
  expect_true(is.na(res$rejections[2]))
  expect_match(res$note[2], "not testable")
  expect_s3_class(attr(res, "tests")$wiener, "gof_test")
})

test_that("compare_models with a seed is reproducible", {
  set.seed(10)
  x <- sim_model(0:10, "gamma", c(a = 2, b = 3), n_traj = 3)
  r1 <- compare_models(x, c("gamma", "ig"), B = 20, n_pool = 100, seed = 3)
  r2 <- compare_models(x, c("ig", "gamma"), B = 20, n_pool = 100, seed = 3)
  expect_equal(r1[1, ], r2[2, ], ignore_attr = TRUE)
})
