true_theta <- list(
  wiener = c(mu = 0.5, sigma = 0.8),
  gamma = c(a = 2, b = 3),
  ig = c(mu = 0.5, eta = 2),
  vg = c(theta = 0.5, sigma = 1, nu = 0.5),
  meixner = c(a = 1, b = 0.3, d = 1.5, m = 0.2),
  wiener_power = c(mu = 1, sigma = 0.3, r = 0.5),
  gamma_power = c(a = 1, b = 2, r = 1.5)
)

test_that("gof_models lists the seven models", {
  g <- gof_models()
  expect_s3_class(g, "data.frame")
  expect_equal(g$model, names(true_theta))
  expect_equal(g$increasing, c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE))
})

test_that("sim_model returns levysim matrices starting at 0", {
  set.seed(1)
  for (m in names(true_theta)) {
    x <- sim_model(c(0, 1, 2.5, 4), m, true_theta[[m]], n_traj = 3)
    expect_s3_class(x, "levysim")
    expect_equal(dim(x), c(4, 4))
    expect_equal(unname(x[1, -1]), c(0, 0, 0))
    expect_true(all(is.finite(x)))
  }
  expect_equal(dim(sim_model(1:5, "wiener", true_theta$wiener)), c(6, 2))
  expect_error(sim_model(1:5, "gamma", c(a = 2)), "missing parameter")
  expect_error(sim_model(1:5, "lognormal", c(a = 2)), "unknown model")
})

test_that("the increasing models give increasing trajectories", {
  set.seed(1)
  for (m in gof_models()$model[gof_models()$increasing]) {
    X <- gof_data(sim_model(0:20, m, true_theta[[m]], n_traj = 20))$X
    expect_true(all(diff(t(X)) >= 0))
  }
})

test_that("fit_model recovers the parameters on large samples", {
  set.seed(2)
  t <- 0:30
  tol <- c(wiener = 0.05, gamma = 0.05, ig = 0.05, vg = 0.2, wiener_power = 0.05, gamma_power = 0.05)
  for (m in names(tol)) {
    x <- sim_model(t, m, true_theta[[m]], n_traj = 300)
    e <- fit_model(x, m)
    expect_equal(attr(e, "model"), m)
    expect_equal(as.vector(e), unname(true_theta[[m]]), tolerance = tol[[m]], info = m)
  }
})

test_that("fit_model works with irregular inspection times", {
  set.seed(3)
  t <- c(0, round(cumsum(runif(40, 0.5, 1.5)), 2))
  for (m in c("wiener", "gamma", "ig", "wiener_power", "gamma_power")) {
    x <- sim_model(t, m, true_theta[[m]], n_traj = 300)
    expect_equal(as.vector(fit_model(x, m)), unname(true_theta[[m]]), tolerance = 0.06, info = m)
  }
})

test_that("the Meixner simulation has the moments of the Meixner process", {
  set.seed(4)
  th <- true_theta$meixner
  t <- c(0.5, 2, 2.3, 5)
  X <- gof_data(sim_model(t, "meixner", th, n_traj = 20000))$X
  k1 <- th[["m"]] + th[["a"]] * th[["d"]] * tan(th[["b"]] / 2)
  k2 <- th[["a"]]^2 * th[["d"]] / (2 * cos(th[["b"]] / 2)^2)
  expect_equal(unname(colMeans(X)), k1 * t, tolerance = 0.03)
  expect_equal(unname(apply(X, 2, stats::var)), k2 * t, tolerance = 0.05)
})

test_that("the Meixner model can be fitted", {
  set.seed(5)
  x <- sim_model(0:200, "meixner", true_theta$meixner, n_traj = 5)
  e <- suppressWarnings(fit_model(x, "meixner"))
  expect_named(e, c("a", "b", "d", "m"))
})

test_that("increasing models: negative increments give an error, zero increments a warning", {
  t <- 1:5
  expect_error(fit_model(c(1, 2, 1.5, 3, 4), "gamma", t = t), "negative increments")
  expect_error(fit_model(c(1, 2, 1.5, 3, 4), "ig", t = t), "negative increments")
  expect_error(fit_model(c(1, 2, 1.5, 3, 4), "gamma_power", t = t), "negative increments")
  expect_warning(fit_model(c(1, 2, 2, 3, 4), "gamma", t = t), "1 zero increment")
  expect_warning(fit_model(c(1, 2, 2, 3, 4), "ig", t = t), "zero increment")
})

test_that("a variance-gamma model fitted to increasing data gives NA with a warning", {
  set.seed(6)
  x <- sim_model(0:50, "gamma", true_theta$gamma, n_traj = 3)
  expect_warning(e <- fit_model(x, "vg"))
  expect_true(anyNA(e))
})

test_that("custom_model can be fitted and simulated", {
  lin <- custom_model("linear path", c("mu", "s"),
                      fit = function(X, t) c(mu = mean(X / matrix(t, nrow(X), length(t), byrow = TRUE)), s = 1),
                      simulate = function(theta, t, n) outer(rep(theta[["mu"]], n), t))
  expect_s3_class(lin, "gof_model")
  x <- sim_model(1:5, lin, c(mu = 2, s = 1), n_traj = 2)
  expect_equal(unname(x[, "X1"]), 2 * (0:5))
  expect_equal(as.vector(fit_model(x, lin)), c(2, 1))
  bad <- custom_model("bad", "mu", fit = function(X, t) c(nu = 1), simulate = function(theta, t, n) 0)
  expect_error(fit_model(x, bad), "named numeric vector")
})

test_that("sim_path simulates the three general path models", {
  set.seed(7)
  x <- sim_path(0:10, "linear", 0.2, 0, sigma_eps = 0, n_traj = 3)
  expect_equal(unname(x[, "X1"]) / 0.2, 0:10)
  x <- sim_path(0:10, "power", c(1, 0.5), c(0, 0), sigma_eps = 0, n_traj = 2)
  expect_equal(unname(x[, "X2"]), sqrt(0:10))
  x <- sim_path(0:10, "exponential", 0.5, 0, sigma_eps = 0, omega = 0.1)
  expect_equal(unname(x[, "X1"]), 0.5 * (exp(0.1 * (0:10)) - 1))
  expect_error(sim_path(0:10, "power", 1, 0.1, sigma_eps = 0), "length 2")
  expect_error(sim_path(0:10, "exponential", 1, 0.1, sigma_eps = 0), "omega")
})
