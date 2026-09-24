test_that("sim_proc Wiener respecte E et Var du Tableau 1 (homogene)", {
  set.seed(1)
  T <- c(1, 3, 7, 12, 20)
  mu_W <- 0.35; sigma <- 0.6
  X <- sim_proc(T, list(mu_W = mu_W, sigma = sigma), "Wiener", mp = 100000)
  expect_equal(mean(X[, 5]), mu_W * T[5], tolerance = 0.02)
  expect_equal(stats::var(X[, 5]), sigma^2 * T[5], tolerance = 0.05)
})

test_that("sim_proc Wiener respecte E et Var du Tableau 1 (non-homogene)", {
  set.seed(1)
  T <- c(1, 3, 7, 12, 20)
  mu_W <- 0.35; sigma <- 0.6; r <- 1.5
  X <- sim_proc(T, list(mu_W = mu_W, sigma = sigma), "Wiener", mp = 100000,
                Lam = function(t) t^r)
  expect_equal(mean(X[, 5]), mu_W * T[5]^r, tolerance = 0.02)
  expect_equal(stats::var(X[, 5]), sigma^2 * T[5], tolerance = 0.05)  # Var ~ t, PAS Lambda(t)
})

test_that("sim_proc Gamma respecte E et Var du Tableau 1", {
  set.seed(1)
  T <- c(1, 3, 7, 12, 20)
  alpha <- 1.7; beta <- 0.45
  X <- sim_proc(T, list(alpha = alpha, beta = beta), "Gamma", mp = 100000)
  expect_equal(mean(X[, 5]), alpha * beta * T[5], tolerance = 0.02)
  expect_equal(stats::var(X[, 5]), alpha * beta^2 * T[5], tolerance = 0.05)
})

test_that("sim_proc IG respecte E et Var du Tableau 1", {
  set.seed(1)
  T <- c(1, 3, 7, 12, 20)
  mu_IG <- 0.9; eta <- 1.4
  X <- sim_proc(T, list(mu_IG = mu_IG, eta = eta), "IG", mp = 100000)
  expect_equal(mean(X[, 5]), mu_IG * T[5], tolerance = 0.02)
  expect_equal(stats::var(X[, 5]), mu_IG^3 * T[5] / eta, tolerance = 0.05)
})

test_that("sim_proc produit des increments independants", {
  set.seed(1)
  T <- c(1, 3, 7, 12, 20)
  X <- sim_proc(T, list(mu_W = 0.35, sigma = 0.6), "Wiener", mp = 50000)
  dX <- X - cbind(0, X[, -ncol(X), drop = FALSE])
  cors <- stats::cor(dX)
  diag(cors) <- 0
  expect_true(all(abs(cors) < 0.02))
})

test_that("sim_traj lineaire respecte E et Var du Tableau 2", {
  set.seed(1)
  T <- c(1, 3, 7, 12, 20)
  mu <- 0.2; sd <- 0.05; s_eps <- 0.03
  X <- sim_traj(T, list(s_eps = s_eps), "lineaire", mp = 100000, gen_phi_norm(mu, sd))
  expect_equal(mean(X[, 5]), mu * T[5], tolerance = 0.02)
  expect_equal(stats::var(X[, 5]), sd^2 * T[5]^2 + s_eps^2, tolerance = 0.05)
})

test_that("sim_traj exponentiel respecte E et Var du Tableau 2", {
  set.seed(1)
  T <- c(1, 3, 7, 12, 20)
  mu <- 0.2; sd <- 0.05; s_eps <- 0.03; omega <- 0.05
  X <- sim_traj(T, list(omega = omega, s_eps = s_eps), "exponentiel", mp = 100000,
                gen_phi_norm(mu, sd))
  g <- exp(omega * T[5]) - 1
  expect_equal(mean(X[, 5]), mu * g, tolerance = 0.02)
  expect_equal(stats::var(X[, 5]), sd^2 * g^2 + s_eps^2, tolerance = 0.05)
})

test_that("sim_traj : phi_k constant le long d'une trajectoire", {
  set.seed(1)
  T <- c(1, 3, 7, 12, 20)
  X <- sim_traj(T, list(s_eps = 1e-9), "lineaire", mp = 1000, gen_phi_norm(0.2, 0.05))
  pentes <- X / matrix(T, 1000, length(T), byrow = TRUE)
  expect_true(all(apply(pentes, 1, stats::sd) < 1e-6))
})
