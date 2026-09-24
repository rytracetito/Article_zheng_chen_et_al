grille_irreguliere <- function() {
  set.seed(1)
  round(cumsum(runif(50, 0.5, 1.5)), 2)
}

test_that("sim_proc respecte la variance sigma^2 * dt sur une grille irreguliere", {
  T <- grille_irreguliere()
  dt <- diff(c(0, T))
  set.seed(2)
  X <- sim_proc(T, list(mu_W = 0.3, sigma = 0.5), "Wiener", mp = 20000)
  dX <- X - cbind(0, X[, -ncol(X)])
  expect_equal(unname(apply(dX, 2, var)) / dt, rep(0.25, length(T)), tolerance = 0.05)
  expect_equal(unname(colMeans(dX)) / dt, rep(0.3, length(T)), tolerance = 0.05)
})

test_that("les estimateurs recuperent les parametres sur une grille irreguliere", {
  T <- grille_irreguliere()
  set.seed(3)

  est <- estimer_parametres(sim_proc(T, list(mu_W = 0.3, sigma = 0.08), "Wiener", 5000), T, "Wiener")
  expect_equal(unlist(est), c(mu_W = 0.3, sigma = 0.08), tolerance = 0.02)

  est <- estimer_parametres(sim_proc(T, list(alpha = 1.5, beta = 0.7), "Gamma", 5000), T, "Gamma")
  expect_equal(unlist(est), c(alpha = 1.5, beta = 0.7), tolerance = 0.03)

  est <- estimer_parametres(sim_proc(T, list(mu_IG = 0.8, eta = 1.2), "IG", 5000), T, "IG")
  expect_equal(unlist(est), c(mu_IG = 0.8, eta = 1.2), tolerance = 0.03)

  X <- sim_proc(T, list(mu_W = 0.15, sigma = 0.4), "Wiener", 5000, Lam = function(t) t^1.7)
  est <- estimer_parametres(X, T, "Wiener_nonhomog")
  expect_equal(unlist(est), c(mu_W = 0.15, sigma = 0.4, r = 1.7), tolerance = 0.03)

  X <- sim_proc(T, list(alpha = 0.8, beta = 0.5), "Gamma", 5000, Lam = function(t) t^1.4)
  est <- estimer_parametres(X, T, "Gamma_nonhomog")
  expect_equal(unlist(est), c(alpha = 0.8, beta = 0.5, r = 1.4), tolerance = 0.03)
})

test_that("un increment nul est remplace avec un avertissement, un negatif reste une erreur", {
  T <- 1:5
  X_nul <- matrix(c(1, 2, 2, 3, 4), nrow = 1)
  expect_warning(estimer_gamma_homogene(X_nul, T), "1 increment\\(s\\) nul\\(s\\)")
  expect_warning(estimer_ig_homogene(X_nul, T), "increment\\(s\\) nul\\(s\\)")
  expect_warning(estimer_gamma_nonhomog(X_nul, T), "increment\\(s\\) nul\\(s\\)")

  X_neg <- matrix(c(1, 2, 1.5, 3, 4), nrow = 1)
  expect_error(estimer_gamma_homogene(X_neg, T), "increments negatifs")
  expect_error(estimer_gamma_nonhomog(X_neg, T), "increments negatifs")
})

test_that("test_gof aboutit quand le bootstrap Gamma produit des increments nuls", {
  T <- grille_irreguliere()
  set.seed(7)
  X <- sim_proc(T, list(mu_IG = 1.5, eta = 0.1), "IG", mp = 6)
  res <- test_gof(X, T, est_modele("Gamma"), sim_modele("Gamma"), mp = 600, L = 50)
  expect_true(all(is.finite(res$S_obs)))
  expect_true(all(res$p_valeur >= 0 & res$p_valeur <= 1))
})
