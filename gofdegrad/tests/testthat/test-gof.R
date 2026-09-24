test_that("dist_mat reproduit les equations 3 et 4 (transcription naive)", {
  set.seed(1)
  m <- 5; mp <- 13; n <- 8
  X  <- matrix(stats::rnorm(m * n), m, n)
  Xs <- matrix(stats::rnorm(mp * n), mp, n)

  D_KS_naif  <- matrix(0, m, mp)
  D_CvM_naif <- matrix(0, m, mp)
  for (i in seq_len(m)) for (k in seq_len(mp)) {
    D_KS_naif[i, k]  <- max(abs(X[i, ] - Xs[k, ]))
    D_CvM_naif[i, k] <- sum((X[i, ] - Xs[k, ])^2)
  }

  expect_equal(gofdegrad:::dist_mat(X, Xs, "KS"),  D_KS_naif)
  expect_equal(gofdegrad:::dist_mat(X, Xs, "CvM"), D_CvM_naif)
})

test_that("stats_sim : S_sup == S_sum quand m == 1 (eq. 8 et 11)", {
  set.seed(1)
  X1 <- matrix(stats::rnorm(8), 1, 8)
  Xs <- matrix(stats::rnorm(13 * 8), 13, 8)
  s <- stats_sim(X1, Xs)
  expect_equal(unname(s["KS_sup"]),  unname(s["KS_sum"]))
  expect_equal(unname(s["CvM_sup"]), unname(s["CvM_sum"]))
})

test_that("test_gof : decision coherente avec la p-valeur", {
  set.seed(1)
  T <- 1:30
  sim_w <- function(theta, T, m_traj) sim_proc(T, theta, "Wiener", m_traj)
  X <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
  res <- test_gof(X, T, estimer_wiener_homogene, sim_w, mp = 300, L = 100)
  expect_equal(res$rejet, res$p_valeur <= 0.05)
})

test_that("test_gof : niveau empirique proche de alpha sous H0 (Wiener)", {
  skip_on_cran()
  set.seed(42)
  T <- 1:50
  sim_w <- function(theta, T, m_traj) sim_proc(T, theta, "Wiener", m_traj)
  gen_h0 <- function() sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
  niveau <- etude_puissance(h = 80, gen_h0, T, estimer_wiener_homogene, sim_w, mp = 300, L = 100)
  expect_true(all(niveau < 0.15))  # large marge : h modeste, juste garde-fou contre une derive grossiere
})

test_that("test_gof rejette une alternative nette (lineaire vs H0 Wiener)", {
  set.seed(1)
  T <- 1:50
  sim_w <- function(theta, T, m_traj) sim_proc(T, theta, "Wiener", m_traj)
  phi <- gen_phi_norm(mu = 0.2, sd = 0.05)
  X <- sim_traj(T, list(s_eps = 0.05), "lineaire", mp = 6, phi)
  res <- test_gof(X, T, estimer_wiener_homogene, sim_w, mp = 600, L = 100)
  expect_true(all(res$rejet))
})
