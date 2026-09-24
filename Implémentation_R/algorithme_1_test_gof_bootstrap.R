# Algorithme 1 - test d'adequation par similarite et bootstrap parametrique

# X  : trajectoires reelles m x n
# Xs : trajectoires simulees m' x n
# D  : "KS" ou "CvM"

dist_mat <- function(X, Xs, D) {
  if (D == "CvM") {
    outer(rowSums(X^2), rowSums(Xs^2), "+") - 2 * X %*% t(Xs)
  } else {
    n <- ncol(X)
    Dm <- matrix(0, nrow(X), nrow(Xs))
    for (j in seq_len(n)) Dm <- pmax(Dm, abs(outer(X[, j], Xs[, j], "-")))
    Dm
  }
}

# X  : trajectoires reelles (ou bootstrap)
# Xs : pool de trajectoires simulees

stats_sim <- function(X, Xs) {
  d_KS  <- apply(dist_mat(X, Xs, "KS"),  1, min)
  d_CvM <- apply(dist_mat(X, Xs, "CvM"), 1, min)
  c(KS_sup  = max(d_KS),  KS_sum  = sum(d_KS),
    CvM_sup = max(d_CvM), CvM_sum = sum(d_CvM))
}

# X     : donnees reelles m x n
# T     : instants d'observation
# est   : function(X, T) -> theta
# sim   : function(theta, T, m_traj) -> matrice
# mp    : taille du pool simule (m')
# L     : nb de replications bootstrap
# alpha : niveau du test

test_gof <- function(X, T, est, sim, mp = 100 * nrow(X), L = 1000, alpha = 0.05) {
  m <- nrow(X)
  th <- est(X, T)
  S_obs <- stats_sim(X, sim(th, T, mp))

  S_star <- matrix(NA_real_, L, 4, dimnames = list(NULL, names(S_obs)))
  for (v in seq_len(L)) {
    Xb <- sim(th, T, m)
    th_b <- est(Xb, T)
    S_star[v, ] <- stats_sim(Xb, sim(th_b, T, mp))
  }
  s_crit <- apply(S_star, 2, quantile, probs = 1 - alpha, names = FALSE)
  list(rejet = S_obs > s_crit, S_obs = S_obs, s_crit = s_crit,
       p_valeur = colMeans(sweep(S_star, 2, S_obs, ">=")))
}
