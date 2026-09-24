library(SMRD.data)

source("algorithme_1_test_gof_bootstrap.R")
source("algorithme_2_simulateur_processus_stochastique.R")
source("estimateurs_mle_homogenes.R")

# Un processus de Wiener homogène ajuste-t-il ces données ?

d <- gaaslaser
unites <- unique(d$unit)
heures <- sort(unique(d$hours))

X_full <- t(sapply(unites, function(k) {
  sous <- d[d$unit == k, ]
  sous$increase[order(sous$hours)]
}))

T_gaas <- heures[-1]
X_gaas <- X_full[, -1]

dim(X_gaas)
matplot(T_gaas, t(X_gaas), type = "b", pch = 1, lty = 1,
        xlab = "temps (h)", ylab = "augmentation du courant (%)")
abline(h = 10, lty = 2)


# ---- Estimation du Wiener homogene ----

th_gaas <- estimer_wiener_homogene(X_gaas, T_gaas)
th_gaas                         # article : mu_W = 2.0372e-03, sigma = 1.2657e-02


# ---- Test GoF ----

sim_wiener <- function(theta, T, m_traj){
  sim_proc(T, theta, "Wiener", m_traj)
}

set.seed(30)
res_gaas <- test_gof(X_gaas, T_gaas, estimer_wiener_homogene, sim_wiener,
                     mp = 100 * nrow(X_gaas), L = 1000)

res_gaas$rejet                  # article : rejet pour les 4 statistiques
round(res_gaas$S_obs, 4)        # article : 1.8752  6.4109  12.7108  19.9976
round(res_gaas$s_crit, 4)       # article : 0.5661  5.0354   1.4075   7.7271
round(res_gaas$p_valeur, 4)     # article : 0.001   0.022    0.001    0.006


# Un processus de Wiener non homogène ajuste-t-il ces données ?

estimer_wiener_nonhomog <- function(X, T) {
  m <- nrow(X); n <- ncol(X)
  dt <- diff(c(0, T))
  dX <- X - cbind(0, X[, -n, drop = FALSE])
  Sj <- colSums(dX)
  
  sigma2_hat_r <- function(r) {
    dLam <- diff(c(0, T)^r)
    w <- dLam / dt
    mu_W_hat <- sum(w * Sj) / (m * sum(dLam^2 / dt))
    sum((dX - matrix(mu_W_hat * dLam, m, n, byrow = TRUE))^2 /
          matrix(dt, m, n, byrow = TRUE)) / (m * n)
  }
  
  opt <- optimize(sigma2_hat_r, interval = c(0.05, 10))
  r_hat <- opt$minimum
  dLam <- diff(c(0, T)^r_hat)
  w <- dLam / dt
  mu_W_hat <- sum(w * Sj) / (m * sum(dLam^2 / dt))
  list(mu_W = mu_W_hat, sigma = sqrt(sigma2_hat_r(r_hat)), r = r_hat)
}

th_gaas_nh <- estimer_wiener_nonhomog(X_gaas, T_gaas)
th_gaas_nh


#Estimation des paramètres

set.seed(30)
T_verif <- 1:50
th_vrai <- list(mu_W = 0.15, sigma = 0.4, r = 1.7)
X_verif <- sim_proc(T_verif, th_vrai, "Wiener", mp = 20000, Lam = function(t) t^th_vrai$r)

th_estime <- estimer_wiener_nonhomog(X_verif, T_verif)
th_estime


# ---- Gamma et IG homogenes sur les donnees GaAs ----

sim_gamma <- function(th, T, m_traj) sim_proc(T, th, "Gamma", m_traj)
res_gamma_gaas <- test_gof(X_gaas, T_gaas, estimer_gamma_homogene,
                           sim_gamma, mp = 100 * nrow(X_gaas), L = 1000)

sim_ig <- function(th, T, m_traj) sim_proc(T, th, "IG", m_traj)
res_ig_gaas <- test_gof(X_gaas, T_gaas, estimer_ig_homogene, sim_ig,
                        mp = 100 * nrow(X_gaas), L = 1000)

res_gamma_gaas$rejet
res_ig_gaas$rejet
