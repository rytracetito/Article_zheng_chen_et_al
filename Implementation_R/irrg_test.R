library(gofdegrad)
library(ggplot2)

set.seed(1)
T <- round(cumsum(runif(50, 0.5, 1.5)), 2)
range(diff(c(0, T)))

h  <- 1000
mp <- 600
L  <- 200

tracer <- function(X, titre) {
  df <- data.frame(t = rep(T, each = nrow(X)),
                   x = as.vector(X),
                   traj = factor(rep(seq_len(nrow(X)), times = length(T))))
  ggplot(df, aes(t, x, group = traj)) +
    geom_line(colour = "#2563EB", linewidth = 0.5, alpha = 0.8) +
    geom_rug(data = data.frame(t = T), aes(x = t), inherit.aes = FALSE,
             sides = "b", colour = "grey50") +
    labs(title = titre, x = "Temps", y = "Dégradation") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
}


# ---- Algorithme 2 : sim_proc ----

X_proc <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 10)
dim(X_proc)
tracer(X_proc, "Wiener homogène")

X_proc_gamma <- sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 10)
dim(X_proc_gamma)
tracer(X_proc_gamma, "Gamma homogène")

X_proc_ig <- sim_proc(T, list(mu_IG = 0.5, eta = 0.5), "IG", mp = 10)
dim(X_proc_ig)
tracer(X_proc_ig, "Inverse Gaussien homogène")

X_proc_nhwp <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 10, Lam = function(t) t^1.5)
dim(X_proc_nhwp)
tracer(X_proc_nhwp, "Wiener non homogène (r = 1.5)")

X_proc_nhgp <- sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 10, Lam = function(t) t^1.5)
dim(X_proc_nhgp)
tracer(X_proc_nhgp, "Gamma non homogène (r = 1.5)")


# ---- Algorithme 3 : sim_traj ----

X_gen <- sim_traj(T, list(omega = 0.05, s_eps = 0.05), "exponentiel", mp = 10,
                  gen_phi_norm(mu = 0.2, sd = 0.05))
dim(X_gen)
tracer(X_gen, "Trajectoire générale exponentielle")

X_gen_lin <- sim_traj(T, list(s_eps = 0.05), "lineaire", mp = 10,
                      gen_phi_norm(mu = 0.2, sd = 0.05))
dim(X_gen_lin)
tracer(X_gen_lin, "Trajectoire générale linéaire")

X_gen_puiss <- sim_traj(T, list(s_eps = 0.05), "loi_puissance", mp = 10,
                        gen_phi_norm(mu = c(0.5, 0.5), sd = c(0.1, 0.1)))
dim(X_gen_puiss)
tracer(X_gen_puiss, "Trajectoire générale loi de puissance")


# ---- Estimateurs MLE ----

th1 <- list(mu_W = 0.3, sigma = 0.08)
X1 <- sim_proc(T, th1, "Wiener", mp = 20000)
estimer_parametres(X1, T, "Wiener")

th2 <- list(alpha = 1.5, beta = 0.7)
X2 <- sim_proc(T, th2, "Gamma", mp = 20000)
estimer_parametres(X2, T, "Gamma")

th3 <- list(mu_IG = 0.8, eta = 1.2)
X3 <- sim_proc(T, th3, "IG", mp = 20000)
estimer_parametres(X3, T, "IG")

th4 <- list(mu_W = 0.3, sigma = 0.08, r = 1.5)
X4 <- sim_proc(T, th4, "Wiener", mp = 20000, Lam = function(t) t^th4$r)
estimer_parametres(X4, T, "Wiener_nonhomog")

th5 <- list(alpha = 1.5, beta = 0.7, r = 1.5)
X5 <- sim_proc(T, th5, "Gamma", mp = 20000, Lam = function(t) t^th5$r)
estimer_parametres(X5, T, "Gamma_nonhomog")


# ---- Algorithme 1 : test_gof ----

## Application 1 - H0 : HWP

set.seed(123)
X_wiener <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
resultat_wiener <- test_gof(X_wiener, T, est_modele("Wiener"), sim_modele("Wiener"), mp = mp, L = L)
resultat_wiener$rejet
resultat_wiener$p_valeur

X_h1_nhwp <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6, Lam = function(t) t^1.5)
resultat_h1_nhwp <- test_gof(X_h1_nhwp, T, est_modele("Wiener"), sim_modele("Wiener"), mp = mp, L = L)
resultat_h1_nhwp$rejet
resultat_h1_nhwp$p_valeur

X_h1_hgp <- sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 6) +
  matrix(rnorm(6 * length(T), 0, 0.05), 6, length(T))
resultat_h1_hgp <- test_gof(X_h1_hgp, T, est_modele("Wiener"), sim_modele("Wiener"), mp = mp, L = L)
resultat_h1_hgp$rejet
resultat_h1_hgp$p_valeur

X_h1_ig <- sim_proc(T, list(mu_IG = 0.5, eta = 0.5), "IG", mp = 6) +
  matrix(rnorm(6 * length(T), 0, 0.05), 6, length(T))
resultat_h1_ig <- test_gof(X_h1_ig, T, est_modele("Wiener"), sim_modele("Wiener"), mp = mp, L = L)
resultat_h1_ig$rejet
resultat_h1_ig$p_valeur

X_lineaire <- sim_traj(T, list(s_eps = 0.05), "lineaire", mp = 6, gen_phi_norm(mu = 0.2, sd = 0.05))
resultat_lineaire <- test_gof(X_lineaire, T, est_modele("Wiener"), sim_modele("Wiener"), mp = mp, L = L)
resultat_lineaire$rejet


## Application 2 - H0 : HGP

X_gamma <- sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 6)
resultat_gamma <- test_gof(X_gamma, T, est_modele("Gamma"), sim_modele("Gamma"), mp = mp, L = L)
resultat_gamma$rejet
resultat_gamma$p_valeur

X_h1_hwp <- sim_proc(T, list(mu_W = 0.2, sigma = 5), "Wiener", mp = 6)
resultat_h1_hwp <- tryCatch(
  test_gof(X_h1_hwp, T, est_modele("Gamma"), sim_modele("Gamma"), mp = mp, L = L),
  error = function(e) conditionMessage(e))
resultat_h1_hwp

X_h1_nhwp <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6, Lam = function(t) t^1.5)
resultat_h1_nhwp <- test_gof(X_h1_nhwp, T, est_modele("Gamma"), sim_modele("Gamma"), mp = mp, L = L)
resultat_h1_nhwp$rejet

X_h1_ig <- sim_proc(T, list(mu_IG = 1.5, eta = 0.1), "IG", mp = 6)
resultat_h1_ig <- test_gof(X_h1_ig, T, est_modele("Gamma"), sim_modele("Gamma"), mp = mp, L = L)
resultat_h1_ig$rejet


## Application 3 - H0 : IGP

X_ig <- sim_proc(T, list(mu_IG = 0.5, eta = 0.5), "IG", mp = 6)
resultat_ig <- test_gof(X_ig, T, est_modele("IG"), sim_modele("IG"), mp = mp, L = L)
resultat_ig$rejet
resultat_ig$p_valeur

X_h1_hwp <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
resultat_h1_hwp <- test_gof(X_h1_hwp, T, est_modele("IG"), sim_modele("IG"), mp = mp, L = L)
resultat_h1_hwp$rejet

X_h1_hgp <- sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 6)
resultat_h1_hgp <- test_gof(X_h1_hgp, T, est_modele("IG"), sim_modele("IG"), mp = mp, L = L)
resultat_h1_hgp$rejet


## Application 4 - H0 : NHWP

X_nhwp <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6, Lam = function(t) t^2)
resultat_nhwp <- test_gof(X_nhwp, T, est_modele("Wiener_nonhomog"), sim_modele("Wiener_nonhomog"), mp = mp, L = L)
resultat_nhwp$rejet
resultat_nhwp$p_valeur

X_h1_hwp <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
resultat_h1_hwp <- test_gof(X_h1_hwp, T, est_modele("Wiener_nonhomog"), sim_modele("Wiener_nonhomog"), mp = mp, L = L)
resultat_h1_hwp$rejet
resultat_h1_hwp$p_valeur

X_h1_lin <- sim_traj(T, list(s_eps = 0.05), "lineaire", mp = 6, gen_phi_norm(mu = 0.2, sd = 0.05))
resultat_h1_lin <- test_gof(X_h1_lin, T, est_modele("Wiener_nonhomog"), sim_modele("Wiener_nonhomog"), mp = mp, L = L)
resultat_h1_lin$rejet

X_h1_puiss <- sim_traj(T, list(s_eps = 0.05), "loi_puissance", mp = 6,
                       gen_phi_norm(mu = c(0.5, 0.5), sd = c(0.1, 0.1)))
resultat_h1_puiss <- test_gof(X_h1_puiss, T, est_modele("Wiener_nonhomog"), sim_modele("Wiener_nonhomog"), mp = mp, L = L)
resultat_h1_puiss$rejet


## Application 5 - H0 : NHGP

X_nhgp <- sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 6, Lam = function(t) t^4)
resultat_nhgp <- test_gof(X_nhgp, T, est_modele("Gamma_nonhomog"), sim_modele("Gamma_nonhomog"), mp = mp, L = L)
resultat_nhgp$rejet
resultat_nhgp$p_valeur

X_h1_hgp <- sim_proc(T, list(alpha = 0.1, beta = 3), "Gamma", mp = 6)
resultat_h1_hgp <- test_gof(X_h1_hgp, T, est_modele("Gamma_nonhomog"), sim_modele("Gamma_nonhomog"), mp = mp, L = L)
resultat_h1_hgp$rejet


# ---- Niveau et puissance ----

gen_hwp  <- function() sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
gen_hgp  <- function() sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 6)
gen_igp  <- function() sim_proc(T, list(mu_IG = 0.5, eta = 0.5), "IG", mp = 6)
gen_nhwp <- function() sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6, Lam = function(t) t^1.3)
gen_nhgp <- function() sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 6, Lam = function(t) t^1.3)

bruit <- function() matrix(rnorm(6 * length(T), 0, 0.05), 6, length(T))
gen_hgp_bruit <- function() gen_hgp() + bruit()
gen_igp_bruit <- function() gen_igp() + bruit()

gen_lineaire <- function() sim_traj(T, list(s_eps = 0.05), "lineaire", mp = 6,
                                    gen_phi_norm(mu = 0.2, sd = 0.05))
gen_puiss    <- function() sim_traj(T, list(s_eps = 0.05), "loi_puissance", mp = 6,
                                    gen_phi_norm(mu = c(0.5, 0.5), sd = c(0.1, 0.1)))

etude <- function(H0, scenarios) {
  t(sapply(scenarios, function(gen)
    etude_puissance(h, gen, T, est_modele(H0), sim_modele(H0), mp = mp, L = L)))
}

## H0 : HWP
etude("Wiener", list(niveau = gen_hwp, HGP_bruit = gen_hgp_bruit, IGP_bruit = gen_igp_bruit,
                     NHWP = gen_nhwp, lineaire = gen_lineaire))

## H0 : HGP
etude("Gamma", list(niveau = gen_hgp, HWP = gen_hwp, NHWP = gen_nhwp, IGP = gen_igp))

## H0 : IGP
etude("IG", list(niveau = gen_igp, HWP = gen_hwp, HGP = gen_hgp))

## H0 : NHWP
etude("Wiener_nonhomog", list(niveau = gen_nhwp, HGP = gen_hgp, IGP = gen_igp,
                              lineaire = gen_lineaire, puissance = gen_puiss))

## H0 : NHGP
etude("Gamma_nonhomog", list(niveau = gen_nhgp, HWP = gen_hwp, IGP = gen_igp, lineaire = gen_lineaire))
