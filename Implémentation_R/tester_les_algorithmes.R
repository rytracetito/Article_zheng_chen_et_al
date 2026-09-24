source("algorithme_1_test_gof_bootstrap.R")
source("algorithme_2_simulateur_processus_stochastique.R")
source("algorithme_3_simulateur_trajectoire_generale.R")
source("estimateurs_mle_homogenes.R")

set.seed(1)
T <- 1:50


sim_modele <- function(modele) {
  switch(modele,
         Wiener_nonhomog = function(theta, T, m_traj) sim_proc(T, theta, "Wiener", m_traj, Lam = function(t) t^theta$r),
         Gamma_nonhomog  = function(theta, T, m_traj) sim_proc(T, theta, "Gamma",  m_traj, Lam = function(t) t^theta$r),
         function(theta, T, m_traj) sim_proc(T, theta, modele, m_traj)
  )
}

est_modele <- function(modele) function(X, T) estimer_parametres(X, T, modele)


# ---- Algorithme 2 : sim_proc ----

X_proc <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 10)
dim(X_proc)
matplot(t(X_proc), type = "l")

X_proc_gamma <- sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 10)
dim(X_proc_gamma)
matplot(t(X_proc_gamma), type = "l")

X_proc_ig <- sim_proc(T, list(mu_IG = 0.5, eta = 0.5), "IG", mp = 10)
dim(X_proc_ig)
matplot(t(X_proc_ig), type = "l")

X_proc_nhwp <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 10, Lam = function(t) t^1.5)
dim(X_proc_nhwp)
matplot(t(X_proc_nhwp), type = "l")

X_proc_nhgp <- sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 10, Lam = function(t) t^1.5)
dim(X_proc_nhgp)
matplot(t(X_proc_nhgp), type = "l")

# ---- Algorithme 3 : sim_traj ----

th_traj <- list(omega = 0.05, s_eps = 0.05)
phi_ech <- gen_phi_norm(mu = 0.2, sd = 0.05)
X_gen <- sim_traj(T, th_traj, "exponentiel", mp = 10, phi_ech)
dim(X_gen)
matplot(t(X_gen), type = "l")

phi_lin <- gen_phi_norm(mu = 0.2, sd = 0.05)
X_gen_lin <- sim_traj(T, list(s_eps = 0.05), "lineaire", mp = 10, phi_lin)
dim(X_gen_lin)
matplot(t(X_gen_lin), type = "l")

phi_puiss <- gen_phi_norm(mu = c(0.5, 0.5), sd = c(0.1, 0.1))
X_gen_puiss <- sim_traj(T, list(s_eps = 0.05), "loi_puissance", mp = 10, phi_puiss)
dim(X_gen_puiss)
matplot(t(X_gen_puiss), type = "l")


# ---- Estimateurs MLE : retrouver les parametres theoriques ----

# HWP
th1 <- list(mu_W = 0.3, sigma = 0.08)
X1 <- sim_proc(T, th1, "Wiener", mp = 20000)
estimer_parametres(X1, T, "Wiener")   # attendu : proche de mu_W=0.3, sigma=0.08

# HGP
th2 <- list(alpha = 1.5, beta = 0.7)
X2 <- sim_proc(T, th2, "Gamma", mp = 20000)
estimer_parametres(X2, T, "Gamma")  # attendu : proche de alpha=1.5, beta=0.7

# IGP (Inverse Gaussien)
th3 <- list(mu_IG = 0.8, eta = 1.2)
X3 <- sim_proc(T, th3, "IG", mp = 20000)
estimer_parametres(X3, T, "IG")    # attendu : proche de mu_IG=0.8, eta=1.2


# NHWP
th4 <- list(mu_W = 0.3, sigma = 0.08, r = 1.5)
X4 <- sim_proc(T, th4, "Wiener", mp = 20000, Lam = function(t) t^th4$r)
estimer_parametres(X4, T, "Wiener_nonhomog")   # attendu : proche de mu_W=0.3, sigma=0.08, r=1.5

# NHGP
th5 <- list(alpha = 1.5, beta = 0.7, r = 1.5)
X5 <- sim_proc(T, th5, "Gamma", mp = 20000, Lam = function(t) t^th5$r)
estimer_parametres(X5, T, "Gamma_nonhomog")    # attendu : proche de alpha=1.5, beta=0.7, r=1.5

# ---- Algorithme 1 : test_gof ----

## Application 1
## H0 : HWP
set.seed(123)
X_wiener <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
resultat_wiener <- test_gof(X_wiener, T, est_modele("Wiener"), sim_modele("Wiener"),
                            mp = 100 * nrow(X_wiener), L = 200)
resultat_wiener$rejet     # attendu : FALSE (les donnees suivent bien H0 ici)
resultat_wiener$p_valeur


## H1 : NHWP (donnees non-homogenes, testees contre H0 Wiener homogene)

X_h1_nhwp <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6, Lam = function(t) t^1.5)
resultat_h1_nhwp <- test_gof(X_h1_nhwp, T, est_modele("Wiener"), sim_modele("Wiener"),
                             mp = 600, L = 200)
resultat_h1_nhwp$rejet    # attendu : TRUE 
resultat_h1_nhwp$p_valeur

## H1 : HGP (donnees Gamma + bruit, testees contre H0 Wiener)

X_h1_hgp <- sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 6) +
  matrix(rnorm(6 * length(T), 0, 0.05), 6, length(T))
resultat_h1_hgp <- test_gof(X_h1_hgp, T, est_modele("Wiener"), sim_modele("Wiener"),
                            mp = 600, L = 200)
resultat_h1_hgp$rejet     # attendu : TRUE parfois seulement (ecart discret, Table 3 de l'article)
resultat_h1_hgp$p_valeur


## H1 : IG (donnees IG + bruit, testees contre H0 Wiener)

X_h1_ig <- sim_proc(T, list(mu_IG = 0.5, eta = 0.5), "IG", mp = 6) +
  matrix(rnorm(6 * length(T), 0, 0.05), 6, length(T))
resultat_h1_ig <- test_gof(X_h1_ig, T, est_modele("Wiener"), sim_modele("Wiener"),
                           mp = 600, L = 200)
resultat_h1_ig$rejet
resultat_h1_ig$p_valeur

## H1 : modele lineaire general 

phi_lineaire <- gen_phi_norm(mu = 0.2, sd = 0.05)
X_lineaire <- sim_traj(T, list(omega = NA, s_eps = 0.05), "lineaire", mp = 6, phi_lineaire)
resultat_lineaire <- test_gof(X_lineaire, T, est_modele("Wiener"), sim_modele("Wiener"),
                              mp = 600, L = 200)
resultat_lineaire$rejet






## Application 2

## H0 : HGP
X_gamma <- sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 6)
resultat_gamma <- test_gof(X_gamma, T, est_modele("Gamma"), sim_modele("Gamma"),
                           mp = 100 * nrow(X_gamma), L = 200)
resultat_gamma$rejet
resultat_gamma$p_valeur


## H1 : HWP 
X_h1_hwp <- sim_proc(T, list(mu_W = 0.2, sigma = 5), "Wiener", mp = 6)
resultat_h1_hwp <- test_gof(X_h1_hwp, T, est_modele("Gamma"), sim_modele("Gamma"),
                            mp = 600, L = 200)
resultat_h1_hwp$rejet


## H1 : NHWP 
X_h1_nhwp <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6, Lam = function(t) t^1.5)
resultat_h1_nhwp <- test_gof(X_h1_nhwp, T, est_modele("Gamma"), sim_modele("Gamma"),
                             mp = 600, L = 200)
resultat_h1_nhwp$rejet


## H1 : IG 
X_h1_ig <- sim_proc(T, list(mu_IG = 1.5, eta = 0.1), "IG", mp = 6)
resultat_h1_ig <- test_gof(X_h1_ig, T, est_modele("Gamma"), sim_modele("Gamma"),
                           mp = 600, L = 200)
resultat_h1_ig$rejet


# Application 3 

## H0 : IGP
T <- 1:50
X_ig <- sim_proc(T, list(mu_IG = 0.5, eta = 0.5), "IG", mp = 6)
resultat_ig <- test_gof(X_ig, T, est_modele("IG"), sim_modele("IG"),
                        mp = 100 * nrow(X_ig), L = 200)
resultat_ig$rejet
resultat_ig$p_valeur


## H1 : HWP (donnees Wiener, testees contre H0 IG)

X_h1_hwp <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
resultat_h1_hwp <- test_gof(X_h1_hwp, T, est_modele("IG"), sim_modele("IG"),
                            mp = 600, L = 200)
resultat_h1_hwp$rejet


## H1 : HGP (donnees Gamma, testees contre H0 IG)

X_h1_hgp <- sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 6)
resultat_h1_hgp <- test_gof(X_h1_hgp, T, est_modele("IG"), sim_modele("IG"),
                            mp = 600, L = 200)
resultat_h1_hgp$rejet




##Aplication 4

## H0 : NHWP
X_nhwp <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6, Lam = function(t) t^2)
resultat_nhwp <- test_gof(X_nhwp, T, est_modele("Wiener_nonhomog"), sim_modele("Wiener_nonhomog"),
                          mp = 100 * nrow(X_nhwp), L = 200)
resultat_nhwp$rejet
resultat_nhwp$p_valeur

## H1 : HWP 
X_h1_hwp <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
resultat_h1_hwp <- test_gof(X_h1_hwp, T, est_modele("Wiener_nonhomog"), sim_modele("Wiener_nonhomog"),
                            mp = 600, L = 200)
resultat_h1_hwp$rejet     
resultat_h1_hwp$p_valeur


## H1 : lineaire modele a trajectoire generale
phi_lin <- gen_phi_norm(mu = 0.2, sd = 0.05)
X_h1_lin <- sim_traj(T, list(s_eps = 0.05), "lineaire", mp = 6, phi_lin)
resultat_h1_lin <- test_gof(X_h1_lin, T, est_modele("Wiener_nonhomog"), sim_modele("Wiener_nonhomog"),
                            mp = 600, L = 200)
resultat_h1_lin$rejet


## H1 : loi de puissance modele a trajectoire generale
phi_puiss <- gen_phi_norm(mu = c(0.5, 0.5), sd = c(0.1, 0.1))
X_h1_puiss <- sim_traj(T, list(s_eps = 0.05), "loi_puissance", mp = 6, phi_puiss)
resultat_h1_puiss <- test_gof(X_h1_puiss, T, est_modele("Wiener_nonhomog"), sim_modele("Wiener_nonhomog"),
                              mp = 600, L = 200)
resultat_h1_puiss$rejet

## Application 5
## H0 : NHGP

X_nhgp <- sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 6, Lam = function(t) t^4)
resultat_nhgp <- test_gof(X_nhgp, T, est_modele("Gamma_nonhomog"), sim_modele("Gamma_nonhomog"),
                          mp = 100 * nrow(X_nhgp), L = 200)
resultat_nhgp$rejet
resultat_nhgp$p_valeur


## H1 : HGP 
X_h1_hgp <- sim_proc(T, list(alpha = 0.1, beta = 3), "Gamma", mp = 6)
resultat_h1_hgp <- test_gof(X_h1_hgp, T, est_modele("Gamma_nonhomog"), sim_modele("Gamma_nonhomog"),
                            mp = 600, L = 200)
resultat_h1_hgp$rejet     



# ---- Niveau et puissance des tests (plusieurs modeles) ----

# h     : nb de repetitions
# gen_X : function() -> jeu de donnees

etude_puissance <- function(h, gen_X, T, est, sim, mp, L, alpha = 0.05) {
  rejets <- vector("list", h)
  for (r in seq_len(h)) {
    rejets[[r]] <- tryCatch(
      test_gof(gen_X(), T, est, sim, mp = mp, L = L, alpha = alpha)$rejet,
      error = function(e) c(KS_sup = TRUE, KS_sum = TRUE, CvM_sup = TRUE, CvM_sum = TRUE)
    )
  }
  round(colMeans(do.call(rbind, rejets)), 4)   # niveau (si gen_X simule sous H0) ou puissance (sinon)
}


# ---- Generateurs de donnees ----

gen_hwp  <- function() sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
gen_hgp  <- function() sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 6)
gen_igp  <- function() sim_proc(T, list(mu_IG = 0.5, eta = 0.5), "IG", mp = 6)
gen_nhwp <- function() sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6, Lam = function(t) t^1.3)
gen_nhgp <- function() sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 6, Lam = function(t) t^1.3)

gen_hgp_bruit <- function() {
  sim_proc(T, list(alpha = 1, beta = 1/3), "Gamma", mp = 6) +
    matrix(rnorm(6 * length(T), 0, 0.05), 6, length(T))
}
gen_igp_bruit <- function() {
  sim_proc(T, list(mu_IG = 0.5, eta = 0.5), "IG", mp = 6) +
    matrix(rnorm(6 * length(T), 0, 0.05), 6, length(T))
}
gen_lineaire <- function() sim_traj(T, list(s_eps = 0.05), "lineaire", mp = 6, gen_phi_norm(mu = 0.2, sd = 0.05))
gen_puiss    <- function() sim_traj(T, list(s_eps = 0.05), "loi_puissance", mp = 6, gen_phi_norm(mu = c(0.5, 0.5), sd = c(0.1, 0.1)))


# ================= H0 : HWP =================

etude_puissance(1000, gen_hwp,       T, est_modele("Wiener"), sim_modele("Wiener"), mp = 600, L = 200)  # Niveau
etude_puissance(1000, gen_hgp_bruit, T, est_modele("Wiener"), sim_modele("Wiener"), mp = 600, L = 200)  # vs HGP+bruit
etude_puissance(1000, gen_igp_bruit, T, est_modele("Wiener"), sim_modele("Wiener"), mp = 600, L = 200)  # vs IGP+bruit
etude_puissance(1000, gen_nhwp,      T, est_modele("Wiener"), sim_modele("Wiener"), mp = 600, L = 200)  # vs NHWP
etude_puissance(1000, gen_lineaire,  T, est_modele("Wiener"), sim_modele("Wiener"), mp = 600, L = 200)  # vs lineaire


# ================= H0 : HGP =================

etude_puissance(1000, gen_hgp,  T, est_modele("Gamma"), sim_modele("Gamma"), mp = 600, L = 200)  # Niveau
etude_puissance(1000, gen_hwp,  T, est_modele("Gamma"), sim_modele("Gamma"), mp = 600, L = 200)  # vs HWP
etude_puissance(1000, gen_nhwp, T, est_modele("Gamma"), sim_modele("Gamma"), mp = 600, L = 200)  # vs NHWP
etude_puissance(1000, gen_igp,  T, est_modele("Gamma"), sim_modele("Gamma"), mp = 600, L = 200)  # vs IGP


# ================= H0 : IGP =================

etude_puissance(1000, gen_igp, T, est_modele("IG"), sim_modele("IG"), mp = 600, L = 200)  # Niveau
etude_puissance(1000, gen_hwp, T, est_modele("IG"), sim_modele("IG"), mp = 600, L = 200)  # vs HWP
etude_puissance(1000, gen_hgp, T, est_modele("IG"), sim_modele("IG"), mp = 600, L = 200)  # vs HGP


# ================= H0 : NHWP =================

etude_puissance(1000, gen_nhwp,     T, est_modele("Wiener_nonhomog"), sim_modele("Wiener_nonhomog"), mp = 600, L = 200)  # Niveau
etude_puissance(1000, gen_hgp,      T, est_modele("Wiener_nonhomog"), sim_modele("Wiener_nonhomog"), mp = 600, L = 200)  # vs HGP
etude_puissance(1000, gen_igp,      T, est_modele("Wiener_nonhomog"), sim_modele("Wiener_nonhomog"), mp = 600, L = 200)  # vs IGP
etude_puissance(1000, gen_lineaire, T, est_modele("Wiener_nonhomog"), sim_modele("Wiener_nonhomog"), mp = 600, L = 200)  # vs lineaire
etude_puissance(1000, gen_puiss,    T, est_modele("Wiener_nonhomog"), sim_modele("Wiener_nonhomog"), mp = 600, L = 200)  # vs loi de puissance


# ================= H0 : NHGP =================

etude_puissance(1000, gen_nhgp,     T, est_modele("Gamma_nonhomog"), sim_modele("Gamma_nonhomog"), mp = 600, L = 200)  # Niveau
etude_puissance(1000, gen_hwp,      T, est_modele("Gamma_nonhomog"), sim_modele("Gamma_nonhomog"), mp = 600, L = 200)  # vs HWP
etude_puissance(1000, gen_igp,      T, est_modele("Gamma_nonhomog"), sim_modele("Gamma_nonhomog"), mp = 600, L = 200)  # vs IGP
etude_puissance(1000, gen_lineaire, T, est_modele("Gamma_nonhomog"), sim_modele("Gamma_nonhomog"), mp = 600, L = 200)  # vs lineaire