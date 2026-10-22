library(gofdegrad)
library(ggplot2)

set.seed(1)
T <- 1:50

h  <- 1000
mp <- 600
L  <- 200


# titre : Nom du modèle
# x : data (simulation de sim_model ou de sim_path)
tracer <- function(x, titre) {
  df <- levysim::traj_df(x)
  ggplot(df, aes(t, X, group = id)) +
    geom_line(colour = "#2563EB", linewidth = 0.5, alpha = 0.8) +
    geom_rug(data = data.frame(t = unique(df$t)), aes(x = t), inherit.aes = FALSE,
             sides = "b", colour = "grey50") +
    labs(title = titre, x = "Temps", y = "Dégradation") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
}

ajouter_bruit <- function(x, sd = 0.05) {
  x[-1, -1] <- x[-1, -1] + rnorm(length(x[-1, -1]), 0, sd)
  x
}


# ---- Algorithme 2 : processus stochastiques ----

X_proc <- sim_model(T, "wiener", c(mu = 0.2, sigma = 0.05), n_traj = 10)
dim(X_proc)
tracer(X_proc, "Wiener homogène")

X_proc_gamma <- sim_model(T, "gamma", c(a = 1, b = 3), n_traj = 10)
dim(X_proc_gamma)
tracer(X_proc_gamma, "Gamma homogène")

X_proc_ig <- sim_model(T, "ig", c(mu = 0.5, eta = 0.5), n_traj = 10)
dim(X_proc_ig)
tracer(X_proc_ig, "Inverse Gaussien homogène")

X_proc_nhwp <- sim_model(T, "wiener_power", c(mu = 0.2, sigma = 0.05, r = 1.5), n_traj = 10)
dim(X_proc_nhwp)
tracer(X_proc_nhwp, "Wiener non homogène (r = 1.5)")

X_proc_nhgp <- sim_model(T, "gamma_power", c(a = 1, b = 3, r = 1.5), n_traj = 10)
dim(X_proc_nhgp)
tracer(X_proc_nhgp, "Gamma non homogène (r = 1.5)")


# ---- Algorithme 3 : trajectoires générales ----

X_gen <- sim_path(T, "exponential", phi_mean = 0.2, phi_sd = 0.05, sigma_eps = 0.05,
                  omega = 0.05, n_traj = 10)
dim(X_gen)
tracer(X_gen, "Trajectoire générale exponentielle")

X_gen_lin <- sim_path(T, "linear", phi_mean = 0.2, phi_sd = 0.05, sigma_eps = 0.05, n_traj = 10)
dim(X_gen_lin)
tracer(X_gen_lin, "Trajectoire générale linéaire")

X_gen_puiss <- sim_path(T, "power", phi_mean = c(0.5, 0.5), phi_sd = c(0.1, 0.1),
                        sigma_eps = 0.05, n_traj = 10)
dim(X_gen_puiss)
tracer(X_gen_puiss, "Trajectoire générale loi de puissance")


# ---- Estimateurs MLE ----

X1 <- sim_model(T, "wiener", c(mu = 0.3, sigma = 0.08), n_traj = 20000)
fit_model(X1, "wiener")

X2 <- sim_model(T, "gamma", c(a = 1.5, b = 1 / 0.7), n_traj = 20000)
fit_model(X2, "gamma")

X3 <- sim_model(T, "ig", c(mu = 0.8, eta = 1.2), n_traj = 20000)
fit_model(X3, "ig")

X4 <- sim_model(T, "wiener_power", c(mu = 0.3, sigma = 0.08, r = 1.5), n_traj = 20000)
fit_model(X4, "wiener_power")

X5 <- sim_model(T, "gamma_power", c(a = 1.5, b = 1 / 0.7, r = 1.5), n_traj = 20000)
fit_model(X5, "gamma_power")


# ---- Algorithme 1 : test d'adéquation ----

## Application 1 - H0 : HWP

set.seed(123)
X_wiener <- sim_model(T, "wiener", c(mu = 0.2, sigma = 0.05), n_traj = 6)
resultat_wiener <- gof_test(X_wiener, "wiener", B = L, n_pool = mp)
resultat_wiener

X_h1_nhwp <- sim_model(T, "wiener_power", c(mu = 0.2, sigma = 0.05, r = 1.5), n_traj = 6)
resultat_h1_nhwp <- gof_test(X_h1_nhwp, "wiener", B = L, n_pool = mp)
resultat_h1_nhwp

X_h1_hgp <- ajouter_bruit(sim_model(T, "gamma", c(a = 1, b = 3), n_traj = 6))
resultat_h1_hgp <- gof_test(X_h1_hgp, "wiener", B = L, n_pool = mp)
resultat_h1_hgp

X_h1_ig <- ajouter_bruit(sim_model(T, "ig", c(mu = 0.5, eta = 0.5), n_traj = 6))
resultat_h1_ig <- gof_test(X_h1_ig, "wiener", B = L, n_pool = mp)
resultat_h1_ig

X_lineaire <- sim_path(T, "linear", phi_mean = 0.2, phi_sd = 0.05, sigma_eps = 0.05, n_traj = 6)
resultat_lineaire <- gof_test(X_lineaire, "wiener", B = L, n_pool = mp)
resultat_lineaire$reject


## Application 2 - H0 : HGP

X_gamma <- sim_model(T, "gamma", c(a = 1, b = 3), n_traj = 6)
resultat_gamma <- gof_test(X_gamma, "gamma", B = L, n_pool = mp)
resultat_gamma

X_h1_hwp <- sim_model(T, "wiener", c(mu = 0.2, sigma = 5), n_traj = 6)
resultat_h1_hwp <- tryCatch(gof_test(X_h1_hwp, "gamma", B = L, n_pool = mp),
                            error = function(e) conditionMessage(e))
resultat_h1_hwp

X_h1_nhwp <- sim_model(T, "wiener_power", c(mu = 0.2, sigma = 0.05, r = 1.5), n_traj = 6)
resultat_h1_nhwp <- gof_test(X_h1_nhwp, "gamma", B = L, n_pool = mp)
resultat_h1_nhwp$reject

X_h1_ig <- sim_model(T, "ig", c(mu = 1.5, eta = 0.1), n_traj = 6)
resultat_h1_ig <- gof_test(X_h1_ig, "gamma", B = L, n_pool = mp)
resultat_h1_ig$reject


## Application 3 - H0 : IGP

X_ig <- sim_model(T, "ig", c(mu = 0.5, eta = 0.5), n_traj = 6)
resultat_ig <- gof_test(X_ig, "ig", B = L, n_pool = mp)
resultat_ig

X_h1_hwp <- sim_model(T, "wiener", c(mu = 0.2, sigma = 0.05), n_traj = 6)
resultat_h1_hwp <- gof_test(X_h1_hwp, "ig", B = L, n_pool = mp)
resultat_h1_hwp$reject

X_h1_hgp <- sim_model(T, "gamma", c(a = 1, b = 3), n_traj = 6)
resultat_h1_hgp <- gof_test(X_h1_hgp, "ig", B = L, n_pool = mp)
resultat_h1_hgp$reject


## Application 4 - H0 : NHWP

X_nhwp <- sim_model(T, "wiener_power", c(mu = 0.2, sigma = 0.05, r = 2), n_traj = 6)
resultat_nhwp <- gof_test(X_nhwp, "wiener_power", B = L, n_pool = mp)
resultat_nhwp

X_h1_hwp <- sim_model(T, "wiener", c(mu = 0.2, sigma = 0.05), n_traj = 6)
resultat_h1_hwp <- gof_test(X_h1_hwp, "wiener_power", B = L, n_pool = mp)
resultat_h1_hwp

X_h1_lin <- sim_path(T, "linear", phi_mean = 0.2, phi_sd = 0.05, sigma_eps = 0.05, n_traj = 6)
resultat_h1_lin <- gof_test(X_h1_lin, "wiener_power", B = L, n_pool = mp)
resultat_h1_lin$reject

X_h1_puiss <- sim_path(T, "power", phi_mean = c(0.5, 0.5), phi_sd = c(0.1, 0.1),
                       sigma_eps = 0.05, n_traj = 6)
resultat_h1_puiss <- gof_test(X_h1_puiss, "wiener_power", B = L, n_pool = mp)
resultat_h1_puiss$reject


## Application 5 - H0 : NHGP

X_nhgp <- sim_model(T, "gamma_power", c(a = 1, b = 3, r = 4), n_traj = 6)
resultat_nhgp <- gof_test(X_nhgp, "gamma_power", B = L, n_pool = mp)
resultat_nhgp

X_h1_hgp <- sim_model(T, "gamma", c(a = 0.1, b = 1 / 3), n_traj = 6)
resultat_h1_hgp <- gof_test(X_h1_hgp, "gamma_power", B = L, n_pool = mp)
resultat_h1_hgp$reject


# ---- Niveau et puissance ----

gen_hwp  <- function() sim_model(T, "wiener", c(mu = 0.2, sigma = 0.05), n_traj = 6)
gen_hgp  <- function() sim_model(T, "gamma", c(a = 1, b = 3), n_traj = 6)
gen_igp  <- function() sim_model(T, "ig", c(mu = 0.5, eta = 0.5), n_traj = 6)
gen_nhwp <- function() sim_model(T, "wiener_power", c(mu = 0.2, sigma = 0.05, r = 1.3), n_traj = 6)
gen_nhgp <- function() sim_model(T, "gamma_power", c(a = 1, b = 3, r = 1.3), n_traj = 6)

gen_hgp_bruit <- function() ajouter_bruit(gen_hgp())
gen_igp_bruit <- function() ajouter_bruit(gen_igp())

gen_lineaire <- function() sim_path(T, "linear", phi_mean = 0.2, phi_sd = 0.05,
                                    sigma_eps = 0.05, n_traj = 6)
gen_puiss    <- function() sim_path(T, "power", phi_mean = c(0.5, 0.5), phi_sd = c(0.1, 0.1),
                                    sigma_eps = 0.05, n_traj = 6)

etude <- function(H0, scenarios) {
  t(sapply(scenarios, function(gen)
    power_study(gen, H0, n_rep = h, B = L, n_pool = mp)$rate))
}

## H0 : HWP
etude("wiener", list(niveau = gen_hwp, HGP_bruit = gen_hgp_bruit, IGP_bruit = gen_igp_bruit,
                     NHWP = gen_nhwp, lineaire = gen_lineaire))

## H0 : HGP
etude("gamma", list(niveau = gen_hgp, HWP = gen_hwp, NHWP = gen_nhwp, IGP = gen_igp))

## H0 : IGP
etude("ig", list(niveau = gen_igp, HWP = gen_hwp, HGP = gen_hgp))

## H0 : NHWP
etude("wiener_power", list(niveau = gen_nhwp, HGP = gen_hgp, IGP = gen_igp,
                           lineaire = gen_lineaire, puissance = gen_puiss))

## H0 : NHGP
etude("gamma_power", list(niveau = gen_nhgp, HWP = gen_hwp, IGP = gen_igp, lineaire = gen_lineaire))
