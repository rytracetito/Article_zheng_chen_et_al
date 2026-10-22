# =====================================================================
#  Test complet du package gofdegrad sur ce PC
#  A executer dans RStudio : Ctrl + Shift + Entree (ou "Source")
#  Duree totale : environ 2 a 3 minutes
# =====================================================================

ok  <- function(msg) cat("  [OK]    ", msg, "\n")
ko  <- function(msg) cat("  [ECHEC] ", msg, "\n")
titre <- function(msg) cat("\n==", msg, "==\n")

# ---------------------------------------------------------------------
titre("1. Installation depuis GitHub")
# ---------------------------------------------------------------------
# On desinstalle d'abord la version locale pour tester une vraie
# installation, comme le ferait un autre utilisateur.
if ("gofdegrad" %in% rownames(installed.packages())) {
  if ("package:gofdegrad" %in% search()) detach("package:gofdegrad", unload = TRUE)
  remove.packages("gofdegrad")
}

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("rytracetito/Article_zheng_chen_et_al",
                        subdir = "gofdegrad", upgrade = "never", force = TRUE)

library(gofdegrad)
ok(paste("gofdegrad version", packageVersion("gofdegrad"), "installe et charge"))

# ---------------------------------------------------------------------
titre("2. Documentation")
# ---------------------------------------------------------------------
fonctions <- ls("package:gofdegrad")
cat("  Fonctions exportees :", paste(fonctions, collapse = ", "), "\n")
if (length(fonctions) == 10) ok("10 fonctions exportees") else ko("nombre de fonctions inattendu")
?gofdegrad   # doit ouvrir l'aide dans l'onglet Help de RStudio

# ---------------------------------------------------------------------
titre("3. Simulation des 7 modeles et des trajectoires generales")
# ---------------------------------------------------------------------
set.seed(123)
t <- 0:50
vrais <- list(
  wiener = c(mu = 0.3, sigma = 0.08),
  gamma = c(a = 1.5, b = 1.4),
  ig = c(mu = 0.8, eta = 1.2),
  vg = c(theta = 0.5, sigma = 1, nu = 0.5),
  meixner = c(a = 1, b = 0.3, d = 1.5, m = 0.2),
  wiener_power = c(mu = 0.15, sigma = 0.4, r = 1.7),
  gamma_power = c(a = 0.8, b = 2, r = 1.4)
)
for (m in names(vrais)) {
  x <- sim_model(t, m, vrais[[m]], n_traj = 6)
  if (identical(dim(x), c(51L, 7L)) && all(is.finite(x))) ok(paste(m, ": 6 trajectoires x 51 instants"))
  else ko(m)
}
xl <- sim_path(t, "linear", phi_mean = 0.2, phi_sd = 0.05, sigma_eps = 0.05, n_traj = 6)
xp <- sim_path(t, "power", phi_mean = c(0.2, 1.2), phi_sd = c(0.02, 0.05), sigma_eps = 0.05, n_traj = 6)
xe <- sim_path(t, "exponential", phi_mean = 0.5, phi_sd = 0.05, sigma_eps = 0.05, omega = 0.03, n_traj = 6)
ok("trajectoires generales lineaire, puissance et exponentielle simulees")
print(levysim::plot_traj(sim_model(t, "gamma", vrais$gamma, n_traj = 6), "Gamma"))
ok("graphique affiche dans l'onglet Plots")

# ---------------------------------------------------------------------
titre("4. Estimateurs MLE (grand echantillon : on doit retrouver les vrais parametres)")
# ---------------------------------------------------------------------
for (m in setdiff(names(vrais), "meixner")) {
  e <- fit_model(sim_model(t, m, vrais[[m]], n_traj = 2000), m)
  ecart <- max(abs(e - vrais[[m]]) / abs(vrais[[m]]))
  cat(sprintf("  %-13s estime = %-40s vrai = %s\n", m,
              paste(names(e), round(e, 3), sep = "=", collapse = " "),
              paste(names(vrais[[m]]), vrais[[m]], sep = "=", collapse = " ")))
  if (ecart < 0.10) ok(paste(m, ": ecart relatif < 10 %")) else ko(paste(m, ": ecart", round(ecart, 3)))
}

# ---------------------------------------------------------------------
titre("5. Test d'adequation (Algorithme 1)")
# ---------------------------------------------------------------------
# a) Donnees Wiener testees sous H0 = Wiener : on ne doit PAS rejeter
x <- sim_model(1:50, "wiener", c(mu = 0.2, sigma = 0.05), n_traj = 6)
r1 <- gof_test(x, "wiener", B = 200)
print(r1)
if (!any(r1$reject)) ok("H0 vraie (Wiener) : non rejetee") else ko("H0 vraie rejetee (possible par hasard, ~5 %)")

# b) Donnees Wiener non homogenes (r = 2) testees sous H0 = Wiener homogene :
#    on s'attend a un rejet
x2 <- sim_model(1:50, "wiener_power", c(mu = 0.01, sigma = 0.05, r = 2), n_traj = 6)
r2 <- gof_test(x2, "wiener", B = 200)
print(r2)
if (any(r2$reject)) ok("H0 fausse : rejetee par au moins une statistique") else ko("H0 fausse non rejetee")
print(plot(r2))

# ---------------------------------------------------------------------
titre("6. Gestion des erreurs")
# ---------------------------------------------------------------------
x_neg <- sim_model(1:50, "wiener", c(mu = 0.01, sigma = 1), n_traj = 6)
erreur <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))

if (grepl("negative increments", erreur(fit_model(x_neg, "gamma")))) ok("Gamma refuse des increments negatifs") else ko("pas d'erreur Gamma")
if (grepl("negative increments", erreur(fit_model(x_neg, "ig")))) ok("IG refuse des increments negatifs") else ko("pas d'erreur IG")
if (grepl("unknown model", erreur(fit_model(x, "truc")))) ok("modele inconnu refuse") else ko("modele inconnu accepte")
res <- compare_models(x_neg, c("wiener", "gamma"), B = 20, seed = 1)
if (grepl("not testable", res$note[2])) ok("compare_models signale le modele non testable") else ko("compare_models")

cat("\n== Fin du test : si toutes les lignes sont [OK], le package fonctionne. ==\n")
