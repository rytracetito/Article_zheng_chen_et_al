# =====================================================================
#  Test complet du package gofdegrad sur ce PC
#  A executer dans RStudio : Ctrl + Shift + Entree (ou "Source")
#  Duree totale : environ 1 a 2 minutes
# =====================================================================

ok  <- function(msg) cat("  [OK]    ", msg, "\n")
ko  <- function(msg) cat("  [ECHEC] ", msg, "\n")
titre <- function(msg) cat("\n==", msg, "==\n")

# ---------------------------------------------------------------------
titre("1. Installation depuis GitHub (depot prive)")
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
if (length(fonctions) == 14) ok("14 fonctions exportees") else ko("nombre de fonctions inattendu")
?test_gof   # doit ouvrir l'aide dans l'onglet Help de RStudio

# ---------------------------------------------------------------------
titre("3. Simulateurs (Algorithmes 2 et 3)")
# ---------------------------------------------------------------------
set.seed(123)
T <- 1:50

Xw <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
Xg <- sim_proc(T, list(alpha = 1.5, beta = 0.7), "Gamma", mp = 6)
Xi <- sim_proc(T, list(mu_IG = 0.8, eta = 1.2), "IG", mp = 6)
Xl <- sim_traj(T, list(s_eps = 0.05), "lineaire", mp = 6,
               gen_phi_norm(mu = 0.2, sd = 0.05))
Xp <- sim_traj(T, list(s_eps = 0.05), "loi_puissance", mp = 6,
               gen_phi_norm(mu = c(0.2, 1.2), sd = c(0.02, 0.05)))
Xe <- sim_traj(T, list(s_eps = 0.05, omega = 0.03), "exponentiel", mp = 6,
               gen_phi_norm(mu = 0.5, sd = 0.05))

for (nom in c("Xw", "Xg", "Xi", "Xl", "Xp", "Xe")) {
  X <- get(nom)
  if (identical(dim(X), c(6L, 50L)) && all(is.finite(X))) ok(paste(nom, ": matrice 6 x 50"))
  else ko(nom)
}
if (all(apply(Xg, 1, diff) > 0)) ok("Gamma : trajectoires strictement croissantes")

# Graphique : une trajectoire de chaque modele
op <- par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
matplot(T, t(Xw), type = "l", lty = 1, main = "Wiener",       ylab = "X(t)")
matplot(T, t(Xg), type = "l", lty = 1, main = "Gamma",        ylab = "X(t)")
matplot(T, t(Xi), type = "l", lty = 1, main = "IG",           ylab = "X(t)")
matplot(T, t(Xl), type = "l", lty = 1, main = "Lineaire",     ylab = "X(t)")
matplot(T, t(Xp), type = "l", lty = 1, main = "Loi puissance", ylab = "X(t)")
matplot(T, t(Xe), type = "l", lty = 1, main = "Exponentiel",  ylab = "X(t)")
par(op)
ok("graphiques affiches dans l'onglet Plots")

# ---------------------------------------------------------------------
titre("4. Estimateurs MLE (grand echantillon : on doit retrouver les vrais parametres)")
# ---------------------------------------------------------------------
compare <- function(nom, estime, vrai) {
  ecart <- max(abs(unlist(estime[names(vrai)]) - unlist(vrai)) / unlist(vrai))
  cat(sprintf("  %-16s estime = %-35s vrai = %s\n", nom,
              paste(names(estime), round(unlist(estime), 3), sep = "=", collapse = " "),
              paste(names(vrai), unlist(vrai), sep = "=", collapse = " ")))
  if (ecart < 0.10) ok(paste(nom, ": ecart relatif < 10 %")) else ko(paste(nom, ": ecart", round(ecart, 3)))
}

vrai <- list(mu_W = 0.3, sigma = 0.08)
compare("Wiener", estimer_parametres(sim_proc(T, vrai, "Wiener", 2000), T, "Wiener"), vrai)

vrai <- list(alpha = 1.5, beta = 0.7)
compare("Gamma", estimer_parametres(sim_proc(T, vrai, "Gamma", 2000), T, "Gamma"), vrai)

vrai <- list(mu_IG = 0.8, eta = 1.2)
compare("IG", estimer_parametres(sim_proc(T, vrai, "IG", 2000), T, "IG"), vrai)

vrai <- list(mu_W = 0.15, sigma = 0.4, r = 1.7)
compare("Wiener_nonhomog",
        estimer_parametres(sim_proc(T, vrai, "Wiener", 2000, Lam = function(t) t^1.7), T, "Wiener_nonhomog"),
        vrai)

vrai <- list(alpha = 0.8, beta = 0.5, r = 1.4)
compare("Gamma_nonhomog",
        estimer_parametres(sim_proc(T, vrai, "Gamma", 2000, Lam = function(t) t^1.4), T, "Gamma_nonhomog"),
        vrai)

# ---------------------------------------------------------------------
titre("5. Test d'adequation (Algorithme 1)")
# ---------------------------------------------------------------------
# a) Donnees Wiener testees sous H0 = Wiener : on ne doit PAS rejeter
X  <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
r1 <- test_gof(X, T, est_modele("Wiener"), sim_modele("Wiener"), mp = 600, L = 200)
print(r1[c("rejet", "p_valeur")])
if (!any(r1$rejet)) ok("H0 vraie (Wiener) : non rejetee") else ko("H0 vraie rejetee (possible par hasard, ~5 %)")

# b) Donnees Wiener NON-homogenes (r = 2) testees sous H0 = Wiener homogene :
#    on s'attend a un rejet
X2 <- sim_proc(T, list(mu_W = 0.01, sigma = 0.05), "Wiener", mp = 6, Lam = function(t) t^2)
r2 <- test_gof(X2, T, est_modele("Wiener"), sim_modele("Wiener"), mp = 600, L = 200)
print(r2[c("rejet", "p_valeur")])
if (any(r2$rejet)) ok("H0 fausse : rejetee par au moins une statistique") else ko("H0 fausse non rejetee")

# ---------------------------------------------------------------------
titre("6. Gestion des erreurs")
# ---------------------------------------------------------------------
# Wiener tres bruite : contient forcement des increments negatifs
X_neg <- sim_proc(T, list(mu_W = 0.01, sigma = 1), "Wiener", mp = 6)
erreur <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))

if (grepl("increments", erreur(estimer_parametres(X_neg, T, "Gamma")))) ok("Gamma refuse des increments negatifs") else ko("pas d'erreur Gamma")
if (grepl("increments", erreur(estimer_parametres(X_neg, T, "IG"))))    ok("IG refuse des increments negatifs")    else ko("pas d'erreur IG")
if (grepl("modele inconnu", erreur(estimer_parametres(Xw, T, "Truc")))) ok("modele inconnu refuse") else ko("modele inconnu accepte")

cat("\n== Fin du test : si toutes les lignes sont [OK], le package fonctionne. ==\n")
