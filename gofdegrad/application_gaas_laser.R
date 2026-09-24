if (!requireNamespace("remotes", quietly = TRUE))   install.packages("remotes")
if (!requireNamespace("SMRD.data", quietly = TRUE)) install.packages("SMRD.data")
if (!requireNamespace("ggplot2", quietly = TRUE))   install.packages("ggplot2")
if (!requireNamespace("gofdegrad", quietly = TRUE)) {
  remotes::install_github("rytracetito/Article_zheng_chen_et_al", subdir = "gofdegrad")
}

library(gofdegrad)
library(ggplot2)

L <- 1000

d <- SMRD.data::gaaslaser
unites <- unique(d$unit)
heures <- sort(unique(d$hours))

X_full <- t(sapply(unites, function(k) {
  sous <- d[d$unit == k, ]
  sous$increase[order(sous$hours)]
}))

T_gaas <- heures[-1]
X_gaas <- X_full[, -1]
m <- nrow(X_gaas)
cat("Donnees :", m, "trajectoires x", ncol(X_gaas), "instants\n")

seuil <- 10
d$laser <- factor(d$unit)
d$etat  <- ifelse(ave(d$increase, d$unit, FUN = max) >= seuil,
                  "Depasse le seuil", "Reste sous le seuil")

fig_traj <- ggplot(d, aes(hours, increase, group = laser, colour = etat)) +
  geom_hline(yintercept = seuil, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
  annotate("text", x = 0, y = seuil, label = "seuil de defaillance (10 %)",
           hjust = 0, vjust = -0.6, size = 3.3, colour = "grey30") +
  geom_line(linewidth = 0.6, alpha = 0.9) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = c("Depasse le seuil"    = "#C2410C",
                                 "Reste sous le seuil" = "#2563EB")) +
  scale_x_continuous(breaks = seq(0, 4000, 500), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(breaks = seq(0, 14, 2)) +
  labs(title = "GaAs laser : trajectoires de degradation",
       x = "Temps (heures)", y = "Dégradation", colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position      = "top",
        legend.justification = "left",
        panel.grid.minor     = element_blank(),
        plot.title           = element_text(face = "bold"))
print(fig_traj)

modeles <- c("Wiener", "Gamma", "IG", "Wiener_nonhomog", "Gamma_nonhomog")

cat("\n== Parametres estimes (MLE) ==\n")
for (mod in modeles) {
  th <- estimer_parametres(X_gaas, T_gaas, mod)
  cat(sprintf("  %-16s %s\n", mod,
              paste(names(th), signif(unlist(th), 4), sep = " = ", collapse = ", ")))
}

cat("\n== Tests d'adequation (L =", L, ") ==\n")
resultats <- list()
for (mod in modeles) {
  cat("  test du modele", mod, "...\n")
  set.seed(30)
  resultats[[mod]] <- test_gof(X_gaas, T_gaas,
                               est = est_modele(mod),
                               sim = sim_modele(mod),
                               mp  = 100 * m,
                               L   = L)
}

p_valeurs <- t(sapply(resultats, `[[`, "p_valeur"))
rejets    <- t(sapply(resultats, `[[`, "rejet"))

cat("\n== p-valeurs ==\n")
print(round(p_valeurs, 3))

cat("\n== Rejet de H0 ==\n")
print(rejets)

cat("\n== Wiener homogene : comparaison avec l'article ==\n")
print(rbind(
  S_obs_package  = round(resultats$Wiener$S_obs, 4),
  S_obs_article  = c(1.8752, 6.4109, 12.7108, 19.9976),
  s_crit_package = round(resultats$Wiener$s_crit, 4),
  s_crit_article = c(0.5661, 5.0354, 1.4075, 7.7271),
  p_val_package  = round(resultats$Wiener$p_valeur, 4),
  p_val_article  = c(0.001, 0.022, 0.001, 0.006)
))

acceptes <- rownames(rejets)[rowSums(rejets) == 0]
cat("\nModele(s) non rejete(s) :",
    if (length(acceptes)) paste(acceptes, collapse = ", ") else "aucun", "\n")
