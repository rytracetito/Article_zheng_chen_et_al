if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("SMRD.data", quietly = TRUE)) remotes::install_github("Auburngrads/SMRD.data")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("gofdegrad", quietly = TRUE)) {
  remotes::install_github("rytracetito/Article_zheng_chen_et_al", subdir = "gofdegrad", force = TRUE)
}

library(gofdegrad)
library(ggplot2)
L <- 1000
modeles <- c("wiener", "gamma", "ig", "wiener_power", "gamma_power")


# Application 1 : Données complète

d <- SMRD.data::gaaslaser
laser <- data.frame(id = d$unit, t = d$hours, X = d$increase)
donnees <- gof_data(laser)
donnees

# Visualisation des trajectoires
seuil <- 10
d$laser <- factor(d$unit)
d$etat <- ifelse(ave(d$increase, d$unit, FUN = max) >= seuil,
                 "Depasse le seuil", "Reste sous le seuil")

fig_traj <- ggplot(d, aes(hours, increase, group = laser, colour = etat)) +
  geom_hline(yintercept = seuil, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
  annotate("text", x = 0, y = seuil, label = "seuil de defaillance (10 %)",
           hjust = 0, vjust = -0.6, size = 3.3, colour = "grey30") +
  geom_line(linewidth = 0.6, alpha = 0.9) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = c("Depasse le seuil" = "#C2410C",
                                 "Reste sous le seuil" = "#2563EB")) +
  scale_x_continuous(breaks = seq(0, 4000, 500), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(breaks = seq(0, 14, 2)) +
  labs(title = "GaAs laser : trajectoires de degradation",
       x = "Temps (heures)", y = "Degradation", colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        legend.justification = "left",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
print(fig_traj)

# Estimation et test d'adequation pour chaque modele candidat
message("Tests d'adequation, L = ", L)
resultats <- compare_models(laser, modeles, B = L, seed = 30)
resultats

# Comparaison avec les valeurs publiees dans l'article pour le cas Wiener homogene
wiener <- attr(resultats, "tests")$wiener
wiener
plot(wiener)
message("Wiener homogene : comparaison avec l'article")
print(rbind(
  S_obs_package = round(wiener$statistic, 4),
  S_obs_article = c(1.8752, 6.4109, 12.7108, 19.9976),
  s_crit_package = round(wiener$critical, 4),
  s_crit_article = c(0.5661, 5.0354, 1.4075, 7.7271),
  p_val_package = round(wiener$p.value, 4),
  p_val_article = c(0.001, 0.022, 0.001, 0.006)
))

acceptes <- resultats$model[which(resultats$rejections == 0)]
message("Modele(s) non rejete(s) : ", if (length(acceptes)) paste(acceptes, collapse = ", ") else "aucun")


# Application 2 : inspections irregulieres, (frequence decroissante)

temps_conserves <- c(0, 250, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000)
laser_cas1 <- laser[laser$t %in% temps_conserves, ]

# Visualisation des trajectoires
liste_obs <- split(laser_cas1[c("t", "X")], laser_cas1$id)
plot(NULL, xlim = range(temps_conserves),
     ylim = c(0, max(laser_cas1$X)),
     xlab = "Time", ylab = "Performance characteristic")
for (i in seq_along(liste_obs)) {
  lines(liste_obs[[i]]$t, liste_obs[[i]]$X, type = "o", pch = 16, cex = 0.7, lwd = 1)
}
abline(h = 10, lty = 2)

gof_data(laser_cas1)

# Estimation et test d'adequation pour chaque modele candidat
message("Tests d'adequation, L = ", L)
resultats_cas1 <- compare_models(laser_cas1, modeles, B = L, seed = 30)
resultats_cas1

acceptes <- resultats_cas1$model[which(resultats_cas1$rejections == 0)]
message("Modele(s) non rejete(s) : ", if (length(acceptes)) paste(acceptes, collapse = ", ") else "aucun")


# Application 3 : inspections irregulieres, (frequence croissante)

temps_conserves <- c(0, 1000, 1500, 2000, 2500, 3000, 3250, 3500, 3750, 4000)
laser_cas2 <- laser[laser$t %in% temps_conserves, ]

# Visualisation des trajectoires
liste_obs <- split(laser_cas2[c("t", "X")], laser_cas2$id)
plot(NULL, xlim = range(temps_conserves),
     ylim = c(0, max(laser_cas2$X)),
     xlab = "Time", ylab = "Performance characteristic")
for (i in seq_along(liste_obs)) {
  lines(liste_obs[[i]]$t, liste_obs[[i]]$X, type = "o", pch = 16, cex = 0.7, lwd = 1)
}
abline(h = 10, lty = 2)

gof_data(laser_cas2)

# Estimation et test d'adequation pour chaque modele candidat
message("Tests d'adequation, L = ", L)
resultats_cas2 <- compare_models(laser_cas2, modeles, B = L, seed = 30)
resultats_cas2

acceptes <- resultats_cas2$model[which(resultats_cas2$rejections == 0)]
message("Modele(s) non rejete(s) : ", if (length(acceptes)) paste(acceptes, collapse = ", ") else "aucun")
