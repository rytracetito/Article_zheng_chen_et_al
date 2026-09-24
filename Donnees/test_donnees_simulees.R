library(gofdegrad)
library(ggplot2)

dossier <- "C:/Users/ADMIN/Desktop/Ma Thèse/Thèse_1/Code_R/article_zheng_chen/Donnees"
d <- read.csv(file.path(dossier, "degradation_simulee_20traj.csv"))

T <- sort(unique(d$temps))
X <- t(sapply(split(d, d$trajectoire), function(s) s$degradation[order(s$temps)]))
dim(X)

L <- 1000

fig <- ggplot(d, aes(temps, degradation, group = trajectoire)) +
  geom_line(colour = "#2563EB", linewidth = 0.5, alpha = 0.8) +
  geom_point(colour = "#2563EB", size = 0.9) +
  labs(title = "Données simulées : 20 trajectoires", x = "Temps", y = "Dégradation") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
print(fig)

modeles <- c("Wiener", "Gamma", "IG", "Wiener_nonhomog", "Gamma_nonhomog")

estimations <- lapply(modeles, function(mod)
  tryCatch(estimer_parametres(X, T, mod), error = function(e) conditionMessage(e)))
names(estimations) <- modeles
estimations

resultats <- list()
for (mod in modeles) {
  set.seed(30)
  resultats[[mod]] <- tryCatch(
    test_gof(X, T, est_modele(mod), sim_modele(mod), mp = 100 * nrow(X), L = L),
    error = function(e) NULL)
}

testes <- Filter(Negate(is.null), resultats)
setdiff(modeles, names(testes))

round(t(sapply(testes, `[[`, "p_valeur")), 3)
t(sapply(testes, `[[`, "rejet"))
