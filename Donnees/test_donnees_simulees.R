library(gofdegrad)
library(ggplot2)

dossier <- "C:/Users/ADMIN/Desktop/Ma Thèse/Thèse_1/Recherche_1/article_zheng_chen/Donnees"
d <- read.csv(file.path(dossier, "degradation_simulee_20traj.csv"))

donnees <- data.frame(id = d$trajectoire, t = d$temps, X = d$degradation)
gof_data(donnees)

L <- 1000

fig <- ggplot(d, aes(temps, degradation, group = trajectoire)) +
  geom_line(colour = "#2563EB", linewidth = 0.5, alpha = 0.8) +
  geom_point(colour = "#2563EB", size = 0.9) +
  labs(title = "Données simulées : 20 trajectoires", x = "Temps", y = "Dégradation") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
print(fig)

modeles <- c("wiener", "gamma", "ig", "wiener_power", "gamma_power")

resultats <- compare_models(donnees, modeles, B = L, seed = 30)
resultats
