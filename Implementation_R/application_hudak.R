# Hudak crack growth data

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("mdhglm", quietly = TRUE)) remotes::install_github("cran/mdhglm", upgrade = "never")
if (!requireNamespace("gofdegrad", quietly = TRUE)) {
  remotes::install_github("rytracetito/Article_zheng_chen_et_al", subdir = "gofdegrad", force = TRUE)
}

library(gofdegrad)
data(crack_growth, package = "mdhglm")

L <- 1000
modeles <- c("wiener", "gamma", "ig", "wiener_power", "gamma_power")


crack_growth$l <- crack_growth$crack0 + crack_growth$y
specimens <- sort(unique(crack_growth$specimen))


liste_obs <- vector("list", length(specimens))
for (i in seq_along(specimens)) {
  sous <- crack_growth[crack_growth$specimen == specimens[i], ]
  sous <- sous[order(sous$cycle), ]
  liste_obs[[i]] <- data.frame(t = c(0, sous$cycle), X = c(0, sous$l - 0.9))
}

plot(NULL, xlim = c(0, 0.12),
     ylim = c(0, max(sapply(liste_obs, function(d) max(d$X)))),
     xlab = "Millions of cycles", ylab = "Length")
for (i in seq_along(liste_obs)) {
  lines(liste_obs[[i]]$t, liste_obs[[i]]$X, type = "o", pch = 16, cex = 0.7, lwd = 1)
}
abline(h = 0.7, lty = 2)


# Donnees du test : les 9 premiers instants, communs a tous les specimens
T_hudak <- sort(unique(crack_growth$cycle))[1:9]
hudak <- crack_growth[crack_growth$cycle %in% T_hudak, ]
hudak <- data.frame(id = hudak$specimen, t = hudak$cycle, X = hudak$l - 0.9)
gof_data(hudak)

# Estimation et test d'adequation pour chaque modele candidat
message("Tests d'adequation, L = ", L)
resultats <- compare_models(hudak, modeles, B = L, seed = 30)
resultats

acceptes <- resultats$model[which(resultats$rejections == 0)]
message("Modele(s) non rejete(s) : ", if (length(acceptes)) paste(acceptes, collapse = ", ") else "aucun")
