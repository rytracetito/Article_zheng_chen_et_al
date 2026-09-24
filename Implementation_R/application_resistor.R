library(gofdegrad)
library(ggplot2)

L <- 1000
modeles <- c("Wiener", "Gamma", "IG", "Wiener_nonhomog", "Gamma_nonhomog")

analyser <- function(d, y, u, t) {
  d <- d[d[[t]] > 0, ]
  T <- sort(unique(d[[t]]))
  X <- t(sapply(split(d, d[[u]]), function(s) s[[y]][order(s[[t]])]))
  lignes <- lapply(modeles, function(mod) {
    set.seed(30)
    r <- tryCatch(test_gof(X, T, est_modele(mod), sim_modele(mod), mp = 100 * nrow(X), L = L),
                  error = function(e) NULL)
    if (is.null(r)) return(data.frame(modele = mod, KS_sup = NA, KS_sum = NA,
                                      CvM_sup = NA, CvM_sum = NA, rejets = "non testable"))
    data.frame(modele = mod, t(round(r$p_valeur, 3)), rejets = paste0(sum(r$rejet), "/4"))
  })
  do.call(rbind, lignes)
}

par_condition <- function(d, cond, ...) {
  do.call(rbind, lapply(sort(unique(d[[cond]])), function(g)
    cbind(condition = g, analyser(d[d[[cond]] == g, ], ...))))
}


# ---- derive relative (R - R0) / R0 ----

r1 <- SMRD.data::resistor

origine <- unique(r1[, c("resistor", "celsius")])
origine$percent <- 0
origine$hours <- 0
r1_graph <- rbind(origine[, names(r1)], r1)
figures <- lapply(sort(unique(r1$celsius)), function(g)
  ggplot(r1_graph[r1_graph$celsius == g, ], aes(hours, percent, group = resistor)) +
    geom_line(colour = "#2563EB", linewidth = 0.6) +
    geom_point(colour = "#2563EB", size = 1.6) +
    labs(title = paste0("Résistances à ", g, " °C : dérive relative"),
         x = "Temps (milliers d'heures)", y = "(R - R0) / R0") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
names(figures) <- paste0(sort(unique(r1$celsius)), "C")

print(figures$`83C`)
print(figures$`133C`)
print(figures$`173C`)

res_resistor <- par_condition(r1, "celsius", y = "percent", u = "resistor", t = "hours")
res_resistor
