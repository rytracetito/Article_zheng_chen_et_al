library(gofdegrad)
library(ggplot2)

L <- 1000
modeles <- c("wiener", "gamma", "ig", "wiener_power", "gamma_power")

par_condition <- function(d, cond, y, u, t) {
  do.call(rbind, lapply(sort(unique(d[[cond]])), function(g) {
    s <- d[d[[cond]] == g, ]
    cbind(condition = g,
          compare_models(data.frame(id = s[[u]], t = s[[t]], X = s[[y]]), modeles, B = L, seed = 30))
  }))
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
