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

mw <- SMRD.data::metalwear
mw$charge <- factor(paste(mw$grams, "g"), levels = paste(sort(unique(mw$grams)), "g"))

ggplot(mw, aes(cycles, microns, group = unit)) +
  geom_line(colour = "#2563EB", linewidth = 0.5) +
  geom_point(colour = "#2563EB", size = 1.2) +
  scale_x_log10(breaks = sort(unique(mw$cycles))) +
  facet_wrap(~ charge, scales = "free_y") +
  labs(title = "Usure d'un alliage métallique", x = "Cycles (échelle log)", y = "Usure (microns)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

res_metalwear <- par_condition(mw, "grams", y = "microns", u = "unit", t = "cycles")
res_metalwear
