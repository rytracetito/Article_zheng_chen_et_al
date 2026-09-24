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
