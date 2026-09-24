# Algorithme 3 - simulation de trajectoires generales (lineaire / puissance / exponentiel)

# mu : moyenne(s) des effets aleatoires phi
# sd : ecart(s)-type des effets aleatoires phi
#      (longueur 1 pour lineaire/exponentiel, longueur 2 pour loi de puissance)

gen_phi_norm <- function(mu, sd) {
  function(mp) {
    if (length(mu) == 1L) rnorm(mp, mu, sd)
    else sapply(seq_along(mu), function(k) rnorm(mp, mu[k], sd[k]))
  }
}

# T       : instants d'observation
# th      : list(omega, s_eps) - effets fixes et ecart-type de l'erreur de mesure
# mod     : "lineaire", "loi_puissance" ou "exponentiel"
# mp      : nombre m' de trajectoires a simuler
# gen_phi : fonction de tirage des effets aleatoires phi_k

sim_traj <- function(T, th, mod, mp, gen_phi) {
  n <- length(T)
  Tm <- matrix(T, nrow = mp, ncol = n, byrow = TRUE)
  phi <- gen_phi(mp)
  
  g <- switch(mod,
              lineaire      = phi * Tm,
              loi_puissance = phi[, 1] * Tm^phi[, 2],
              exponentiel   = phi * (exp(th$omega * Tm) - 1)
  )
  g + matrix(rnorm(mp * n, 0, th$s_eps), mp, n)
}
