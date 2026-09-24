# Algorithme 2 - simulation de trajectoires (Wiener / Gamma / IG)

library(statmod)

# T   : instants d'observation
# th  : parametres du modele
#       (mu_W, sigma) Wiener / (alpha, beta) Gamma / (mu_IG, eta) IG
# mod : "Wiener", "Gamma" ou "IG"
# mp  : nombre m' de trajectoires a simuler
# Lam : fonction de tendance Lambda(t), identity = cas homogene

sim_proc <- function(T, th, mod, mp, Lam = identity) {
  n <- length(T)
  t0T <- c(0, T)
  dt <- diff(t0T)
  dLam <- diff(Lam(t0T))
  
  dX <- matrix(0, mp, n)
  for (j in seq_len(n)) {
    dX[, j] <- switch(mod,
                      Wiener = rnorm(mp, th$mu_W * dLam[j], th$sigma * sqrt(dt[j])),
                      Gamma  = rgamma(mp, shape = th$alpha * dLam[j], scale = th$beta),
                      IG     = rinvgauss(mp, mean = th$mu_IG * dLam[j], shape = th$eta * dLam[j]^2)
    )
  }
  t(apply(dX, 1, cumsum))
}
