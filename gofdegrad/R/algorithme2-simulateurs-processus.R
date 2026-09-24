#' Simulation de trajectoires sous un modele a processus stochastique
#'
#' Genere \code{mp} trajectoires de degradation simulees sous un modele
#' Wiener, Gamma ou Inverse Gaussien, homogene ou non-homogene, aux memes
#' instants d'observation que des donnees reelles (Algorithme 2 de
#' Zheng, Wang, Chen & Xie, 2026).
#'
#' @param T Vecteur numerique des \code{n} instants d'observation
#'   \eqn{(t_1, \dots, t_n)}. \eqn{t_0 = 0} est ajoute automatiquement.
#' @param th Liste des parametres du modele : \code{list(mu_W, sigma)} pour
#'   \code{"Wiener"}, \code{list(alpha, beta)} pour \code{"Gamma"},
#'   \code{list(mu_IG, eta)} pour \code{"IG"}.
#' @param mod Chaine de caracteres : \code{"Wiener"}, \code{"Gamma"} ou
#'   \code{"IG"}.
#' @param mp Nombre entier de trajectoires \eqn{m'} a simuler.
#' @param Lam Fonction de tendance \eqn{\Lambda(t)}. Par defaut
#'   \code{identity} (cas homogene, \eqn{\Lambda(t) = t}).
#'
#' @return Matrice numerique \eqn{m' \times n} des trajectoires simulees ;
#'   la ligne \code{k} est la trajectoire \eqn{X'_k = [x'_{k,1}, \dots, x'_{k,n}]}.
#'
#' @examples
#' T <- 1:50
#' X <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 5)
#' dim(X)
#'
#' @export
sim_proc <- function(T, th, mod, mp, Lam = identity) {
  n <- length(T)
  t0T <- c(0, T)
  dt <- diff(t0T)
  dLam <- diff(Lam(t0T))

  dX <- matrix(0, mp, n)
  for (j in seq_len(n)) {
    dX[, j] <- switch(mod,
      Wiener = stats::rnorm(mp, th$mu_W * dLam[j], th$sigma * sqrt(dt[j])),
      Gamma  = stats::rgamma(mp, shape = th$alpha * dLam[j], scale = th$beta),
      IG     = statmod::rinvgauss(mp, mean = th$mu_IG * dLam[j], shape = th$eta * dLam[j]^2)
    )
  }
  t(apply(dX, 1, cumsum))
}
