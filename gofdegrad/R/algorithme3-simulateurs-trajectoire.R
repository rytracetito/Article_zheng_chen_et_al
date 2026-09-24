#' Generateur gaussien d'effets aleatoires
#'
#' Construit une fonction de tirage des effets aleatoires \eqn{\phi_k} sous
#' hypothese de normalite, a utiliser comme argument \code{gen_phi} de
#' \code{\link{sim_traj}}.
#'
#' @param mu Moyenne(s) des composantes de \eqn{\phi} : un seul nombre pour
#'   les modeles \code{"lineaire"}/\code{"exponentiel"}, un vecteur de
#'   longueur 2 pour \code{"loi_puissance"}.
#' @param sd Ecart(s)-type associe(s), meme longueur que \code{mu}.
#'
#' @return Une fonction \code{function(mp)} qui retourne \code{mp} tirages
#'   de \eqn{\phi} : un vecteur si \code{length(mu) == 1}, une matrice
#'   \code{mp x 2} sinon.
#'
#' @examples
#' gen_phi <- gen_phi_norm(mu = 0.2, sd = 0.05)
#' gen_phi(5)
#'
#' @export
gen_phi_norm <- function(mu, sd) {
  function(mp) {
    if (length(mu) == 1L) stats::rnorm(mp, mu, sd)
    else sapply(seq_along(mu), function(k) stats::rnorm(mp, mu[k], sd[k]))
  }
}

#' Simulation de trajectoires sous un modele a trajectoire generale
#'
#' Genere \code{mp} trajectoires de degradation simulees sous un modele a
#' trajectoire generale (lineaire, loi de puissance ou exponentiel),
#' Algorithme 3 de Zheng, Wang, Chen & Xie (2026).
#'
#' @param T Vecteur numerique des \code{n} instants d'observation.
#' @param th Liste \code{list(omega, s_eps)} : effets fixes et ecart-type de
#'   l'erreur de mesure. \code{omega} n'est utilise que pour le modele
#'   \code{"exponentiel"}.
#' @param mod Chaine de caracteres : \code{"lineaire"}, \code{"loi_puissance"}
#'   ou \code{"exponentiel"}.
#' @param mp Nombre entier de trajectoires \eqn{m'} a simuler.
#' @param gen_phi Fonction de tirage des effets aleatoires \eqn{\phi_k},
#'   typiquement produite par \code{\link{gen_phi_norm}}.
#'
#' @return Matrice numerique \eqn{m' \times n} des trajectoires simulees.
#'
#' @examples
#' T <- 1:50
#' phi <- gen_phi_norm(mu = 0.2, sd = 0.05)
#' X <- sim_traj(T, list(s_eps = 0.05), "lineaire", mp = 5, phi)
#' dim(X)
#'
#' @export
sim_traj <- function(T, th, mod, mp, gen_phi) {
  n <- length(T)
  Tm <- matrix(T, nrow = mp, ncol = n, byrow = TRUE)
  phi <- gen_phi(mp)

  g <- switch(mod,
    lineaire      = phi * Tm,
    loi_puissance = phi[, 1] * Tm^phi[, 2],
    exponentiel   = phi * (exp(th$omega * Tm) - 1)
  )
  g + matrix(stats::rnorm(mp * n, 0, th$s_eps), mp, n)
}
