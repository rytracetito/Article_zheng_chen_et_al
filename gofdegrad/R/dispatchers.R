#' Fermeture de simulation liee a un modele
#'
#' Construit une fonction de simulation \code{function(theta, T, m_traj)}
#' liee a un nom de modele, utilisable directement comme argument \code{sim}
#' de \code{\link{test_gof}}. Gere en interne la traduction des modeles
#' non-homogenes (\code{"Wiener_nonhomog"}, \code{"Gamma_nonhomog"}) vers
#' \code{\link{sim_proc}} avec la fonction de tendance \eqn{\Lambda(t) = t^r}
#' correspondante, \code{r} etant lu dans \code{theta$r}.
#'
#' @param modele Chaine de caracteres : \code{"Wiener"}, \code{"Gamma"},
#'   \code{"IG"}, \code{"Wiener_nonhomog"} ou \code{"Gamma_nonhomog"}.
#'
#' @return Une fonction \code{function(theta, T, m_traj) -> matrice}.
#'
#' @examples
#' sim_w <- sim_modele("Wiener")
#' T <- 1:50
#' X <- sim_w(list(mu_W = 0.2, sigma = 0.05), T, 6)
#' dim(X)
#'
#' @export
sim_modele <- function(modele) {
  switch(modele,
    Wiener_nonhomog = function(theta, T, m_traj) sim_proc(T, theta, "Wiener", m_traj, Lam = function(t) t^theta$r),
    Gamma_nonhomog  = function(theta, T, m_traj) sim_proc(T, theta, "Gamma",  m_traj, Lam = function(t) t^theta$r),
    function(theta, T, m_traj) sim_proc(T, theta, modele, m_traj)
  )
}

#' Fermeture d'estimation liee a un modele
#'
#' Construit une fonction d'estimation \code{function(X, T)} liee a un nom
#' de modele, utilisable directement comme argument \code{est} de
#' \code{\link{test_gof}}.
#'
#' @param modele Chaine de caracteres : \code{"Wiener"}, \code{"Gamma"},
#'   \code{"IG"}, \code{"Wiener_nonhomog"} ou \code{"Gamma_nonhomog"}.
#'
#' @return Une fonction \code{function(X, T) -> theta}.
#'
#' @examples
#' est_w <- est_modele("Wiener")
#' T <- 1:50
#' X <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 2000)
#' est_w(X, T)
#'
#' @export
est_modele <- function(modele) function(X, T) estimer_parametres(X, T, modele)
