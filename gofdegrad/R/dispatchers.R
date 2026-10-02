#' Simulation closure bound to a model
#'
#' Builds a simulation function \code{function(theta, T, m_traj)} bound to
#' a model name, usable directly as the \code{sim} argument of
#' \code{\link{test_gof}}. Internally handles the translation of the
#' non-homogeneous models (\code{"Wiener_nonhomog"}, \code{"Gamma_nonhomog"})
#' into \code{\link{sim_proc}} with the corresponding trend function
#' \eqn{\Lambda(t) = t^r}, \code{r} being read from \code{theta$r}.
#'
#' @param modele Character string: \code{"Wiener"}, \code{"Gamma"},
#'   \code{"IG"}, \code{"Wiener_nonhomog"} or \code{"Gamma_nonhomog"}.
#'
#' @return A function \code{function(theta, T, m_traj) -> matrix}.
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

                                                                     #' Estimation closure bound to a model
                                                                     #'
                                                                     #' Builds an estimation function \code{function(X, T)} bound to a model
                                                                     #' name, usable directly as the \code{est} argument of
                                                                     #' \code{\link{test_gof}}.
                                                                     #'
                                                                     #' @param modele Character string: \code{"Wiener"}, \code{"Gamma"},
                                                                     #'   \code{"IG"}, \code{"Wiener_nonhomog"} or \code{"Gamma_nonhomog"}.
                                                                     #'
                                                                     #' @return A function \code{function(X, T) -> theta}.
                                                                     #'
                                                                     #' @examples
                                                                     #' est_w <- est_modele("Wiener")
                                                                     #' T <- 1:50
                                                                     #' X <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 2000)
                                                                     #' est_w(X, T)
                                                                     #'
                                                                     #' @export
                                                                     est_modele <- function(modele) function(X, T) estimer_parametres(X, T, modele)
                                                                       
