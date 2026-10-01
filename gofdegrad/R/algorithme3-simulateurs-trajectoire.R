#' Gaussian generator of random effects
#'
#' Builds a sampling function for the random effects \eqn{\phi_k} under a
#' normality assumption, for use as the \code{gen_phi} argument of
#' \code{\link{sim_traj}}.
#'
#' @param mu Mean(s) of the \eqn{\phi} components: a single number for the
#'   \code{"lineaire"}/\code{"exponentiel"} models, a vector of length 2
#'   for \code{"loi_puissance"}.
#' @param sd Standard deviation(s), same length as \code{mu}.
#'
#' @return A function \code{function(mp)} that returns \code{mp} draws of
#'   \eqn{\phi}: a vector if \code{length(mu) == 1}, an \code{mp x 2}
#'   matrix otherwise.
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

                      #' Simulate trajectories under a general path model
                      #'
                      #' Generates \code{mp} simulated degradation trajectories under a general
                      #' path model (linear, power law or exponential), Algorithm 3 of
                      #' Zheng, Wang, Chen & Xie (2026).
                      #'
                      #' @param T Numeric vector of the \code{n} observation times.
                      #' @param th List \code{list(omega, s_eps)}: fixed effect and residual
                      #'   standard deviation. \code{omega} is only used for the
                      #'   \code{"exponentiel"} model.
                      #' @param mod Character string: \code{"lineaire"}, \code{"loi_puissance"}
                      #'   or \code{"exponentiel"}.
                      #' @param mp Integer number of trajectories \eqn{m'} to simulate.
                      #' @param gen_phi Sampling function for the random effects \eqn{\phi_k},
                      #'   typically produced by \code{\link{gen_phi_norm}}.
                      #'
                      #' @return Numeric \eqn{m' \times n} matrix of simulated trajectories.
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
                      
