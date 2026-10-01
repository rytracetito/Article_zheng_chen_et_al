#' Simulate trajectories under a stochastic process model
#'
#' Generates \code{mp} simulated degradation trajectories under a Wiener,
#' Gamma or Inverse Gaussian model, homogeneous or non-homogeneous, at the
#' same observation times as real data (Algorithm 2 of Zheng, Wang, Chen &
#' Xie, 2026).
#'
#' @param T Numeric vector of the \code{n} observation times
#'   \eqn{(t_1, \dots, t_n)}. \eqn{t_0 = 0} is added automatically.
#' @param th List of model parameters: \code{list(mu_W, sigma)} for
#'   \code{"Wiener"}, \code{list(alpha, beta)} for \code{"Gamma"},
#'   \code{list(mu_IG, eta)} for \code{"IG"}.
#' @param mod Character string: \code{"Wiener"}, \code{"Gamma"} or
#'   \code{"IG"}.
#' @param mp Integer number of trajectories \eqn{m'} to simulate.
#' @param Lam Trend function \eqn{\Lambda(t)}. Defaults to \code{identity}
#'   (homogeneous case, \eqn{\Lambda(t) = t}).
#'
#' @return Numeric \eqn{m' \times n} matrix of simulated trajectories; row
#'   \code{k} is the trajectory \eqn{X'_k = [x'_{k,1}, \dots, x'_{k,n}]}.
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
