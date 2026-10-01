#' Distance matrix between trajectories
#'
#' Computes the \eqn{m \times m'} matrix of similarity distances between
#' each real trajectory and each simulated trajectory (eq. 3-4 of
#' Zheng, Wang, Chen & Xie, 2026).
#'
#' @param X Matrix of \eqn{m \times n} real (or bootstrap) trajectories.
#' @param Xs Matrix of \eqn{m' \times n} simulated trajectories.
#' @param D Similarity measure: \code{"KS"} (Chebyshev distance, by
#'   analogy with the Kolmogorov-Smirnov statistic) or \code{"CvM"}
#'   (squared Euclidean distance).
#'
#' @return Numeric \eqn{m \times m'} matrix, term \eqn{(i,k)} equal to
#'   \eqn{D(X_i, X'_k)}.
#'
#' @keywords internal
dist_mat <- function(X, Xs, D) {
    if (D == "CvM") {
          outer(rowSums(X^2), rowSums(Xs^2), "+") - 2 * X %*% t(Xs)
        } else {
          n <- ncol(X)
          Dm <- matrix(0, nrow(X), nrow(Xs))
          for (j in seq_len(n)) Dm <- pmax(Dm, abs(outer(X[, j], Xs[, j], "-")))
          Dm
        }
  }

#' Similarity statistics
#'
#' Computes the 4 test statistics \eqn{S_{KS-sup}}, \eqn{S_{KS-sum}},
#' \eqn{S_{CvM-sup}}, \eqn{S_{CvM-sum}} between a set of trajectories and a
#' simulated pool (eq. 6-10 of the article).
#'
#' @param X Matrix of real or bootstrap trajectories.
#' @param Xs Matrix of the pool of trajectories simulated under the
#'   candidate model.
#'
#' @return Named numeric vector of length 4:
#'   \code{c(KS_sup, KS_sum, CvM_sup, CvM_sum)}.
#'
#' @export
stats_sim <- function(X, Xs) {
    d_KS  <- apply(dist_mat(X, Xs, "KS"),  1, min)
    d_CvM <- apply(dist_mat(X, Xs, "CvM"), 1, min)
    c(KS_sup  = max(d_KS),  KS_sum  = sum(d_KS),
          CvM_sup = max(d_CvM), CvM_sum = sum(d_CvM))
  }

#' Goodness-of-fit test by similarity and parametric bootstrap
#'
#' Implements the full Algorithm 1 of Zheng, Wang, Chen & Xie (2026):
#' computes the observed statistic, the critical threshold by parametric
#' bootstrap, and returns the test decision for all 4 statistics at once.
#'
#' @param X Matrix of \eqn{m \times n} real data.
#' @param T Vector of the \code{n} observation times.
#' @param est Estimation function \code{function(X, T) -> theta}, specific
#'   to the candidate model under \eqn{H_0} (see
#'   \code{\link{estimer_parametres}}).
#' @param sim Simulation function \code{function(theta, T, m_traj) ->
#'   matrix}, specific to the candidate model (see \code{\link{sim_proc}},
#'   \code{\link{sim_traj}}).
#' @param mp Size of the simulated pool (\eqn{m'}). Defaults to \eqn{100m},
#'   the recommendation given in the article.
#' @param L Number of bootstrap replications. Defaults to 1000.
#' @param alpha Significance level of the test. Defaults to 0.05.
#'
#' @return A list with elements:
#'   \item{rejet}{logical vector of length 4, \code{TRUE} = \eqn{H_0} rejected}
#'   \item{S_obs}{observed statistics}
#'   \item{s_crit}{critical values (empirical quantile of order \eqn{1-\alpha})}
#'   \item{p_valeur}{empirical p-values}
#'
#' @details
#' Inside the bootstrap loop, warnings raised by \code{est} (for example the
#' zero-increment tolerance warning issued by
#' \code{\link{estimer_gamma_homogene}}, \code{\link{estimer_ig_homogene}}
#' and \code{\link{estimer_gamma_nonhomog}} via their shared helper) are
#' suppressed, since they would otherwise be repeated up to \code{L} times
#' for a single call to \code{test_gof}. A warning raised on the real data
#' (outside the loop) is not suppressed.
#'
#' @examples
#' T <- 1:50
#' X <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
#' sim_w <- function(theta, T, m_traj) sim_proc(T, theta, "Wiener", m_traj)
#' res <- test_gof(X, T, estimer_wiener_homogene, sim_w, mp = 300, L = 50)
#' res$rejet
#'
#' @export
test_gof <- function(X, T, est, sim, mp = 100 * nrow(X), L = 1000, alpha = 0.05) {
    m <- nrow(X)
    th <- est(X, T)
    S_obs <- stats_sim(X, sim(th, T, mp))

    S_star <- matrix(NA_real_, L, length(S_obs), dimnames = list(NULL, names(S_obs)))
    for (v in seq_len(L)) {
          Xb <- sim(th, T, m)
          th_b <- suppressWarnings(est(Xb, T))
          S_star[v, ] <- stats_sim(Xb, sim(th_b, T, mp))
        }
    s_crit <- apply(S_star, 2, stats::quantile, probs = 1 - alpha, names = FALSE)
    list(rejet = S_obs > s_crit, S_obs = S_obs, s_crit = s_crit,
                p_valeur = colMeans(sweep(S_star, 2, S_obs, ">=")))
  }
