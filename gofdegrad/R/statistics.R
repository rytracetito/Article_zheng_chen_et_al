#' Similarity statistics between trajectories and a simulated pool
#'
#' @description
#' Computes the four test statistics of Zheng et al. (2026) between a set of
#' \eqn{m} trajectories and a pool of \eqn{m'} trajectories simulated under
#' the candidate model, all inspected at the same times. Each trajectory is
#' matched with its most similar trajectory of the pool; the statistics
#' summarize these smallest distances by their maximum (`sup`) or their sum
#' (`sum`).
#'
#' @param x Trajectories to test: a numeric matrix with one row per
#'   trajectory and one column per inspection time, or any data form accepted
#'   by [gof_data()] (for instance a `levysim` simulation).
#' @param pool Pool of simulated trajectories, in the same form as `x` and at
#'   the same inspection times.
#'
#' @details
#' For two trajectories \eqn{X = (x_1, \dots, x_n)} and
#' \eqn{Y = (y_1, \dots, y_n)} observed at the same times, the two
#' similarity measures are
#' \deqn{D_{KS}(X, Y) = \max_j |x_j - y_j|, \qquad D_{CvM}(X, Y) = \sum_j (x_j - y_j)^2,}
#' by analogy with the Kolmogorov-Smirnov and Cramér-von Mises statistics.
#' For each trajectory \eqn{X_i} of `x`, let
#' \eqn{d_i = \min_k D(X_i, Y_k)} be its distance to the most similar
#' trajectory of the pool. The statistics are
#' \deqn{S_{sup} = \max_i d_i, \qquad S_{sum} = \sum_i d_i,}
#' for each of the two measures. Large values mean that some (`sup`) or all
#' (`sum`) observed trajectories have no similar trajectory in the pool,
#' that is, are unlikely under the candidate model. With a single
#' trajectory, the `sup` and `sum` statistics are equal.
#'
#' @return A named numeric vector of length 4: `KS_sup`, `KS_sum`,
#'   `CvM_sup`, `CvM_sum`.
#'
#' @references
#' Zheng, J., Wang, Y., Chen, Y. and Xie, M. (2026). Automated verification
#' of degradation models using goodness-of-fit test based on similarity
#' measures and parametric bootstrap. *Quality and Reliability Engineering
#' International*. \doi{10.1002/qre.70382}
#'
#' @seealso [gof_test()], which computes these statistics on the data and
#'   on the bootstrap samples.
#'
#' @examples
#' # 1. Six gamma trajectories against a pool of 600 gamma trajectories
#' set.seed(1)
#' t <- 1:30
#' x <- sim_model(t, "gamma", c(a = 2, b = 3), n_traj = 6)
#' pool <- sim_model(t, "gamma", c(a = 2, b = 3), n_traj = 600)
#' similarity_stats(x, pool)
#'
#' # 2. The same trajectories against a pool from another model: larger values
#' pool_ig <- sim_model(t, "ig", c(mu = 2 / 3, eta = 0.5), n_traj = 600)
#' similarity_stats(x, pool_ig)
#'
#' # 3. Small matrices, one row per trajectory
#' X <- matrix(c(1, 2, 3,
#'               1, 3, 4), nrow = 2, byrow = TRUE)
#' Y <- matrix(c(1, 2, 2,
#'               2, 3, 5,
#'               0, 1, 3), nrow = 3, byrow = TRUE)
#' similarity_stats(X, Y)
#' @export
similarity_stats <- function(x, pool) {
  X <- as_level_matrix(x)
  Y <- as_level_matrix(pool)
  if (ncol(X) != ncol(Y)) stop("x and pool must have the same number of inspection times")
  d_ks <- apply(distance_matrix(X, Y, "KS"), 1, min)
  d_cvm <- apply(distance_matrix(X, Y, "CvM"), 1, min)
  c(KS_sup = max(d_ks), KS_sum = sum(d_ks), CvM_sup = max(d_cvm), CvM_sum = sum(d_cvm))
}

as_level_matrix <- function(x) {
  if (is.matrix(x) && !inherits(x, "levysim")) {
    if (!is.numeric(x)) stop("the trajectories must be numeric")
    return(x)
  }
  gof_data(x)$X
}

distance_matrix <- function(X, Y, measure) {
  if (measure == "CvM") {
    D <- outer(rowSums(X^2), rowSums(Y^2), "+") - 2 * X %*% t(Y)
    return(pmax(D, 0))
  }
  D <- matrix(0, nrow(X), nrow(Y))
  for (j in seq_len(ncol(X))) D <- pmax(D, abs(outer(X[, j], Y[, j], "-")))
  D
}
