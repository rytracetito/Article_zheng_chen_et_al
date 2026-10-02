#' Strictly positive increments for the Gamma and IG models
#'
#' Computes the trajectory increments. A strictly negative increment
#' raises an error (data incompatible with a monotone process). An
#' increment of exactly zero is treated as machine-precision rounding: a
#' Gamma or IG increment that is very small compared to the level of the
#' trajectory can vanish during the cumulative sum. It is replaced by
#' \eqn{\epsilon_{machine} \max|X|}, an upper bound on its true value,
#' with a warning.
#'
#' @param X Matrix of \eqn{m \times n} trajectories.
#' @param nom Name of the calling estimator, reused in the messages.
#' @param modele Name of the model, reused in the messages.
#'
#' @return Matrix of \eqn{m \times n} increments, strictly positive.
#'
#' @keywords internal
increments_positifs <- function(X, nom, modele) {
    n <- ncol(X)
    dX <- X - cbind(0, X[, -n, drop = FALSE])
    if (any(dX < 0)) {
          stop(nom, ": negative increments present, incompatible with a ", modele, " model", call. = FALSE)
        }
    nuls <- (dX == 0)
    if (any(nuls)) {
          plancher <- max(abs(X)) * .Machine$double.eps
          if (plancher <= 0) stop(nom, ": no strictly positive increment", call. = FALSE)
          dX[nuls] <- plancher
          warning(sprintf("%s: %d zero increment(s) replaced by %.3g (machine precision)",
                                               nom, sum(nuls), plancher), call. = FALSE)
        }
    dX
  }

#' MLE estimation - homogeneous Wiener
#'
#' Maximum likelihood estimation of the parameters of a homogeneous Wiener
#' process (\eqn{\Lambda(t) = t}).
#'
#' @param X Matrix of \eqn{m \times n} observed trajectories.
#' @param T Vector of the \code{n} observation times.
#'
#' @return List \code{list(mu_W, sigma)}.
#'
#' @examples
#' T <- 1:50
#' X <- sim_proc(T, list(mu_W = 0.3, sigma = 0.08), "Wiener", mp = 2000)
#' estimer_wiener_homogene(X, T)
#'
#' @export
estimer_wiener_homogene <- function(X, T) {
    m <- nrow(X); n <- ncol(X)
    dt <- diff(c(0, T)); t_n <- T[n]
    dX <- X - cbind(0, X[, -n, drop = FALSE])
    dt_mat <- matrix(dt, m, n, byrow = TRUE)
    mu_W_hat <- sum(dX) / (m * t_n)
    sigma2_hat <- sum((dX - mu_W_hat * dt_mat)^2 / dt_mat) / (m * n)
    list(mu_W = mu_W_hat, sigma = sqrt(sigma2_hat))
  }

#' MLE estimation - homogeneous Gamma
#'
#' Maximum likelihood estimation of the parameters of a homogeneous Gamma
#' process (\eqn{\Lambda(t) = t}). There is no closed form for \eqn{\alpha}
#' (digamma function): \eqn{\beta} is profiled analytically, \eqn{\alpha}
#' is searched numerically in 1D.
#'
#' @param X Matrix of \eqn{m \times n} observed trajectories.
#' @param T Vector of the \code{n} observation times.
#'
#' @return List \code{list(alpha, beta)}.
#'
#' @examples
#' T <- 1:50
#' X <- sim_proc(T, list(alpha = 1.5, beta = 0.7), "Gamma", mp = 2000)
#' estimer_gamma_homogene(X, T)
#'
#' @export
estimer_gamma_homogene <- function(X, T) {
    m <- nrow(X); n <- ncol(X)
    dt <- diff(c(0, T)); t_n <- T[n]
    dX <- increments_positifs(X, "estimer_gamma_homogene", "Gamma")
    dt_mat <- matrix(dt, m, n, byrow = TRUE)
    S1 <- sum(dX)
    beta_de_alpha <- function(alpha) S1 / (alpha * m * t_n)
    neg_profile_loglik <- function(log_alpha) {
          alpha <- exp(log_alpha); beta <- beta_de_alpha(alpha)
          -sum(stats::dgamma(dX, shape = alpha * dt_mat, scale = beta, log = TRUE))
        }
    opt <- stats::optimize(neg_profile_loglik, interval = log(c(1e-3, 1e3)))
    alpha_hat <- exp(opt$minimum)
    list(alpha = alpha_hat, beta = beta_de_alpha(alpha_hat))
  }

#' MLE estimation - homogeneous Inverse Gaussian
#'
#' Maximum likelihood estimation of the parameters of a homogeneous
#' Inverse Gaussian process (\eqn{\Lambda(t) = t}).
#'
#' @param X Matrix of \eqn{m \times n} observed trajectories.
#' @param T Vector of the \code{n} observation times.
#'
#' @return List \code{list(mu_IG, eta)}.
#'
#' @examples
#' T <- 1:50
#' X <- sim_proc(T, list(mu_IG = 0.8, eta = 1.2), "IG", mp = 2000)
#' estimer_ig_homogene(X, T)
#'
#' @export
estimer_ig_homogene <- function(X, T) {
    m <- nrow(X); n <- ncol(X)
    dt <- diff(c(0, T)); t_n <- T[n]
    dX <- increments_positifs(X, "estimer_ig_homogene", "IG")
    dt_mat <- matrix(dt, m, n, byrow = TRUE)
    A <- sum(dX); B <- m * t_n; C <- sum(dt_mat^2 / dX)
    mu_IG_hat <- A / B
    Q <- A - 2 * mu_IG_hat * B + mu_IG_hat^2 * C
    eta_hat <- (m * n) * mu_IG_hat^2 / Q
    list(mu_IG = mu_IG_hat, eta = eta_hat)
  }

#' MLE estimation - non-homogeneous Wiener
#'
#' Maximum likelihood estimation of the parameters of a non-homogeneous
#' Wiener process, \eqn{\Lambda(t) = t^r}. \eqn{\hat\mu_W(r)} and
#' \eqn{\hat\sigma^2(r)} have a closed form; \eqn{r} is searched
#' numerically in 1D (maximizing the profile likelihood is equivalent to
#' minimizing \eqn{\hat\sigma^2(r)}).
#'
#' @param X Matrix of \eqn{m \times n} observed trajectories.
#' @param T Vector of the \code{n} observation times.
#'
#' @return List \code{list(mu_W, sigma, r)}.
#'
#' @examples
#' T <- 1:50
#' th0 <- list(mu_W = 0.15, sigma = 0.4)
#' X <- sim_proc(T, th0, "Wiener", mp = 2000, Lam = function(t) t^1.7)
#' estimer_wiener_nonhomog(X, T)
#'
#' @export
estimer_wiener_nonhomog <- function(X, T) {
    m <- nrow(X); n <- ncol(X)
    dt <- diff(c(0, T))
    dX <- X - cbind(0, X[, -n, drop = FALSE])
    Sj <- colSums(dX)

    sigma2_hat_r <- function(r) {
          dLam <- diff(c(0, T)^r)
          w <- dLam / dt
          mu_W_hat <- sum(w * Sj) / (m * sum(dLam^2 / dt))
          sum((dX - matrix(mu_W_hat * dLam, m, n, byrow = TRUE))^2 /
                        matrix(dt, m, n, byrow = TRUE)) / (m * n)
        }

    opt <- stats::optimize(sigma2_hat_r, interval = c(0.05, 10))
    r_hat <- opt$minimum
    dLam <- diff(c(0, T)^r_hat)
    w <- dLam / dt
    mu_W_hat <- sum(w * Sj) / (m * sum(dLam^2 / dt))
    list(mu_W = mu_W_hat, sigma = sqrt(sigma2_hat_r(r_hat)), r = r_hat)
  }

#' MLE estimation - non-homogeneous Gamma
#'
#' Maximum likelihood estimation of the parameters of a non-homogeneous
#' Gamma process, \eqn{\Lambda(t) = t^r}. \eqn{\beta} has a closed form
#' for fixed \eqn{(\alpha, r)}; \eqn{(\alpha, r)} are searched numerically
#' in 2D.
#'
#' @param X Matrix of \eqn{m \times n} observed trajectories.
#' @param T Vector of the \code{n} observation times.
#'
#' @return List \code{list(alpha, beta, r)}.
#'
#' @examples
#' T <- 1:50
#' th0 <- list(alpha = 0.8, beta = 0.5)
#' X <- sim_proc(T, th0, "Gamma", mp = 2000, Lam = function(t) t^1.4)
#' estimer_gamma_nonhomog(X, T)
#'
#' @export
estimer_gamma_nonhomog <- function(X, T) {
    m <- nrow(X); n <- ncol(X)
    dt <- diff(c(0, T))
    dX <- increments_positifs(X, "estimer_gamma_nonhomog", "Gamma")
    S1 <- sum(dX)

    neg_profile_loglik <- function(par) {
          alpha <- exp(par[1]); r <- exp(par[2])
          dLam <- diff(c(0, T)^r)
          beta <- S1 / (alpha * m * sum(dLam))
          dLam_mat <- matrix(dLam, m, n, byrow = TRUE)
          -sum(stats::dgamma(dX, shape = alpha * dLam_mat, scale = beta, log = TRUE))
        }

    opt <- stats::optim(c(0, 0), neg_profile_loglik, method = "Nelder-Mead")
    alpha_hat <- exp(opt$par[1]); r_hat <- exp(opt$par[2])
    dLam <- diff(c(0, T)^r_hat)
    beta_hat <- S1 / (alpha_hat * m * sum(dLam))
    list(alpha = alpha_hat, beta = beta_hat, r = r_hat)
  }

#' MLE estimation dispatcher
#'
#' Selects and applies the MLE estimator corresponding to the requested
#' candidate model.
#'
#' @param X Matrix of \eqn{m \times n} observed trajectories.
#' @param T Vector of the \code{n} observation times.
#' @param modele Character string: \code{"Wiener"}, \code{"Gamma"},
#'   \code{"IG"}, \code{"Wiener_nonhomog"} or \code{"Gamma_nonhomog"}.
#'
#' @return The parameter list returned by the corresponding estimator.
#'
#' @examples
#' T <- 1:50
#' X <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 2000)
#' estimer_parametres(X, T, "Wiener")
#'
#' @export
estimer_parametres <- function(X, T, modele) {
    switch(modele,
               Wiener          = estimer_wiener_homogene(X, T),
               Gamma           = estimer_gamma_homogene(X, T),
               IG              = estimer_ig_homogene(X, T),
               Wiener_nonhomog = estimer_wiener_nonhomog(X, T),
               Gamma_nonhomog  = estimer_gamma_nonhomog(X, T),
               stop("unknown model: ", modele)
             )
  }
