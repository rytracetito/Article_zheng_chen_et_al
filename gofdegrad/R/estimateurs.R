#' Estimation MLE - Wiener homogene
#'
#' Estime par maximum de vraisemblance les parametres d'un processus de
#' Wiener homogene (\eqn{\Lambda(t) = t}).
#'
#' @param X Matrice \eqn{m \times n} des trajectoires observees.
#' @param T Vecteur des \code{n} instants d'observation.
#'
#' @return Liste \code{list(mu_W, sigma)}.
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

#' Estimation MLE - Gamma homogene
#'
#' Estime par maximum de vraisemblance les parametres d'un processus Gamma
#' homogene (\eqn{\Lambda(t) = t}). Pas de forme close pour \eqn{\alpha}
#' (fonction digamma) : \eqn{\beta} est profile analytiquement, \eqn{\alpha}
#' recherche numeriquement en 1D.
#'
#' @param X Matrice \eqn{m \times n} des trajectoires observees.
#' @param T Vecteur des \code{n} instants d'observation.
#'
#' @return Liste \code{list(alpha, beta)}.
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
  dX <- X - cbind(0, X[, -n, drop = FALSE])
  if (any(dX <= 0)) stop("estimer_gamma_homogene : increments negatifs ou nuls presents, incompatible avec un modele Gamma")
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

#' Estimation MLE - Inverse Gaussien homogene
#'
#' Estime par maximum de vraisemblance les parametres d'un processus
#' Inverse Gaussien homogene (\eqn{\Lambda(t) = t}).
#'
#' @param X Matrice \eqn{m \times n} des trajectoires observees.
#' @param T Vecteur des \code{n} instants d'observation.
#'
#' @return Liste \code{list(mu_IG, eta)}.
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
  dX <- X - cbind(0, X[, -n, drop = FALSE])
  if (any(dX <= 0)) stop("estimer_ig_homogene : increments negatifs ou nuls presents, incompatible avec un modele IG")
  dt_mat <- matrix(dt, m, n, byrow = TRUE)
  A <- sum(dX); B <- m * t_n; C <- sum(dt_mat^2 / dX)
  mu_IG_hat <- A / B
  Q <- A - 2 * mu_IG_hat * B + mu_IG_hat^2 * C
  eta_hat <- (m * n) * mu_IG_hat^2 / Q
  list(mu_IG = mu_IG_hat, eta = eta_hat)
}

#' Estimation MLE - Wiener non-homogene
#'
#' Estime par maximum de vraisemblance les parametres d'un processus de
#' Wiener non-homogene, \eqn{\Lambda(t) = t^r}. \eqn{\hat\mu_W(r)} et
#' \eqn{\hat\sigma^2(r)} ont une forme close ; \eqn{r} est recherche
#' numeriquement en 1D (maximiser la vraisemblance profilee equivaut a
#' minimiser \eqn{\hat\sigma^2(r)}).
#'
#' @param X Matrice \eqn{m \times n} des trajectoires observees.
#' @param T Vecteur des \code{n} instants d'observation.
#'
#' @return Liste \code{list(mu_W, sigma, r)}.
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

#' Estimation MLE - Gamma non-homogene
#'
#' Estime par maximum de vraisemblance les parametres d'un processus Gamma
#' non-homogene, \eqn{\Lambda(t) = t^r}. \eqn{\beta} a une forme close a
#' \eqn{(\alpha, r)} fixes ; \eqn{(\alpha, r)} recherches numeriquement en 2D.
#'
#' @param X Matrice \eqn{m \times n} des trajectoires observees.
#' @param T Vecteur des \code{n} instants d'observation.
#'
#' @return Liste \code{list(alpha, beta, r)}.
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
  dX <- X - cbind(0, X[, -n, drop = FALSE])
  if (any(dX <= 0)) stop("estimer_gamma_nonhomog : increments negatifs ou nuls presents, incompatible avec un modele Gamma")
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

#' Dispatcher d'estimation MLE
#'
#' Selectionne et applique l'estimateur MLE correspondant au modele candidat
#' demande.
#'
#' @param X Matrice \eqn{m \times n} des trajectoires observees.
#' @param T Vecteur des \code{n} instants d'observation.
#' @param modele Chaine de caracteres : \code{"Wiener"}, \code{"Gamma"},
#'   \code{"IG"}, \code{"Wiener_nonhomog"} ou \code{"Gamma_nonhomog"}.
#'
#' @return La liste de parametres retournee par l'estimateur correspondant.
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
    stop("modele inconnu : ", modele)
  )
}
