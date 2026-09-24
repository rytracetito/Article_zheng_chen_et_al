estimer_wiener_homogene <- function(X, T) {
  m <- nrow(X); n <- ncol(X)
  dt <- diff(c(0, T)); t_n <- T[n]
  dX <- X - cbind(0, X[, -n, drop = FALSE])
  dt_mat <- matrix(dt, m, n, byrow = TRUE)
  mu_W_hat <- sum(dX) / (m * t_n)
  sigma2_hat <- sum((dX - mu_W_hat * dt_mat)^2 / dt_mat) / (m * n)
  list(mu_W = mu_W_hat, sigma = sqrt(sigma2_hat))
}
estimer_gamma_homogene <- function(X, T) {
  m <- nrow(X); n <- ncol(X)
  dt <- diff(c(0, T)); t_n <- T[n]
  dX <- X - cbind(0, X[, -n, drop = FALSE])
  dt_mat <- matrix(dt, m, n, byrow = TRUE)
  S1 <- sum(dX)
  beta_de_alpha <- function(alpha) S1 / (alpha * m * t_n)
  neg_profile_loglik <- function(log_alpha) {
    alpha <- exp(log_alpha); beta <- beta_de_alpha(alpha)
    -sum(dgamma(dX, shape = alpha * dt_mat, scale = beta, log = TRUE))
  }
  opt <- optimize(neg_profile_loglik, interval = log(c(1e-3, 1e3)))
  alpha_hat <- exp(opt$minimum)
  list(alpha = alpha_hat, beta = beta_de_alpha(alpha_hat))
}


estimer_ig_homogene <- function(X, T) {
  m <- nrow(X); n <- ncol(X)
  dt <- diff(c(0, T)); t_n <- T[n]
  dX <- X - cbind(0, X[, -n, drop = FALSE])
  dt_mat <- matrix(dt, m, n, byrow = TRUE)
  A <- sum(dX); B <- m * t_n; C <- sum(dt_mat^2 / dX)
  mu_IG_hat <- A / B
  Q <- A - 2 * mu_IG_hat * B + mu_IG_hat^2 * C
  eta_hat <- (m * n) * mu_IG_hat^2 / Q
  list(mu_IG = mu_IG_hat, eta = eta_hat)
}
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
  opt <- optimize(sigma2_hat_r, interval = c(0.05, 10))
  r_hat <- opt$minimum
  dLam <- diff(c(0, T)^r_hat)
  w <- dLam / dt
  mu_W_hat <- sum(w * Sj) / (m * sum(dLam^2 / dt))
  list(mu_W = mu_W_hat, sigma = sqrt(sigma2_hat_r(r_hat)), r = r_hat)
}
estimer_gamma_nonhomog <- function(X, T) {
  m <- nrow(X); n <- ncol(X)
  dt <- diff(c(0, T))
  dX <- X - cbind(0, X[, -n, drop = FALSE])
  S1 <- sum(dX)
  neg_profile_loglik <- function(par) {
    alpha <- exp(par[1]); r <- exp(par[2])
    dLam <- diff(c(0, T)^r)
    beta <- S1 / (alpha * m * sum(dLam))
    dLam_mat <- matrix(dLam, m, n, byrow = TRUE)
    -sum(dgamma(dX, shape = alpha * dLam_mat, scale = beta, log = TRUE))
  }
  opt <- optim(c(0, 0), neg_profile_loglik, method = "Nelder-Mead")
  alpha_hat <- exp(opt$par[1]); r_hat <- exp(opt$par[2])
  dLam <- diff(c(0, T)^r_hat)
  beta_hat <- S1 / (alpha_hat * m * sum(dLam))
  list(alpha = alpha_hat, beta = beta_hat, r = r_hat)
}
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
