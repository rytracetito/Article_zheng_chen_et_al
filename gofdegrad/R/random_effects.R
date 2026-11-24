ig_re_loglik <- function(par, inc, T) {
  xi <- par[1]
  s <- exp(par[2])
  eta <- exp(par[3])
  S <- rowSums(inc$dX)
  Q <- rowSums(inc$dt^2 / inc$dX)
  A <- eta * S + 1 / s^2
  B <- eta * T + xi / s^2
  C <- xi^2 / s^2
  sum(0.5 * log(eta * inc$dt^2 / (2 * pi * inc$dX^3))) +
    sum(-eta * Q / 2 - log(s) - 0.5 * log(A) + 0.5 * (B^2 / A - C) +
          stats::pnorm(B / sqrt(A), log.p = TRUE) - stats::pnorm(xi / s, log.p = TRUE))
}

fit_ig_re <- function(d) {
  inc <- positive_increments(d, "inverse Gaussian")
  m <- nrow(inc$dX)
  if (m < 2) stop("the random-effects model needs at least two trajectories", call. = FALSE)
  T <- d$t[length(d$t)]
  S <- rowSums(inc$dX)
  Q <- rowSums(inc$dt^2 / inc$dX)

  delta <- T / S
  eta0 <- length(inc$dX) / sum(Q - T^2 / S)
  within <- mean(delta / (eta0 * T) + 2 / (eta0^2 * T^2))
  s0 <- sqrt(max(stats::var(delta) - within, 0.01 * stats::var(delta)))
  start <- c(mean(delta), log(s0), log(eta0))

  opt <- stats::optim(start, function(p) -ig_re_loglik(p, inc, T), method = "BFGS",
                      control = list(parscale = c(abs(start[1]), 1, 1), maxit = 500, reltol = 1e-12))
  if (opt$convergence != 0 || !is.finite(opt$value)) {
    stop("the maximization of the likelihood of the random-effects model did not converge", call. = FALSE)
  }
  c(xi = opt$par[1], sd_delta = exp(opt$par[2]), eta = exp(opt$par[3]))
}

sim_ig_re_levels <- function(theta, t, n) {
  xi <- theta[["xi"]]
  s <- theta[["sd_delta"]]
  eta <- theta[["eta"]]
  if (s <= 0 || eta <= 0) stop("sd_delta and eta must be positive")
  u <- stats::runif(n, stats::pnorm(-xi / s), 1)
  delta <- xi + s * stats::qnorm(u)
  dt <- matrix(diff(c(0, t)), n, length(t), byrow = TRUE)
  cumulate(rig(dt / delta, eta * dt^2))
}

rig <- function(mean, shape) {
  w <- mean * stats::rnorm(length(mean))^2
  s <- sqrt(w^2 + 4 * shape * w)
  x1 <- mean * 4 * shape * w / (s + w)^2
  x1 <- ifelse(w == 0, mean, x1)
  out <- ifelse(stats::runif(length(mean)) <= mean / (mean + x1), x1, mean^2 / x1)
  matrix(out, nrow = nrow(mean))
}
