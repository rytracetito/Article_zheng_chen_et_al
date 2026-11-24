r_bounds <- c(0.01, 20)

check_r <- function(r) {
  if (r < r_bounds[1] * 1.01 || r > r_bounds[2] / 1.01) {
    warning(sprintf("the estimate of r (%.4g) is at the boundary of the search interval [%g, %g]",
                    r, r_bounds[1], r_bounds[2]), call. = FALSE)
  }
}

wiener_given_trend <- function(inc, dL) {
  m <- nrow(inc$dX)
  n <- ncol(inc$dX)
  dt <- inc$dt[1, ]
  S <- colSums(inc$dX)
  mu <- sum(dL / dt * S) / (m * sum(dL^2 / dt))
  res <- inc$dX - matrix(mu * dL, m, n, byrow = TRUE)
  list(mu = mu, sigma2 = sum(res^2 / inc$dt) / (m * n))
}

fit_wiener_power <- function(d) {
  inc <- increments_of(d)
  n <- ncol(inc$dX)
  if (n < 2) stop("the power-law models need at least two inspections after time 0", call. = FALSE)
  tn <- d$t[n]
  u <- c(0, d$t / tn)
  given_r <- function(r) wiener_given_trend(inc, diff(u^r))
  opt <- stats::optimize(function(lr) log(given_r(exp(lr))$sigma2), log(r_bounds))
  r <- exp(opt$minimum)
  check_r(r)
  g <- given_r(r)
  c(mu = g$mu / tn^r, sigma = sqrt(g$sigma2), r = r)
}

kappa_bounds <- c(1e-3, 50)

fit_wiener_exp <- function(d) {
  inc <- increments_of(d)
  n <- ncol(inc$dX)
  if (n < 2) stop("the exponential-trend model needs at least two inspections after time 0", call. = FALSE)
  tn <- d$t[n]
  u <- c(0, d$t / tn)
  given_k <- function(k) wiener_given_trend(inc, diff(expm1(k * u)))
  opt <- stats::optimize(function(lk) log(given_k(exp(lk))$sigma2), log(kappa_bounds))
  k <- exp(opt$minimum)
  if (k < kappa_bounds[1] * 1.01) {
    warning(sprintf("the estimate of lambda (%.4g) is at the lower boundary of the search interval: ",
                    k / tn), "the trend is almost linear", call. = FALSE)
  } else if (k > kappa_bounds[2] / 1.01) {
    warning(sprintf("the estimate of lambda (%.4g) is at the upper boundary of the search interval",
                    k / tn), call. = FALSE)
  }
  g <- given_k(k)
  c(mu = g$mu, sigma = sqrt(g$sigma2), lambda = k / tn)
}

fit_gamma_power <- function(d) {
  inc <- positive_increments(d, "gamma")
  n <- ncol(inc$dX)
  if (n < 2) stop("the power-law models need at least two inspections after time 0", call. = FALSE)
  tn <- d$t[n]
  u <- c(0, d$t / tn)
  dx <- as.vector(inc$dX)
  col <- as.vector(col(inc$dX))

  given_r <- function(r) {
    dL <- diff(u^r)[col]
    e <- suppressWarnings(levyfit::gamma_mle(levyfit::increments_data(dx, dL)))
    if (anyNA(e)) return(list(a = NA, b = NA, loglik = -Inf))
    list(a = e[["a"]], b = e[["b"]],
         loglik = sum(stats::dgamma(dx, shape = e[["a"]] * dL, rate = e[["b"]], log = TRUE)))
  }
  opt <- stats::optimize(function(lr) {
    ll <- given_r(exp(lr))$loglik
    if (is.finite(ll)) -ll else .Machine$double.xmax
  }, log(r_bounds))
  r <- exp(opt$minimum)
  g <- given_r(r)
  if (is.na(g$a)) stop("the gamma model with power-law shape function cannot be fitted to these data", call. = FALSE)
  check_r(r)
  c(a = g$a / tn^r, b = g$b, r = r)
}
