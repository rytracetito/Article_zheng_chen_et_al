#' Simulation of general path degradation models
#'
#' @description
#' Simulates trajectories of a general path model (Algorithm 3 of Zheng et
#' al., 2026): a deterministic degradation path with random effects, specific
#' to each unit, plus independent measurement errors. These models are the
#' usual alternatives to the stochastic process models in a power study (see
#' [power_study()]); they can be tested themselves with [custom_model()].
#'
#' @param t Inspection times, strictly increasing, with or without time 0.
#' @param shape Shape of the path: `"linear"`, `"power"` or `"exponential"`
#'   (see Details).
#' @param phi_mean,phi_sd Mean and standard deviation of the normal random
#'   effects: one number each for the linear and exponential paths, a vector
#'   of length 2 each for the power path (\eqn{\phi_1} then \eqn{\phi_2}).
#' @param sigma_eps Standard deviation of the measurement errors,
#'   \eqn{\sigma_\varepsilon \ge 0}.
#' @param omega Rate of the exponential path, a positive number; used only
#'   by `shape = "exponential"`.
#' @param n_traj Number of trajectories (default 1).
#'
#' @details
#' Each unit \eqn{k} has its own random effect \eqn{\phi_k}, drawn from a
#' normal distribution, and is observed with independent measurement errors
#' \eqn{\varepsilon_{kj} \sim N(0, \sigma_\varepsilon^2)}:
#' \deqn{X_k(t_j) = g(t_j; \phi_k) + \varepsilon_{kj},}
#' with
#' * `"linear"`: \eqn{g(t; \phi) = \phi t};
#' * `"power"`: \eqn{g(t; \phi) = \phi_1 t^{\phi_2}}, with independent
#'   \eqn{\phi_1} and \eqn{\phi_2};
#' * `"exponential"`: \eqn{g(t; \phi) = \phi (e^{\omega t} - 1)}.
#'
#' The level at time 0 is 0, without measurement error. The normality of the
#' random effects is an assumption of the package: the article does not
#' specify their distribution.
#'
#' @return An object of class `levysim`: a matrix whose first column `t`
#'   contains the inspection times (time 0 included) and the columns `X1`,
#'   ..., `Xn` the simulated levels. It can be plotted with
#'   [levysim::plot_traj()] and given to [fit_model()] and [gof_test()].
#'
#' @references
#' Zheng, J., Wang, Y., Chen, Y. and Xie, M. (2026). Automated verification
#' of degradation models using goodness-of-fit test based on similarity
#' measures and parametric bootstrap. *Quality and Reliability Engineering
#' International*. \doi{10.1002/qre.70382}
#'
#' Lu, C. J. and Meeker, W. Q. (1993). Using degradation measures to
#' estimate a time-to-failure distribution. *Technometrics*, 35(2), 161-174.
#' \doi{10.1080/00401706.1993.10485038}
#'
#' @seealso [sim_model()] for the stochastic process models,
#'   [power_study()], [custom_model()].
#'
#' @examples
#' t <- 0:50
#' set.seed(1)
#'
#' # 1. Linear paths with random slopes
#' x <- sim_path(t, "linear", phi_mean = 0.2, phi_sd = 0.05, sigma_eps = 0.1, n_traj = 6)
#' levysim::plot_traj(x)
#'
#' # 2. Power paths: concave when phi_2 < 1
#' x <- sim_path(t, "power", phi_mean = c(1, 0.5), phi_sd = c(0.2, 0.05),
#'               sigma_eps = 0.1, n_traj = 6)
#' levysim::plot_traj(x, "Power paths")
#'
#' # 3. Exponential paths
#' x <- sim_path(t, "exponential", phi_mean = 0.2, phi_sd = 0.05, sigma_eps = 0.1,
#'               omega = 0.05, n_traj = 6)
#' levysim::plot_traj(x, "Exponential paths")
#'
#' # 4. Irregular inspection times
#' t <- levysim::inspection_times("dec-inc", n = 20, tmax = 50)$time
#' x <- sim_path(t, "linear", phi_mean = 0.2, phi_sd = 0.05, sigma_eps = 0.1, n_traj = 6)
#' head(levysim::traj_df(x))
#'
#' # 5. Is a linear path model with random slopes taken for a Wiener process?
#' gof_test(x, "wiener", B = 50, n_pool = 200)
#' @export
sim_path <- function(t, shape = c("linear", "power", "exponential"), phi_mean, phi_sd,
                     sigma_eps, omega = NULL, n_traj = 1) {
  shape <- match.arg(shape)
  if (!is.numeric(t) || anyNA(t) || any(t < 0)) stop("t must be non-negative inspection times")
  if (t[1] == 0) t <- t[-1]
  if (length(t) == 0 || any(diff(t) <= 0)) stop("inspection times must be strictly increasing")
  k <- if (shape == "power") 2 else 1
  if (length(phi_mean) != k || length(phi_sd) != k) {
    stop("phi_mean and phi_sd must have length ", k, " for the ", shape, " path")
  }
  if (any(phi_sd < 0) || sigma_eps < 0) stop("phi_sd and sigma_eps must be non-negative")
  if (shape == "exponential" && (is.null(omega) || omega <= 0)) stop("omega must be a positive number")

  n <- n_traj
  tm <- matrix(t, n, length(t), byrow = TRUE)
  phi <- matrix(stats::rnorm(n * k, phi_mean, phi_sd), n, k, byrow = TRUE)
  g <- switch(shape,
    linear = phi[, 1] * tm,
    power = phi[, 1] * tm^phi[, 2],
    exponential = phi[, 1] * (exp(omega * tm) - 1)
  )
  X <- g + matrix(stats::rnorm(n * length(t), 0, sigma_eps), n)
  par <- c(phi_mean = phi_mean, phi_sd = phi_sd, sigma_eps = sigma_eps)
  if (shape == "exponential") par <- c(par, omega = omega)
  as_levysim(X, t, paste(shape, "path"), par)
}
