#' Maximum likelihood estimation of a candidate model
#'
#' @description
#' Estimates the parameters of one of the candidate models of
#' [gof_models()] (or of a model built with [custom_model()]) from one or
#' several degradation trajectories inspected at the same, possibly
#' irregular, times. This is the estimation step of [gof_test()], applied
#' once to the observed data and once to each bootstrap sample.
#'
#' @param x Degradation data, in any form accepted by [gof_data()]: result
#'   of a `levysim` simulation, matrix with a column of times, matrix with
#'   one row per trajectory (with `t`), data frame (`id`, `t`, `X`) or list
#'   of data frames.
#' @param model Name of the model, one of `"wiener"`, `"gamma"`, `"ig"`,
#'   `"vg"`, `"meixner"`, `"wiener_power"`, `"gamma_power"` (see
#'   [gof_models()]), or a model built with [custom_model()].
#' @param t Inspection times, only when `x` is a matrix with one row per
#'   trajectory or a vector (see [gof_data()]).
#'
#' @details
#' **Estimators.** For the five Lévy models the estimators of the package
#' `levyfit` are used, applied to the increments of all the trajectories
#' pooled together:
#' * Wiener: closed form, [levyfit::wiener_mle()]; the returned `sigma` is
#'   the square root of the maximum likelihood estimate of \eqn{\sigma^2};
#' * gamma: root of one equation, [levyfit::gamma_mle()];
#' * inverse Gaussian: closed form, [levyfit::ig_mle()];
#' * variance-gamma: EM algorithm then quasi-Newton, [levyfit::vg_mle()];
#' * Meixner: profile likelihood, [levyfit::meixner_mle()].
#'
#' For the two power-law models, the likelihood is profiled in \eqn{r}:
#' * `"wiener_power"`: for a given \eqn{r}, with
#'   \eqn{\Delta\Lambda_j = t_j^r - t_{j-1}^r}, the estimators are explicit,
#'   \deqn{\hat\mu(r) = \frac{\sum_j (\Delta\Lambda_j/\Delta t_j) S_j}{m \sum_j \Delta\Lambda_j^2/\Delta t_j},
#'   \qquad \hat\sigma^2(r) = \frac{1}{mn}\sum_{i,j}\frac{(\Delta x_{ij} - \hat\mu(r)\Delta\Lambda_j)^2}{\Delta t_j},}
#'   where \eqn{S_j} is the sum of the increments of the \eqn{j}-th
#'   interval; the profile likelihood is maximal where \eqn{\hat\sigma^2(r)}
#'   is minimal;
#' * `"gamma_power"`: for a given \eqn{r}, the process is a homogeneous gamma
#'   process in the time \eqn{t^r}, estimated by [levyfit::gamma_mle()].
#'
#' \eqn{r} is searched by [stats::optimize()] in \eqn{[0.01, 20]}, on the
#' time scale \eqn{t/t_n} for numerical stability; a warning is given when
#' the estimate is at a boundary of this interval (for instance with a
#' single trajectory, whose trend is poorly identified).
#'
#' **Increasing models.** The gamma, inverse Gaussian and power-law gamma
#' models only increase: a negative increment in the data gives an error. An
#' increment exactly equal to 0 is usually a rounding effect: a very small
#' increment, added to a large level, disappears in the computer. It is
#' replaced by the machine precision times the largest level, an upper bound
#' of its true value, with a warning.
#'
#' **Estimation failures.** When the maximum likelihood estimator does not
#' exist for the data (for instance a variance-gamma model fitted to
#' increasing data, whose likelihood is maximal at the gamma limit),
#' `levyfit` gives a warning and the estimates are `NA`.
#'
#' @return A named numeric vector with the estimates of the parameters of
#'   the model (see [gof_models()]), with the attribute `model` (name of the
#'   model).
#'
#' @references
#' Zheng, J., Wang, Y., Chen, Y. and Xie, M. (2026). Automated verification
#' of degradation models using goodness-of-fit test based on similarity
#' measures and parametric bootstrap. *Quality and Reliability Engineering
#' International*. \doi{10.1002/qre.70382}
#'
#' @seealso [gof_models()], [sim_model()], [gof_test()], and the package
#'   `levyfit` for confidence intervals and more estimation options.
#'
#' @examples
#' # 1. Six trajectories of a gamma process, inspected every time unit
#' set.seed(1)
#' x <- sim_model(0:50, "gamma", c(a = 2, b = 3), n_traj = 6)
#' fit_model(x, "gamma")
#'
#' # 2. Several models fitted to the same data
#' fit_model(x, "ig")
#' fit_model(x, "wiener")
#'
#' # 3. Irregular inspection times (common to all the trajectories)
#' t <- levysim::inspection_times("increasing", n = 30, tmax = 50)$time
#' x <- sim_model(t, "wiener_power", c(mu = 0.5, sigma = 0.3, r = 0.6), n_traj = 10)
#' fit_model(x, "wiener_power")
#'
#' # 4. Gamma process with power-law shape function
#' x <- sim_model(0:40, "gamma_power", c(a = 1, b = 2, r = 1.5), n_traj = 10)
#' fit_model(x, "gamma_power")
#'
#' # 5. A matrix with one row per trajectory and the times given in t
#' X <- matrix(c(0.6, 1.3, 2.2, 3.1,
#'               0.4, 1.1, 1.9, 2.4,
#'               0.7, 1.2, 2.5, 3.3), nrow = 3, byrow = TRUE)
#' fit_model(X, "gamma", t = c(1, 2, 3.5, 5))
#'
#' # 6. A decrease makes the increasing models impossible
#' try(fit_model(c(0, 1, 0.8, 2), "gamma", t = c(0, 1, 2, 3)))
#'
#' # 7. Real data: GaAs laser (package SMRD.data)
#' if (requireNamespace("SMRD.data", quietly = TRUE)) {
#'   laser <- SMRD.data::gaaslaser
#'   laser <- data.frame(id = laser$unit, t = laser$hours, X = laser$increase)
#'   fit_model(laser, "wiener")
#'   fit_model(laser, "wiener_power")
#' }
#' @export
fit_model <- function(x, model, t = NULL) {
  mod <- get_model(model)
  d <- gof_data(x, t)
  est <- mod$fit(d)
  est <- check_estimate(est, mod)
  structure(est, model = mod$name)
}

check_estimate <- function(est, mod) {
  est <- unlist(est)
  if (!is.numeric(est) || is.null(names(est)) || !all(mod$parameters %in% names(est))) {
    stop("the estimation function of the ", mod$name, " model must return a named numeric vector with ",
         paste(mod$parameters, collapse = ", "), call. = FALSE)
  }
  est[mod$parameters]
}

#' Simulation of a candidate model
#'
#' @description
#' Simulates exactly one or several degradation trajectories of one of the
#' candidate models of [gof_models()] (or of a model built with
#' [custom_model()]) at given inspection times, common to all the
#' trajectories. This is the simulation step of [gof_test()], used to
#' generate the bootstrap samples and the pools of simulated trajectories.
#'
#' @param t Inspection times, strictly increasing, with or without time 0
#'   (the trajectories start from 0 at time 0).
#' @param model Name of the model (see [gof_models()]) or a model built with
#'   [custom_model()].
#' @param theta Parameters of the model: a named numeric vector or list, with
#'   the names given in [gof_models()] (for instance `c(a = 2, b = 3)` for the
#'   gamma model). The result of [fit_model()] can be given directly.
#' @param n_traj Number of trajectories (default 1).
#'
#' @details
#' The five Lévy models and the two power-law models are simulated by the
#' package `levysim` (exact simulation of the independent increments), except
#' the Meixner model: [levysim::sim_meixner()] draws the increments one at a
#' time, which is too slow for the thousands of trajectories of a bootstrap.
#' The Meixner increments of a same time step are therefore drawn together,
#' with the same law (numerical inversion of the distribution function with
#' `Runuran`, parameters \eqn{a}, \eqn{b}, \eqn{d \Delta t},
#' \eqn{m \Delta t}).
#'
#' The parametrizations are those of [gof_models()]: in particular the gamma
#' model uses the **rate** \eqn{b}, as in `levysim` and `levyfit`.
#'
#' @return An object of class `levysim`: a matrix whose first column `t`
#'   contains the inspection times (time 0 included) and the columns `X1`,
#'   ..., `Xn` the simulated levels, equal to 0 at time 0. It can be plotted
#'   with [levysim::plot_traj()], converted to a data frame with
#'   [levysim::traj_df()], and given directly to [fit_model()] and
#'   [gof_test()].
#'
#' @seealso [gof_models()], [fit_model()], [sim_path()] for general path
#'   models, and the package `levysim` for more simulation options
#'   (different inspection times per trajectory, other time
#'   transformations).
#'
#' @examples
#' # 1. Ten trajectories of a Wiener process, inspected every 0.5 time unit
#' set.seed(1)
#' x <- sim_model(seq(0, 10, by = 0.5), "wiener", c(mu = 0.5, sigma = 0.8), n_traj = 10)
#' x
#' levysim::plot_traj(x)
#'
#' # 2. The seven models, with the same inspection times
#' t <- 0:20
#' sim_model(t, "gamma", c(a = 2, b = 3), n_traj = 3)
#' sim_model(t, "ig", c(mu = 0.5, eta = 2), n_traj = 3)
#' sim_model(t, "vg", c(theta = 0.5, sigma = 1, nu = 0.5), n_traj = 3)
#' sim_model(t, "meixner", c(a = 1, b = 0.3, d = 1.5, m = 0.2), n_traj = 3)
#' sim_model(t, "wiener_power", c(mu = 1, sigma = 0.3, r = 0.5), n_traj = 3)
#' sim_model(t, "gamma_power", c(a = 1, b = 2, r = 1.5), n_traj = 3)
#'
#' # 3. Irregular inspection times, more and more frequent
#' t <- levysim::inspection_times("increasing", n = 25, tmax = 50)$time
#' x <- sim_model(t, "gamma", c(a = 2, b = 3), n_traj = 8)
#' levysim::plot_traj(x, "Gamma process, irregular inspections")
#'
#' # 4. Simulate from estimated parameters (parametric bootstrap)
#' theta <- fit_model(x, "gamma")
#' y <- sim_model(t, "gamma", theta, n_traj = 8)
#'
#' # 5. All the data in one data frame, for instance to save them
#' head(levysim::traj_df(y))
#' @export
sim_model <- function(t, model, theta, n_traj = 1) {
  mod <- get_model(model)
  theta <- check_theta(theta, mod)
  if (!is.numeric(t) || anyNA(t) || any(t < 0)) stop("t must be non-negative inspection times")
  if (t[1] == 0) t <- t[-1]
  if (length(t) == 0 || any(diff(t) <= 0)) stop("inspection times must be strictly increasing")
  if (!is.numeric(n_traj) || length(n_traj) != 1 || n_traj < 1) stop("n_traj must be a positive integer")
  X <- simulate_levels(mod, theta, t, n_traj)
  as_levysim(X, t, mod$name, theta)
}

simulate_levels <- function(mod, theta, t, n) {
  X <- mod$simulate(theta, t, n)
  X <- matrix(X, nrow = n)
  if (ncol(X) != length(t)) {
    stop("the simulation function of the ", mod$name, " model must return a matrix with ", length(t),
         " columns, one per inspection time", call. = FALSE)
  }
  X
}

as_levysim <- function(X, t, model, parameters) {
  out <- cbind(t = c(0, t), rbind(0, t(X)))
  colnames(out) <- c("t", paste0("X", seq_len(nrow(X))))
  class(out) <- c("levysim", "matrix", "array")
  attr(out, "model") <- model
  attr(out, "parameters") <- unlist(parameters)
  out
}

#' Custom candidate model
#'
#' @description
#' Builds a candidate model from its own estimation and simulation functions,
#' so that any parametric degradation model can be tested by [gof_test()],
#' fitted by [fit_model()], simulated by [sim_model()] and studied by
#' [power_study()]: for instance a general path model, a model with random
#' effects, or a Lévy model with another time transformation.
#'
#' @param name Name of the model, a character string used in the printed
#'   results.
#' @param parameters Names of the parameters, a character vector.
#' @param fit Estimation function `function(X, t)`, where `X` is the
#'   \eqn{m \times n} matrix of degradation levels (one row per trajectory,
#'   measured from time 0) and `t` the \eqn{n} inspection times after time
#'   0. It returns a named numeric vector of estimates, with the names given
#'   in `parameters`.
#' @param simulate Simulation function `function(theta, t, n)`, where `theta`
#'   is a named numeric vector of parameters, `t` the inspection times after
#'   time 0 and `n` the number of trajectories. It returns an
#'   \eqn{n \times} `length(t)` matrix of degradation levels, one row per
#'   trajectory.
#' @param increasing `TRUE` if the trajectories of the model only increase
#'   (for information only).
#'
#' @details
#' The validity of the test relies on the estimation function: it must be
#' consistent, and it should be the same as the one used in practice for the
#' model, because the bootstrap reproduces it on every bootstrap sample. The
#' two functions are called about \eqn{2B} times by [gof_test()], so they
#' should be fast.
#'
#' @return An object of class `gof_model`, to be given as argument `model`
#'   of [gof_test()], [fit_model()], [sim_model()] and [power_study()].
#'
#' @seealso [gof_models()] for the built-in models, [sim_path()] for general
#'   path models.
#'
#' @examples
#' # A linear general path model: X_i(t) = phi_i t + error, with
#' # phi_i ~ N(mu, s_phi^2) and measurement errors N(0, s_eps^2)
#' linear <- custom_model(
#'   name = "linear path",
#'   parameters = c("mu", "s_phi", "s_eps"),
#'   fit = function(X, t) {
#'     slopes <- apply(X, 1, function(x) sum(x * t) / sum(t^2))
#'     res <- X - outer(slopes, t)
#'     c(mu = mean(slopes), s_phi = stats::sd(slopes),
#'       s_eps = sqrt(sum(res^2) / (length(res) - nrow(X))))
#'   },
#'   simulate = function(theta, t, n) {
#'     phi <- stats::rnorm(n, theta[["mu"]], theta[["s_phi"]])
#'     outer(phi, t) + matrix(stats::rnorm(n * length(t), 0, theta[["s_eps"]]), n)
#'   }
#' )
#' linear
#'
#' # Data simulated from this model, fitted, then tested
#' set.seed(1)
#' x <- sim_model(1:20, linear, c(mu = 0.5, s_phi = 0.1, s_eps = 0.2), n_traj = 6)
#' fit_model(x, linear)
#' gof_test(x, linear, B = 50, n_pool = 200)
#'
#' # The same data tested against a Wiener process
#' gof_test(x, "wiener", B = 50, n_pool = 200)
#'
#' # A Wiener process whose trend and variance follow different time scales,
#' # X(t) = mu t^0.5 + sigma B(t^0.25), which is not among the built-in
#' # models: given the two time functions, the estimators are explicit
#' Lambda <- function(t) t^0.5
#' tau <- function(t) t^0.25
#' wiener_tt <- custom_model(
#'   name = "Wiener, trend t^0.5, variance t^0.25",
#'   parameters = c("mu", "sigma"),
#'   fit = function(X, t) {
#'     dX <- X - cbind(0, X[, -ncol(X), drop = FALSE])
#'     dL <- matrix(diff(Lambda(c(0, t))), nrow(X), length(t), byrow = TRUE)
#'     dT <- matrix(diff(tau(c(0, t))), nrow(X), length(t), byrow = TRUE)
#'     mu <- sum(dX * dL / dT) / sum(dL^2 / dT)
#'     c(mu = mu, sigma = sqrt(mean((dX - mu * dL)^2 / dT)))
#'   },
#'   simulate = function(theta, t, n) {
#'     x <- levysim::sim_wiener(c(0, t), mu = theta[["mu"]], sigma = theta[["sigma"]],
#'                              Lambda = Lambda, tau = tau, n_traj = n)
#'     gof_data(x)$X
#'   }
#' )
#' t <- levysim::inspection_times("increasing", n = 29, tmax = 35, digits = 2)$time
#' y <- sim_model(t, wiener_tt, c(mu = 2, sigma = 0.78), n_traj = 6)
#' fit_model(y, wiener_tt)
#' gof_test(y, wiener_tt, B = 50, n_pool = 300)
#' @export
custom_model <- function(name, parameters, fit, simulate, increasing = FALSE) {
  if (!is.character(name) || length(name) != 1) stop("name must be a character string")
  if (!is.character(parameters) || length(parameters) == 0) stop("parameters must be a character vector")
  if (!is.function(fit) || !is.function(simulate)) stop("fit and simulate must be functions")
  structure(list(
    name = name, process = name, parameters = parameters, increasing = isTRUE(increasing),
    law = "user-defined", estimation = "user-defined", simulation = "user-defined",
    fit = function(d) fit(d$X, d$t),
    simulate = simulate
  ), class = "gof_model")
}

#' @export
print.gof_model <- function(x, ...) {
  cat("Candidate model:", x$name, "\n")
  cat("process:", x$process, "\n")
  cat("parameters:", paste(x$parameters, collapse = ", "), "\n")
  if (x$law != "user-defined") cat("increments:", x$law, "\n")
  cat("estimation:", x$estimation, "\n")
  cat("simulation:", x$simulation, "\n")
  invisible(x)
}
