model_catalogue <- function() {
  list(
    wiener = list(
      process = "Wiener", parameters = c("mu", "sigma"), increasing = FALSE,
      law = "X(t) - X(s) ~ N(mu (t - s), sigma^2 (t - s))",
      estimation = "levyfit::wiener_mle(), closed form",
      simulation = "DegradSim::sim_wiener()",
      fit = function(d) {
        inc <- increments_of(d)
        e <- levyfit::wiener_mle(levyfit::increments_data(as.vector(inc$dX), as.vector(inc$dt)))
        c(mu = e[["mu"]], sigma = sqrt(e[["sigma2"]]))
      },
      simulate = function(theta, t, n) {
        levels_of(DegradSim::sim_wiener(c(0, t), mu = theta[["mu"]], sigma = theta[["sigma"]], n_traj = n))
      }
    ),
    gamma = list(
      process = "gamma", parameters = c("a", "b"), increasing = TRUE,
      law = "X(t) - X(s) ~ Gamma(shape a (t - s), rate b)",
      estimation = "levyfit::gamma_mle(), root of one equation",
      simulation = "DegradSim::sim_gamma()",
      fit = function(d) {
        inc <- positive_increments(d, "gamma")
        e <- levyfit::gamma_mle(levyfit::increments_data(as.vector(inc$dX), as.vector(inc$dt)))
        c(a = e[["a"]], b = e[["b"]])
      },
      simulate = function(theta, t, n) {
        levels_of(DegradSim::sim_gamma(c(0, t), a = theta[["a"]], b = theta[["b"]], n_traj = n))
      }
    ),
    ig = list(
      process = "inverse Gaussian", parameters = c("mu", "eta"), increasing = TRUE,
      law = "X(t) - X(s) ~ IG(mean mu (t - s), shape eta (t - s)^2)",
      estimation = "levyfit::ig_mle(), closed form",
      simulation = "DegradSim::sim_ig()",
      fit = function(d) {
        inc <- positive_increments(d, "inverse Gaussian")
        e <- levyfit::ig_mle(levyfit::increments_data(as.vector(inc$dX), as.vector(inc$dt)))
        c(mu = e[["mu"]], eta = e[["eta"]])
      },
      simulate = function(theta, t, n) {
        levels_of(DegradSim::sim_ig(c(0, t), mu = theta[["mu"]], eta = theta[["eta"]], n_traj = n))
      }
    ),
    vg = list(
      process = "variance-gamma", parameters = c("theta", "sigma", "nu"), increasing = FALSE,
      law = "X(s) = theta G(s) + sigma W(G(s)), G gamma process of mean s and variance nu s",
      estimation = "levyfit::vg_mle(), EM algorithm then L-BFGS-B",
      simulation = "DegradSim::sim_vg()",
      fit = function(d) {
        inc <- increments_of(d)
        e <- levyfit::vg_mle(levyfit::increments_data(as.vector(inc$dX), as.vector(inc$dt)))
        c(theta = e[["theta"]], sigma = e[["sigma"]], nu = e[["nu"]])
      },
      simulate = function(theta, t, n) {
        levels_of(DegradSim::sim_vg(c(0, t), theta = theta[["theta"]], sigma = theta[["sigma"]],
                                  nu = theta[["nu"]], n_traj = n))
      }
    ),
    meixner = list(
      process = "Meixner", parameters = c("a", "b", "d", "m"), increasing = FALSE,
      law = "X(t) - X(s) ~ Meixner(a, b, d (t - s), m (t - s))",
      estimation = "levyfit::meixner_mle(), profile likelihood",
      simulation = "numerical inversion (Runuran), same law as DegradSim::sim_meixner()",
      fit = function(d) {
        inc <- increments_of(d)
        e <- levyfit::meixner_mle(levyfit::increments_data(as.vector(inc$dX), as.vector(inc$dt)))
        c(a = e[["a"]], b = e[["b"]], d = e[["d"]], m = e[["m"]])
      },
      simulate = function(theta, t, n) sim_meixner_levels(theta, t, n)
    ),
    wiener_power = list(
      process = "Wiener with power-law trend", parameters = c("mu", "sigma", "r"), increasing = FALSE,
      law = "X(t) - X(s) ~ N(mu (t^r - s^r), sigma^2 (t - s))",
      estimation = "closed form for (mu, sigma) given r, profile likelihood in r",
      simulation = "exact Gaussian increments, mean mu (t^r - s^r), variance sigma^2 (t - s)",
      fit = function(d) fit_wiener_power(d),
      simulate = function(theta, t, n) sim_wiener_power_levels(theta, t, n)
    ),
    gamma_power = list(
      process = "gamma with power-law shape function", parameters = c("a", "b", "r"), increasing = TRUE,
      law = "X(t) - X(s) ~ Gamma(shape a (t^r - s^r), rate b)",
      estimation = "levyfit::gamma_mle() given r, profile likelihood in r",
      simulation = "DegradSim::sim_gamma() with A(t) = a t^r",
      fit = function(d) fit_gamma_power(d),
      simulate = function(theta, t, n) {
        a <- theta[["a"]]
        r <- theta[["r"]]
        levels_of(DegradSim::sim_gamma(c(0, t), b = theta[["b"]], A = function(s) a * s^r, n_traj = n))
      }
    ),
    ig_re = list(
      process = "inverse Gaussian with random drift", parameters = c("xi", "sd_delta", "eta"),
      increasing = TRUE,
      law = paste("delta_i ~ N(xi, sd_delta^2) truncated to delta > 0, then",
                  "X_i(t) - X_i(s) ~ IG(mean (t - s) / delta_i, shape eta (t - s)^2)"),
      estimation = "exact marginal likelihood (model M3 of Peng, 2015), quasi-Newton",
      simulation = "truncated normal delta_i, then inverse Gaussian increments (Michael, Schucany and Haas)",
      fit = function(d) fit_ig_re(d),
      simulate = function(theta, t, n) sim_ig_re_levels(theta, t, n)
    )
  )
}

#' Candidate models available in the package
#'
#' @description
#' Lists the eight candidate models that can be tested by [gof_test()],
#' fitted by [fit_model()] and simulated by [sim_model()], with their
#' parameters, the law of their increments and the methods used for the
#' estimation and the simulation.
#'
#' @details
#' All the models start from \eqn{X(0) = 0}. The first five are Lévy
#' processes (independent and stationary increments): all the units follow
#' the same process. The next two have a power-law time transformation,
#' \eqn{t^r}, which gives a concave (\eqn{r < 1}) or convex (\eqn{r > 1})
#' mean degradation path. The last one is a random-effects model: each unit
#' degrades at its own mean rate, drawn from a common distribution.
#'
#' * `"wiener"`: Wiener process with drift,
#'   \eqn{X(t) = \mu t + \sigma B(t)}; increments can be negative.
#' * `"gamma"`: gamma process with shape \eqn{a} and **rate** \eqn{b}
#'   (scale \eqn{1/b}), \eqn{E[X(t)] = a t / b}; increasing.
#' * `"ig"`: inverse Gaussian process with the parametrization of Wang and
#'   Xu (2010), \eqn{E[X(t)] = \mu t},
#'   \eqn{Var(X(t)) = \mu^3 t / \eta}; increasing.
#' * `"vg"`: variance-gamma process, a Brownian motion with drift
#'   \eqn{\theta} and diffusion \eqn{\sigma} evaluated at a gamma time of
#'   variance \eqn{\nu} per unit of time; increments can be negative.
#' * `"meixner"`: Meixner process with parameters \eqn{a > 0},
#'   \eqn{-\pi < b < \pi}, \eqn{d > 0} and real \eqn{m}; increments can be
#'   negative.
#' * `"wiener_power"`: Wiener process with power-law trend,
#'   \eqn{X(t) = \mu t^r + \sigma B(t)}, the non-homogeneous Wiener model
#'   of Zheng et al. (2026): the trend is transformed, not the variance.
#' * `"gamma_power"`: gamma process with shape function \eqn{a t^r} and
#'   rate \eqn{b}; increasing.
#' * `"ig_re"`: inverse Gaussian process with random drift, the model
#'   \eqn{M_3} of Peng (2015). Unit \eqn{i} has its own mean degradation
#'   rate \eqn{\mu_i = 1/\delta_i}, where \eqn{\delta_i} follows the normal
#'   distribution \eqn{N(\xi, \sigma_\delta^2)} truncated to
#'   \eqn{\delta_i > 0}; given \eqn{\delta_i}, the trajectory is an inverse
#'   Gaussian process, \eqn{X_i(t) - X_i(s) \sim IG((t - s)/\delta_i, \eta (t - s)^2)}.
#'   \eqn{1/\xi} is a typical rate and \eqn{\sigma_\delta} measures the
#'   differences of rate between units; \eqn{\sigma_\delta \to 0} gives back
#'   the model `"ig"` with \eqn{\mu = 1/\xi}. Increasing.
#'
#'   Peng (2015) writes \eqn{\delta_i \sim N(\xi, \sigma_\mu^2/\lambda)}
#'   with \eqn{\lambda = \eta}, so that \eqn{\sigma_\mu = \sigma_\delta \sqrt{\eta}},
#'   and does not truncate \eqn{\delta_i}, assuming
#'   \eqn{P(\delta_i \le 0)} negligible; the truncation used here keeps the
#'   model valid in all cases and changes nothing when this probability is
#'   negligible. The random-effects inverse Gaussian model of Wang and Xu
#'   (2010) is different: it lets the shape parameter vary between units
#'   (model \eqn{M_2} of Peng, 2015), not the rate.
#'
#' The model `"ig"` assumes that all the units degrade at the same mean rate,
#' the differences between trajectories coming only from the randomness of
#' the increments: the spread of the trajectories then grows like
#' \eqn{\sqrt{t}}. When some units are consistently faster than others, the
#' spread grows like \eqn{t}, and `"ig_re"` is the suitable model.
#'
#' The increasing models (gamma, inverse Gaussian, power-law gamma,
#' random-effects inverse Gaussian) cannot be fitted to data with a negative
#' increment. An increment exactly equal to
#' 0, due to the rounding of the levels, is replaced by a tiny positive value
#' with a warning (see [fit_model()]).
#'
#' @return A data frame with one row per model and the columns `model` (name
#'   to give to the other functions), `process`, `parameters`, `increasing`
#'   (`TRUE` for the models with increasing trajectories), `law` (law of the
#'   increment between times \eqn{s < t}), `estimation` and `simulation`.
#'
#' @references
#' Zheng, J., Wang, Y., Chen, Y. and Xie, M. (2026). Automated verification
#' of degradation models using goodness-of-fit test based on similarity
#' measures and parametric bootstrap. *Quality and Reliability Engineering
#' International*. \doi{10.1002/qre.70382}
#'
#' Wang, X. and Xu, D. (2010). An inverse Gaussian process model for
#' degradation data. *Technometrics*, 52(2), 188-197.
#' \doi{10.1198/TECH.2009.08197}
#'
#' Peng, C.-Y. (2015). Inverse Gaussian processes with random effects and
#' explanatory variables for degradation data. *Technometrics*, 57(1),
#' 100-111. \doi{10.1080/00401706.2013.879077}
#'
#' Michael, J. R., Schucany, W. R. and Haas, R. W. (1976). Generating random
#' variates using transformations with multiple roots. *The American
#' Statistician*, 30(2), 88-90. \doi{10.1080/00031305.1976.10479147}
#'
#' @seealso [gof_test()], [fit_model()], [sim_model()], [custom_model()] to
#'   test another model.
#'
#' @examples
#' # All the models
#' gof_models()
#'
#' # Names and parameters only
#' gof_models()[, c("model", "parameters")]
#'
#' # Models with increasing trajectories
#' gof_models()$model[gof_models()$increasing]
#' @export
gof_models <- function() {
  cat_ <- model_catalogue()
  data.frame(
    model = names(cat_),
    process = vapply(cat_, `[[`, "", "process"),
    parameters = vapply(cat_, function(m) paste(m$parameters, collapse = ", "), ""),
    increasing = vapply(cat_, `[[`, TRUE, "increasing"),
    law = vapply(cat_, `[[`, "", "law"),
    estimation = vapply(cat_, `[[`, "", "estimation"),
    simulation = vapply(cat_, `[[`, "", "simulation"),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

get_model <- function(model) {
  if (inherits(model, "gof_model")) return(model)
  if (!is.character(model) || length(model) != 1) {
    stop("model must be the name of a model of gof_models() or the result of custom_model()")
  }
  cat_ <- model_catalogue()
  if (!model %in% names(cat_)) {
    stop("unknown model '", model, "'; available models: ", paste(names(cat_), collapse = ", "),
         "; use custom_model() for another model")
  }
  structure(c(list(name = model), cat_[[model]]), class = "gof_model")
}

levels_of <- function(x) {
  x <- unclass(x)
  t(x[-1, -1, drop = FALSE])
}

check_theta <- function(theta, model) {
  theta <- unlist(theta)
  if (!is.numeric(theta) || is.null(names(theta))) stop("theta must be a named numeric vector or list")
  missing <- setdiff(model$parameters, names(theta))
  if (length(missing)) {
    stop("missing parameter(s) for the ", model$name, " model: ", paste(missing, collapse = ", "))
  }
  theta[model$parameters]
}

sim_meixner_levels <- function(theta, t, n) {
  a <- theta[["a"]]
  b <- theta[["b"]]
  d <- theta[["d"]]
  m <- theta[["m"]]
  if (a <= 0 || d <= 0) stop("a and d must be positive")
  if (abs(b) >= pi) stop("b must be in (-pi, pi)")
  dt <- diff(c(0, t))
  dX <- matrix(0, n, length(t))
  for (s in unique(dt)) {
    cols <- which(dt == s)
    gen <- Runuran::pinvd.new(Runuran::udmeixner(alpha = a, beta = b, delta = d * s, mu = m * s))
    dX[, cols] <- Runuran::ur(gen, n * length(cols))
  }
  cumulate(dX)
}

sim_wiener_power_levels <- function(theta, t, n) {
  if (theta[["sigma"]] <= 0) stop("sigma must be positive")
  if (theta[["r"]] <= 0) stop("r must be positive")
  dL <- diff(c(0, t)^theta[["r"]])
  dt <- diff(c(0, t))
  dX <- matrix(stats::rnorm(n * length(t)), n, length(t)) * matrix(theta[["sigma"]] * sqrt(dt), n, length(t), byrow = TRUE) +
    matrix(theta[["mu"]] * dL, n, length(t), byrow = TRUE)
  cumulate(dX)
}

cumulate <- function(dX) {
  X <- dX
  for (j in seq_len(ncol(dX))[-1]) X[, j] <- X[, j - 1] + dX[, j]
  X
}
