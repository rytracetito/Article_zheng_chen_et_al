#' Goodness-of-fit test by similarity measures and parametric bootstrap
#'
#' @description
#' Tests whether degradation trajectories follow a candidate model, with the
#' method of Zheng, Wang, Chen and Xie (2026) (their Algorithm 1). The
#' observed trajectories are compared with a large pool of trajectories
#' simulated under the fitted model through four similarity statistics, and
#' the distribution of these statistics under the null hypothesis is
#' obtained by parametric bootstrap.
#'
#' @param x Degradation data, in any form accepted by [gof_data()]: result
#'   of a `levysim` simulation, matrix with a column of times, matrix with
#'   one row per trajectory (with `t`), data frame (`id`, `t`, `X`) or list
#'   of data frames. All the trajectories must share the same inspection
#'   times, which can be irregular.
#' @param model Candidate model under the null hypothesis: the name of a
#'   model of [gof_models()] (`"wiener"`, `"gamma"`, `"ig"`, `"vg"`,
#'   `"meixner"`, `"wiener_power"`, `"gamma_power"`) or a model built with
#'   [custom_model()].
#' @param t Inspection times, only when `x` is a matrix with one row per
#'   trajectory or a vector (see [gof_data()]).
#' @param B Number of bootstrap samples, \eqn{L} in the article (default
#'   1000).
#' @param n_pool Number of trajectories of each simulated pool, \eqn{m'} in
#'   the article; by default \eqn{100 m}, the value recommended by the
#'   authors, \eqn{m} being the number of observed trajectories.
#' @param alpha Significance level of the test (default 0.05).
#' @param progress `TRUE` to display a progress bar during the bootstrap.
#'
#' @details
#' **Null hypothesis.** \eqn{H_0}: the \eqn{m} trajectories are independent
#' realizations of the candidate model \eqn{M(\theta)}, for some unknown
#' parameter \eqn{\theta}.
#'
#' **Algorithm** (Zheng et al., 2026, Algorithm 1).
#' 1. Estimate \eqn{\hat\theta} by maximum likelihood from the data (see
#'    [fit_model()]).
#' 2. Simulate a pool of \eqn{m'} trajectories under \eqn{M(\hat\theta)} at
#'    the inspection times of the data, and compute the four statistics
#'    \eqn{S} between the data and the pool (see [similarity_stats()]).
#' 3. For \eqn{b = 1, \dots, B}: simulate a bootstrap sample of \eqn{m}
#'    trajectories under \eqn{M(\hat\theta)}, estimate \eqn{\hat\theta^*_b}
#'    from it, simulate a pool of \eqn{m'} trajectories under
#'    \eqn{M(\hat\theta^*_b)} and compute the statistics \eqn{S^*_b}.
#'    Re-estimating the parameters on each bootstrap sample reproduces the
#'    effect of the estimation on the statistics.
#' 4. For each statistic, the critical value is the empirical quantile of
#'    order \eqn{1 - \alpha} of the \eqn{S^*_b}, the p-value is the
#'    proportion of \eqn{S^*_b \ge S}, and \eqn{H_0} is rejected when
#'    \eqn{S} exceeds the critical value.
#'
#' **Choice of the statistics.** The `sup` statistics react to one
#' trajectory far from the model, the `sum` statistics to a misfit shared by
#' all the trajectories. The authors recommend the `sup` statistics for an
#' automated use; looking at the four is more informative.
#'
#' **Bootstrap failures.** When the estimation fails on a bootstrap sample
#' (for instance no maximum likelihood estimate for a variance-gamma or
#' Meixner sample), this sample is discarded and counted in `n_failed`; the
#' critical values and p-values use the other samples. Warnings raised on the
#' bootstrap samples are not shown. A failure on the observed data stops the
#' test with an error explaining why the model cannot be fitted.
#'
#' **Computing time.** The test calls the estimation and the simulation of
#' the model about \eqn{2B} times. With the default values it takes a few
#' seconds for the Wiener, gamma and inverse Gaussian models, and one or a
#' few minutes for the variance-gamma, Meixner and power-law models. Use a
#' small `B` (50 to 200) to explore, then `B = 1000` for the final results,
#' after `set.seed()` for reproducibility.
#'
#' **Power.** With few trajectories or few inspections, the test has little
#' power: a model that does not fit the data well, for instance a
#' homogeneous model fitted to concave paths, may not be rejected because
#' its estimated variance is inflated. A non-rejection means that the data
#' are compatible with the model, not that the model is validated. Use
#' [power_study()] to assess the power of the test in a given design.
#'
#' @return An object of class `gof_test`, a list with
#'   * `model`, `process`: name of the candidate model and of its process;
#'   * `estimate`: estimates of the parameters from the data;
#'   * `statistic`: observed values of the four statistics;
#'   * `critical`: critical values;
#'   * `p.value`: p-values;
#'   * `reject`: `TRUE` when \eqn{H_0} is rejected, for each statistic;
#'   * `alpha`, `B`, `n_pool`: settings of the test;
#'   * `n_failed`: number of discarded bootstrap samples;
#'   * `m`, `t`: number of trajectories and inspection times of the data;
#'   * `boot`: matrix of the bootstrap statistics, one row per sample.
#'
#' Printing the object shows a summary table; [as.data.frame()] gives this
#' table as a data frame and [plot()] draws the bootstrap distributions.
#'
#' @references
#' Zheng, J., Wang, Y., Chen, Y. and Xie, M. (2026). Automated verification
#' of degradation models using goodness-of-fit test based on similarity
#' measures and parametric bootstrap. *Quality and Reliability Engineering
#' International*. \doi{10.1002/qre.70382}
#'
#' @seealso [compare_models()] to test several models on the same data,
#'   [power_study()] for the level and the power of the test,
#'   [gof_models()], [fit_model()], [similarity_stats()],
#'   [plot.gof_test()].
#'
#' @examples
#' # Small values of B and n_pool to keep the examples fast;
#' # use the default values (B = 1000, n_pool = 100 m) in practice.
#'
#' # 1. Gamma data tested against the gamma model: H0 is true
#' set.seed(1)
#' x <- sim_model(0:30, "gamma", c(a = 2, b = 3), n_traj = 6)
#' res <- gof_test(x, "gamma", B = 100, n_pool = 300)
#' res
#'
#' # 2. The elements of the result
#' res$estimate
#' res$p.value
#' res$reject
#' as.data.frame(res)
#'
#' # 3. Bootstrap distributions with the observed values (solid lines)
#' #    and the critical values (dashed lines)
#' plot(res)
#'
#' # 4. Concave Wiener paths tested against the homogeneous Wiener model
#' x <- sim_model(1:30, "wiener_power", c(mu = 2, sigma = 0.3, r = 0.5), n_traj = 6)
#' gof_test(x, "wiener", B = 100, n_pool = 300)
#' gof_test(x, "wiener_power", B = 50, n_pool = 300)
#'
#' # 5. Irregular inspection times common to all the trajectories
#' t <- levysim::inspection_times("increasing", n = 20, tmax = 40)$time
#' x <- sim_model(t, "ig", c(mu = 0.5, eta = 2), n_traj = 8)
#' gof_test(x, "ig", B = 100, n_pool = 400)
#'
#' # 6. Data given as a matrix, one row per trajectory, with the times
#' X <- gof_data(x)$X
#' gof_test(X, "ig", t = gof_data(x)$t, B = 50, n_pool = 400)
#'
#' \donttest{
#' # 7. Real data: GaAs laser (package SMRD.data); the homogeneous Wiener
#' #    model is rejected, as in Zheng et al. (2026)
#' if (requireNamespace("SMRD.data", quietly = TRUE)) {
#'   laser <- SMRD.data::gaaslaser
#'   laser <- data.frame(id = laser$unit, t = laser$hours, X = laser$increase)
#'   set.seed(30)
#'   gof_test(laser, "wiener", B = 200)
#' }
#' }
#' @export
gof_test <- function(x, model, t = NULL, B = 1000, n_pool = NULL, alpha = 0.05, progress = FALSE) {
  mod <- get_model(model)
  d <- gof_data(x, t)
  m <- nrow(d$X)
  if (is.null(n_pool)) n_pool <- 100 * m
  check_settings(B, n_pool, alpha)

  theta <- tryCatch(check_estimate(mod$fit(d), mod), error = function(e) {
    stop("the ", mod$name, " model cannot be fitted to the data: ", conditionMessage(e), call. = FALSE)
  })
  if (anyNA(theta)) {
    stop("the ", mod$name, " model cannot be fitted to the data: no maximum likelihood estimate",
         call. = FALSE)
  }
  S <- similarity_stats(d$X, simulate_levels(mod, theta, d$t, n_pool))

  boot <- matrix(NA_real_, B, 4, dimnames = list(NULL, names(S)))
  bar <- if (progress) utils::txtProgressBar(max = B, style = 3)
  for (b in seq_len(B)) {
    boot[b, ] <- tryCatch(suppressWarnings({
      Xb <- simulate_levels(mod, theta, d$t, m)
      th_b <- check_estimate(mod$fit(structure(list(t = d$t, X = Xb), class = "gof_data")), mod)
      if (anyNA(th_b)) NA_real_ else similarity_stats(Xb, simulate_levels(mod, th_b, d$t, n_pool))
    }), error = function(e) NA_real_)
    if (progress) utils::setTxtProgressBar(bar, b)
  }
  if (progress) close(bar)
  ok <- !is.na(boot[, 1])
  if (!any(ok)) stop("the estimation failed on all the bootstrap samples", call. = FALSE)
  boot <- boot[ok, , drop = FALSE]

  critical <- apply(boot, 2, stats::quantile, probs = 1 - alpha, names = FALSE)
  structure(list(
    model = mod$name, process = mod$process, estimate = theta,
    statistic = S, critical = critical,
    p.value = colMeans(sweep(boot, 2, S, ">=")),
    reject = S > critical,
    alpha = alpha, B = B, n_pool = n_pool, n_failed = sum(!ok),
    m = m, t = d$t, boot = boot
  ), class = "gof_test")
}

check_settings <- function(B, n_pool, alpha) {
  if (!is.numeric(B) || length(B) != 1 || B < 1) stop("B must be a positive integer")
  if (!is.numeric(n_pool) || length(n_pool) != 1 || n_pool < 1) stop("n_pool must be a positive integer")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
}

#' @export
print.gof_test <- function(x, digits = 4, ...) {
  cat("\nGoodness-of-fit test by similarity measures and parametric bootstrap\n")
  cat("(Zheng, Wang, Chen and Xie, 2026)\n\n")
  cat("H0: the data follow the ", x$model, " model (", x$process, " process)\n", sep = "")
  cat("data:", x$m, "trajectory(ies),", length(x$t), "inspection times after time 0\n")
  cat("estimates:", paste(names(x$estimate), signif(x$estimate, digits), sep = " = ", collapse = ", "), "\n")
  cat("bootstrap:", x$B, "samples", if (x$n_failed) paste0("(", x$n_failed, " failed and discarded)"),
      "| pools of", x$n_pool, "trajectories | alpha =", x$alpha, "\n\n")
  tab <- as.data.frame(x)
  tab$statistic <- signif(tab$statistic, digits)
  tab$critical <- signif(tab$critical, digits)
  print(tab, row.names = FALSE)
  cat("\nH0 rejected by", sum(x$reject), "of the 4 statistics.\n")
  invisible(x)
}

#' Results of a goodness-of-fit test as a data frame
#'
#' @param x Result of [gof_test()].
#' @param row.names,optional,... Not used, for compatibility with the generic
#'   function.
#'
#' @return A data frame with one row per statistic and the columns `stat`
#'   (name of the statistic), `statistic` (observed value), `critical`
#'   (critical value), `p.value` and `decision` (`"rejected"` or
#'   `"not rejected"`).
#'
#' @examples
#' set.seed(1)
#' x <- sim_model(0:20, "wiener", c(mu = 1, sigma = 0.5), n_traj = 5)
#' as.data.frame(gof_test(x, "wiener", B = 50, n_pool = 200))
#' @export
as.data.frame.gof_test <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    stat = names(x$statistic),
    statistic = unname(x$statistic),
    critical = unname(x$critical),
    p.value = unname(x$p.value),
    decision = ifelse(x$reject, "rejected", "not rejected"),
    stringsAsFactors = FALSE
  )
}

#' Plot the bootstrap distributions of a goodness-of-fit test
#'
#' @description
#' Draws, for each of the four statistics, the histogram of the bootstrap
#' statistics with the observed value (solid line) and the critical value
#' (dashed line). \eqn{H_0} is rejected by a statistic when the solid line is
#' to the right of the dashed line.
#'
#' @param x Result of [gof_test()].
#' @param ... Not used.
#'
#' @return A `ggplot` object, which can be modified with the usual
#'   `ggplot2` functions.
#'
#' @examples
#' set.seed(1)
#' x <- sim_model(0:20, "gamma", c(a = 2, b = 3), n_traj = 5)
#' res <- gof_test(x, "gamma", B = 100, n_pool = 250)
#' plot(res)
#'
#' # Change the title
#' plot(res) + ggplot2::labs(title = "Gamma model, bootstrap distributions")
#' @export
plot.gof_test <- function(x, ...) {
  stats_names <- colnames(x$boot)
  boot <- data.frame(stat = factor(rep(stats_names, each = nrow(x$boot)), levels = stats_names),
                     value = as.vector(x$boot))
  marks <- data.frame(stat = factor(stats_names, levels = stats_names),
                      observed = unname(x$statistic), critical = unname(x$critical))
  ggplot2::ggplot(boot, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 30, fill = "grey75", colour = "white") +
    ggplot2::geom_vline(data = marks, ggplot2::aes(xintercept = critical), linetype = "dashed") +
    ggplot2::geom_vline(data = marks, ggplot2::aes(xintercept = observed), colour = "#C2410C",
                        linewidth = 0.8) +
    ggplot2::facet_wrap(~ stat, scales = "free") +
    ggplot2::labs(title = paste("Bootstrap distributions,", x$model, "model"),
                  x = "Statistic", y = "Number of bootstrap samples") +
    ggplot2::theme_minimal()
}

#' Test several candidate models on the same data
#'
#' @description
#' Runs [gof_test()] for several candidate models on the same data and
#' gathers the estimates, the p-values and the decisions in one table. A
#' model that cannot be fitted to the data (for instance an increasing model
#' when the data contain a decrease) is reported as not testable, with the
#' reason, instead of stopping.
#'
#' @inheritParams gof_test
#' @param models Candidate models: names of models of [gof_models()] (by
#'   default all of them) or a list mixing names and models built with
#'   [custom_model()].
#' @param seed If not `NULL`, the random seed set before each test, so that
#'   each result does not depend on the other models tested.
#'
#' @return A data frame with one row per model and the columns `model`,
#'   `estimate` (estimates as text), the p-values `KS_sup`, `KS_sum`,
#'   `CvM_sup`, `CvM_sum`, `rejections` (number of statistics rejecting
#'   \eqn{H_0}, out of 4) and `note` (reason why a model is not testable, or
#'   number of discarded bootstrap samples). The list of the complete test
#'   results is in the attribute `tests`.
#'
#' @seealso [gof_test()]
#'
#' @examples
#' # Inverse Gaussian data tested against three increasing models
#' set.seed(1)
#' x <- sim_model(0:30, "ig", c(mu = 0.5, eta = 2), n_traj = 6)
#' res <- compare_models(x, c("gamma", "ig", "gamma_power"), B = 50, n_pool = 300, seed = 1)
#' res
#'
#' # The complete result of one test
#' attr(res, "tests")$ig
#'
#' # Data with negative increments: the increasing models are not testable
#' y <- sim_model(0:30, "wiener", c(mu = 0.5, sigma = 1), n_traj = 6)
#' compare_models(y, c("wiener", "gamma"), B = 50, n_pool = 300, seed = 1)
#' @export
compare_models <- function(x, models = gof_models()$model, t = NULL, B = 1000, n_pool = NULL,
                           alpha = 0.05, seed = NULL) {
  d <- gof_data(x, t)
  if (!is.list(models)) models <- as.list(models)
  labels <- vapply(models, function(m) get_model(m)$name, "")
  tests <- stats::setNames(vector("list", length(models)), labels)
  rows <- lapply(seq_along(models), function(k) {
    if (!is.null(seed)) set.seed(seed)
    res <- tryCatch(suppressWarnings(gof_test(d, models[[k]], B = B, n_pool = n_pool, alpha = alpha)),
                    error = function(e) e)
    if (inherits(res, "error")) {
      return(data.frame(model = labels[k], estimate = NA_character_, KS_sup = NA_real_, KS_sum = NA_real_,
                        CvM_sup = NA_real_, CvM_sum = NA_real_, rejections = NA_integer_,
                        note = sub("^the .* model cannot be fitted to the data: ", "not testable: ",
                                   conditionMessage(res)), stringsAsFactors = FALSE))
    }
    tests[[k]] <<- res
    data.frame(model = labels[k],
               estimate = paste(names(res$estimate), signif(res$estimate, 4), sep = " = ", collapse = ", "),
               t(res$p.value), rejections = sum(res$reject),
               note = if (res$n_failed) paste(res$n_failed, "bootstrap samples discarded") else "",
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "tests") <- tests
  out
}
