#' gofdegrad: goodness-of-fit tests for degradation models
#'
#' @description
#' `gofdegrad` implements the goodness-of-fit test for degradation models of
#' Zheng, Wang, Chen and Xie (2026). Is a gamma process, a Wiener process or
#' another model adequate for a set of degradation trajectories? The test
#' compares the observed trajectories with trajectories simulated under the
#' fitted model, through similarity measures, and calibrates the decision by
#' parametric bootstrap. It applies to trajectories inspected at the same,
#' possibly irregular, times.
#'
#' @section Functions of the package:
#' Test (click on a function name for its help page):
#'
#' | Function | Purpose |
#' |:---|:---|
#' | [gof_test()] | goodness-of-fit test of one candidate model |
#' | [compare_models()] | tests of several candidate models on the same data |
#' | [power_study()] | empirical level or power of the test |
#'
#' Models and data:
#'
#' | Function | Purpose |
#' |:---|:---|
#' | [gof_models()] | the seven available candidate models |
#' | [fit_model()] | maximum likelihood estimation of a candidate model |
#' | [sim_model()] | simulation of a candidate model |
#' | [sim_path()] | simulation of general path models |
#' | [custom_model()] | any other candidate model, from its own functions |
#' | [gof_data()] | conversion and check of the data |
#' | [similarity_stats()] | the four similarity statistics |
#'
#' The result of [gof_test()] is printed as a table, converted to a data
#' frame by [as.data.frame()][as.data.frame.gof_test()] and plotted by
#' [plot()][plot.gof_test()].
#'
#' @section Candidate models:
#'
#' | Model | Process | Parameters |
#' |:---|:---|:---|
#' | `"wiener"` | Wiener with drift | `mu`, `sigma` |
#' | `"gamma"` | gamma | shape `a`, rate `b` |
#' | `"ig"` | inverse Gaussian | `mu`, `eta` |
#' | `"vg"` | variance-gamma | `theta`, `sigma`, `nu` |
#' | `"meixner"` | Meixner | `a`, `b`, `d`, `m` |
#' | `"wiener_power"` | Wiener, trend \eqn{\mu t^r} | `mu`, `sigma`, `r` |
#' | `"gamma_power"` | gamma, shape function \eqn{a t^r} | `a`, `b`, `r` |
#'
#' The first five are Lévy processes, estimated with the package `levyfit`
#' and simulated with the package `levysim`; [gof_models()] gives the law
#' of their increments. Any other model, for instance a general path model,
#' can be tested with [custom_model()].
#'
#' @section Getting started:
#' A test always follows the same three steps.
#' 1. **Prepare the data**: all the trajectories must be inspected at the
#'    same times. Simulated data from `levysim` or [sim_model()] can be used
#'    directly; real data can be given as a data frame with the columns
#'    `id` (trajectory), `t` (time) and `X` (degradation level), or as a
#'    matrix with one row per trajectory and the times in `t` (see
#'    [gof_data()]). The degradation is measured from time 0.
#' 2. **Test** each candidate model with [gof_test()], or several at once
#'    with [compare_models()], after `set.seed()` for reproducible results.
#'    Use a small number of bootstrap samples `B` to explore, then
#'    `B = 1000`.
#' 3. **Read the result**: the p-values and decisions of the four statistics,
#'    with [plot()][plot.gof_test()] for the bootstrap distributions. A
#'    non-rejection with few trajectories or few inspections only means that
#'    the data are compatible with the model: check the power with
#'    [power_study()].
#'
#' @section The test in brief:
#' Let \eqn{X_1, \dots, X_m} be the observed trajectories, inspected at the
#' times \eqn{t_1 < \dots < t_n}. The parameters of the candidate model are
#' estimated by maximum likelihood, and \eqn{m'} trajectories are simulated
#' under the fitted model (the pool, \eqn{m' = 100m} by default). Each
#' observed trajectory is matched with the closest trajectory of the pool,
#' for the distance \eqn{\max_j |x_j - y_j|} (KS type) or
#' \eqn{\sum_j (x_j - y_j)^2} (CvM type), and these smallest distances are
#' summarized by their maximum (`sup`) or their sum (`sum`): four statistics.
#' Their distribution under the null hypothesis is obtained by parametric
#' bootstrap, the parameters being estimated again on each bootstrap sample.
#' See [gof_test()] for the details.
#'
#' @references
#' Zheng, J., Wang, Y., Chen, Y. and Xie, M. (2026). Automated verification
#' of degradation models using goodness-of-fit test based on similarity
#' measures and parametric bootstrap. *Quality and Reliability Engineering
#' International*. \doi{10.1002/qre.70382}
#'
#' Ye, Z.-S. and Xie, M. (2015). Stochastic modelling and analysis of
#' degradation for highly reliable products. *Applied Stochastic Models in
#' Business and Industry*, 31(1), 16-32. \doi{10.1002/asmb.2063}
#'
#' Meeker, W. Q., Escobar, L. A. and Pascual, F. G. (2021). *Statistical
#' Methods for Reliability Data*, 2nd edition. Wiley.
#'
#' @seealso
#' Useful links:
#' * <https://github.com/rytracetito/Article_zheng_chen_et_al>
#' * Report bugs at <https://github.com/rytracetito/Article_zheng_chen_et_al/issues>
#' * The packages `levysim` (simulation) and `levyfit` (estimation):
#'   <https://github.com/rytracetito/levy_packages>
#'
#' @examples
#' # 1. Data: six trajectories of a gamma process, inspected every time unit
#' set.seed(1)
#' x <- sim_model(0:30, "gamma", c(a = 2, b = 3), n_traj = 6)
#' levysim::plot_traj(x)
#'
#' # 2. Test of the gamma model (small B to keep the example fast)
#' res <- gof_test(x, "gamma", B = 100, n_pool = 300)
#' res
#' plot(res)
#'
#' # 3. Several candidate models at once
#' compare_models(x, c("gamma", "ig", "wiener"), B = 50, n_pool = 300, seed = 1)
"_PACKAGE"

utils::globalVariables(c("value", "critical", "observed"))
