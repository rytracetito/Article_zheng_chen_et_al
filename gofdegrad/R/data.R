#' Degradation data in the form used by the package
#'
#' @description
#' Converts degradation data given in any of the forms accepted by the
#' package into a single internal form: the inspection times after time 0
#' and the matrix of degradation levels, one row per trajectory. This
#' function is called by [gof_test()], [fit_model()] and
#' [similarity_stats()]; it is exported so that the data can be checked
#' before a test.
#'
#' @param x Degradation data, in one of the following forms:
#'   * the result of a simulation function of the package `levysim` (for
#'     instance `levysim::sim_gamma()`), of [sim_model()] or of [sim_path()];
#'   * a numeric matrix whose first column contains the inspection times and
#'     the other columns the degradation levels of the trajectories (the
#'     format of `levysim` and `levyfit`), when `t` is not given;
#'   * a numeric matrix with one row per trajectory and one column per
#'     inspection time, when `t` is given;
#'   * a numeric vector of degradation levels of one trajectory, with `t`;
#'   * a data frame with the columns `id` (trajectory), `t` (inspection time)
#'     and `X` (degradation level), in any row order, for instance read from
#'     a csv file or returned by `levysim::traj_df()`; the trajectories are
#'     sorted by `id`;
#'   * a list of data frames, one per trajectory, whose first column contains
#'     the inspection times and second column the degradation levels.
#' @param t Inspection times, only when `x` is a matrix with one row per
#'   trajectory or a vector: one time per column of `x` (or per element).
#'
#' @details
#' **Common inspection times.** The test of Zheng et al. (2026) compares
#' trajectories inspected at the same times, so all the trajectories must
#' share the same inspection times. These times can be irregular.
#'
#' **Origin of the degradation.** The models start from \eqn{X(0) = 0}.
#' When the data contain time 0, the level at time 0 of each trajectory is
#' subtracted from its levels, so that the degradation is measured from the
#' initial state (for instance from the initial length of a crack). When the
#' first inspection is after time 0, the degradation is assumed to be 0 at
#' time 0: the first increment goes from 0 at time 0 to the first observed
#' level, and it weighs on the estimation of the models.
#'
#' @return An object of class `gof_data`: a list with
#'   * `t`: the \eqn{n} inspection times after time 0, strictly increasing;
#'   * `X`: the \eqn{m \times n} matrix of degradation levels, measured from
#'     the level at time 0, one row per trajectory.
#'
#' @seealso [gof_test()], [fit_model()]
#'
#' @examples
#' # 1. Simulated data from levysim: a matrix with the column t
#' set.seed(1)
#' x <- levysim::sim_gamma(seq(0, 10, by = 1), a = 2, b = 3, n_traj = 4)
#' d <- gof_data(x)
#' d
#' d$t
#' dim(d$X)
#'
#' # 2. One row per trajectory, with the times given separately
#' X <- matrix(c(0.5, 1.1, 1.8,
#'               0.4, 0.9, 1.5), nrow = 2, byrow = TRUE)
#' gof_data(X, t = c(1, 2, 4))
#'
#' # 3. A data frame (id, t, X), as read from a csv file: the levels at
#' #    time 0 (here 0.9) are subtracted
#' df <- data.frame(id = rep(1:2, each = 3), t = rep(c(0, 1, 2), 2),
#'                  X = c(0.9, 1.0, 1.2, 0.9, 1.1, 1.4))
#' gof_data(df)$X
#'
#' # 4. A list of data frames (t, X), one per trajectory
#' liste <- list(data.frame(t = c(0, 5, 10), X = c(0, 1.2, 2.9)),
#'               data.frame(t = c(0, 5, 10), X = c(0, 0.8, 2.1)))
#' gof_data(liste)
#'
#' # 5. Real data in long format: GaAs laser data (package SMRD.data)
#' if (requireNamespace("SMRD.data", quietly = TRUE)) {
#'   laser <- SMRD.data::gaaslaser
#'   d <- gof_data(data.frame(id = laser$unit, t = laser$hours, X = laser$increase))
#'   dim(d$X)
#' }
#' @export
gof_data <- function(x, t = NULL) {
  if (inherits(x, "gof_data")) return(x)
  paths <- as_paths(x, t)
  if (length(paths) == 0) stop("no trajectory in the data")
  times <- lapply(paths, function(p) p$t)
  t0 <- times[[1]]
  if (!all(vapply(times, function(tt) length(tt) == length(t0) && isTRUE(all.equal(tt, t0)), TRUE))) {
    stop("all the trajectories must be inspected at the same times")
  }
  if (!is.numeric(t0) || anyNA(t0) || any(t0 < 0)) stop("inspection times must be non-negative numbers")
  if (any(diff(t0) <= 0)) stop("inspection times must be strictly increasing")
  X <- matrix(unlist(lapply(paths, function(p) p$X)), nrow = length(paths), byrow = TRUE)
  if (!is.numeric(X) || anyNA(X) || any(!is.finite(X))) stop("degradation levels must be finite numbers")
  if (t0[1] == 0) {
    X <- X[, -1, drop = FALSE] - X[, 1]
    t0 <- t0[-1]
  }
  if (length(t0) == 0) stop("at least one inspection after time 0 is needed")
  structure(list(t = t0, X = X), class = "gof_data")
}

as_paths <- function(x, t) {
  if (is.data.frame(x)) {
    if (!all(c("id", "t", "X") %in% names(x))) {
      stop("a data frame must have the columns id, t and X; use a list of data frames for one data frame per trajectory")
    }
    ids <- sort(unique(x$id))
    return(lapply(ids, function(i) {
      s <- x[x$id == i, ]
      s <- s[order(s$t), ]
      list(t = s$t, X = s$X)
    }))
  }
  if (is.matrix(x)) {
    x <- unclass(x)
    if (is.null(t)) {
      if (ncol(x) < 2) stop("the matrix must contain the column of times and at least one trajectory")
      return(lapply(seq_len(ncol(x) - 1) + 1, function(k) list(t = x[, 1], X = x[, k])))
    }
    if (length(t) != ncol(x)) stop("t must give one inspection time per column of x")
    return(lapply(seq_len(nrow(x)), function(i) list(t = t, X = x[i, ])))
  }
  if (is.list(x)) {
    return(lapply(x, function(d) {
      d <- as.data.frame(d)
      if (ncol(d) < 2) stop("each trajectory must contain the times and the levels")
      o <- order(d[[1]])
      list(t = d[[1]][o], X = d[[2]][o])
    }))
  }
  if (is.numeric(x)) {
    if (is.null(t) || length(t) != length(x)) stop("for one trajectory, give t with one time per level")
    return(list(list(t = t, X = x)))
  }
  stop("unsupported data format, see ?gof_data")
}

#' @export
print.gof_data <- function(x, ...) {
  cat("Degradation data:", nrow(x$X), "trajectory(ies),", length(x$t), "inspection times after time 0\n")
  cat("times:", format(utils::head(x$t, 8), digits = 4), if (length(x$t) > 8) "...", "\n")
  cat("levels at the last inspection:", format(utils::head(x$X[, ncol(x$X)], 8), digits = 4),
      if (nrow(x$X) > 8) "...", "\n")
  invisible(x)
}

increments_of <- function(d) {
  n <- ncol(d$X)
  dX <- d$X - cbind(0, d$X[, -n, drop = FALSE])
  dt <- matrix(diff(c(0, d$t)), nrow(d$X), n, byrow = TRUE)
  list(dX = dX, dt = dt)
}

positive_increments <- function(d, process) {
  inc <- increments_of(d)
  if (any(inc$dX < 0)) {
    stop("negative increments in the data: a ", process, " process only increases", call. = FALSE)
  }
  zero <- inc$dX == 0
  if (any(zero)) {
    floor_value <- max(abs(d$X)) * .Machine$double.eps
    if (floor_value <= 0) stop("no strictly positive increment in the data", call. = FALSE)
    inc$dX[zero] <- floor_value
    warning(sprintf("%d zero increment(s) replaced by %.3g (machine precision)", sum(zero), floor_value),
            call. = FALSE)
  }
  inc
}
