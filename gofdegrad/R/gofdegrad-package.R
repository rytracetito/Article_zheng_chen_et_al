#' gofdegrad: Goodness-of-Fit Tests for Degradation Processes
#'
#' Implementation of the goodness-of-fit testing method by similarity
#' measures and parametric bootstrap of Zheng, Wang, Chen & Xie (2026,
#' Quality and Reliability Engineering International, \doi{10.1002/qre.70382}).
#'
#' @section Available models:
#' \describe{
#'   \item{\code{"Wiener"}}{Homogeneous Wiener}
#'   \item{\code{"Gamma"}}{Homogeneous Gamma}
#'   \item{\code{"IG"}}{Homogeneous Inverse Gaussian}
#'   \item{\code{"Wiener_nonhomog"}}{Non-homogeneous Wiener, \eqn{\Lambda(t) = t^r}}
#'   \item{\code{"Gamma_nonhomog"}}{Non-homogeneous Gamma, \eqn{\Lambda(t) = t^r}}
#' }
#'
#' @section Typical entry points:
#' \code{\link{sim_proc}} / \code{\link{sim_traj}} to simulate,
#' \code{\link{estimer_parametres}} to estimate, \code{\link{test_gof}}
#' to test, \code{\link{etude_puissance}} for the level/power study.
#'
#' @keywords internal
"_PACKAGE"
