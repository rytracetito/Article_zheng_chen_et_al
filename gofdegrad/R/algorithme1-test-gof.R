#' Matrice des distances entre trajectoires
#'
#' Calcule la matrice \eqn{m \times m'} des distances de similarite entre
#' chaque trajectoire reelle et chaque trajectoire simulee (eq. 3-4 de
#' Zheng, Wang, Chen & Xie, 2026).
#'
#' @param X Matrice \eqn{m \times n} de trajectoires reelles (ou bootstrap).
#' @param Xs Matrice \eqn{m' \times n} de trajectoires simulees.
#' @param D Mesure de similarite : \code{"KS"} (distance de Chebyshev, par
#'   analogie avec la statistique de Kolmogorov-Smirnov) ou \code{"CvM"}
#'   (distance euclidienne au carre).
#'
#' @return Matrice numerique \eqn{m \times m'}, terme \eqn{(i,k)} egal a
#'   \eqn{D(X_i, X'_k)}.
#'
#' @keywords internal
dist_mat <- function(X, Xs, D) {
  if (D == "CvM") {
    outer(rowSums(X^2), rowSums(Xs^2), "+") - 2 * X %*% t(Xs)
  } else {
    n <- ncol(X)
    Dm <- matrix(0, nrow(X), nrow(Xs))
    for (j in seq_len(n)) Dm <- pmax(Dm, abs(outer(X[, j], Xs[, j], "-")))
    Dm
  }
}

#' Statistiques de similarite
#'
#' Calcule les 4 statistiques de test \eqn{S_{KS-sup}}, \eqn{S_{KS-sum}},
#' \eqn{S_{CvM-sup}}, \eqn{S_{CvM-sum}} entre des trajectoires et un pool
#' simule (eq. 6-10 de l'article).
#'
#' @param X Matrice de trajectoires reelles ou bootstrap.
#' @param Xs Matrice du pool de trajectoires simulees sous le modele candidat.
#'
#' @return Vecteur numerique nomme de longueur 4 :
#'   \code{c(KS_sup, KS_sum, CvM_sup, CvM_sum)}.
#'
#' @export
stats_sim <- function(X, Xs) {
  d_KS  <- apply(dist_mat(X, Xs, "KS"),  1, min)
  d_CvM <- apply(dist_mat(X, Xs, "CvM"), 1, min)
  c(KS_sup  = max(d_KS),  KS_sum  = sum(d_KS),
    CvM_sup = max(d_CvM), CvM_sum = sum(d_CvM))
}

#' Test d'adequation par similarite et bootstrap parametrique
#'
#' Implemente l'Algorithme 1 complet de Zheng, Wang, Chen & Xie (2026) :
#' calcule la statistique observee, le seuil critique par bootstrap
#' parametrique, et rend la decision de test pour les 4 statistiques
#' simultanement.
#'
#' @param X Matrice \eqn{m \times n} des donnees reelles.
#' @param T Vecteur des \code{n} instants d'observation.
#' @param est Fonction d'estimation \code{function(X, T) -> theta}, propre
#'   au modele candidat sous \eqn{H_0} (voir \code{\link{estimer_parametres}}).
#' @param sim Fonction de simulation \code{function(theta, T, m_traj) ->
#'   matrice}, propre au modele candidat (voir \code{\link{sim_proc}},
#'   \code{\link{sim_traj}}).
#' @param mp Taille du pool simule (\eqn{m'}). Par defaut \eqn{100m},
#'   recommandation de l'article.
#' @param L Nombre de replications bootstrap. Par defaut 1000.
#' @param alpha Niveau de signification du test. Par defaut 0.05.
#'
#' @return Une liste avec les elements :
#'   \item{rejet}{vecteur logique de longueur 4, \code{TRUE} = \eqn{H_0} rejetee}
#'   \item{S_obs}{statistiques observees}
#'   \item{s_crit}{valeurs critiques (quantile empirique d'ordre \eqn{1-\alpha})}
#'   \item{p_valeur}{p-valeurs empiriques}
#'
#' @examples
#' T <- 1:50
#' X <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
#' sim_w <- function(theta, T, m_traj) sim_proc(T, theta, "Wiener", m_traj)
#' res <- test_gof(X, T, estimer_wiener_homogene, sim_w, mp = 300, L = 50)
#' res$rejet
#'
#' @export
test_gof <- function(X, T, est, sim, mp = 100 * nrow(X), L = 1000, alpha = 0.05) {
  m <- nrow(X)
  th <- est(X, T)
  S_obs <- stats_sim(X, sim(th, T, mp))

  S_star <- matrix(NA_real_, L, length(S_obs), dimnames = list(NULL, names(S_obs)))
  for (v in seq_len(L)) {
    Xb <- sim(th, T, m)
    th_b <- suppressWarnings(est(Xb, T))
    S_star[v, ] <- stats_sim(Xb, sim(th_b, T, mp))
  }
  s_crit <- apply(S_star, 2, stats::quantile, probs = 1 - alpha, names = FALSE)
  list(rejet = S_obs > s_crit, S_obs = S_obs, s_crit = s_crit,
       p_valeur = colMeans(sweep(S_star, 2, S_obs, ">=")))
}
