#' Niveau empirique ou puissance d'un test d'adequation
#'
#' Repete \code{\link{test_gof}} sur \code{h} jeux de donnees independants
#' generes par \code{gen_X}, et retourne le taux de rejet moyen pour chacune
#' des 4 statistiques. C'est le niveau empirique si \code{gen_X} simule sous
#' \eqn{H_0}, la puissance sinon. Un tirage pour lequel \code{test_gof}
#' echoue (par exemple des increments negatifs incompatibles avec un modele
#' Gamma ou IG) est compte comme un rejet sur les 4 statistiques : des
#' donnees qui ne peuvent meme pas etre ajustees par le modele candidat ne
#' lui sont, par nature, pas conformes.
#'
#' @param h Nombre entier de repetitions.
#' @param gen_X Fonction \code{function() -> matrice}, generatrice du jeu de
#'   donnees a chaque repetition.
#' @param T Vecteur des instants d'observation.
#' @param est Fonction d'estimation, voir \code{\link{test_gof}}.
#' @param sim Fonction de simulation, voir \code{\link{test_gof}}.
#' @param mp Taille du pool simule (\eqn{m'}).
#' @param L Nombre de replications bootstrap.
#' @param alpha Niveau de signification du test. Par defaut 0.05.
#'
#' @return Vecteur numerique nomme de longueur 4, taux de rejet arrondi a
#'   4 decimales pour chaque statistique.
#'
#' @examples
#' T <- 1:50
#' gen_h0 <- function() sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
#' etude_puissance(h = 10, gen_h0, T, est_modele("Wiener"), sim_modele("Wiener"),
#'                 mp = 300, L = 50)
#'
#' @export
etude_puissance <- function(h, gen_X, T, est, sim, mp, L, alpha = 0.05) {
  rejets <- vector("list", h)
  for (r in seq_len(h)) {
    rejets[[r]] <- tryCatch(
      test_gof(gen_X(), T, est, sim, mp = mp, L = L, alpha = alpha)$rejet,
      error = function(e) c(KS_sup = TRUE, KS_sum = TRUE, CvM_sup = TRUE, CvM_sum = TRUE)
    )
  }
  round(colMeans(do.call(rbind, rejets)), 4)
}
