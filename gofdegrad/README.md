# gofdegrad

Test d'adéquation par similarité de trajectoires et bootstrap paramétrique
pour les modèles de dégradation, d'après :

Zheng J., Wang Y., Chen Y., Xie M. (2026). *Automated Verification of
Degradation Models Using Goodness-of-Fit Test Based on Similarity Measures
and Parametric Bootstrap*. Quality and Reliability Engineering International.
https://doi.org/10.1002/qre.70382

## Installation

Depuis GitHub (depot prive : un jeton personnel GitHub est necessaire,
voir `usethis::create_github_token()` puis `gitcreds::gitcreds_set()`) :

```r
# install.packages("remotes")
remotes::install_github("rytracetito/Article_zheng_chen_et_al", subdir = "gofdegrad")
```

Depuis le dossier source local :

```r
devtools::install("chemin/vers/gofdegrad")
```

## Modeles disponibles

| `modele` | Description |
|---|---|
| `"Wiener"` | Wiener homogene |
| `"Gamma"` | Gamma homogene |
| `"IG"` | Inverse Gaussien homogene |
| `"Wiener_nonhomog"` | Wiener non-homogene, Lambda(t) = t^r |
| `"Gamma_nonhomog"` | Gamma non-homogene, Lambda(t) = t^r |

## Utilisation minimale

```r
library(gofdegrad)

T <- 1:50
X <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)

resultat <- test_gof(X, T, est_modele("Wiener"), sim_modele("Wiener"),
                      mp = 100 * nrow(X), L = 1000)

resultat$rejet
resultat$p_valeur
```

## Etude de niveau/puissance

```r
gen_H0 <- function() sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
etude_puissance(h = 1000, gen_H0, T, est_modele("Wiener"), sim_modele("Wiener"),
                mp = 600, L = 200)
```

## Structure

- `R/algorithme1-test-gof.R` : Algorithme 1 (distances, statistiques, test complet)
- `R/algorithme2-simulateurs-processus.R` : Algorithme 2 (Wiener, Gamma, IG)
- `R/algorithme3-simulateurs-trajectoire.R` : Algorithme 3 (lineaire, loi de puissance, exponentiel)
- `R/estimateurs.R` : estimateurs MLE, un par modele, plus le dispatcher `estimer_parametres()`
- `R/dispatchers.R` : `sim_modele()` et `est_modele()`, fermetures generiques pour `test_gof()`
- `R/etude-puissance.R` : etude de niveau et de puissance par repetition

## Limites connues

- Les estimateurs Gamma/IG (homogenes et non-homogenes) refusent les increments
  negatifs (erreur explicite). Un increment exactement nul, du a l'arrondi machine
  quand la forme `alpha * dt` est tres petite, est remplace par
  `.Machine$double.eps * max(abs(X))` avec un avertissement ; cet avertissement
  n'est pas affiche pour les repliques bootstrap de `test_gof()`.
- Le pas de temps peut etre irregulier, mais toutes les trajectoires doivent
  partager les memes instants `T`.
- `estimer_wiener_nonhomog()` et `estimer_gamma_nonhomog()` recherchent `r`
  sur un intervalle borne (`[0.05, 10]` pour Wiener, recherche libre pour
  Gamma) ; aucune garantie mathematique d'optimalite globale n'est etablie
  pour cette recherche sur `r` (contrairement aux autres parametres, qui ont
  une forme close).
- La forme de la loi des effets aleatoires `Phi` dans les modeles a
  trajectoire generale (`gen_phi_norm`) est une hypothese de normalite non
  explicitement donnee par l'article source.
