# Fiche d'implémentation — Génération de trajectoires simulées (modèles à processus stochastique)

## Référence de l'article
Zheng J., Wang Y., Chen Y., Xie M. (2026). *Automated Verification of Degradation Models Using
Goodness-of-Fit Test Based on Similarity Measures and Parametric Bootstrap*. Quality and
Reliability Engineering International. https://doi.org/10.1002/qre.70382

## Emplacement dans l'article
Section 3, **Algorithme 2** (p.6) : "Simulated degradation path generation for stochastic process
degradation models". Bloc de base réutilisé par le test GoF (fiche `_algorithm_01`) pour construire
le pool de trajectoires simulées, aussi bien pour la statistique observée que pour le bootstrap.

## Objectif
Générer $m'$ trajectoires simulées sous un modèle candidat de type processus stochastique
(Wiener, Gamma ou IG, homogène ou non-homogène via $\Lambda(t)$), aux mêmes instants d'observation
que les données réelles.

## Hypothèses
- Incréments indépendants sur intervalles non chevauchants, valeur initiale nulle p.s.
  (définition d'un modèle de dégradation stochastique, Section 2.1 de l'article).
- La loi de l'incrément $\Delta SX(t,t+\Delta t)$ suit exactement la spécification du Tableau 1
  de l'article (Wiener = Normale, Gamma = Gamma, IG = Inverse Gaussienne), paramétrée par
  $\hat{\boldsymbol\theta}$ et par la fonction de tendance $\Lambda(t)$.
- $\Lambda(t)$ est connue/estimée en amont (ex. $\Lambda(t)=t$ pour le cas homogène,
  $\Lambda(t)=t^{\hat r}$ pour le cas non-homogène en loi de puissance testé dans l'article).

## Entrées
- `T` : vecteur des instants d'observation $[t_1,\dots,t_n]$ (avec $t_0=0$ par convention).
- Modèle estimé $M(\hat{\boldsymbol\theta})$ : type (Wiener / Gamma / IG) + valeurs des paramètres.
- Fonction de tendance $\Lambda(t;\hat r)$ (identité si modèle homogène).

## Paramètres
- `m_prime` : nombre de trajectoires à simuler. Valeur recommandée par l'article : `100 * m`
  (où `m` est le nombre d'unités réellement observées).

## Notations
- $x'_{k,j}$ : valeur simulée de la $k$-ième trajectoire à l'instant $t_j$.
- $\Delta x'_{k,j}$ : incrément simulé sur $[t_{j-1},t_j]$.
- $\Delta\Lambda(t_{j-1},t_j) = \Lambda(t_j)-\Lambda(t_{j-1})$.
- $p(\Delta x; t_{j-1},t_j;\hat{\boldsymbol\theta})$ : densité de l'incrément sur cet intervalle,
  selon le modèle.

## Étapes détaillées

1. Poser $t_0 = 0$.
2. Pour $k = 1,\dots,m'$ :
   1. Poser $x'_{k,0} = 0$.
   2. Pour $j = 1,\dots,n$ :
      1. Calculer $\Delta\Lambda(t_{j-1},t_j) = \Lambda(t_j) - \Lambda(t_{j-1})$.
      2. Tirer l'incrément $\Delta x'_{k,j}$ selon le modèle :
         - **Wiener** : $\Delta x'_{k,j} \sim \mathcal N\big(\mu_W\,\Delta\Lambda(t_{j-1},t_j),\ \sigma^2 (t_j-t_{j-1})\big)$
           — attention : la variance dépend de $\Delta t = t_j - t_{j-1}$, **pas** de $\Delta\Lambda$
           (cf. Tableau 1 de l'article).
         - **Gamma** : $\Delta x'_{k,j} \sim \mathcal G\big(\text{shape}=\alpha\,\Delta\Lambda(t_{j-1},t_j),\ \text{scale}=\beta\big)$
         - **IG** : $\Delta x'_{k,j} \sim \mathcal{IG}\big(\text{mean}=\mu_{IG}\,\Delta\Lambda(t_{j-1},t_j),\ \text{shape}=\eta\,\Delta\Lambda^2(t_{j-1},t_j)\big)$
      3. Poser $x'_{k,j} = x'_{k,j-1} + \Delta x'_{k,j}$.
3. Retourner $\{X'_k = [x'_{k,1},\dots,x'_{k,n}]\}_{k=1}^{m'}$.

## Pseudo-code
```
fonction simuler_processus_stochastique(T, theta_hat, modele, m_prime):
    n <- longueur(T)
    t0 <- 0
    X_prime <- matrice(m_prime, n)
    pour k de 1 à m_prime:
        x_prec <- 0
        pour j de 1 à n:
            dLambda <- Lambda(T[j], theta_hat) - Lambda(T[j-1] ou t0, theta_hat)
            dt <- T[j] - (T[j-1] ou t0)
            si modele == "Wiener":
                dx <- rnorm(1, mean = mu_W * dLambda, sd = sigma * sqrt(dt))
            sinon si modele == "Gamma":
                dx <- rgamma(1, shape = alpha * dLambda, scale = beta)
            sinon si modele == "IG":
                dx <- rinvgauss(1, mean = mu_IG * dLambda, shape = eta * dLambda^2)
            x_courant <- x_prec + dx
            X_prime[k, j] <- x_courant
            x_prec <- x_courant
    retourner X_prime
```

## Sortie
Matrice $m' \times n$ de trajectoires simulées $X'_k = [x'_{k,1},\dots,x'_{k,n}]$, $k=1,\dots,m'$.

## Informations nécessaires à l'implémentation
- Générateurs aléatoires disponibles nativement en R (`rnorm`, `rgamma`) et Julia
  (`Distributions.jl` : `Normal`, `Gamma`). Pour la loi Inverse Gaussienne, package `statmod`
  (R : `rinvgauss`) ou `Distributions.jl` (Julia : `InverseGaussian`).

## Points nécessitant une vérification
- **Paramétrisation Gamma** : le texte affirme explicitement "$\beta>0$ refers to the scale
  parameter of the Gamma process" (Tableau 1), ce qui suggère la convention *shape-scale*
  ($\mathbb E[X]=\text{shape}\times\text{scale}$, convention native de `rgamma(shape, scale=...)`
  en R). À **vérifier par cohérence dimensionnelle** avant implémentation définitive (ex. en
  recalculant une moyenne de trajectoire simulée et en la comparant à une valeur attendue), car
  une confusion scale/rate ($\text{rate}=1/\text{scale}$) est une source d'erreur fréquente.
- **Paramétrisation IG** : la notation $\mathcal{IG}(\mu_{IG}\Delta\Lambda,\ \eta\Delta\Lambda^2)$
  correspond à la paramétrisation usuelle (moyenne, paramètre de forme) de la loi de Wald,
  $f(x;\mu,\lambda) = \sqrt{\lambda/(2\pi x^3)}\exp\big(-\lambda(x-\mu)^2/(2\mu^2 x)\big)$, avec ici
  $\mu = \mu_{IG}\Delta\Lambda$, $\lambda=\eta\Delta\Lambda^2$. Correspond à
  `statmod::rinvgauss(n, mean = mu, shape = lambda)` en R — à vérifier également par cohérence
  numérique avant utilisation en production.

## Ambiguïtés éventuelles
- L'article ne détaille pas explicitement la forme de $\Lambda(t)$ pour des cas non-homogènes
  autres que ceux testés dans l'étude de simulation ($\Lambda(t)=t$, $t^r$, $\sqrt t$, $t^2$) —
  $\Lambda$ doit être fournie ou estimée en dehors de cet algorithme, qui se contente de l'évaluer.

## Remarques importantes pour la reproduction
- Algorithme purement génératif, sans étape d'optimisation ni de bootstrap interne. Réutilisé tel
  quel par la fiche `_algorithm_01` (test GoF complet), aussi bien pour construire le pool initial
  de $m'$ trajectoires que pour chaque réplication bootstrap.
- Pour une implémentation efficace en R/Julia : vectoriser la génération sur les $m'$ trajectoires
  simultanément à chaque pas de temps $j$ (tirage vectorisé de $m'$ incréments) plutôt que de
  boucler unité par unité, afin de limiter le coût pour $m' = 100m$.
