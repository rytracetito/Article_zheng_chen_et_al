# Fiche d'implémentation — Génération de trajectoires simulées (modèles à trajectoire générale)

## Référence de l'article
Zheng J., Wang Y., Chen Y., Xie M. (2026). *Automated Verification of Degradation Models Using
Goodness-of-Fit Test Based on Similarity Measures and Parametric Bootstrap*. Quality and
Reliability Engineering International. https://doi.org/10.1002/qre.70382

## Emplacement dans l'article
Section 3, **Algorithme 3** (p.6) : "Simulated degradation path generation for general path
degradation models". Bloc de base réutilisé par le test GoF (fiche `_algorithm_01`) lorsque le
modèle candidat appartient à la classe des modèles à trajectoire générale (Tableau 2 :
linéaire, loi de puissance, exponentiel).

## Objectif
Générer $m'$ trajectoires simulées sous un modèle à trajectoire générale ajusté
$M(\hat{\boldsymbol\theta})$, avec $\hat{\boldsymbol\theta}=(\hat{\boldsymbol\omega},\hat{\boldsymbol\varphi},\hat\sigma_\varepsilon)$,
aux mêmes instants d'observation que les données réelles.

## Hypothèses
- Modèle mixte à effets aléatoires : $x_{i,j} = g(t_j;\boldsymbol\omega,\boldsymbol\phi_i) + \varepsilon_{i,j}$,
  $\varepsilon_{i,j}\sim\mathcal N(0,\sigma_\varepsilon^2)$ i.i.d., $\boldsymbol\phi_i \sim \boldsymbol\Phi(\boldsymbol\varphi)$
  i.i.d. entre unités simulées (Section 2.1 de l'article).
- $\boldsymbol\omega$ (effets fixes) est commun à toutes les unités simulées ; $\boldsymbol\phi_k$
  (effets aléatoires) varie d'une unité simulée à l'autre.
- La forme fonctionnelle $g$ est fixée (linéaire, loi de puissance ou exponentielle, Tableau 2).

## Entrées
- `T` : vecteur des instants d'observation $[t_1,\dots,t_n]$ (avec $t_0=0$).
- $\hat{\boldsymbol\theta} = (\hat{\boldsymbol\omega}, \hat{\boldsymbol\varphi}, \hat\sigma_\varepsilon)$ :
  effets fixes estimés, hyperparamètres estimés de la loi des effets aléatoires, écart-type
  résiduel estimé.
- Forme fonctionnelle `g` (linéaire / loi de puissance / exponentiel, Tableau 2).

## Paramètres
- `m_prime` : nombre de trajectoires à simuler. Valeur recommandée par l'article : `100 * m`.

## Notations
- $\boldsymbol\phi_k$ : vecteur d'effets aléatoires tiré pour la $k$-ième unité simulée.
- $\varepsilon_{k,j}$ : erreur de mesure simulée pour l'unité $k$ à l'instant $t_j$.
- $g(t;\boldsymbol\omega,\boldsymbol\phi_k)$ : valeur moyenne de dégradation selon le Tableau 2 :
  - Linéaire : $g(t;\phi_L) = \phi_L\, t$
  - Loi de puissance : $g(t;\phi_{P1},\phi_{P2}) = \phi_{P1}\, t^{\phi_{P2}}$
  - Exponentiel : $g(t;\omega_E,\phi_E) = \phi_E\big[\exp(\omega_E t)-1\big]$

## Étapes détaillées

1. Fixer $\boldsymbol\omega = \hat{\boldsymbol\omega}$.
2. Pour $k = 1,\dots,m'$ :
   1. Tirer $\boldsymbol\phi_k \sim \boldsymbol\Phi(\hat{\boldsymbol\varphi})$.
   2. Pour $j = 1,\dots,n$ :
      1. Tirer $\varepsilon_{k,j} \sim \mathcal N(0,\hat\sigma_\varepsilon^2)$.
      2. Poser $x'_{k,j} = g(t_j;\boldsymbol\omega,\boldsymbol\phi_k) + \varepsilon_{k,j}$.
3. Retourner $\{X'_k = [x'_{k,1},\dots,x'_{k,n}]\}_{k=1}^{m'}$.

## Pseudo-code
```
fonction simuler_trajectoire_generale(T, theta_hat, modele, m_prime):
    n <- longueur(T)
    omega <- theta_hat.omega
    X_prime <- matrice(m_prime, n)
    pour k de 1 à m_prime:
        phi_k <- tirer(Phi_hat)               # loi des effets aléatoires, hyperparamètres estimés
        pour j de 1 à n:
            eps <- rnorm(1, mean = 0, sd = sigma_eps_hat)
            X_prime[k, j] <- g(T[j], omega, phi_k, modele) + eps
    retourner X_prime

fonction g(t, omega, phi, modele):
    si modele == "lineaire":       retourner phi.L * t
    si modele == "loi_puissance":  retourner phi.P1 * t^phi.P2
    si modele == "exponentiel":    retourner phi.E * (exp(omega.E * t) - 1)
```

## Sortie
Matrice $m' \times n$ de trajectoires simulées $X'_k$, $k=1,\dots,m'$.

## Informations nécessaires à l'implémentation
- Générateur pour $\varepsilon$ : `rnorm` (R) / `Normal` (Julia, `Distributions.jl`).
- Générateur pour $\boldsymbol\phi_k$ : dépend de la forme retenue pour $\boldsymbol\Phi$
  (voir ambiguïté ci-dessous) — typiquement une loi normale (multivariée si plusieurs
  composantes, ex. loi de puissance à deux paramètres $\phi_{P1},\phi_{P2}$).

## Points nécessitant une vérification
- Pour le modèle exponentiel, l'attribution effet fixe/effet aléatoire n'est **pas donnée
  explicitement** dans la formule générale de la Section 2.1, mais elle est **déduite du contexte
  des tableaux de simulation** (Tables 7–8) : $\phi_E$ y est traité comme aléatoire
  (ex. $\phi_E \sim \mathcal N(0.2, 0.05^2)$) et $\omega_E$ comme fixe (ex. $\omega_E = 0.05$,
  valeur unique). Cette lecture est cohérente avec la notation $\boldsymbol\omega$ (vecteur des
  effets fixes *communs à toutes les unités*) et $\boldsymbol\phi_i$ (effets *individuels*) de la
  Section 2.1. **Interprétation déduite du contexte, non énoncé explicite** — à confirmer si une
  version plus détaillée de l'article (annexe) est disponible.

## Ambiguïtés éventuelles
- **La forme exacte de $\boldsymbol\Phi$ (distribution des effets aléatoires) n'est pas précisée**
  dans le corps de l'article fourni : le texte se limite à "$\boldsymbol\phi_i$ suit une
  distribution multivariée $\boldsymbol\Phi$ avec hyperparamètres $\boldsymbol\varphi$" (Section 2.1),
  sans donner la famille (normale multivariée, log-normale, etc.). Dans les tableaux de simulation
  (Tables 7–8), les effets aléatoires apparaissent notés $\mathcal N(\mu,\sigma^2)$ (ex.
  $\phi_L \sim \mathcal N(0.2, 0.05^2)$), ce qui suggère fortement une hypothèse de normalité
  (éventuellement multivariée pour le modèle en loi de puissance à deux paramètres
  $(\phi_{P1},\phi_{P2})$, avec corrélation non précisée). **Ne pas supposer une corrélation nulle
  entre $\phi_{P1}$ et $\phi_{P2}$ sans justification** — si l'implémentation nécessite une
  matrice de covariance, la poser comme diagonale est une **extension/interprétation
  personnelle**, pas une donnée de l'article, et doit être documentée comme telle.

## Remarques importantes pour la reproduction
- Algorithme purement génératif, sans étape d'optimisation ni de bootstrap interne. Réutilisé tel
  quel par la fiche `_algorithm_01` (test GoF complet) pour les modèles à trajectoire générale.
- Contrairement à `_algorithm_02` (processus stochastiques), il n'y a pas de dépendance
  séquentielle entre les instants $t_1,\dots,t_n$ : chaque $x'_{k,j}$ ne dépend que de $t_j$, de
  $\boldsymbol\omega$, de $\boldsymbol\phi_k$ et de $\varepsilon_{k,j}$ — la génération peut donc
  être entièrement vectorisée sur $(k,j)$ simultanément, ce qui la rend nettement moins coûteuse
  que l'Algorithme 2.
