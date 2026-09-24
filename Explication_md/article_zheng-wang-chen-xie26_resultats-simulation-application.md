# Synthèse des résultats — simulation et application à des données réelles

## Référence de l'article
Zheng J., Wang Y., Chen Y., Xie M. (2026). *Automated Verification of Degradation Models Using
Goodness-of-Fit Test Based on Similarity Measures and Parametric Bootstrap*. Quality and
Reliability Engineering International. https://doi.org/10.1002/qre.70382

---

# Section 4 — Étude de simulation

## 4.1 Protocole général

**Ce qui a été fait** : avant d'étudier la puissance du test sur des cas particuliers, les auteurs
vérifient deux choses systématiquement pour chaque scénario : (i) que le **niveau empirique**
(taux de rejet quand $H_0$ est vraie) est proche du niveau théorique $\alpha_S$, et (ii) que la
**puissance** (taux de rejet quand $H_0$ est fausse) est la plus élevée possible.

$$1-\beta_p = \frac{\text{nombre de rejets de } H_0}{h}\times 100\%$$

**Paramètres communs à toutes les expériences** : $\alpha_S = 5\%$, estimation par MLE,
$L=1000$ réplications bootstrap, $m'=100m$.

## 4.2 Résultats de simulation — cinq expériences distinctes

### Expérience A — $H_0$ : Wiener homogène (éq. 15), Table 3

**Ce qui a été fait** : ajuster un Wiener homogène et simuler des données sous trois processus
différents pour voir si le test (i) accepte bien quand c'est vrai, (ii) rejette quand c'est faux.

**Paramètres étudiés** : $h=1000$, $m=6$, $n=50$.
- Wiener (pour le niveau) : $\mu_W\in\{0.2,0.4\}$, $\sigma\in\{0.05,0.1\}$
- Gamma bruité (alternative) : $\alpha=1$, $\beta\in\{1/3,3\}$, $\sigma_\varepsilon\in\{0.05,0.1\}$
- IG bruité (alternative) : $\mu_{IG}\in\{0.5,1\}$, $\eta\in\{0.5,1\}$, $\sigma_\varepsilon\in\{0.05,0.1\}$

**Résultats** :
- Sous Wiener (H0 vraie) : niveaux empiriques ≈ 3,7 % à 5,5 % — proches de 5 %, le test est bien
  calibré.
- Sous Gamma bruité : **puissance faible**, entre 3,8 % et 12,8 % ; $S_{KS-sup}$ domine
  (jusqu'à 12,8 %).
- Sous IG bruité : puissance faible à modérée, entre 4,2 % et 32,5 % ; $S_{KS-sup}$ domine encore
  (jusqu'à 32,5 %).
- $S_{KS-sup}$ et $S_{CvM-sup}$ surpassent systématiquement $S_{KS-sum}$ et $S_{CvM-sum}$.
- Quelques puissances **sous** le niveau nominal → signe de biais potentiel dans certains cas.

**Figure 3** illustre pourquoi la puissance est si faible : en ajustant un Wiener sur des données
Gamma bruitées ($\alpha=1,\beta=1/3,\sigma_\varepsilon=0.05$ → $\hat\mu_W=0.338,\hat\sigma=0.355$)
ou IG bruitées ($\mu_{IG}=0.5,\eta=0.5,\sigma_\varepsilon=0.05$ → $\hat\mu_W=0.492,\hat\sigma=0.498$),
les trajectoires simulées et réelles se ressemblent visuellement — le test peine donc à les
distinguer.

### Expérience B — même $H_0$, mais $n=100$ (Table 4)

**Ce qui a été fait** : refaire les cas Gamma/IG de l'expérience A avec deux fois plus de points
d'observation, pour voir l'effet de $n$ sur la puissance.

**Résultats** : la puissance de $S_{KS-sup}$ et $S_{CvM-sup}$ augmente nettement — Gamma :
21,6–22,4 % ; IG : 36,8–45,8 % pour $S_{KS-sup}$. En revanche $S_{CvM-sup}$ reste faible pour IG
(7,2–7,9 %) : l'amélioration n'est pas uniforme entre statistiques.

### Expérience C — une seule unité de test, $m=1$ (Table 5)

**Ce qui a été fait** : même protocole que l'expérience A mais avec $m=1$ (donc seulement
$S_{KS}$ et $S_{CvM}$, qui coïncident avec les versions sup/sum).

**Résultats** : puissance systématiquement plus basse qu'avec 6 unités — Wiener (H0) : 3,1–5,4 %
(niveau correct) ; Gamma : 4,2–5,3 % ; IG : 3,3–5,4 %. Confirme que **le nombre d'unités testées
influence fortement la puissance**.

### Expérience D — $H_0$ : Wiener non-homogène à tendance loi de puissance (éq. 16), Table 6

**Ce qui a été fait** : ajuster $SX(t)=\mu_W t^r+\sigma B(t)$ ($r$ estimé) et simuler sous des
versions **non-homogènes** de Wiener, Gamma, IG.

**Paramètres étudiés** : $h=1000$, $m=6$, $n=100$ ; deux formes de tendance testées,
$\Lambda(t)=\sqrt t$ et $\Lambda(t)=t^2$, avec des jeux de paramètres différents pour chaque
processus (ex. Gamma : $\alpha\in\{2,0.5,0.002,0.0005\}$, $\beta\in\{5,20\}$, $\sigma_\varepsilon=1$).

**Résultats** :
- Niveaux empiriques sous Wiener (H0) proches de 5 % → validité confirmée en non-homogène
  également.
- $S_{KS-sup}$ reste la statistique la plus puissante (jusqu'à 62,6 % pour l'IG avec
  $\Lambda(t)=\sqrt t$), suivi de $S_{CvM-sup}$.
- La puissance **augmente avec la variance entre trajectoires générées**.

### Expérience E — distinguer processus stochastique vs trajectoire générale (éq. 17), Tables 7 et 8

**Ce qui a été fait** : ajuster un Wiener généralisé (tendance $\Lambda(t)$ libre) sur des données
simulées sous trois modèles à trajectoire générale (linéaire, loi de puissance, exponentiel),
pour voir si le test arrive à rejeter le Wiener quand le vrai modèle est d'une nature
complètement différente.

**Paramètres étudiés** : $\sigma_\varepsilon\in\{0.05,0.1\}$ (bruit de mesure), $m\in\{1,6\}$
(Tables 7 et 8), $n=50$.

**Résultats** (avec $m=6$, Table 8) :

| Modèle générateur | $\sigma_\varepsilon=0.05$ | $\sigma_\varepsilon=0.1$ |
|---|---|---|
| Linéaire | 92–94 % | 34–50 % |
| Loi de puissance | 94–99 % | 31–43 % |
| Exponentiel | 99 % | 51–72 % |

**Constat central** : puissance proche de 100 % à faible bruit, mais **chute significative** dès
que $\sigma_\varepsilon$ double — le bruit de mesure dilue le signal discriminant. Avec $m=1$
(Table 7), les puissances sont systématiquement plus basses qu'avec $m=6$.

## 4.3 Synthèse des auteurs

- Recommandation : privilégier $S_{KS-sup}$ et $S_{CvM-sup}$ dans le cadre automatisé.
- Pour $H_0$=Wiener contre Gamma/IG : puissance globalement faible, $S_{KS-sup}$ le plus efficace
  — attribué à la ressemblance graphique subtile entre ces processus une fois ajustés.
- Pour distinguer processus stochastique vs trajectoire générale : test efficace, mais puissance
  qui diminue avec $\sigma_\varepsilon$.
- Facteurs qui améliorent la puissance : plus d'observations $n$, plus d'unités $m$ ; un $m'$
  plus grand améliore aussi la puissance mais avec rendements décroissants et coût croissant.

## 4.4 Discussion (positionnement + limites pratiques)

**Ce qui a été fait** : replacer la méthode par rapport à la littérature et discuter ses limites
opérationnelles.

- Pour les modèles homogènes à observation périodique, les tests classiques (KS, χ²) restent
  utilisables ; la méthode proposée les complète pour les cas non-périodiques/non-homogènes.
- Pour les processus non-stationnaires (Gamma non-stationnaire, Wiener à dérive/diffusion
  variable, transformations d'échelle de temps), les tests classiques sont inadéquats — la
  méthode proposée s'applique.
- Comparaison à Li et al. : leur méthode (fondée sur la normalité multivariée) ne s'applique
  qu'au Wiener, pas au Gamma/IG à erreurs non-normales.
- **Limite reconnue** : hypothèse de synchronisation temporelle de toutes les unités — extension
  aux plans asynchrones proposée comme travail futur.
- **Coût computationnel** : ~80 s par test unique, ~20 h pour l'étude de simulation complète
  (1000 jeux) — dominé par $m'$.
- **Compromis sur $m'$** : plus grand → estimation plus stable, mais (a) coût de calcul
  explosif et (b) la statistique tend vers 0 quand $m'\to\infty$, ce qui **réduit** le pouvoir
  discriminant du test. Un $m'$ mal calibré peut donc nuire à la puissance autant que l'améliorer.

---

# Section 5 — Application à trois jeux de données réels

## 5.1 Compteurs électriques intelligents (SEM) — erreur chronométrique (CE)

**Ce qui a été fait** : Xu et al. avaient ajusté un Wiener homogène sur ces données et obtenu
une durée de vie fiable théorique nettement inférieure à celle observée — signal d'une possible
mauvaise spécification du modèle. Les auteurs testent formellement cette hypothèse.

**Paramètres étudiés** : 64 échantillons SEM (points aberrants exclus),
$\hat\mu_W=1.9601\times10^{-4}$, $\hat\sigma=7.1187\times10^{-3}$.

**Résultat** : les 4 statistiques **rejettent** le Wiener homogène (Table 9, $p$-valeurs de
0,004 à 0,032). Fait intéressant : le test de Ye et al. (2013), qui ne compare que la trajectoire
moyenne, ne rejetterait *pas* ce modèle (les moyennes coïncident) — mais il ignore les
fluctuations dues à $\sigma$, que la méthode proposée capture et qui suffisent à rejeter le
modèle. Un Wiener non-homogène à tendance loi de puissance a aussi été testé
($\hat\mu_W=1.8015\times10^{-5}$, $\hat\sigma=7.5161\times10^{-3}$, $\hat r=1.4087$) : également
rejeté (Table 10, $p$-valeurs de 0,001 à 0,028) — aucun des deux Wiener n'est adéquat pour ces
données.

## 5.2 Résistance en courant alternatif (ACR) d'une batterie

**Ce qui a été fait** : ajuster un Wiener **non-homogène** à tendance loi de puissance sur des
données de batterie à trois températures, avec 4 trajectoires par température.

**Paramètres étudiés/estimés** (Table 11) :

| Température | $\Lambda(t;\hat r)$ | $\hat\mu_W$ | $\hat\sigma$ |
|---|---|---|---|
| 25 °C | $t^{0.6092}$ | 0,0031 | 0,0010 |
| 35 °C | $t^{0.4543}$ | 0,0115 | 0,0020 |
| 45 °C | $t^{0.7256}$ | 0,0041 | 0,0036 |

**Résultat** : aux trois températures, $S(\hat{\boldsymbol\theta})$ reste **sous** le seuil
critique pour les 4 statistiques ($p$-valeurs entre 0,089 et 0,523, Table 12) → **non-rejet** de
$H_0$. Le Wiener non-homogène à tendance loi de puissance est jugé adéquat pour modéliser la
dégradation ACR à ces trois températures — cas où le modèle candidat est validé.

## 5.3 Courant de fonctionnement d'un laser GaAs

**Ce qui a été fait** : des travaux antérieurs (Omshi & Shemehsavar) avaient utilisé un Wiener
homogène pour ces données **sans effectuer de test d'adéquation**. Les auteurs y remédient.

**Paramètres étudiés** : 15 lasers testés à 80 °C, durée totale 4000 h, intervalles de 250 h. Le
Wiener homogène semble adapté visuellement (Fig. 7a-b). $\hat\mu_W=2.0372\times10^{-3}$,
$\hat\sigma=1.2657\times10^{-2}$.

**Résultat** : les 4 statistiques **rejettent catégoriquement** le Wiener homogène (Table 13,
$p$-valeurs de 0,001 à 0,022 — nettement plus tranché que pour le cas ACR). Point notable :
**malgré une bonne apparence visuelle**, le test formel rejette le modèle — cela illustre
concrètement l'intérêt d'un test rigoureux plutôt que d'une inspection graphique, argument déjà
mis en avant pour le cas SEM.

---

# Section 6 — Conclusion (synthèse finale de l'article)

- Cadre développé : deux mesures de similarité (type KS, type CvM) + quatre statistiques +
  calibration par bootstrap paramétrique.
- Validité confirmée par la correspondance niveau empirique/théorique dans l'étude de simulation.
- Le test rejette efficacement le Wiener quand les données proviennent de modèles à trajectoire
  générale.
- La puissance croît avec le nombre d'observations et d'unités testées.
- $S_{KS-sup}$ et $S_{CvM-sup}$ recommandés pour un usage automatisé.
- Travaux futurs annoncés : modèles non-paramétriques et processus corrélés ; extension aux
  plans d'observation asynchrones.
