# Déroulé pas à pas de la méthode — exemple concret et version théorique générale

> **Statut de ce document** : reconstruction pédagogique produite pour accompagner la
> compréhension de l'article ci-dessous. Le jeu de données numérique, les valeurs estimées et
> les résultats du bootstrap sont **calculés par une exécution de code originale** (non tirés de
> l'article), à des fins d'illustration uniquement — ce ne sont ni des données ni des résultats
> publiés par les auteurs.

## Référence de l'article

Zheng J., Wang Y., Chen Y., Xie M. (2026). *Automated Verification of Degradation Models Using
Goodness-of-Fit Test Based on Similarity Measures and Parametric Bootstrap*. Quality and
Reliability Engineering International. https://doi.org/10.1002/qre.70382

---

## Partie 1 — Exemple numérique concret (jouet)

Jeu de données observé, $m=2$ unités, $n=3$ instants communs, modèle candidat Wiener homogène
($\Lambda(t)=t$) :

$$
X = \begin{pmatrix} x_{1,1} & x_{1,2} & x_{1,3} \\ x_{2,1} & x_{2,2} & x_{2,3} \end{pmatrix}
= \begin{pmatrix} 1.1 & 2.3 & 2.9 \\ 0.9 & 1.8 & 3.2 \end{pmatrix}, \qquad T = [1,2,3]
$$

### Étape 1 — Hypothèse
$H_0$ : "$M_T \in \mathcal M = \{M(\mu_W,\sigma)\}$" (Wiener homogène) contre $H_1$ : "non".

### Étape 2 — Estimation
Incréments (avec $x_{i,0}=0$) :

$$
\Delta X = \begin{pmatrix} 1.1 & 1.2 & 0.6 \\ 0.9 & 0.9 & 1.4 \end{pmatrix}, \qquad \Delta t = [1,1,1]
$$

MLE (incréments $\sim \mathcal N(\mu_W \Delta t,\ \sigma^2 \Delta t)$) :

$$\hat\mu_W = 1.0167, \qquad \hat\sigma = 0.2544$$

### Étape 3 — Génération du pool de comparaison
$m' = 100 \times m = 200$ trajectoires simulées sous $M(\hat\mu_W,\hat\sigma)$, aux instants $T$.

### Étape 4 — Statistique observée
Pour chaque trajectoire réelle, distance à la plus proche parmi les 200 simulées (KS = écart max,
CvM = écart quadratique cumulé), puis agrégation (max ou somme sur les 2 unités) :

| Statistique | Valeur observée $S(\hat{\boldsymbol\theta})$ |
|---|---|
| $S_{KS-sup}$  | 0.0899 |
| $S_{KS-sum}$  | 0.1722 |
| $S_{CvM-sup}$ | 0.0132 |
| $S_{CvM-sum}$ | 0.0251 |

### Étape 5 — Seuil critique (bootstrap paramétrique, $L=200$ réplications)
Pour chaque réplication $v$ : simuler $X^*$ ($m=2$ trajectoires) sous $M(\hat\mu_W,\hat\sigma)$,
ré-estimer $\hat{\boldsymbol\theta}^*_v$, regénérer un pool de $m'=200$ trajectoires sous
$M(\hat{\boldsymbol\theta}^*_v)$ et recalculer $S^*_v$. Quantile empirique à $95\%$ des
$L=200$ valeurs obtenues :

| Statistique | Seuil critique $s^*$ | $S(\hat{\boldsymbol\theta})$ | Rejet ? | $p$-valeur empirique |
|---|---|---|---|---|
| $S_{KS-sup}$  | 0.1822 | 0.0899 | Non | 0.480 |
| $S_{KS-sum}$  | 0.2694 | 0.1722 | Non | 0.330 |
| $S_{CvM-sup}$ | 0.0497 | 0.0132 | Non | 0.445 |
| $S_{CvM-sum}$ | 0.0598 | 0.0251 | Non | 0.360 |

### Étape 6 — Décision
Dans les 4 cas, $S(\hat{\boldsymbol\theta})$ est sous le seuil critique $\Rightarrow$ **on ne
rejette pas** $H_0$ au niveau $\alpha_S=5\%$ : rien ne permet de rejeter l'hypothèse Wiener
homogène sur ce jeu de données jouet.

**Remarque** : avec seulement $m=2$ et $n=3$, le seuil critique est très large (peu
d'information disponible), ce qui illustre concrètement pourquoi la puissance du test croît
avec $m$ et $n$ (constat empirique de l'article, §4.2–4.3).

---

## Partie 2 — Déroulé théorique général

Jeu de données observé, en toute généralité :

$$
\mathcal X = [x_{i,j}]_{m\times n}, \qquad X_i = [x_{i,1},\dots,x_{i,n}], \qquad T=[t_1,\dots,t_n]
$$

$m$ trajectoires, $n$ instants **communs** à toutes les unités (hypothèse de synchronisation de
l'article). Modèle candidat $M(\boldsymbol\theta)$, $\boldsymbol\theta\in\boldsymbol\Theta$,
dont on sait (a) estimer les paramètres et (b) simuler des trajectoires complètes.

### Étape 1 — Formulation de l'hypothèse
$$H_0 : \text{"}M_T \in \mathcal M = \{M(\boldsymbol\theta)\mid\boldsymbol\theta\in\boldsymbol\Theta\}\text{"} \qquad\text{vs}\qquad H_1 : \text{"}M_T \notin \mathcal M\text{"}$$

### Étape 2 — Estimation ponctuelle
$$\hat{\boldsymbol\theta} = \hat{\boldsymbol\theta}(\mathcal X)$$
Seul moment où $\mathcal X$ sert d'entrée à une procédure d'estimation (méthode externe au
schéma, MLE dans les expériences du papier).

### Étape 3 — Génération d'un pool de trajectoires simulées
$$X'_1,\dots,X'_{m'} \sim M(\hat{\boldsymbol\theta}), \qquad m' = 100\,m$$
- Processus stochastique (Wiener/Gamma/IG) : tirage séquentiel des incréments,
  $x'_{k,j} = x'_{k,j-1} + \Delta x'_{k,j}$ (Algorithme 2).
- Trajectoire générale (linéaire/loi de puissance/exponentiel) : tirage indépendant par point,
  $x'_{k,j} = g(t_j;\boldsymbol\omega,\boldsymbol\phi_k) + \varepsilon_{k,j}$ (Algorithme 3).

### Étape 4 — Mesure de similarité et construction de la statistique
$$D_{KS}(X_i,X'_k) = \sup_{1\le j\le n} |x_{i,j}-x'_{k,j}| \qquad D_{CvM}(X_i,X'_k) = \sum_{j=1}^n |x_{i,j}-x'_{k,j}|^2$$
$$\tilde X_i' = \arg\min_{X'_k,\ k=1,\dots,m'} D(X_i,X'_k)$$
$$S_{sup}(\hat{\boldsymbol\theta}) = \max_{1\le i\le m} D(X_i,\tilde X_i') \qquad S_{sum}(\hat{\boldsymbol\theta}) = \sum_{i=1}^m D(X_i,\tilde X_i')$$
(chacune déclinée en version $D_{KS}$ ou $D_{CvM}$, soit quatre statistiques).

### Étape 5 — Calibration du seuil par bootstrap paramétrique
Pour $v=1,\dots,L$ :
1. Générer $\mathcal X^*$ ($m\times n$) sous $M(\hat{\boldsymbol\theta})$ ;
2. Ré-estimer $\hat{\boldsymbol\theta}^*_v = \hat{\boldsymbol\theta}(\mathcal X^*)$ ;
3. Refaire intégralement l'étape 3 (nouveau pool de $m'$ sous $M(\hat{\boldsymbol\theta}^*_v)$)
   puis l'étape 4 sur $\mathcal X^*$ : obtenir $S^*_v$.

$$s^* = \text{quantile}_{1-\alpha_S}\big(S^*_1,\dots,S^*_L\big)$$

*Pourquoi ré-estimer à chaque réplication* : la variabilité de $S$ cumule le bruit
d'échantillonnage et l'incertitude d'estimation ; le test doit être calibré sur cette variabilité
totale, identique à celle subie par $S(\hat{\boldsymbol\theta})$.

### Étape 6 — Décision
$$\text{Rejeter } H_0 \text{ au niveau } \alpha_S \iff S(\hat{\boldsymbol\theta}) > s^*$$
$$p\text{-valeur empirique} = \frac{1}{L}\sum_{v=1}^L \mathbb 1\{S^*_v \ge S(\hat{\boldsymbol\theta})\}$$

### Ce qui rend la méthode unifiée
Seules les étapes 2 (estimation) et 3 (simulation) dépendent du type de modèle $M(\boldsymbol\theta)$.
Les étapes 4, 5, 6 — mesure de similarité, bootstrap, décision — sont **identiques** quel que
soit le modèle testé : toute la spécificité du modèle est encapsulée dans deux modules
interchangeables (estimer / simuler), le squelette du test restant fixe.
