# Fiche d'implémentation — Test d'adéquation par similarité de trajectoires et bootstrap paramétrique

## Référence de l'article
Zheng J., Wang Y., Chen Y., Xie M. (2026). *Automated Verification of Degradation Models Using
Goodness-of-Fit Test Based on Similarity Measures and Parametric Bootstrap*. Quality and
Reliability Engineering International. https://doi.org/10.1002/qre.70382

## Emplacement dans l'article
Section 2.2 (construction des statistiques, éq. 3–11) + Section 2.3 (Algorithme 1, calcul du seuil
critique par bootstrap paramétrique). Combine deux blocs de l'article en une procédure unique et
exécutable, puisque l'Algorithme 1 du papier renvoie lui-même à la construction des statistiques
de la Section 2.2.

## Objectif
Tester $H_0$: "$M_T \in \mathcal M = \{M(\boldsymbol\theta) \mid \boldsymbol\theta \in \boldsymbol\Theta\}$"
contre $H_1$: "$M_T \notin \mathcal M$", à partir de $m$ trajectoires de dégradation observées
$X_1,\dots,X_m$ mesurées aux mêmes $n$ instants $t_1,\dots,t_n$.

## Hypothèses
- Les $m$ trajectoires observées sont mutuellement indépendantes et issues du même modèle vrai $M_T$.
- **Toutes les unités sont mesurées simultanément aux mêmes instants** $t_1,\dots,t_n$ (hypothèse
  explicite de l'article, reconnue comme limitante en pratique — §4.4 de l'article).
- $\boldsymbol\theta$ est estimable à partir des données (méthode au choix ; MLE dans les
  simulations du papier).
- Le modèle candidat $M(\boldsymbol\theta)$ dispose d'un simulateur permettant de générer des
  trajectoires complètes sous ses paramètres (cf. fiches `_algorithm_02` pour les processus
  stochastiques et `_algorithm_03` pour les trajectoires générales).

## Entrées
- `X` : matrice $m \times n$ des trajectoires réelles observées, $x_{i,j}$ = valeur de dégradation
  de l'unité $i$ à l'instant $t_j$.
- `T` : vecteur des $n$ instants d'observation $t_1,\dots,t_n$ (communs à toutes les unités).
- Spécification du modèle candidat $M(\boldsymbol\theta)$ (famille paramétrique + fonction de
  tendance $\Lambda(t)$ le cas échéant).
- Fonction d'estimation de $\boldsymbol\theta$ à partir d'un jeu de données $m \times n$ (ex. MLE).
- Fonction de génération de trajectoires sous $M(\boldsymbol\theta)$ donné (Algorithme 2 ou 3).

## Paramètres
- `D` : mesure de similarité, `"KS"` ($D_{KS}$) ou `"CvM"` ($D_{CvM}$).
- `agg` : agrégation, `"sup"` ou `"sum"`.
- `m_prime` : nombre de trajectoires simulées pour la recherche du plus proche voisin.
  Valeur recommandée par l'article : `m_prime = 100 * m`.
- `L` : nombre de réplications bootstrap pour le calcul du quantile critique.
  Valeur utilisée dans le papier : `L = 1000`.
- `alpha_S` : niveau de signification. Valeur par défaut de l'article : `0.05`.

## Notations
- $D_{KS}(X_i, X_k') = \sup_{1\le j\le n} |x_{i,j} - x'_{k,j}|$
- $D_{CvM}(X_i, X_k') = \sum_{j=1}^n |x_{i,j} - x'_{k,j}|^2$
- $\tilde X_i' = \arg\min_{X_k',\,k=1,\dots,m'} D(X_i, X_k')$ : trajectoire simulée la plus proche de $X_i$
- $S_{KS-sup}(\hat{\boldsymbol\theta}) = \max_{1\le i\le m} D_{KS}(X_i,\tilde X_i')$,
  $S_{KS-sum}(\hat{\boldsymbol\theta}) = \sum_{i=1}^m D_{KS}(X_i,\tilde X_i')$ (idem pour CvM)
- $S^*_v(\hat{\boldsymbol\theta}^*)$ : statistique recalculée sur la $v$-ième réplication bootstrap

## Étapes détaillées

### Étape A — Statistique observée
1. Estimer $\hat{\boldsymbol\theta}$ à partir de `X` (fonction d'estimation fournie en entrée).
2. Générer $m'$ trajectoires $X'_1,\dots,X'_{m'}$ sous $M(\hat{\boldsymbol\theta})$, aux instants `T`.
3. Pour $i=1,\dots,m$ : calculer $D(X_i, X'_k)$ pour $k=1,\dots,m'$ ; retenir
   $\tilde X_i' = \arg\min_k D(X_i,X'_k)$ et $d_i = D(X_i,\tilde X_i')$.
4. Calculer $S(\hat{\boldsymbol\theta}) = \max_i d_i$ (variante `sup`) ou $\sum_i d_i$ (variante `sum`).

### Étape B — Seuil critique (bootstrap paramétrique, Algorithme 1 du papier)
Pour $v = 1,\dots,L$ :
1. Générer un jeu de données bootstrap $X^* = [x^*_{i,j}]_{m\times n}$ sous $M(\hat{\boldsymbol\theta})$
   (mêmes $m$, mêmes instants `T`).
2. Ré-estimer $\boldsymbol\theta$ à partir de $X^*$ pour obtenir $\hat{\boldsymbol\theta}^*_v$
   (même méthode d'estimation qu'à l'étape A.1).
3. Répéter les étapes A.2–A.4 en substituant `X` par $X^*$ et $\hat{\boldsymbol\theta}$ par
   $\hat{\boldsymbol\theta}^*_v$ (nouveau pool de $m'$ trajectoires simulé sous
   $M(\hat{\boldsymbol\theta}^*_v)$) pour obtenir $S^*_v$.

Trier $S^*_1,\dots,S^*_L$. Le seuil critique $s^*$ est le quantile empirique d'ordre $(1-\alpha_S)$
de cet échantillon.

### Étape C — Décision
- Si $S(\hat{\boldsymbol\theta}) > s^*$ : rejeter $H_0$ au niveau $\alpha_S$.
- Sinon : ne pas rejeter $H_0$.
- $p$-valeur empirique optionnelle : proportion de $\{S^*_v\} \ge S(\hat{\boldsymbol\theta})$ parmi
  les $L$ réplications (utilisée dans les tableaux d'application du papier, ex. Tables 9–13).

## Sortie
- Décision (rejet / non-rejet de $H_0$)
- Valeur de la statistique observée $S(\hat{\boldsymbol\theta})$
- Seuil critique $s^*$
- $p$-valeur empirique (optionnelle)

## Informations nécessaires à l'implémentation
- Une fonction d'estimation `estimer(X, T, modele) -> theta_hat` (dépend du modèle candidat).
- Une fonction de simulation `simuler(theta, T, m_traj, modele) -> matrice m_traj x n`
  (Algorithme 2 pour processus stochastique, Algorithme 3 pour trajectoire générale — voir fiches
  associées).
- Recommandation de l'article (§4.3, résultats empiriques) : privilégier $S_{KS-sup}$ et
  $S_{CvM-sup}$, qui dominent systématiquement les variantes `-sum` dans les scénarios testés.

## Points nécessitant une vérification
- **Ré-estimation à l'étape B.2** : la méthode de ré-estimation de $\hat{\boldsymbol\theta}^*$ n'est
  pas redéfinie explicitement dans le texte — il faut supposer qu'elle est identique à celle de
  l'étape A.1 (MLE dans les expériences du papier). *Interprétation nécessaire à l'implémentation,
  non une donnée explicite de l'article.*
- **Regénération du pool $m'$ à chaque réplication bootstrap** : la lecture la plus cohérente du
  texte de la Section 2.3 ("we simulate i.i.d. replications $X^*$ from $M(\hat{\boldsymbol\theta})$...
  re-estimate $\boldsymbol\theta$... compute $S^*(\hat{\boldsymbol\theta}^*)$") implique qu'un
  nouveau pool complet de $m'$ trajectoires est simulé à **chaque** réplication $v=1,\dots,L$, sous
  le modèle ré-estimé $M(\hat{\boldsymbol\theta}^*_v)$ — et non sous $M(\hat{\boldsymbol\theta})$
  fixe. Ceci est cohérent avec le coût computationnel annoncé (~80 s pour un test unique). À
  confirmer si possible par correspondance avec les auteurs ou par un test de cohérence numérique
  (reproduction approximative des Tables 3–13).
- Aucune règle formelle de sélection optimale de `m_prime` ou `L` n'est fournie — seulement des
  heuristiques empiriques ($m'=100m$) et un compromis qualitatif coût/puissance discuté en §4.3
  et §4.4.

## Ambiguïtés éventuelles
- Pas de garantie théorique sur le niveau ou la consistance du test — uniquement une validation
  empirique par Monte Carlo (Table 3, comparaison niveau empirique vs $\alpha_S$). Ceci n'est pas
  une ambiguïté d'implémentation mais une limite de fond de la méthode, à garder à l'esprit lors
  de toute comparaison avec une méthode alternative disposant de garanties théoriques.

## Remarques importantes pour la reproduction
- Le test est **agnostique au type de modèle** : il suffit de disposer d'un estimateur et d'un
  simulateur pour n'importe quelle famille paramétrique — d'où le lien direct avec les fiches
  `_algorithm_02` (processus stochastiques : Wiener/Gamma/IG) et `_algorithm_03` (trajectoires
  générales : linéaire/loi de puissance/exponentiel).
- Le coût dominant est la double boucle : $L$ réplications bootstrap $\times$ $m'$ trajectoires
  simulées par réplication (pour construire le pool de plus-proches-voisins) $\times$ ré-estimation
  de $\boldsymbol\theta$ à chaque réplication. Pour une implémentation R/Julia, vectoriser au
  maximum la génération de trajectoires et le calcul des distances (matrice $m \times m'$ de
  distances plutôt qu'une double boucle explicite) est essentiel pour rester praticable.
