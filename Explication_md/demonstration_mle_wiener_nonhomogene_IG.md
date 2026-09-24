# Démonstration rigoureuse — de la vraisemblance à l'estimateur

Je reprends les deux cas en partant strictement de la densité des observations, jusqu'aux estimateurs, avec justification de l'optimalité **globale** (pas seulement locale) à chaque étape où c'est possible — et je signale explicitement l'endroit où cette garantie n'est *pas* établie.

---

## Partie 1 — Wiener non-homogène, $\Lambda(t)=t^r$

### 1.1 Hypothèses et variable observée

$m$ trajectoires indépendantes, observées aux mêmes instants $t_1<\dots<t_n$ ($t_0=0$). Par définition du modèle (§2.1 de l'article) : incréments indépendants sur intervalles non chevauchants, et
$$\Delta x_{i,j} := x_{i,j}-x_{i,j-1} \;\sim\; \mathcal N\big(\mu_W\,\Delta\Lambda_j,\ \sigma^2\Delta t_j\big),\qquad \Delta\Lambda_j = t_j^r-t_{j-1}^r,\quad \Delta t_j = t_j-t_{j-1}$$

Les $mn$ variables $\{\Delta x_{i,j}\}_{i,j}$ sont **mutuellement indépendantes** (indépendance en $j$ par définition du processus ; indépendance en $i$ car les trajectoires le sont).

### 1.2 Construction de la vraisemblance

Chaque $\Delta x_{i,j}$ a pour densité
$$f(\Delta x_{i,j};\mu_W,\sigma^2,r) = \frac{1}{\sqrt{2\pi\sigma^2\Delta t_j}}\exp\left(-\frac{(\Delta x_{i,j}-\mu_W\Delta\Lambda_j)^2}{2\sigma^2\Delta t_j}\right)$$

Par indépendance mutuelle, la vraisemblance est le **produit** de ces $mn$ densités :
$$L(\mu_W,\sigma^2,r) = \prod_{i=1}^m\prod_{j=1}^n \frac{1}{\sqrt{2\pi\sigma^2\Delta t_j}}\exp\left(-\frac{(\Delta x_{i,j}-\mu_W\Delta\Lambda_j)^2}{2\sigma^2\Delta t_j}\right)$$

### 1.3 Log-vraisemblance

$$\ell = \log L = \sum_{i,j}\left[-\frac12\log(2\pi)-\frac12\log\sigma^2-\frac12\log\Delta t_j-\frac{(\Delta x_{i,j}-\mu_W\Delta\Lambda_j)^2}{2\sigma^2\Delta t_j}\right]$$

En séparant les sommes (chaque terme ne dépendant que de $j$ est répété $m$ fois sur $i$) :

$$\boxed{\ell(\mu_W,\sigma^2,r) = -\frac{mn}2\log(2\pi) - \frac{mn}2\log\sigma^2 - \frac m2\sum_j\log\Delta t_j - \frac1{2\sigma^2}\sum_{i,j}\frac{(\Delta x_{i,j}-\mu_W\Delta\Lambda_j)^2}{\Delta t_j}}$$

### 1.4 Score en $\mu_W$, et optimalité globale

$$\frac{\partial\ell}{\partial\mu_W} = \frac1{\sigma^2}\sum_{i,j}\frac{\Delta\Lambda_j}{\Delta t_j}\big(\Delta x_{i,j}-\mu_W\Delta\Lambda_j\big)$$

Posons $w_j=\Delta\Lambda_j/\Delta t_j$ et $S_j=\sum_i\Delta x_{i,j}$. L'équation $\partial\ell/\partial\mu_W=0$ donne :
$$\sum_j w_j S_j = \mu_W\, m\sum_j w_j\Delta\Lambda_j = \mu_W\, m\sum_j\frac{\Delta\Lambda_j^2}{\Delta t_j} \quad\Longrightarrow\quad \hat\mu_W(r) = \frac{\sum_j \Delta\Lambda_j S_j/\Delta t_j}{m\sum_j\Delta\Lambda_j^2/\Delta t_j}$$

**Optimalité globale (pas seulement locale)** : à $(\sigma^2,r)$ fixés, $\ell$ est une fonction **exactement quadratique** de $\mu_W$, de coefficient dominant
$$\frac{\partial^2\ell}{\partial\mu_W^2} = -\frac m{\sigma^2}\sum_j\frac{\Delta\Lambda_j^2}{\Delta t_j} < 0 \quad\text{(constante, ne dépend pas de }\mu_W\text{)}$$
Une parabole à coefficient dominant strictement négatif est **globalement** concave : $\hat\mu_W(r)$ est donc le maximum global en $\mu_W$, quel que soit $\sigma^2>0$ — pas seulement un point critique local.

### 1.5 Score en $\sigma^2$, et optimalité globale

$$\frac{\partial\ell}{\partial\sigma^2} = -\frac{mn}{2\sigma^2}+\frac1{2\sigma^4}\sum_{i,j}\frac{(\Delta x_{i,j}-\mu_W\Delta\Lambda_j)^2}{\Delta t_j}$$

Notons $Q=\sum_{i,j}(\Delta x_{i,j}-\mu_W\Delta\Lambda_j)^2/\Delta t_j$ (à $\mu_W=\hat\mu_W(r)$ fixé). En posant $u=\sigma^2$ :
$$\ell'(u) = -\frac{mn}{2u}+\frac{Q}{2u^2} = \frac{1}{2u^2}\big(Q-mn\,u\big)$$

Ceci s'annule **une seule fois**, en $u^*=Q/(mn)$, et change de signe de **positif** (pour $u<u^*$) à **négatif** (pour $u>u^*$). Une fonction dont la dérivée change de signe une seule fois, de $+$ à $-$, admet en ce point un **maximum global** sur $(0,\infty)$ — argument de monotonicité, plus fort qu'un simple test de dérivée seconde locale (qui, ici, ne serait même pas concluant globalement : $\ell''(u)=mn/(2u^2)-Q/u^3$ change de signe selon $u$). D'où :
$$\hat\sigma^2(r) = \frac1{mn}\sum_{i,j}\frac{(\Delta x_{i,j}-\hat\mu_W(r)\Delta\Lambda_j)^2}{\Delta t_j} \quad\text{— maximum global confirmé.}$$

### 1.6 Profilage sur $r$

En substituant $\hat\mu_W(r)$ et $\hat\sigma^2(r)$ — et en utilisant que, par construction, $\sum_{i,j}(\Delta x_{i,j}-\hat\mu_W(r)\Delta\Lambda_j)^2/\Delta t_j = mn\,\hat\sigma^2(r)$ exactement (c'est la définition même de $\hat\sigma^2(r)$) — le dernier terme de $\ell$ devient une constante $-mn/2$ :

$$\ell_p(r) = \underbrace{-\frac{mn}2\log(2\pi)-\frac m2\sum_j\log\Delta t_j-\frac{mn}2}_{\text{ne dépend pas de }r} \;-\; \frac{mn}2\log\hat\sigma^2(r)$$

$x\mapsto-\log x$ étant strictement décroissante :
$$\boxed{\hat r = \arg\max_r \ell_p(r) = \arg\min_r \hat\sigma^2(r)}$$

**Limite honnête de cette étape** : contrairement à $\hat\mu_W(r)$ et $\hat\sigma^2(r)$ (optimalité globale démontrée ci-dessus), je n'ai **pas** de preuve que $r\mapsto\hat\sigma^2(r)$ soit unimodale sur $(0,10)$. `optimize()` (recherche par section dorée) suppose l'unimodalité et peut converger vers un minimum local si le profil en a plusieurs. Ce que j'ai fait pour compenser cette absence de preuve, c'est une vérification **empirique** (cross-check contre `optim()` sur les 3 paramètres simultanément, convergence identique) — ce n'est pas une garantie mathématique, seulement une évidence numérique sur les cas testés.

---

## Partie 2 — Inverse Gaussienne homogène

### 2.1 Hypothèses

$\Lambda(t)=t$ (homogène), donc $\Delta\Lambda_j=\Delta t_j$. Par le Tableau 1 :
$$\Delta x_{i,j} \sim \mathcal{IG}(\mu_j,\lambda_j),\qquad \mu_j=\mu_{IG}\Delta t_j,\quad \lambda_j=\eta\,\Delta t_j^2$$

### 2.2 Densité et vraisemblance

Densité IG standard (paramétrisation moyenne/forme) :
$$f(x;\mu,\lambda) = \sqrt{\frac{\lambda}{2\pi x^3}}\,\exp\left(-\frac{\lambda(x-\mu)^2}{2\mu^2x}\right),\quad x>0$$

$$L(\mu_{IG},\eta) = \prod_{i,j}\sqrt{\frac{\lambda_j}{2\pi(\Delta x_{i,j})^3}}\exp\left(-\frac{\lambda_j(\Delta x_{i,j}-\mu_j)^2}{2\mu_j^2\Delta x_{i,j}}\right)$$

### 2.3 Log-vraisemblance et simplification décisive

$$\ell = \sum_{i,j}\left[\frac12\log\lambda_j - \frac12\log\big(2\pi(\Delta x_{i,j})^3\big) - \frac{\lambda_j(\Delta x_{i,j}-\mu_j)^2}{2\mu_j^2\Delta x_{i,j}}\right]$$

**Point clé de la dérivation** — calculons le rapport $\lambda_j/\mu_j^2$ :
$$\frac{\lambda_j}{\mu_j^2} = \frac{\eta\,\Delta t_j^2}{\mu_{IG}^2\,\Delta t_j^2} = \frac{\eta}{\mu_{IG}^2} \quad\text{— ne dépend pas de }j\text{ !}$$

C'est cette annulation exacte de $\Delta t_j$ qui rend le problème traitable analytiquement. En substituant $\lambda_j=\eta\Delta t_j^2$, $\log\lambda_j=\log\eta+2\log\Delta t_j$ :

$$\ell(\mu_{IG},\eta) = \frac{mn}2\log\eta + \underbrace{m\sum_j\log\Delta t_j - \frac12\sum_{i,j}\log\big(2\pi(\Delta x_{i,j})^3\big)}_{=:K,\ \text{constante (ne dépend ni de }\mu_{IG}\text{ ni de }\eta)} - \frac{\eta}{2\mu_{IG}^2}\,Q(\mu_{IG})$$

où $Q(\mu_{IG})=\sum_{i,j}\dfrac{(\Delta x_{i,j}-\mu_{IG}\Delta t_j)^2}{\Delta x_{i,j}}$. En développant le carré :
$$Q(\mu_{IG}) = \underbrace{\sum_{i,j}\Delta x_{i,j}}_{A} - 2\mu_{IG}\underbrace{\sum_{i,j}\Delta t_j}_{B} + \mu_{IG}^2\underbrace{\sum_{i,j}\frac{\Delta t_j^2}{\Delta x_{i,j}}}_{C} = A-2\mu_{IG}B+\mu_{IG}^2C$$

($B=m\sum_j\Delta t_j=m\,t_n$ par télescopage.)

### 2.4 Score en $\mu_{IG}$ — et optimalité globale, indépendante de $\eta$

$$\frac{\partial\ell}{\partial\mu_{IG}} = -\frac\eta2\,\frac{d}{d\mu_{IG}}\left[\frac{Q(\mu_{IG})}{\mu_{IG}^2}\right]$$

Posons $R(\mu)=Q(\mu)/\mu^2 = A\mu^{-2}-2B\mu^{-1}+C$. Alors :
$$R'(\mu) = -2A\mu^{-3}+2B\mu^{-2} = \frac{2}{\mu^3}\big(B\mu-A\big)$$

**Racine unique**, en $\mu=A/B$, avec changement de signe de $-$ (pour $\mu<A/B$) à $+$ (pour $\mu>A/B$) : $R$ est donc **strictement décroissante puis strictement croissante**, avec un **minimum global unique** en $\mu=A/B$ sur $(0,\infty)$ — confirmé aussi par $R''(A/B)=2B^4/A^3>0$, mais l'argument de monotonicité établit en plus qu'il n'y a **aucun autre** point critique.

Comme $\eta>0$, $\ell(\mu_{IG})=-\frac\eta2R(\mu_{IG})+\text{cste}$ est maximisée exactement là où $R$ est minimisée :
$$\boxed{\hat\mu_{IG} = \frac AB = \frac{\sum_{i,j}\Delta x_{i,j}}{m\,t_n}} \qquad \text{— maximum global, \textbf{et ce quel que soit }}\eta>0$$

C'est ce dernier point (indépendance de $\hat\mu_{IG}$ par rapport à $\eta$) qui justifiera la validité de la procédure séquentielle à l'étape 2.6.

### 2.5 Score en $\eta$ — et optimalité globale

À $\mu_{IG}=\hat\mu_{IG}$ fixé, notons $Q_0=Q(\hat\mu_{IG})>0$ (strictement positif car c'est une somme de carrés, nulle seulement dans le cas dégénéré où $\Delta x_{i,j}=\hat\mu_{IG}\Delta t_j$ pour tout $(i,j)$ — exclu par la présence de bruit dans des données réelles). Alors :
$$\ell(\eta) = \frac{mn}2\log\eta - \frac{Q_0}{2\hat\mu_{IG}^2}\,\eta + K$$

C'est de la forme $a\log\eta - b\,\eta$ avec $a=mn/2>0$, $b=Q_0/(2\hat\mu_{IG}^2)>0$. Alors :
$$\ell'(\eta) = \frac a\eta - b,\qquad \ell''(\eta) = -\frac a{\eta^2} < 0 \quad\text{pour tout } \eta>0$$

$\ell$ est donc **strictement concave sur tout son domaine** $(0,\infty)$ : le point critique est nécessairement l'unique maximum **global**.
$$\ell'(\eta)=0 \iff \eta = \frac ab \quad\Longrightarrow\quad \boxed{\hat\eta = \frac{mn\,\hat\mu_{IG}^2}{Q(\hat\mu_{IG})}}$$

### 2.6 Pourquoi la démarche séquentielle donne bien le maximum global conjoint

C'est le point qui mérite d'être justifié explicitement, et qui ne l'avait pas été avant. On a montré à l'étape 2.4 que $\hat\mu_{IG}=A/B$ maximise $\ell(\cdot,\eta)$ pour **n'importe quel** $\eta>0$ — la valeur optimale de $\mu_{IG}$ ne dépend pas de $\eta$. Formellement :
$$\max_{\mu_{IG}>0}\ \ell(\mu_{IG},\eta) = \ell(\hat\mu_{IG},\eta) \quad \text{pour tout } \eta>0$$

Donc :
$$\max_{\mu_{IG}>0,\,\eta>0}\ell(\mu_{IG},\eta) = \max_{\eta>0}\Big[\max_{\mu_{IG}>0}\ell(\mu_{IG},\eta)\Big] = \max_{\eta>0}\ \ell(\hat\mu_{IG},\eta)$$

et ce dernier maximum est atteint en $\hat\eta$ par l'étape 2.5 (concavité stricte globale). La paire $(\hat\mu_{IG},\hat\eta)$ est donc bien le **maximum global conjoint** de $\ell$ sur $(0,\infty)^2$ — pas seulement un point stationnaire trouvé par une heuristique en deux étapes.

---

## Ce que cette version rigoureuse ajoute à la précédente

| | Version précédente | Ici |
|---|---|---|
| $\hat\mu_W(r)$ | dérivée seconde négative (local) | quadratique exacte → concavité **globale** |
| $\hat\sigma^2(r)$ | dérivée seconde au point critique | changement de signe de $\ell'$ → maximum **global** par monotonicité |
| $\hat r$ | affirmé comme équivalent à $\arg\min\hat\sigma^2(r)$ | même résultat, **plus** la limite honnête : pas de preuve d'unimodalité du profil |
| $\hat\mu_{IG}$ | dérivée seconde positive de $R$ au point critique | comportement de $R'$ sur tout $(0,\infty)$ → minimum **global** unique, indépendant de $\eta$ |
| $\hat\eta$ | dérivée seconde négative au point critique | concavité stricte de $\ell(\eta)$ sur **tout** le domaine → maximum **global** |
| Validité de la procédure séquentielle | implicite | démontrée explicitement (§2.6) via l'indépendance de $\hat\mu_{IG}$ vis-à-vis de $\eta$ |
