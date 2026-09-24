# Plan de validation par simulation — reproduction des Tables 3 à 8

## Référence de l'article
Zheng J., Wang Y., Chen Y., Xie M. (2026). *Automated Verification of Degradation Models Using
Goodness-of-Fit Test Based on Similarity Measures and Parametric Bootstrap*. Quality and
Reliability Engineering International. https://doi.org/10.1002/qre.70382

## Remarque préalable
L'article ne précise pas les instants exacts $t_1,\dots,t_n$ utilisés en simulation (seulement
$n$). Convention retenue ici, faute d'indication contraire dans le texte : $t_j = j$ pour
$j=1,\dots,n$ (pas de temps unitaire) — **hypothèse d'implémentation, pas une donnée de
l'article**. Procédure d'exécution commune à toutes les expériences ci-dessous : répéter $h$
fois {générer $m$ trajectoires sous le processus générateur → ajuster $H_0$ par MLE
(paramètres libres estimés à chaque répétition) → calculer les statistiques via l'Étape A de la
fiche `_algorithm_01` → calibrer le seuil par bootstrap ($L$ réplications, Étape B) → décision},
puis rapporter le pourcentage de rejets. Marge Monte Carlo à garder en tête pour juger un écart
entre une reproduction et les valeurs cibles : $\pm 2\sqrt{\hat p(1-\hat p)/h}$ autour de chaque
valeur.

---

## Expérience 1 — Table 3 ($n=50$, $m=6$, $h=1000$, $L=1000$, $\alpha_S=5\%$, $m'=600$)

### Processus générateur
Incréments, $\Lambda(t)=t$ (cas homogène), $\Delta t = t_j - t_{j-1}$ :

$$
\Delta SX \sim \begin{cases}
\mathcal N(\mu_W\Delta t,\ \sigma^2\Delta t) & \text{Wiener}\\[4pt]
\mathcal{Gamma}(\text{shape}=\alpha\Delta t,\ \text{scale}=\beta),\ \ x_{i,j}=SX(t_j)+\varepsilon_{i,j} & \text{Gamma bruité}\\[4pt]
\mathcal{IG}(\text{mean}=\mu_{IG}\Delta t,\ \text{shape}=\eta\Delta t^2),\ \ x_{i,j}=SX(t_j)+\varepsilon_{i,j} & \text{IG bruité}
\end{cases}
$$
avec $\varepsilon_{i,j}\overset{iid}{\sim}\mathcal N(0,\sigma_\varepsilon^2)$ pour les deux derniers cas.

### $H_0$ testé
$SX(t;\mu_W,\sigma) = \mu_W t + \sigma B(t)$ (Wiener homogène), pour toutes les lignes du tableau.

### Résultats attendus

| Processus générateur | Paramètres | $S_{KS-sup}$ | $S_{KS-sum}$ | $S_{CvM-sup}$ | $S_{CvM-sum}$ |
|---|---|---|---|---|---|
| Wiener | $\mu_W=0.2,\ \sigma=0.05$ | 5,2 | 3,1 | 3,8 | 4,2 |
| Wiener | $\mu_W=0.2,\ \sigma=0.1$ | 4,8 | 5,3 | 4,1 | 3,5 |
| Wiener | $\mu_W=0.4,\ \sigma=0.05$ | 3,7 | 5,1 | 4,5 | 5,2 |
| Wiener | $\mu_W=0.4,\ \sigma=0.1$ | 3,9 | 4,1 | 5,5 | 5,7 |
| Gamma | $\alpha=1,\ \beta=1/3,\ \sigma_\varepsilon=0.05$ | 11,5 | 3,9 | 10,2 | 4,6 |
| Gamma | $\alpha=1,\ \beta=3,\ \sigma_\varepsilon=0.1$ | 12,8 | 4,8 | 10,3 | 3,8 |
| IG | $\mu_{IG}=0.5,\ \eta=0.5,\ \sigma_\varepsilon=0.05$ | 21,0 | 4,4 | 4,8 | 5,0 |
| IG | $\mu_{IG}=0.5,\ \eta=1,\ \sigma_\varepsilon=0.05$ | 12,9 | 7,6 | 4,9 | 4,2 |
| IG | $\mu_{IG}=1,\ \eta=0.5,\ \sigma_\varepsilon=0.1$ | 32,5 | 4,5 | 5,1 | 4,8 |
| IG | $\mu_{IG}=1,\ \eta=1,\ \sigma_\varepsilon=0.1$ | 20,7 | 4,5 | 5,1 | 4,2 |

*(4 premières lignes = niveau empirique ; 6 suivantes = puissance ; valeurs en %)*

### Conclusion
Le niveau empirique reste dans $[3{,}7\%,5{,}5\%]$ sous $H_0$ vraie — le test est correctement
calibré. Contre Gamma/IG bruités, la puissance reste **faible** (3,8 % à 32,5 %) : conséquence
directe de la ressemblance géométrique des trajectoires une fois le Wiener ajusté (cf. Figure 3
de l'article, déjà discutée), le test étant structurellement peu sensible à ce type d'écart à
cette taille d'échantillon. $S_{KS-sup}$ et $S_{CvM-sup}$ dominent systématiquement les variantes
`-sum`, confirmant le motif "pire cas local" vs "moyenné" identifié en §2.2 : un écart de
spécification concentré sur quelques unités est capté par les statistiques de type `sup`,
noyé par celles de type `sum`.

---

## Expérience 2 — Table 4 (identique à l'Expérience 1, mais $n=100$)

### Processus générateur
Mêmes formules que l'Expérience 1 (Gamma et IG uniquement), $\Delta t$ recalculé pour $n=100$.

### $H_0$ testé
Identique : $\mu_W t + \sigma B(t)$.

### Résultats attendus

| Processus générateur | Paramètres | $S_{KS-sup}$ | $S_{KS-sum}$ | $S_{CvM-sup}$ | $S_{CvM-sum}$ |
|---|---|---|---|---|---|
| Gamma | $\alpha=1,\ \beta=1/3,\ \sigma_\varepsilon=0.05$ | 21,6 | 5,2 | 20,1 | 5,9 |
| Gamma | $\alpha=1,\ \beta=3,\ \sigma_\varepsilon=0.1$ | 22,4 | 6,2 | 16,7 | 4,9 |
| IG | $\mu_{IG}=1,\ \eta=0.5,\ \sigma_\varepsilon=0.1$ | 45,8 | 5,3 | 7,2 | 6,0 |
| IG | $\mu_{IG}=1,\ \eta=1,\ \sigma_\varepsilon=0.1$ | 36,8 | 5,9 | 7,9 | 5,3 |

### Conclusion
Doubler $n$ (50→100) augmente nettement la puissance de $S_{KS-sup}$ et $S_{CvM-sup}$ (ex. IG :
20,7 %→36,8 %) — davantage de points d'observation par trajectoire offre davantage de chances de
détecter un écart local entre trajectoire réelle et son plus proche voisin simulé. L'amélioration
**n'est pas uniforme** entre statistiques : $S_{CvM-sup}$ reste faible pour l'IG (7,2–7,9 %),
montrant que l'effet de $n$ interagit avec le choix de la mesure de similarité — augmenter $n$
n'est pas un levier universel, il profite surtout aux statistiques déjà les plus sensibles.

---

## Expérience 3 — Table 5 (identique à l'Expérience 1, mais $m=1$)

### Processus générateur
Mêmes formules que l'Expérience 1, avec $m=1$ (donc $m'=100$). Seules $S_{KS}$ et $S_{CvM}$
sont calculées (coïncidence sup/sum quand $m=1$, cf. éq. 8, 11).

### $H_0$ testé
Identique.

### Résultats attendus

| Processus générateur | Paramètres | $S_{KS}$ | $S_{CvM}$ |
|---|---|---|---|
| Wiener | $\mu_W=0.2,\ \sigma=0.05$ | 3,1 | 3,4 |
| Wiener | $\mu_W=0.2,\ \sigma=0.1$ | 3,9 | 3,5 |
| Wiener | $\mu_W=0.4,\ \sigma=0.05$ | 4,1 | 3,5 |
| Wiener | $\mu_W=0.4,\ \sigma=0.1$ | 3,1 | 4,2 |
| Gamma | $\alpha=1,\ \beta=1/3,\ \sigma_\varepsilon=0.05$ | 4,2 | 3,1 |
| Gamma | $\alpha=1,\ \beta=3,\ \sigma_\varepsilon=0.1$ | 5,3 | 5,1 |
| IG | $\mu_{IG}=0.5,\ \eta=0.5,\ \sigma_\varepsilon=0.05$ | 3,3 | 4,5 |
| IG | $\mu_{IG}=0.5,\ \eta=1,\ \sigma_\varepsilon=0.05$ | 5,2 | 3,6 |
| IG | $\mu_{IG}=1,\ \eta=0.5,\ \sigma_\varepsilon=0.1$ | 3,9 | 3,6 |
| IG | $\mu_{IG}=1,\ \eta=1,\ \sigma_\varepsilon=0.1$ | 5,4 | 4,7 |

### Conclusion
Avec une seule unité testée, la puissance sous alternative reste quasiment **au niveau du bruit
Monte Carlo autour de $\alpha_S$** (3,1 % à 5,4 %) — le test devient presque aveugle aux
alternatives Gamma/IG. Ceci confirme que $m$ (le nombre d'unités agrégées, indépendamment de $n$)
est un second levier de puissance à part entière : agréger sur plusieurs unités multiplie les
occasions de détecter un écart, même si celui-ci n'est marqué que sur une partie des unités —
un résultat cohérent avec, et distinct de, l'effet de $n$ observé en Expérience 2.

---

## Expérience 4 — Table 6 (Wiener généralisé, $n=100$, $m=6$)

### Processus générateur
Forme non-homogène, $\Delta\Lambda(t_{j-1},t_j)=\Lambda(t_j)-\Lambda(t_{j-1})$ :

$$
\Delta SX \sim \begin{cases}
\mathcal N(\mu_W\Delta\Lambda,\ \sigma^2\Delta t) & \text{Wiener non-homogène}\\[4pt]
\mathcal{Gamma}(\text{shape}=\alpha\Delta\Lambda,\ \text{scale}=\beta),\ \ x_{i,j}=SX(t_j)+\varepsilon_{i,j} & \text{Gamma non-homogène bruité}\\[4pt]
\mathcal{IG}(\text{mean}=\mu_{IG}\Delta\Lambda,\ \text{shape}=\eta\Delta\Lambda^2),\ \ x_{i,j}=SX(t_j)+\varepsilon_{i,j} & \text{IG non-homogène bruité}
\end{cases}
$$

### $H_0$ testé
$SX(t;\mu_W,\sigma,r)=\mu_W t^r+\sigma B(t)$, avec $r$ **estimé** conjointement avec $\mu_W,\sigma$.

### Résultats attendus

| $\Lambda(t)$ générateur | Processus | Paramètres | $S_{KS-sup}$ | $S_{KS-sum}$ | $S_{CvM-sup}$ | $S_{CvM-sum}$ |
|---|---|---|---|---|---|---|
| $\sqrt t$ | Wiener | $\mu_W=10,\ \sigma=1$ | 4,7 | 4,8 | 5,8 | 5,6 |
| $t^2$ | Wiener | $\mu_W=0.01,\ \sigma=1$ | 4,8 | 4,3 | 5,1 | 5,1 |
| $\sqrt t$ | Gamma | $\alpha=2,\ \beta=5,\ \sigma_\varepsilon=1$ | 32,2 | 22,4 | 23,5 | 14,8 |
| $\sqrt t$ | Gamma | $\alpha=0.5,\ \beta=20,\ \sigma_\varepsilon=1$ | 44,9 | 18,8 | 41,8 | 30,6 |
| $t^2$ | Gamma | $\alpha=0.002,\ \beta=5,\ \sigma_\varepsilon=1$ | 20,0 | 6,4 | 15,3 | 8,4 |
| $t^2$ | Gamma | $\alpha=0.0005,\ \beta=20,\ \sigma_\varepsilon=1$ | 39,7 | 4,8 | 30,6 | 8,8 |
| $\sqrt t$ | IG | $\mu_{IG}=10,\ \eta=50,\ \sigma_\varepsilon=1$ | 31,3 | 18,3 | 27,5 | 20,3 |
| $\sqrt t$ | IG | $\mu_{IG}=10,\ \eta=5,\ \sigma_\varepsilon=1$ | 62,6 | 25,7 | 53,3 | 27,6 |
| $t^2$ | IG | $\mu_{IG}=0.01,\ \eta=5\times10^{-5},\ \sigma_\varepsilon=1$ | 20,4 | 7,5 | 13,3 | 6,3 |
| $t^2$ | IG | $\mu_{IG}=0.01,\ \eta=5\times10^{-6},\ \sigma_\varepsilon=1$ | 48,3 | 6,3 | 41,7 | 6,9 |

*(2 premières lignes = niveau ; 8 suivantes = puissance)*

### Conclusion
Le niveau empirique reste proche de $5\%$ (4,3–5,8 %) même lorsque $r$ doit être estimé en plus
de $\mu_W,\sigma$ : la validité du test s'étend donc au cas non-homogène sans dégradation de
calibration. La puissance suit le même motif de domination des statistiques `-sup` qu'en
Expérience 1, mais atteint des valeurs nettement plus élevées (jusqu'à 62,6 % pour l'IG) —
la non-homogénéité en tant que telle n'est pas un obstacle au test ; dans certains régimes de
paramètres, elle s'accompagne même d'une meilleure discrimination Wiener/Gamma-IG que dans le
cas homogène de l'Expérience 1.

---

## Expérience 5 — Tables 7 et 8 (Wiener vs trajectoire générale, $n=50$)

### Processus générateur
$$
x_{i,j} = g(t_j;\boldsymbol\omega,\boldsymbol\phi_i) + \varepsilon_{i,j}, \qquad \varepsilon_{i,j}\overset{iid}{\sim}\mathcal N(0,\sigma_\varepsilon^2)
$$
- **Linéaire** : $g=\phi_{L,i}\,t_j$, $\phi_{L,i}\overset{iid}{\sim}\mathcal N(0.2,0.05^2)$
- **Loi de puissance** : $g=\phi_{P1,i}\,t_j^{\phi_{P2,i}}$, $\phi_{P1,i},\phi_{P2,i}\overset{iid}{\sim}\mathcal N(0.5,0.1^2)$
- **Exponentiel** : $g=\phi_{E,i}[\exp(\omega_E t_j)-1]$, $\phi_{E,i}\overset{iid}{\sim}\mathcal N(0.2,0.05^2)$, $\omega_E=0.05$ fixe

### $H_0$ testé
$SX(t;\mu_W,\sigma,\Lambda)=\mu_W\Lambda(t)+\sigma B(t)$ ; trois formes de $\Lambda$ testées
séparément selon le générateur : $\Lambda(t)=t$ (vs linéaire), $\Lambda(t)=t^r$, $r$ estimé
(vs loi de puissance), $\Lambda(t)=\exp(\lambda t)-1$, $\lambda$ estimé (vs exponentiel).

### Résultats attendus — Table 7 ($m=1$)

| Modèle générateur | Paramètres | $\sigma_\varepsilon$ | $S_{KS}$ | $S_{CvM}$ |
|---|---|---|---|---|
| Linéaire | $\phi_L\sim\mathcal N(0.2,0.05^2)$ | 0,05 | 91,1 | 90,3 |
| Linéaire | idem | 0,1 | 31,6 | 25,7 |
| Loi de puissance | $\phi_{P1},\phi_{P2}\sim\mathcal N(0.5,0.1^2)$ | 0,05 | 94,6 | 93,1 |
| Loi de puissance | idem | 0,1 | 32,7 | 35,6 |
| Exponentiel | $\phi_E\sim\mathcal N(0.2,0.05^2),\ \omega_E=0.05$ | 0,05 | 56,7 | 62,3 |
| Exponentiel | idem | 0,1 | 20,3 | 18,2 |

### Résultats attendus — Table 8 ($m=6$)

| Modèle générateur | $\sigma_\varepsilon$ | $S_{KS-sup}$ | $S_{KS-sum}$ | $S_{CvM-sup}$ | $S_{CvM-sum}$ |
|---|---|---|---|---|---|
| Linéaire | 0,05 | 94,1 | 92,3 | 94,4 | 93,1 |
| Linéaire | 0,1 | 49,6 | 33,6 | 45,6 | 33,5 |
| Loi de puissance | 0,05 | 97,4 | 94,1 | 99,1 | 98,3 |
| Loi de puissance | 0,1 | 35,6 | 31,3 | 43,0 | 38,7 |
| Exponentiel | 0,05 | 99,2 | 98,7 | 99,3 | 98,9 |
| Exponentiel | 0,1 | 71,7 | 56,2 | 63,3 | 50,8 |

### Conclusion
C'est le scénario le plus favorable de toute l'étude : puissance proche de 100 % à faible bruit
($\sigma_\varepsilon=0.05$) pour les trois formes de générateur — le test discrimine très
nettement un Wiener d'un modèle à trajectoire générale, confirmant empiriquement la Contribution
3 de l'introduction. Mais la puissance **chute fortement** dès que $\sigma_\varepsilon$ double
(ex. exponentiel : 99,2 %→71,7 %) : la performance du test dépend donc autant du rapport
signal/bruit que de la distance structurelle entre les deux classes de modèles comparées — un
bruit de mesure important peut masquer une différence de nature autrement très détectable.
Comme en Expérience 3, $m=6$ (Table 8) donne des puissances systématiquement supérieures à
$m=1$ (Table 7).

---

## Synthèse globale des cinq expériences

- **Validité** : dans toutes les expériences, le niveau empirique sous $H_0$ vraie reste proche
  de $\alpha_S=5\%$ — y compris en présence de paramètres de tendance supplémentaires à estimer
  ($r$, $\lambda$).
- **Statistiques recommandées** : $S_{KS-sup}$ et $S_{CvM-sup}$ dominent systématiquement
  $S_{KS-sum}$ et $S_{CvM-sum}$ — confirmé dans les cinq expériences sans exception notable.
- **Facteurs qui augmentent la puissance** : davantage d'observations $n$ (Expérience 2),
  davantage d'unités testées $m$ (Expériences 3 et 5), une différence structurelle marquée entre
  modèle testé et modèle générateur (Expérience 5 vs Expérience 1).
- **Facteur qui réduit la puissance** : un bruit de mesure $\sigma_\varepsilon$ élevé, qui dilue
  le signal discriminant même dans le cas le plus favorable (Expérience 5).
- **Cas le plus difficile pour le test** : distinguer un Wiener d'un Gamma/IG bruité à faible
  $n$ et $m$ (Expérience 1) — la ressemblance géométrique des trajectoires, une fois les deux
  modèles ajustés, limite structurellement la puissance atteignable par une méthode fondée
  uniquement sur la similarité de trajectoires.
