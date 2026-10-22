# Article Zheng, Wang, Chen & Xie (2026)

Implémentation R et notes d'étude de :

Zheng J., Wang Y., Chen Y., Xie M. (2026). *Automated Verification of
Degradation Models Using Goodness-of-Fit Test Based on Similarity Measures
and Parametric Bootstrap*. Quality and Reliability Engineering International.
https://doi.org/10.1002/qre.70382

## Contenu

- `gofdegrad/` : package R du test d'adéquation (version 1.0.0), fondé sur les
  packages `levysim` (simulation) et `levyfit` (estimation)
- `manuel/gofdegrad.pdf` : manuel de référence du package
- `Implementation_R/` : scripts d'utilisation du package (grilles régulière et
  irrégulière, applications aux données GaAs laser, résistances, usure, Hudak)
  et premiers scripts des algorithmes
- `Donnees/` : données simulées et script de test
- `Explication_md/` : explications détaillées des algorithmes et démonstrations
  (`.md` et `.html`)

## Installation du package

```r
# install.packages("remotes")
remotes::install_github("rytracetito/Article_zheng_chen_et_al", subdir = "gofdegrad")
```

Les dépendances `levysim` et `levyfit` sont installées en même temps depuis
GitHub. Si les dépôts sont privés, un jeton GitHub est nécessaire
(`usethis::create_github_token()` puis `gitcreds::gitcreds_set()`).
