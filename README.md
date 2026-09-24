# Article Zheng, Wang, Chen & Xie (2026)

Implémentation R et notes d'étude de :

Zheng J., Wang Y., Chen Y., Xie M. (2026). *Automated Verification of
Degradation Models Using Goodness-of-Fit Test Based on Similarity Measures
and Parametric Bootstrap*. Quality and Reliability Engineering International.
https://doi.org/10.1002/qre.70382

## Contenu

- `Implémentation_R/` : scripts R (algorithmes 1 à 3, estimateurs MLE, tests, application aux données Laser)
- `Explication_md/` : explications détaillées des algorithmes et démonstrations (`.md` et `.html`)
- `gofdegrad/` : package R issu de ces scripts

## Installation du package

Le dépôt étant privé, un jeton GitHub (PAT) est nécessaire
(`usethis::create_github_token()` puis `gitcreds::gitcreds_set()`) :

```r
# install.packages("remotes")
remotes::install_github("rytracetito/Article_zheng_chen_et_al", subdir = "gofdegrad")
```
