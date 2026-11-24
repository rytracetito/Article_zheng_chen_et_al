# gofdegrad

Goodness-of-fit tests for degradation models by similarity measures and
parametric bootstrap, after:

Zheng J., Wang Y., Chen Y., Xie M. (2026). *Automated Verification of
Degradation Models Using Goodness-of-Fit Test Based on Similarity Measures
and Parametric Bootstrap*. Quality and Reliability Engineering International.
<https://doi.org/10.1002/qre.70382>

Is a gamma process, a Wiener process or another model adequate for a set of
degradation trajectories? The test compares the observed trajectories with
trajectories simulated under the fitted model, through Kolmogorov-Smirnov and
Cramér-von Mises type similarity measures, and calibrates the decision by
parametric bootstrap. The trajectories must share the same inspection times,
which can be irregular.

The simulations rely on the package `DegradSim` and the maximum likelihood
estimators on the package `levyfit`
(<https://github.com/rytracetito/levy_packages>).

## Installation

```r
# install.packages("remotes")
remotes::install_github("rytracetito/Article_zheng_chen_et_al", subdir = "gofdegrad")
```

The dependencies `DegradSim` and `levyfit` are installed from GitHub at the same
time (field `Remotes` of the `DESCRIPTION` file). If the repositories are
private, a GitHub token is needed: `usethis::create_github_token()`, then
`gitcreds::gitcreds_set()`.

## Candidate models

| `model` | Process | Parameters |
|:---|:---|:---|
| `"wiener"` | Wiener with drift | `mu`, `sigma` |
| `"gamma"` | gamma (shape, **rate**) | `a`, `b` |
| `"ig"` | inverse Gaussian | `mu`, `eta` |
| `"vg"` | variance-gamma | `theta`, `sigma`, `nu` |
| `"meixner"` | Meixner | `a`, `b`, `d`, `m` |
| `"wiener_power"` | Wiener, trend `mu t^r` | `mu`, `sigma`, `r` |
| `"wiener_exp"` | Wiener, trend `mu (exp(lambda t) - 1)` | `mu`, `sigma`, `lambda` |
| `"gamma_power"` | gamma, shape function `a t^r` | `a`, `b`, `r` |
| `"ig_re"` | inverse Gaussian with random drift (Peng, 2015, model M3) | `xi`, `sd_delta`, `eta` |

`gof_models()` gives the law of the increments of each model. In the first
eight models all the units follow the same process, with the same parameters.

The model `"ig_re"` lets the mean degradation rate vary from unit to unit:
unit `i` degrades at its own rate `mu_i = 1/delta_i`, with

- `xi`: mean of `delta_i = 1/mu_i`, so that `1/xi` is a typical rate;
- `sd_delta`: standard deviation of the `delta_i`, the unit-to-unit
  variability (`delta_i ~ N(xi, sd_delta^2)`, truncated to `delta_i > 0`);
- `eta`: shape parameter common to all the units, the variability within a
  trajectory (given `delta_i`, `X_i(t) - X_i(s) ~ IG((t - s)/delta_i, eta (t - s)^2)`).

With `sd_delta` close to 0 it reduces to the model `"ig"`. It suits data where
some units are consistently faster than others, such as the GaAs laser data.
Any other model (general path model, other random effects, ...) can be tested
with `custom_model()`.

## Getting started

```r
library(gofdegrad)

# 1. Data: six gamma trajectories inspected every time unit
set.seed(1)
x <- sim_model(0:50, "gamma", c(a = 2, b = 3), n_traj = 6)
DegradSim::plot_traj(x)

# 2. Test of the gamma model (B = 1000 bootstrap samples, pools of 600 trajectories)
res <- gof_test(x, "gamma")
res
plot(res)

# 3. Several candidate models at once
compare_models(x, c("gamma", "ig", "wiener", "gamma_power"), seed = 1)
```

Real data are given as a data frame with the columns `id`, `t` and `X`, or as a
matrix with one row per trajectory and the times in `t`:

```r
laser <- SMRD.data::gaaslaser
laser <- data.frame(id = laser$unit, t = laser$hours, X = laser$increase)
set.seed(30)
gof_test(laser, "wiener")   # rejected, as in Zheng et al. (2026)

# The lasers degrade at different rates: one rate per laser
set.seed(30)
gof_test(laser, "ig_re")    # not rejected
```

## Level and power

```r
# Empirical level: data generated under the candidate model
power_study(function() sim_model(0:50, "gamma", c(a = 2, b = 3), n_traj = 6),
            model = "gamma", n_rep = 1000, B = 1000)

# Power against linear general paths
power_study(function() sim_path(0:50, "linear", phi_mean = 0.2, phi_sd = 0.05,
                                sigma_eps = 0.05, n_traj = 6),
            model = "wiener", n_rep = 1000, B = 1000)
```

## Functions

| Function | Purpose |
|:---|:---|
| `gof_test()` | goodness-of-fit test of one candidate model |
| `compare_models()` | tests of several candidate models on the same data |
| `power_study()` | empirical level or power of the test |
| `gof_models()` | the nine candidate models |
| `fit_model()` | maximum likelihood estimation of a candidate model |
| `sim_model()` | simulation of a candidate model |
| `sim_path()` | simulation of general path models |
| `custom_model()` | any other candidate model, from its own functions |
| `gof_data()` | conversion and check of the data |
| `similarity_stats()` | the four similarity statistics |

The reference manual is in `manuel/gofdegrad.pdf`, at the root of the
repository.

## Known limitations

- All the trajectories must share the same inspection times (the article
  announces the extension to asynchronous inspections as future work).
- The models start from 0 at time 0. When the data contain time 0, the level
  at time 0 is subtracted; otherwise the degradation is assumed to be 0 at
  time 0, and the first increment weighs on the estimation.
- With few trajectories or few inspections the test has little power: a
  non-rejection means that the data are compatible with the model, not that
  the model is validated.
- The Meixner parameters `a` and `d` are strongly correlated and need many
  increments to be estimated precisely.
- The normality of the random effects in `sim_path()` is an assumption of the
  package: the article does not specify their distribution.
