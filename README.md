# irtident

<!-- badges: start -->
[![R-CMD-check](https://github.com/kimcourage/irtident/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kimcourage/irtident/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Identification diagnostics for item response models fitted with
[mirt](https://cran.r-project.org/package=mirt).

Convergence does not establish that every item parameter is well determined.
irtident examines parameter uncertainty using variance decomposition, profile
likelihoods and repeated fits from different starting values. It also simulates
parameter precision across sample sizes.

## Installation

```r
# install.packages("remotes")
remotes::install_github("kimcourage/irtident", build_vignettes = TRUE)
```

## Example

```r
library(mirt)
library(irtident)

dat <- expand.table(LSAT7)
fit <- mirt(dat, 1, itemtype = "3PL", verbose = FALSE)

set.seed(1)
dx <- id_diagnose(fit)
dx
```

<!-- dx-output:start -->

```text
irtident identification diagnostics
Model: 3PL, 5 items, N = 1000, logLik = -2658.05, converged: yes
mirt information-matrix condition number: 3.84e+06

Items:
  Item.1  3PL  boundary
      g   at bound (0.000), flat profile likelihood
  Item.2  3PL  weak
      a1  profile interval open on one side, large SE
      d   profile interval open on one side, confounded with Item.2:g (VIF 21.5)
      g   profile interval reaches the bound (g = 0.00), confounded with Item.2:d (VIF 23.6)
  Item.3  3PL  boundary
      g   at bound (0.000), flat profile likelihood
  Item.4  3PL  boundary
      g   at bound (0.000)
  Item.5  3PL  boundary
      g   at bound (0.000), flat profile likelihood

Advice:
  Item.1:g: consider a simpler model (fix g = 0), or assess sensitivity to a prior on g
  Item.2:a1: slope may be unbounded (near-Guttman item); consider a prior on a1
  Item.2:d: estimate d and Item.2:g; inspect their joint uncertainty and sensitivity to priors
  Item.2:g: the data do not rule out g = 0.00; report the profile interval rather than the Wald interval
  Item.3:g: consider a simpler model (fix g = 0), or assess sensitivity to a prior on g
  Item.4:g: consider a simpler model (fix g = 0), or assess sensitivity to a prior on g
  Item.5:g: consider a simpler model (fix g = 0), or assess sensitivity to a prior on g

Multi-start: 20 starts, 1 equivalent-fit cluster(s), 0 inferior optimum/optima; the original fit is among the best solutions.

Details: $params, $local, $profile, $multistart; plot(x) shows the profiles.
```

<!-- dx-output:end -->

In this example, four items have guessing estimates near zero. Item 2 has a
positive guessing estimate, but its profile interval also reaches zero. A
positive point estimate alone therefore does not establish a nonzero guessing
parameter.

```r
plot(dx, params = c("Item.1:g", "Item.2:g", "Item.4:g", "Item.2:a1"))
```

![Profile deviances](man/figures/README-profile.png)

The next simulation shows how the variability of item 2's guessing and slope
estimates changes as N increases, using the fitted model to generate data.

```r
set.seed(1)
ss <- id_samplesize(fit, params = c("Item.2:g", "Item.2:a1"),
                    N = c(1000, 2000, 4000, 8000), reps = 20)
ss
```

<!-- ss-output:start -->

```text
<irtident_samplesize> 20 replications per sample size; current N = 1000

     param    N current  mean    sd mirt_se n_used success conv se_failure
  Item.2:g 1000   0.295 0.153 0.162   0.278     15       1 0.85      0.118
  Item.2:g 2000   0.295 0.243 0.157   0.171     18       1 0.95      0.053
  Item.2:g 4000   0.295 0.214 0.112   0.211     16       1 1.00      0.200
  Item.2:g 8000   0.295 0.265 0.071   0.088     17       1 0.90      0.056
 Item.2:a1 1000   1.928 2.411 2.890   5.154     15       1 0.85      0.118
 Item.2:a1 2000   1.928 2.118 1.147   1.748     18       1 0.95      0.053
 Item.2:a1 4000   1.928 1.712 0.509   0.662     16       1 1.00      0.200
 Item.2:a1 8000   1.928 1.807 0.333   0.469     17       1 0.90      0.056

mean, sd and mirt_se use the same n_used converged fits with usable SEs.
mirt_se = sqrt(mean(SE^2)); se_failure is conditional on convergence.
Estimates and SEs are on the reporting scale.
```

<!-- ss-output:end -->

`current` is the estimate from the original fit. For each N, `mean` and `sd`
describe the repeated estimates; `mirt_se` summarises the SEs reported by mirt
as `sqrt(mean(SE^2))`. The plot compares the two measures of uncertainty.

```r
plot(ss)
```

![Repeated-estimate SD and mirt SE by sample size](man/figures/README-samplesize.png)

These figures were generated with R 4.5.2 and mirt 1.45.1. Twenty replications
illustrate the procedure; the simulated summaries will vary between runs.
Results also depend on the generating model and the fitting settings.

## Functions

| Function | Output |
|---|---|
| `id_local()` | Conditional standard errors, variance inflation within and between items, and parameter combinations flagged by the eigenvalue diagnostic |
| `id_profile()` | Profile likelihood intervals; each side is `closed`, `bound` or `open` |
| `id_multistart()` | Clusters of estimates, their log-likelihoods and convergence status |
| `id_samplesize()` | Current estimate, repeated-estimate mean and SD, mirt SE, fitting success, convergence and SE failure across sample sizes |
| `id_diagnose()` | Local, profile and multi-start results combined by item, with suggested follow-up checks |

The local variance decomposition is

$$SE_j^2 = \underbrace{1/I_{jj}}_{\text{conditional variance}} \times VIF^{within}_j \times VIF^{between}_j.$$

The conditional variance assumes all other parameters are known. The inflation
factors describe uncertainty associated with estimating the other parameters
of the same item and of the rest of the model. These quantities use mirt's
internal scale: `g` and `u` are logits. With `SE.type = "Fisher"`, the
calculation can use planned parameter values and sample size.

`id_samplesize()` helps assess how much precision a larger sample may provide.
The mean, SD and mirt SE use the same `n_used` converged fits with a usable
estimate and SE for that parameter. `success` and `conv` count successful
and converged fits out of all attempts; `se_failure` is the fraction of
converged fits without a usable estimate/SE pair. A small SD or SE based on
few usable fits can be misleading, so inspect these counts and rates together.

Estimates and SEs use the same reporting scale. For `g` and `u`, mirt's logit
SEs are converted to the probability scale by the delta method. This
approximation can be unreliable near a boundary. The function does not set
a target precision or choose a required N.

## Scope

Supported models are unidimensional, single-group fits from `mirt::mirt()`,
with itemtypes `Rasch`, `2PL`, `3PL`, `3PLu`, `4PL` and `graded`, without equality
constraints.

Diagnostic labels use numerical thresholds and a finite profile grid. They
are screening results, not proofs of structural identification. Profile
intervals use a chi-squared approximation, whose nominal coverage need not
hold at parameter boundaries or under weak identification.

Refitting 3PL and 4PL models can be slow. Functions that refit accept
`parallel = TRUE` through `future.apply`:

```r
future::plan(future::multisession, workers = 4)
dx <- id_diagnose(fit, parallel = TRUE)
```

See `vignette("irtident")` for the variance decomposition, profile diagnostics
and sample-size examples.

## References

- Chalmers, R. P. (2012). mirt: A multidimensional item response theory package for the R environment. *Journal of Statistical Software, 48*(6), 1–29.
- Raue, A., Kreutz, C., Maiwald, T., Bachmann, J., Schilling, M., Klingmüller, U., & Timmer, J. (2009). Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood. *Bioinformatics, 25*(15), 1923–1929.

## License

MIT
