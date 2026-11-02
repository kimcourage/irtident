# irtident

<!-- badges: start -->
[![R-CMD-check](https://github.com/kimcourage/irtident/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kimcourage/irtident/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Identification diagnostics for item response models fitted with
[mirt](https://cran.r-project.org/package=mirt).

Convergence does not establish that every item parameter is well determined.
irtident examines parameter uncertainty using variance decomposition, profile
likelihoods and repeated fits from different starting values. It also simulates
parameter recovery across sample sizes.

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
plot(dx, params = c("Item.1:g", "Item.2:g", "Item.4:g", "Item.2:a1"))
```

In this example, four items have guessing estimates near zero. Item 2 has a
positive guessing estimate, but its profile interval also reaches zero. A
positive point estimate alone therefore does not establish a nonzero guessing
parameter.

![Profile deviances](man/figures/README-profile.png)

The next simulation evaluates recovery of item 2's guessing and slope
parameters, using the fitted model as the data-generating model.

```r
set.seed(1)
ss <- id_samplesize(fit, params = c("Item.2:g", "Item.2:a1"),
                    N = c(1000, 2000, 4000, 8000), reps = 20)
ss
plot(ss)
```

![RMSE by sample size](man/figures/README-samplesize.png)

These figures were generated with R 4.5.2 and mirt 1.45.1. Twenty replications
illustrate the procedure; the simulated summaries will vary between runs.
Results also depend on the assumed true parameters and the fitting settings.

## Functions

| Function | Output |
|---|---|
| `id_local()` | Conditional standard errors, variance inflation within and between items, and parameter combinations flagged by the eigenvalue diagnostic |
| `id_profile()` | Profile likelihood intervals; each side is `closed`, `bound` or `open` |
| `id_multistart()` | Clusters of estimates, their log-likelihoods and convergence status |
| `id_samplesize()` | Bias, RMSE, fitting success and convergence across sample sizes |
| `id_diagnose()` | Local, profile and multi-start results combined by item, with suggested follow-up checks |

The local variance decomposition is

$$SE_j^2 = \underbrace{1/I_{jj}}_{\text{conditional variance}} \times VIF^{within}_j \times VIF^{between}_j.$$

The conditional variance assumes all other parameters are known. The inflation
factors describe uncertainty associated with estimating the other parameters
of the same item and of the rest of the model. These quantities use mirt's
internal scale: `g` and `u` are logits. With `SE.type = "Fisher"`, the
calculation can use planned parameter values and sample size.

`id_samplesize()` reports how parameter recovery changes with N. It does not
apply an RMSE cutoff or recommend a required sample size. RMSE includes
non-converged fits that returned estimates; failed fits are excluded and
reported in the fitting success rate. Interpret RMSE alongside convergence
and fitting success.

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
