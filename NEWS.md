# moc.gapbk 0.2.1

## Performance

This release rewrites the hot inner loops of the algorithm using
vectorised R primitives. Numerical results are bit-identical to 0.2.0
under the same random seed; only the runtime changes.

* `generate.groups()`: the original nested for-loop is replaced by a
  single `max.col(-dmatrix[, medoids])` call, which runs in C.
* `calculate.objective.functions()`: the original triple for-loop over
  individuals, clusters and objects is replaced by a single matrix
  slice plus `sum()`. The minimum medoid-to-medoid squared distance
  is now obtained from a (k x k) submatrix with `diag(sub) <- Inf` and
  `min(sub)`, also vectorised.
* `singletons.delete()` and `singletons.repair()` use `tabulate()`
  to count cluster sizes in C instead of `which()` + `length()` in R.
* `singletons.repair()` no longer resets the column index to 1 on
  every replacement; it now retries on the same column up to
  `num_objects` attempts, eliminating an O(n^2) worst case.
* `moc.gapbk()` precomputes the element-wise squared distance matrix
  for each input matrix once and caches it via an attribute, instead
  of squaring entries inside the inner loop on every call.
* Path-Relinking caches the symmetric difference between the
  initial and guide solutions across each `while` iteration (it was
  recomputed three times per cycle in 0.2.x).
* Path-Relinking builds the batch of intermediate solutions in a
  single replicate-and-overwrite step rather than chained
  `rbind()` calls on a growing data frame.
* Pareto Local Search replaces
  `apply(A, 1, function(x) all(x == B))` with vectorised
  `rowSums(A != B) == 0L`.

Observed speedups on the internal benchmark (synthetic data, seed=42):

| Scenario             | 0.2.2   | 0.3.0   | Speedup |
|----------------------|---------|---------|---------|
| NSGA-II, n=200, G=3  | 0.049 s | 0.003 s | 17.8x   |
| NSGA-II, n=300, G=3  | 0.070 s | 0.003 s | 22.3x   |
| NSGA-II, n=500, G=3  | 0.124 s | 0.006 s | 22.3x   |
| Full algo, n=50      | 0.011 s | 0.002 s |  4.3x   |
| Full algo, n=60      | 0.009 s | 0.002 s |  5.3x   |
| Full algo, n=70      | 0.009 s | 0.003 s |  3.6x   |

The asymptotic complexity of the algorithm itself is unchanged
(Pareto Local Search is still O(m^2 * n^2 * k) per generation, where
m is the size of the Pareto front). What changes is the constant
factor: the inner operations now execute in compiled C through R
primitives rather than at the interpreter level. For very large
instances with `local_search = TRUE` a future Rcpp-backed release
will be needed.

# moc.gapbk 0.2.0

## GitHub infrastructure and developer experience

* New `README.Rmd` (dynamic, generates `README.md` via
  `devtools::build_readme()`) with CRAN, R-CMD-check, pkgdown,
  Codecov, lifecycle and license badges.
* Added GitHub Actions workflows under `.github/workflows/`:
  - `R-CMD-check.yaml` on Ubuntu (devel, release, oldrel-1), macOS
    and Windows.
  - `pkgdown.yaml` that deploys the documentation site to
    `gh-pages`.
  - `test-coverage.yaml` that pushes coverage to Codecov.
* Added `_pkgdown.yml` (Bootstrap 5, cosmo theme, navbar, structured
  reference index, OpenGraph metadata, ORCID link).
* Added `codecov.yml` with informational thresholds.
* Added `CODE_OF_CONDUCT.md` (Contributor Covenant 2.1) and
  `.github/CONTRIBUTING.md`, `.github/pull_request_template.md`,
  `.github/ISSUE_TEMPLATE/{bug_report,feature_request}.md`.
* Added a hex-style package logo at `man/figures/logo.png` (and
  `logo.svg` source).
* Added the `Language: en-US` field to `DESCRIPTION`.
* Updated `.Rbuildignore` and `.gitignore` to keep all of the above
  out of the CRAN bundle while preserving them on GitHub.

* Moved `amap` from `Imports` to `Suggests`. The package code itself
  does not use `amap`; it only appeared in the example and in the
  vignette. The example was rewritten to use `stats::dist()` from
  base R, and the vignette wraps the `amap::Dist()` call with
  `requireNamespace()` so it knits in either configuration. This
  removes the `R CMD check` NOTE "Namespace in Imports field not
  imported from: 'amap'".

## Breaking changes

* The main function has been renamed from `moc.gabk()` to `moc.gapbk()`
  to match the package name. The previous name is preserved as a
  deprecated alias that emits a warning and forwards its arguments to
  `moc.gapbk()`. The alias will be removed in a future release.

## CRAN compliance

* Removed `doMPI` from the dependency list. The MPI backend was never
  actually used by the package code and `doMPI` currently fails CRAN
  checks on macOS because it requires `Rmpi`, which depends on a
  system-level MPI installation.
* `nsga2R`, `foreach`, `doParallel`, `parallel`, `utils`, `stats` and
  `amap` are now declared with explicit `importFrom` directives.
* `registerDoParallel()` is now called with its namespace prefix.
* The `Author:` field has been removed from `DESCRIPTION` in favor of a
  single `Authors@R:` declaration. The `Date:` field has been removed.
* A copyright-holder role (`cph`) has been added for the maintainer.
* Examples no longer call `library()` on the package itself or on
  imported packages.

## Code quality

* The single 1100-line `R/main.R` has been split into thematic files:
  `population.R`, `groups.R`, `genetic-ops.R`, `objectives.R`,
  `pareto-local-search.R`, `path-relinking.R`, `results.R` and
  `moc-gapbk.R`.
* Fixed a latent bug in `singletons.repair()` which returned a list
  with duplicated element names. The list now uses distinct,
  meaningful names; callers were updated accordingly.
* Added input validation with informative error messages at the
  beginning of `moc.gapbk()`.
* `parallel::makeCluster()` now uses `on.exit(parallel::stopCluster())`
  to guarantee cluster shutdown even if the parallel section fails.

## Documentation and tests

* Added `tests/testthat/` with smoke tests for `moc.gapbk()` and the
  deprecated alias.
* Added a `vignettes/moc-gapbk-intro.Rmd` introductory vignette.
* Added `README.md`, `NEWS.md` and `cran-comments.md`.
* Minimum R version raised from 3.2.5 to 3.5.0.

 
