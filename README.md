README
================
Amy X. Zhang
12/8/2020

# AXE: (A)pproximate (X)Cross-validated (E)stimates for Bayesian hierarchical regression models

This R package provides code and data to reproduce all examples from the
paper

Zhang, A. X., Bao, L., & Daniels, M. J. (2020). Approximate
Cross-validated Mean Estimates for Bayesian Hierarchical Regression
Models. [*arXiv preprint
arXiv:2011.14238*](https://arxiv.org/abs/2011.14238).

The paper presents a novel method for approximating \(E[Y_j | Y_{-j}]\),
the cross-validated posterior mean for vector of test data \(Y_j\) given
the vector of training data \(Y_{-j}\). The method applies to any CV
scheme; we compared AXE specifically to existing CV approximations
methods under the challenging leave-a-cluster-out CV scheme. In general,
we have found that AXE improves upon existing LCO-CV methods in accuracy
and is typically an order of magnitude faster.

![Figure 2 from paper](data-raw/p_both.png) Figure 2 from the paper
gives point-by-point comparisons of ground truth manual cross-validation
(x-axis) against AXE approximations (panel A, y-axis). In general, AXE approximations lie on or near the 45-degree one-to-one line.
Point-by-point comparisons for other LCO methods integrated importance sampling (iIS-C)
and ghosting (GHOST) are also included, in panel B. (iIS-A is omitted to preserve
the scale of the axes; results summarized in Figure 1 of paper.) For other LCO methods,
the variance of points around the 45-degree line is higher than for AXE.

For a description of how to obtain AXE approximations, see the “Overview” vignette
(`vignettes("overview")`).

# Code description

The R package here bundles together all code used to produce the
examples in the paper. To use, install from github with the following
code:

``` r
library(devtools)
install_github("amytildazhang/AXEexamples")
```

Note that this installs all source code for producing the examples,
including compiled STAN models which can take time and may not be of
interest (and, until the next version of `rstan` hits CRAN, compiling
STAN models within R packages on Windows can be a frustrating endeavor).
Alternatively, manually download and load/source the appropriate files,
i.e. in R do:

``` r
for (file in list.files("R")) {
   source(file.path("R", file))
}  # imports all functions for re-running MCV, AXE, or other LCO approximations

for (obj in list.files("data")) {
   load(file.path("data", obj))
} # loads pre-obtained MCV, AXE, and LCO values
```

Functions that call `rstan::sampling(stanmodels$[model name], ...)` will
need to be replaced with `rstan::stan([filepath], ...)` prior to
sourcing the files.

The main functions in the package are

  - `prep_*()`: Creates a list with data for each example
      - `prep_eight()`: Eight schools (LMM)
      - `prep_radon_full()`: Radon (LMM)
      - `prep_radon_simul()`: Radon subsets (LMM)
      - `prep_lol()`: Esports players (GLMM)
      - `prep_slc()`: Scottish lip cancer (GLMM, CAR)
      - `prep_air()`: Scottish respiratory disease (GLMM,
        spatio-temporal CAR)
  - `pfit_*()`: Fits model to full data and produces LCO approximations,
    based on posterior samples, for iIS-A, iIS-C, and GHOST.
  - `axe_*()`: Uses co/variance posterior means to produce AXE
    estimates.
  - `mcv_*()`: Runs manual cross-validation and saves
    \(E[Y_j | Y_{-j}]\).

Example results are obtained by calling the above functions for each
example in turn, e.g. the Eight schools results are generated using the
following code:

``` r
eight <- prep_eight() 
eight$cv_yhats  <- mcv_eight() 
eight$posteriors <- pfit_eight() 
eight$axe_yhats <- axe_eight() 
```

The object `eight` is available in the package, listed in the `data`
folder. Each of the paper’s examples is saved in `data` under the
following names: - `eight`: Eight schools - `radon_1`: Radon -
`radon_2`: Radon subsets - `lol`: Esports players - `slc`: Scottish lip
cancer - `air`: Scottish respiratory disease Access by loading the
package (`library(AXEexamples)`) and calling the name of the data, e.g.
`eight; str(eight)`.

For code to reproduce paper figures, see the “Overview” vignette
(`vignettes("overview")`).
