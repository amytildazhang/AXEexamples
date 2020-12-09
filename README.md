# AXE: (A)pproximate (X)Cross-validated (E)stimates for Bayesian hierarchical regression models

This R package provides all code and data to reproduce the examples from our paper

Zhang, A. X., Bao, L., & Daniels, M. J. (2020). Approximate Cross-validated Mean Estimates for Bayesian Hierarchical Regression Models. [_arXiv preprint arXiv:2011.14238_](https://arxiv.org/abs/2011.14238).

In short, AXE is a fast and accurate method for approximating <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/091cb68641254b4fe2259ab2b767b445.svg?invert_in_darkmode" align=middle width=69.82668pt height=24.56552999999997pt/>, where <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/59645ee1d242ecd5f3c6e02fdfafc4f0.svg?invert_in_darkmode" align=middle width=15.589530000000002pt height=22.381919999999983pt/> is the vector of test data, and <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/e70dee896df24b6ac8b4b9c6ecbce0c7.svg?invert_in_darkmode" align=middle width=25.82514pt height=22.381919999999983pt/> a vector of training data. 



**AXE method**: Let <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/0edb5f964ed489ea149b4ed86e4cf17a.svg?invert_in_darkmode" align=middle width=60.38505pt height=27.598230000000008pt/> denote a continuous response vector that follows

<p align="center"><img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/2d6863afb385d0d8e213533d65d28557.svg?invert_in_darkmode" align=middle width=641.9060999999999pt height=18.269295pt/></p>

where 

- <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/b804e5c4d96a62a9bda96e61babb8291.svg?invert_in_darkmode" align=middle width=87.40264499999999pt height=27.598230000000008pt/> is positive-definite and typically a diagonal matrix,

- <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/631eae5077367ce1492d07680308f52a.svg?invert_in_darkmode" align=middle width=86.34978000000001pt height=27.598230000000008pt/> is positive-definite, 

- <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/3aebf55562d5d25eb561e37600044391.svg?invert_in_darkmode" align=middle width=54.568470000000005pt height=22.381919999999983pt/>

- <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/81356f1013f18f7a294e948fdd54fece.svg?invert_in_darkmode" align=middle width=180.965895pt height=27.940769999999983pt/> is the design matrix,

- <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/fc597e868bfa2bfc6a26041e84328ed7.svg?invert_in_darkmode" align=middle width=66.17027999999999pt height=27.598230000000008pt/> denote fixed effects, <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/d3af09bf0a49960aad5b866df1c3cabb.svg?invert_in_darkmode" align=middle width=66.17027999999999pt height=27.598230000000008pt/> random effects s.t. <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/e5f466966fecdb34ff036a2ef1fc2cde.svg?invert_in_darkmode" align=middle width=150.36615pt height=32.55549000000001pt/>, 

- WLOG, <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/4de3ab8dcd6b8f58e97e2eccdfcfa0ff.svg?invert_in_darkmode" align=middle width=92.179065pt height=21.10812pt/> 


AXE approximates cross-validated mean estimates using the posterior means for <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/d4261132636819ae7c6f4039dafc4016.svg?invert_in_darkmode" align=middle width=11.827860000000003pt height=31.056300000000004pt/> and <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/bde65f3adefc3c20e7dfdac5b016f5f9.svg?invert_in_darkmode" align=middle width=9.058995000000001pt height=22.745910000000016pt/> as plug-in estimates. Let us refer to <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.681657500000003pt height=21.602129999999985pt/> as the test data indices, so <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/59645ee1d242ecd5f3c6e02fdfafc4f0.svg?invert_in_darkmode" align=middle width=15.589530000000002pt height=22.381919999999983pt/> denotes the vector of test data and <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/e70dee896df24b6ac8b4b9c6ecbce0c7.svg?invert_in_darkmode" align=middle width=25.82514pt height=22.381919999999983pt/> the vector of training data. The AXE approximation for CV mean estimate <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/091cb68641254b4fe2259ab2b767b445.svg?invert_in_darkmode" align=middle width=69.82668pt height=24.56552999999997pt/> is
<p align="center"><img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/3122a703d6ccfe15fe4863dd1be0f56a.svg?invert_in_darkmode" align=middle width=533.511pt height=52.457625pt/></p>

The basic reasoning is that  <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/7f330ed5b7dd756c1109185343b82a48.svg?invert_in_darkmode" align=middle width=297.665445pt height=31.056300000000004pt/>, where <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/4d44968bfd84067b99cd11f3fb43fdcf.svg?invert_in_darkmode" align=middle width=63.107715pt height=31.056300000000004pt/> denote the Empirical Bayes estimates [Kass and Steffey, 1989]. That is, the conditional posterior given variance parameters can approximate the posterior mean for <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode" align=middle width=10.127700000000003pt height=22.745910000000016pt/> when <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/bda33f0358442dd75d7487fa0ba0a279.svg?invert_in_darkmode" align=middle width=17.042355pt height=22.381919999999983pt/> is large enough (whether <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/bda33f0358442dd75d7487fa0ba0a279.svg?invert_in_darkmode" align=middle width=17.042355pt height=22.381919999999983pt/> is large enough can be determined by deriving <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/412226af2aad4e32ff0a61874d8ecb6a.svg?invert_in_darkmode" align=middle width=100.33023pt height=31.056300000000004pt/> and comparing to the posterior mean estimates). Then, so long as the posterior means <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/d4261132636819ae7c6f4039dafc4016.svg?invert_in_darkmode" align=middle width=11.827860000000003pt height=31.056300000000004pt/> and <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/bde65f3adefc3c20e7dfdac5b016f5f9.svg?invert_in_darkmode" align=middle width=9.058995000000001pt height=22.745910000000016pt/> are stable enough across cross-validation folds, they can be used as plug-in estimates in the conditional mean in <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/48d23936f8dd5befe2652a6297bbb7f4.svg?invert_in_darkmode" align=middle width=113.87772pt height=31.056300000000004pt/> to produce approximations of   <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/3bb15193175fa831416f8dac78754980.svg?invert_in_darkmode" align=middle width=78.381105pt height=24.56552999999997pt/>. 

By conditioning on the variance-covariance parameters, we shift the CV problem from probability-based sampling to the same form as maximum likelihood methods for simple linear regression} and is likewise <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/2c18bba38ccecb24edb4382fda296a83.svg?invert_in_darkmode" align=middle width=142.167135pt height=27.940769999999983pt/> in time for each CV fold. It can be used with any CV schema, e.g. K-fold, leave-one-out (LOO), and leave-one-cluster-out (LCO). Our paper focuses on LCO-CV.

**LCO-CV:** LCO-CV is commonly used in models with complex dependency structures, such as spatio-temporal models or models with repeated measures.  In such models, the dependency among the data can cause <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/d6328eaebbcd5c358f426dbea4bdbf70.svg?invert_in_darkmode" align=middle width=15.080505pt height=22.381919999999983pt/>-fold CV to select models which overfit [@arlot2010survey, @opsomer2001nonparametric]. Reducing the amount of correlation between the training and test data typically results in non-random, structured CV, of which LCO-CV is one commoon example. Under LCO-CV, the data are partitioned based on the unique values of one or more random intercepts, e.g. for a mixed-effects model with repeated measures, all repeated measures for a test unit are set aside from the training data. This provides a more realistic estimate of a model's predictive capability for a new unit. 

In general, we have found that AXE improves upon existing LCO-CV methods in accuracy and is typically an order of magnitude faster. 


# Generalized Linear Mixed Models (GLMMs)

Essentially, we use expectation propagation with a Gaussian approximating density, which allows us to use the same equations as the LMM case, just with the transformed response and variance parameters. 

Take as example a one-way GLMM with clustered response data <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/81cd5d7f4642dc17bd1cbd70e9b10b9d.svg?invert_in_darkmode" align=middle width=15.589530000000002pt height=22.381919999999983pt/>, where <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.681657500000003pt height=21.602129999999985pt/> denotes the <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/95291b39ba5d9dba052b40bf07b12cd2.svg?invert_in_darkmode" align=middle width=20.29962pt height=27.852989999999977pt/> cluster, <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/54158e2c605c3ecf783cdc13e7235676.svg?invert_in_darkmode" align=middle width=15.911775pt height=14.102549999999994pt/> the size of the cluster, and <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/18782f74837a5ea1c78b0477c82e634c.svg?invert_in_darkmode" align=middle width=88.37928pt height=21.10812pt/> indexes values within the <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/95291b39ba5d9dba052b40bf07b12cd2.svg?invert_in_darkmode" align=middle width=20.29962pt height=27.852989999999977pt/> cluster. The data <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/5915d85d514588f02f78c0d62a42c632.svg?invert_in_darkmode" align=middle width=20.5392pt height=22.381919999999983pt/> have some probability density function <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/f30fdded685c83b0e7b446aa9c9aa120.svg?invert_in_darkmode" align=middle width=9.922935000000003pt height=14.102549999999994pt/> and are modeled as a regression through link function <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode" align=middle width=8.398995000000005pt height=14.102549999999994pt/> so that
<p align="center"><img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/73a24cd48fb6548bc1a720d2794039cd.svg?invert_in_darkmode" align=middle width=235.52100000000002pt height=42.926895pt/></p>
where <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/391778de429bd69950e52d640c84082d.svg?invert_in_darkmode" align=middle width=169.99174499999998pt height=26.70657pt/> and <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/8217d9e1e73d921529eaa09fd2639ae8.svg?invert_in_darkmode" align=middle width=249.30229500000002pt height=26.70657pt/>. Then taking the normal approximation with equivalent moments converges as <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/a8f31f4eaabba25bc3f6b9603cbdc63f.svg?invert_in_darkmode" align=middle width=15.911775pt height=14.102549999999994pt/> becomes large:
<p align="center"><img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/44053ec472a47caecfc622d32bb8ed22.svg?invert_in_darkmode" align=middle width=531.6134999999999pt height=66.47487pt/></p>

Denote in \ereft{eqn:glmm_approx} the transformed response <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/44849e8252c5150da94f717b14de782f.svg?invert_in_darkmode" align=middle width=42.57429pt height=24.56552999999997pt/> as <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/cabd60a15cecb9c03cd7992e8de12089.svg?invert_in_darkmode" align=middle width=20.5392pt height=25.670699999999986pt/>, the variance as <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/fb95b7a172e28cc065861be926a2b187.svg?invert_in_darkmode" align=middle width=15.548115000000001pt height=26.70657pt/>, the diagonal matrix of <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/fb95b7a172e28cc065861be926a2b187.svg?invert_in_darkmode" align=middle width=15.548115000000001pt height=26.70657pt/>'s as <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/28661a9e9fcbc4ed7e83a5b94d96f5a4.svg?invert_in_darkmode" align=middle width=35.665905pt height=24.56552999999997pt/>, and the normal approximation as <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/a985530bd7b6ca7c39b5c776712987ba.svg?invert_in_darkmode" align=middle width=132.7524pt height=25.670699999999986pt/>. When <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/813cd865c037c89fcdc609b25c465a05.svg?invert_in_darkmode" align=middle width=11.827860000000003pt height=22.381919999999983pt/> and <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode" align=middle width=10.127700000000003pt height=22.745910000000016pt/> are unknown, we plug in the marginal posterior mean <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/5458790cb780f379ee7fdeb47d5c72c2.svg?invert_in_darkmode" align=middle width=85.41489pt height=31.056300000000004pt/> and <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/a71e9245f4eabee2f68c1ce37c9cada3.svg?invert_in_darkmode" align=middle width=82.01127pt height=31.421609999999983pt/>. Then,
<p align="center"><img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/5585182c754f786a889bfdce13ccbcda.svg?invert_in_darkmode" align=middle width=377.4309pt height=99.19337999999999pt/></p>


# Examples

![Figure 2 from paper](data-raw/p_both.png)
Figure 2 from the paper gives point-by-point comparisons of ground truth manual cross-validation (x-axis) against AXE approximations (panel A, y-axis). Point-by-point comparisons for other LCO methods iIS-C and GHOST are also included. (iIS-A is omitted to preserve the scale of the axes; results summarized in Figure 1 of paper.)

# Code
The R package here bundles together all code used to produce the examples in the paper. To use, install from github 

```{r}
library(devtools)
install_github("amytildazhang/AXEexamples")

```

Alternatively, download and load/source the appropriate files, i.e. in R do:

```{r}
for (file in list.files("R")) {
   source(file.path("R", file))
}  # imports all functions for re-running MCV, AXE, or other LCO approximations

for (obj in list.files("data")) {
   load(file.path("data", obj))
} # loads pre-obtained MCV, AXE, and LCO values

```

The main functions in the package are

- `prep_*()`: Creates a list with data for each example
   + `prep_eight()`: Eight schools (LMM)
   + `prep_radon_full()`: Radon (LMM)
   + `prep_radon_simul()`: Radon subsets (LMM)
   + `prep_lol()`: Esports players (GLMM)
   + `prep_slc()`: Scottish lip cancer (GLMM, CAR)
   + `prep_air()`: Scottish respiratory disease (GLMM, spatio-temporal CAR)
- `pfit_*()`: Fits model to full data and produces LCO approximations, based on posterior samples, for iIS-A, iIS-C, and GHOST.
- `axe_*()`: Uses co/variance posterior means to produce AXE estimates.
- `mcv_*()`: Runs manual cross-validation and saves <img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/091cb68641254b4fe2259ab2b767b445.svg?invert_in_darkmode" align=middle width=69.82668pt height=24.56552999999997pt/>.


Example results are obtained by calling the above functions for each example in turn, e.g. the Eight schools results are generated using the following code:
```{r eightex, eval=FALSE}
eight <- prep_eight() 
eight<img src="https://rawgit.com/amytildazhang/AXEexamples (fetch/None/svgs/742add498f9675f93dd62958b12be573.svg?invert_in_darkmode" align=middle width=204.829845pt height=24.56552999999997pt/>posteriors <- pfit_eight() 
eight$axe_yhats <- axe_eight() 

```

The object `eight` is available in the package, listed in the `data` folder. Each of the paper's examples is saved in  `data` under the following names:
- `eight`: Eight schools
- `radon_1`: Radon
- `radon_2`: Radon subsets
- `lol`: Esports players
- `slc`: Scottish lip cancer
- `air`: Scottish respiratory disease
Access simply by loading the package (`library(AXEexamples)`) and calling the name of the data, e.g. `eight; str(eight)`.

For code to produce paper figures, see the vignette in `doc/overview.html`


# References

Robert E Kass and Duane Steffey. Approximate Bayesian inference in conditionally independent hierarchical models (parametric empirical Bayes models). Journal of the American Statistical Association, 84(407):717–726, 1989.




