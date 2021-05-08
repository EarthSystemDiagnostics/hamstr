hamstr: Hierarchical Accumulation Modelling with Stan and R.
================
Andrew M. Dolman
2021-05-08

------------------------------------------------------------------------

**hamstr** implements a *Bacon-like* (Blaauw and Christen, 2011)
sediment accumulation or age-depth model with hierarchically structured
multi-resolution sediment sections. The Bayesian model is implemented in
the Stan probabilistic programming language (<https://mc-stan.org/>).

We aim to develop this package following the [Guidelines for R packages
providing interfaces to
Stan](https://cran.r-project.org/web/packages/rstantools/vignettes/developer-guidelines.html)

### Installation

**hamstr** can be installed directly from Github

``` r
if (!require("remotes")) {
  install.packages("remotes")
}

remotes::install_github("earthsystemdiagnostics/hamstr", args = "--preclean", build_vignettes = FALSE)
```

### Using hamstr

Examples using the example core “MSB2K” from the
[rbacon](https://cran.r-project.org/web/packages/rbacon/index.html)
package.

``` r
library(hamstr)
library(rstan)
library(tidyverse)
#> Warning: package 'tibble' was built under R version 4.0.3
#> Warning: package 'readr' was built under R version 4.0.3

set.seed(20200827)
```

#### Converting radiocarbon ages to calendar ages.

Unlike Bacon, hamstr does not do the conversion of radiocarbon dates to
calendar ages as part of the model fitting process. This must be done in
advance. hamstr includes the helper function `calibrate_14C_age` to do
this, which in turn uses the function `BchronCalibrate` from the
[Bchron](https://cran.r-project.org/web/packages/Bchron/index.html)
package.

Additionally, unlike Bacon, hamstr summarises the complex empirical
calendar age PDF that results from calibration into a single point
estimate and 1-sigma uncertainty. This is a necessary compromise in
order to be able to use the power of the Stan platform. Viewed in
context with the many other uncertainties in radiocarbon dates and the
resulting age-models this will not usually be a major issue.

The function `calibrate_14C_age` will append columns to a data.frame
with the calendar ages and 1-sigma uncertainties.

``` r
MSB2K_cal <- calibrate_14C_age(MSB2K, age.14C = "age", age.14C.se = "error")
```

#### Fitting age-models with hamstr

Age-depth (sediment accumulation) models are fit with the function
`hamstr`.

By default hamstr runs with three Markov chains and these can be run in
parallel. This code will assign 3 (processor) cores as long as the
machine has at least 3.

``` r
if (parallel::detectCores() >= 3) options(mc.cores = 3)
```

The key structural difference between Bacon and hamstr models is that
with hamstr the sediment profile is modelled at multiple resolutions
with a hierarchical structure. This removes the need to trade-off
smoothness and flexibility. One of the more critical tuning parameters
in the Bacon model is the parameter `thick`, which determines the
thickness and number of discrete core sections. Finding a good or
optimal value for a given core is often critical to getting a good
age-depth model. Too few sections and the resulting age-model is very
“blocky” and can miss changes in sedimentation rate; however,
counter-intuitively, too many very thin sections can also often result
in an age-model that “under-fits” the data - ploughing a straight line
through the age-control points when a lower resolution model shows
variation in accumulation rate.

The parameter `K` controls the number and structure of the hierarchical
sections. It is specified as a vector, where each value indicates the
number of new child sections for each parent section at each finer
hierarchical level. E.g. `c(10, 10)` specifies 10 sections at the
coarsest level, with 10 new sections for each coarse section at the next
finer level, for a total of 100 sections at the highest / finest
resolution level. `c(10, 10, 10)` would specify 1000 sections at the
finest level and 3 hierarchical levels of 10, 100 and 1000 sections.

The structure is hierarchical in the sense that the accumulation rates
for the parent sections act as priors for their child sections,
specifically the mean accumulation rate for a given parent is the mean
of the gamma prior for it’s child sections. In turn, the overall mean
accumulation rate for the whole core is itself a parameter estimated by
the fitting process. The hierarchical structure of increasing resolution
allows the model to adapt to low-frequency changes in the accumulation
rate, that is changes between “regimes” of high or low accumulation that
persist for long periods, perhaps due to glacial-interglacial cycles.

By default, hamstr uses robust linear regression (`MASS::rlm`) to
estimate the mean accumulation rate from the data, and then uses this to
parametrise a prior distribution for the overall mean accumulation rate.
This prior is a half-normal with zero mean and standard deviation equal
to 10 times the estimated mean. Although this does introduce an element
of “data-peaking” or using the data twice, for both the prior and
likelihood, the resulting prior is only weakly-informative and unless
the data set is very poor, the data exert very strong control over this
parameter anyway. The advantage is that the prior is automatically
scaled appropriately regardless of the units of depth or age.

Here we fit a model with a 10, 100 section structure.

``` r
hamstr_fit_1 <- hamstr(depth = MSB2K_cal$depth,
                   obs_age = MSB2K_cal$age.14C.cal,
                   obs_err = MSB2K_cal$age.14C.cal.se,
                   K = c(10, 10))
```

The default plotting method shows the fitted age models together with
some diagnostic plots: a traceplot of the log-posterior to assess
convergence of the overall model; a plot of accumulation rate against
depth at each hierarchical level; the prior and posterior of the memory
parameter. By default the age-models are summarised to show the mean,
median, 25% and 95% posterior intervals. The data are shown as points
with their 1-sigma uncertainties. The structure of the sections is shown
along the top of the age-model plot.

``` r
plot(hamstr_fit_1)
#> Joining, by = "idx"
#> Joining, by = "alpha_idx"
```

![](readme_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

A “spaghetti” plot can be created instead of shaded regions. This shows
a random sample of iterations from the posterior distribution
(realisation of the age-depth model). This can be slow if lots of
iterations are plotted, the default is to plot 1000 iterations.
Additionally, plotting of the diagnostic plots can be switched off.

``` r
plot(hamstr_fit_1, summarise = FALSE, plot_diagnostics = FALSE)
#> Joining, by = "idx"
```

![](readme_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Additional diagnostic plots are available. A plot of the prior and
posterior distributions of the overall mean accumulation rate for
example. See ?plot.hamstr\_fit for options.

#### Prior and posterior of overall mean accumulation rate

``` r
plot(hamstr_fit_1, type = "acc_mean_prior_post")
```

![](readme_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

#### Plot modelled accumulation rates at each hierarchical level

``` r
plot(hamstr_fit_1, type = "hier_acc")
#> Joining, by = "alpha_idx"
```

![](readme_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

#### Plot memory prior and posterior

As for this example the highest resolution sections are approximately 1
cm thick, there is not much difference between R and w.

``` r
plot(hamstr_fit_1, type = "mem")
```

![](readme_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

#### Examine calendar age PDFs

As mentioned earlier, unlike Bacon, hamstr approximates the complicated
empirical calendar age PDFs with t-distributions. These can be examined
with the function `compare_14C_PDF`

A sample of six dates are plotted here.

``` r
i <- seq(1, 40, by = floor(40/6))[1:6]
compare_14C_PDF(MSB2K$age[i], MSB2K$error[i], cal_curve = "intcal20")
```

![](readme_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

### Varying the number and structure of the modelled sections.

To get to approximately 100 sections at the finest level, we could have
a single level `c(100)` or at the other extreme c(2,2,2,2,2,2,2), which
would result in 2^7 = 128 sections.

There is a trade-off between these two approaches. Having the maximum
number of hierarchical levels maximises the ultimate flexibility of the
model, but comes at a cost of computation time and model stability. More
hierarchical levels means that more parameters are ultimately being
estimated to end up with the same number at the finest level. For
example, `c(10, 10)` results in 110 fine sections but 110 accumulation
rates are estimated. `c(2, 2, 2, 2, 2, 2, 2)` results in 128 sections
but 254 parameters are estimated. However, more importantly, having
fewer child sections means that less information is being used “pooled”
to estimate the parent section mean. Each parameter is estimated with
less certainty and therefore is a less informative prior for its child
sections. This makes sampling the posterior distribution more difficult,
slower, and less stable. The optimal structure will likely depend on the
given dataset and optimising this is under investigation, but
hierarchical levels of about 10 seem to work well.

Here we compare a 2^7 structure with a single level of 128 sections.

``` r
hamstr_fit_2pow7 <- hamstr(depth = MSB2K_cal$depth,
                   obs_age = MSB2K_cal$age.14C.cal,
                   obs_err = MSB2K_cal$age.14C.cal.se,
                   K = c(2, 2, 2, 2, 2, 2, 2))
#> Warning: There were 20 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#> Warning: Examine the pairs() plot to diagnose sampling problems

hamstr_fit_5pow3 <- hamstr(depth = MSB2K_cal$depth,
                   obs_age = MSB2K_cal$age.14C.cal,
                   obs_err = MSB2K_cal$age.14C.cal.se,
                   K = c(5,5,5))

hamstr_fit_128 <- hamstr(depth = MSB2K_cal$depth,
                   obs_age = MSB2K_cal$age.14C.cal,
                   obs_err = MSB2K_cal$age.14C.cal.se,
                   K = c(128))
```

``` r
rstan::get_elapsed_time(hamstr_fit_2pow7$fit)
#>         warmup sample
#> chain:1 12.887 28.445
#> chain:2 14.395 26.182
#> chain:3 13.857 11.897
rstan::get_elapsed_time(hamstr_fit_5pow3$fit)
#>         warmup sample
#> chain:1  4.518  2.899
#> chain:2  4.483  3.143
#> chain:3  4.377  2.893
rstan::get_elapsed_time(hamstr_fit_128$fit)
#>         warmup sample
#> chain:1  5.678  5.591
#> chain:2  5.408  5.388
#> chain:3  5.432  5.343
```

``` r
plot(hamstr_fit_2pow7)
#> Joining, by = "idx"
#> Joining, by = "alpha_idx"
```

![](readme_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
plot(hamstr_fit_5pow3)
#> Joining, by = "idx"
#> Joining, by = "alpha_idx"
```

![](readme_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

``` r
plot(hamstr_fit_128)
#> Joining, by = "idx"
#> Joining, by = "alpha_idx"
```

![](readme_files/figure-gfm/unnamed-chunk-14-3.png)<!-- -->

#### Refit with 21 non-hierarchical sections for comparison with figures in Blaauw and Christen (2011)

``` r
hamstr_fit_21 <- hamstr(depth = MSB2K_cal$depth,
                   obs_age = MSB2K_cal$age.14C.cal,
                   obs_err = MSB2K_cal$age.14C.cal.se,
                   K = c(21))
```

``` r
plot(hamstr_fit_21, summarise = FALSE)
#> Joining, by = "idx"
#> Joining, by = "alpha_idx"
```

![](readme_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

### Convenience functions

#### Extracting the age models

`get_posterior_ages` extracts all iterations (realisations) of the
age-depth model

``` r
age.mods <- get_posterior_ages(hamstr_fit_5pow3)
#> Joining, by = "idx"
age.mods
#> # A tibble: 378,000 x 5
#>     iter par      age   idx depth
#>    <int> <chr>  <dbl> <dbl> <dbl>
#>  1     1 c_ages 4403.     1  1.5 
#>  2     1 c_ages 4409.     2  2.28
#>  3     1 c_ages 4424.     3  3.07
#>  4     1 c_ages 4440.     4  3.85
#>  5     1 c_ages 4458.     5  4.64
#>  6     1 c_ages 4472.     6  5.42
#>  7     1 c_ages 4484.     7  6.20
#>  8     1 c_ages 4492.     8  6.99
#>  9     1 c_ages 4503.     9  7.77
#> 10     1 c_ages 4515.    10  8.56
#> # ... with 377,990 more rows
```

The hierarchical structure of the sections makes it difficult to specify
the exact depth resolution that you want for your resulting age-depth
model. The function `interpolate_age_models` makes it easier to
interpolate to a specific set of depths. The function returns NA for
depths that are outside the modelled depths.

``` r
age.mods.interp <- predict(hamstr_fit_5pow3,
                           new_depth = seq(0, 100, by = 1))
#> Joining, by = "idx"
```

These interpolated age models can summarised with the same function as
the original fitted objects, but the n\_eff and Rhat information is
lost.

``` r
summary(hamstr_fit_5pow3)
#> Joining, by = "idx"
#> # A tibble: 126 x 13
#>    depth   idx par    mean se_mean    sd `2.5%` `25%` `50%` `75%` `97.5%` n_eff
#>    <dbl> <dbl> <chr> <dbl>   <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>
#>  1  1.5      1 c_ag~ 4510.    3.64  73.0  4350. 4466. 4515. 4561.   4640.  402.
#>  2  2.28     2 c_ag~ 4522.    3.41  69.1  4370. 4479. 4525. 4570.   4645.  411.
#>  3  3.07     3 c_ag~ 4533.    3.19  65.6  4391. 4493. 4536. 4579.   4653.  422.
#>  4  3.85     4 c_ag~ 4545.    2.98  62.4  4412. 4507. 4548. 4588.   4661.  438.
#>  5  4.64     5 c_ag~ 4557.    2.77  59.6  4428. 4520. 4559. 4597.   4668.  461.
#>  6  5.42     6 c_ag~ 4568.    2.56  57.2  4447. 4532. 4570. 4607.   4675.  498.
#>  7  6.20     7 c_ag~ 4580.    2.36  54.7  4463. 4545. 4582. 4618.   4682.  537.
#>  8  6.99     8 c_ag~ 4591.    2.17  52.2  4482. 4559. 4593. 4628.   4689.  579.
#>  9  7.77     9 c_ag~ 4603.    1.97  49.8  4499. 4571. 4604. 4637.   4697.  637.
#> 10  8.56    10 c_ag~ 4614.    1.79  47.7  4515. 4583. 4615. 4647.   4706.  707.
#> # ... with 116 more rows, and 1 more variable: Rhat <dbl>
```

``` r
summary(age.mods.interp)
#> # A tibble: 101 x 8
#>    depth  mean    sd `2.5%` `25%` `50%` `75%` `97.5%`
#>    <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl>   <dbl>
#>  1     0   NA   NA      NA    NA    NA    NA      NA 
#>  2     1   NA   NA      NA    NA    NA    NA      NA 
#>  3     2 4518.  70.5  4362. 4474. 4521. 4566.   4644.
#>  4     3 4532.  65.9  4390. 4492. 4535. 4578.   4652.
#>  5     4 4547.  61.8  4417. 4510. 4550. 4590.   4662.
#>  6     5 4562.  58.4  4437. 4526. 4565. 4602.   4671.
#>  7     6 4577.  55.3  4458. 4541. 4579. 4615.   4680.
#>  8     7 4592.  52.1  4482. 4559. 4593. 4628.   4690.
#>  9     8 4606.  49.1  4504. 4575. 4608. 4640.   4700.
#> 10     9 4621.  46.7  4522. 4591. 4622. 4653.   4710.
#> # ... with 91 more rows
```

Within the hamstr\_fit object is an rstan object on which all the
standard rstan functions should operate correctly.

``` r
rstan::check_divergences(hamstr_fit_5pow3$fit)
#> 0 of 3000 iterations ended with a divergence.

rstan::stan_rhat(hamstr_fit_5pow3$fit)
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
#> Warning: Removed 3 rows containing non-finite values (stat_bin).
```

![](readme_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

The first `alpha` parameter is the overall mean accumulation rate.

``` r
rstan::traceplot(hamstr_fit_5pow3$fit, par = c("alpha[1]"),
                 inc_warmup = TRUE)
```

![](readme_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

`alpha[2]` to `alpha[6]` will be the first hierarchical level of 5 rates
for this model.

``` r
plot(hamstr_fit_5pow3, type = "hier")
#> Joining, by = "alpha_idx"
```

![](readme_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
rstan::traceplot(hamstr_fit_5pow3$fit, par = paste0("alpha[", 2:6, "]"),
                 inc_warmup = TRUE)
```

![](readme_files/figure-gfm/unnamed-chunk-23-2.png)<!-- -->

### References

-   Blaauw, Maarten, and J. Andrés Christen. 2011. Flexible Paleoclimate
    Age-Depth Models Using an Autoregressive Gamma Process. Bayesian
    Analysis 6 (3): 457-74. <doi:10.1214/ba/1339616472>.

-   Parnell, Andrew. 2016. Bchron: Radiocarbon Dating, Age-Depth
    Modelling, Relative Sea Level Rate Estimation, and Non-Parametric
    Phase Modelling. R package version 4.2.6.
    <https://CRAN.R-project.org/package=Bchron>

-   Stan Development Team (2020). RStan: the R interface to Stan. R
    package version 2.21.2. <http://mc-stan.org/>.
