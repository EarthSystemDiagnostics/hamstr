hamstr: Hierarchical Accumulation Modelling with Stan and R.
================
Andrew M. Dolman
2021-10-04

------------------------------------------------------------------------

**hamstr** implements a *Bacon-like* (Blaauw and Christen, 2011)
sediment accumulation or age-depth model with hierarchically structured
multi-resolution sediment sections. The Bayesian model is implemented in
the Stan probabilistic programming language (<https://mc-stan.org/>).

## Installation

**hamstr** can be installed directly from Github

``` r
if (!require("remotes")) {
  install.packages("remotes")
}

remotes::install_github("earthsystemdiagnostics/hamstr", args = "--preclean", build_vignettes = FALSE)
```

## Using **hamstr**

Examples using the example core “MSB2K” from the
[rbacon](https://cran.r-project.org/web/packages/rbacon/index.html)
package.

``` r
library(hamstr)
library(rstan)
library(tidyverse)

set.seed(20200827)
```

### Converting radiocarbon ages to calendar ages.

Unlike Bacon, **hamstr** does not do the conversion of radiocarbon dates
to calendar ages as part of the model fitting process. This must be done
in advance. **hamstr** includes the helper function `calibrate_14C_age`
to do this, which in turn uses the function `BchronCalibrate` from the
[Bchron](https://cran.r-project.org/web/packages/Bchron/index.html)
package.

Additionally, unlike Bacon, **hamstr** approximates the complex
empirical calendar age PDF that results from calibration into a single
point estimate and 1-sigma uncertainty. This is a necessary compromise
in order to be able to use the power of the Stan platform. Viewed in
context with the many other uncertainties in radiocarbon dates and the
resulting age-models this will not usually be a major issue.

The function `calibrate_14C_age` will append columns to a data.frame
with the calendar ages and 1-sigma uncertainties.

``` r
MSB2K_cal <- calibrate_14C_age(MSB2K, age.14C = "age", age.14C.se = "error")
```

The approximated calendar age PDFs can be compared with the empirical
PDFs with the function `compare_14C_PDF`

A sample of six dates are plotted here for the IntCal20 and Marine20
calibrations. This approximation is much less of an issue for marine
radiocarbon dates, as the cosmogenic radiocarbon signal has been
smoothed by mixing in the ocean.

``` r
i <- seq(1, 40, by = floor(40/6))[1:6]
compare_14C_PDF(MSB2K$age[i], MSB2K$error[i], cal_curve = "intcal20")+
  labs(title = "Intcal20")
```

![](readme_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
compare_14C_PDF(MSB2K$age[i], MSB2K$error[i], cal_curve = "marine20") +
  labs(title = "Marine20")
```

![](readme_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

### Fitting age-models with **hamstr**

By default **hamstr** runs with three Markov chains and these can be run
in parallel. This code will assign 3 processor cores as long as the
machine has at least 3. The number of cores can also be set for specific
calls of the `hamstr` function using the `cores` argument.

``` r
if (parallel::detectCores() >= 3) options(mc.cores = 3)
```

Age-depth (sediment accumulation) models are fit with the function
`hamstr`. A vectors of depth, observed age and age uncertainty are
passed as arguments to the function.

``` r
hamstr_fit_1 <- hamstr(depth = MSB2K_cal$depth,
                       obs_age = MSB2K_cal$age.14C.cal,
                       obs_err = MSB2K_cal$age.14C.cal.se)
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
#> Joining, by = "alpha_idx"
```

![](readme_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

A “spaghetti” plot can be created instead of shaded regions. This shows
a random sample of iterations from the posterior distribution
(realisation of the age-depth model). This can be slow if lots of
iterations are plotted, the default is to plot 1000 iterations.
Additionally, plotting of the diagnostic plots can be switched off.

``` r
plot(hamstr_fit_1, summarise = FALSE, plot_diagnostics = FALSE)
```

![](readme_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

#### Mean accumulation rate

There is no need to specify a prior value for the mean accumulation rate
(parameter `acc.mean` in Bacon) as in **hamstr**, this overall mean
accumulation rate is a full parameter estimated from the data.

By default, **hamstr** uses robust linear regression (`MASS::rlm`) to
estimate the mean accumulation rate from the data, and then uses this to
parametrise a prior distribution for the overall mean accumulation rate.
This prior is a half-normal with zero mean and standard deviation equal
to 10 times the estimated mean. Although this does introduce a slight
element of “double-dipping”, using the data twice (for both the prior
and likelihood), the resulting prior is only weakly-informative. The
advantage of this approach is that the prior is automatically scaled
appropriately regardless of the units of depth or age.

This prior can be checked visually against the posterior. The posterior
distribution should be much narrower than the weakly informative prior.

``` r
plot(hamstr_fit_1, type = "acc_mean_prior_post")
```

![](readme_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

#### Other hyperparameters

Default parameter values for the shape of the gamma distributed
accumulation rates `acc_shape = 1.5`, the memory mean `mem_mean = 0.5`
and memory strength `mem_strength = 10`, are the same as for Bacon &gt;=
2.5.1.

### Setting the thickness, number, and hierarchical structure of the discrete sections

One of the more critical tuning parameters in the **Bacon** model is the
parameter `thick`, which determines the thickness and number of discrete
down-core sediment sections modelled. Finding a good or optimal value
for a given core is often critical to getting a good age-depth model.
Too few sections and the resulting age-model is very “blocky” and can
miss changes in sedimentation rate; however, counter-intuitively, too
many very thin sections can also often result in an age-model that
“under-fits” the data - a straight line through the age-control points
when a lower resolution model shows variation in accumulation rate.

The key structural difference between **Bacon** and **hamstr** models is
that with **hamstr** the sediment core is modelled at multiple
resolutions simultaneously with a hierarchical structure. This removes
the need to trade-off smoothness and flexibility.

The parameter `K` controls the number and structure of the hierarchical
sections. It is specified as a vector, where each value indicates the
number of new child sections for each parent section at each finer
hierarchical level. E.g. `c(10, 10)` would specify 10 sections at the
coarsest level, with 10 new sections at the next finer level for each
coarse section, giving a total of 100 sections at the highest / finest
resolution level. `c(10, 10, 10)` would specify 1000 sections at the
finest level and 3 hierarchical levels of 10, 100 and 1000 sections.

The structure is hierarchical in the sense that the modelled
accumulation rates for the parent sections act as priors for their child
sections; specifically, the mean accumulation rate for a given parent is
the mean of the gamma prior for it’s child sections. In turn, the
overall mean accumulation rate for the whole core is itself a parameter
estimated by the fitting process. The hierarchical structure of
increasing resolution allows the model to adapt to low-frequency changes
in the accumulation rate, that is changes between “regimes” of high or
low accumulation that persist for long periods.

By default `K` is chosen so that the number of hierarchical levels, and
the number of new child sections per level, are approximately equal,
e.g. c(4, 4, 4, 4). The total number of sections at the finest level is
set so that the resolution is 1 cm per section, up to a total length of
900 cm, above which the default remains 900 sections and a coarser
resolution is used. This can be changed from the default via the
parameter `K`.

For a given shape parameter `acc_shape`, increasing the number of
modelled hierarchical levels increases the total variance in the
accumulation rates at the highest / finest resolution level. From
**hamstr** version 0.5.0 and onwards, the total variance is controlled
by modifying the shape parameter according to the number of hierarchical
levels.

### Getting the fitted age models

The fitted age models can be obtained with the `predict` and `summary`
methods. *iter* is the iteration of the sampler, or “realisation” of the
age model.

``` r
predict(hamstr_fit_1)
#> # A tibble: 243,000 x 3
#>     iter depth   age
#>    <int> <dbl> <dbl>
#>  1     1  1.5  4467.
#>  2     1  2.72 4514.
#>  3     1  3.95 4549.
#>  4     1  5.18 4581.
#>  5     1  6.4  4606.
#>  6     1  7.62 4628.
#>  7     1  8.85 4646.
#>  8     1 10.1  4662.
#>  9     1 11.3  4676.
#> 10     1 12.5  4688.
#> # ... with 242,990 more rows
```

`summary` returns the age model summarised over the realisations.

``` r
summary(hamstr_fit_1)
#> # A tibble: 81 x 13
#>    depth   idx par     mean se_mean    sd `2.5%` `25%` `50%` `75%` `97.5%` n_eff
#>    <dbl> <dbl> <chr>  <dbl>   <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>
#>  1  1.5      1 c_age~ 4512.   1.72   68.2  4362. 4471. 4517. 4562.   4626. 1575.
#>  2  2.72     2 c_age~ 4529.   1.54   62.6  4396. 4490. 4532. 4574.   4636. 1646.
#>  3  3.95     3 c_age~ 4545.   1.37   58.1  4422. 4508. 4549. 4586.   4646. 1789.
#>  4  5.18     4 c_age~ 4562.   1.24   54.5  4448. 4526. 4565. 4600.   4659. 1942.
#>  5  6.4      5 c_age~ 4578.   1.11   51.9  4471. 4544. 4581. 4614.   4672. 2168.
#>  6  7.62     6 c_age~ 4595.   1.04   50.3  4492. 4563. 4598. 4630.   4685. 2322.
#>  7  8.85     7 c_age~ 4612.   0.955  47.5  4515. 4581. 4614. 4645.   4698. 2474.
#>  8 10.1      8 c_age~ 4630.   0.858  44.3  4541. 4601. 4631. 4660.   4713. 2666.
#>  9 11.3      9 c_age~ 4648.   0.807  41.8  4564. 4620. 4649. 4676.   4728. 2686.
#> 10 12.5     10 c_age~ 4666.   0.751  40.3  4587. 4639. 4666. 4693.   4745. 2881.
#> # ... with 71 more rows, and 1 more variable: Rhat <dbl>
```

The hierarchical structure of the sections makes it difficult to specify
the exact depth resolution that you want for your resulting age-depth
model. The `predict` method takes an additional argument `depth` to
interpolate to a specific set of depths. The function returns NA for
depths that are outside the modelled depths.

``` r
age.mods.interp <- predict(hamstr_fit_1, depth = seq(0, 100, by = 1))
```

These interpolated age models can summarised with the same function as
the original fitted objects, but the n\_eff and Rhat information is
lost.

``` r
summary(age.mods.interp)
#> # A tibble: 101 x 8
#>    depth  mean    sd `2.5%` `25%` `50%` `75%` `97.5%`
#>    <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl>   <dbl>
#>  1     0   NA   NA      NA    NA    NA    NA      NA 
#>  2     1   NA   NA      NA    NA    NA    NA      NA 
#>  3     2 4519.  65.7  4375. 4478. 4523. 4566.   4629.
#>  4     3 4532.  61.5  4402. 4494. 4536. 4577.   4639.
#>  5     4 4546.  57.9  4424. 4509. 4549. 4587.   4647.
#>  6     5 4559.  54.9  4445. 4523. 4563. 4598.   4657.
#>  7     6 4573.  52.6  4464. 4538. 4576. 4610.   4668.
#>  8     7 4586.  50.9  4481. 4553. 4589. 4622.   4678.
#>  9     8 4600.  49.3  4499. 4568. 4602. 4634.   4690.
#> 10     9 4614.  47.0  4517. 4584. 4616. 4646.   4700.
#> # ... with 91 more rows
```

### Getting and plotting the accumulation rate

The down-core accumulation rates are returned and plotted in both
depth-per-time, and time-per-depth units. If the input data are in years
and cm then the units will be cm/kyr and yrs/cm respectively. Note that
the acc\_mean parameter in both **hamstr** and Bacon is parametrised in
terms of time per depth.

``` r
plot(hamstr_fit_1, type = "acc_rates")
#> Joining, by = "idx"
```

![](readme_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
summary(hamstr_fit_1, type = "acc_rates") 
#> Joining, by = "idx"
#> # A tibble: 160 x 12
#>    depth c_depth_top c_depth_bottom acc_rate_unit   idx  mean    sd `2.5%` `25%`
#>    <dbl>       <dbl>          <dbl> <chr>         <dbl> <dbl> <dbl>  <dbl> <dbl>
#>  1  1.5         1.5            2.72 depth_per_ti~     1 122.  114.    27.3  57.7
#>  2  2.72        2.72           3.95 depth_per_ti~     2 108.   94.4   29.6  57.8
#>  3  3.95        3.95           5.18 depth_per_ti~     3 105.   78.9   31.2  58.2
#>  4  5.18        5.18           6.4  depth_per_ti~     4 104.   76.9   30.6  57.8
#>  5  6.4         6.4            7.62 depth_per_ti~     5 104.   82.6   31.5  57.7
#>  6  7.62        7.62           8.85 depth_per_ti~     6  86.3  44.6   33.8  57.0
#>  7  8.85        8.85          10.1  depth_per_ti~     7  86.6  48.8   32.1  55.4
#>  8 10.1        10.1           11.3  depth_per_ti~     8  88.6  56.3   31.7  54.0
#>  9 11.3        11.3           12.5  depth_per_ti~     9  90.0  61.5   31.1  54.4
#> 10 12.5        12.5           13.8  depth_per_ti~    10  90.2  63.6   30.2  53.1
#> # ... with 150 more rows, and 3 more variables: 50% <dbl>, 75% <dbl>,
#> #   97.5% <dbl>
```

### Diagnostic plots

Additional diagnostic plots are available. See ?plot.hamstr\_fit for
options.

#### Plot modelled accumulation rates at each hierarchical level

``` r
plot(hamstr_fit_1, type = "hier_acc")
#> Joining, by = "alpha_idx"
```

![](readme_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

#### Plot memory prior and posterior

As for this example the highest resolution sections are approximately 1
cm thick, there is not much difference between R and w.

``` r
plot(hamstr_fit_1, type = "mem")
```

![](readme_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

### Other `rstan` functions

Within the hamstr\_fit object is an *rstan* object on which all the
standard rstan functions should operate correctly.

For example:

``` r
rstan::check_divergences(hamstr_fit_1$fit)
#> 0 of 3000 iterations ended with a divergence.

rstan::stan_rhat(hamstr_fit_1$fit)
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
#> Warning: Removed 3 rows containing non-finite values (stat_bin).
```

![](readme_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

The first `alpha` parameter is the overall mean accumulation rate.

``` r
rstan::traceplot(hamstr_fit_1$fit, par = c("alpha[1]"),
                 inc_warmup = TRUE)
```

![](readme_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

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
