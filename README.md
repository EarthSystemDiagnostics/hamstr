# baconr: An rstan implementation of the Bayesian Age-Depth model *Bacon* (Blaauw and Christen, 2011).
Andrew M. Dolman  
`r Sys.Date()`  

-------------------------------

**baconr** implements the Bayesian Age~Depth model, **Bacon**, in the **Stan** probabilistic programming language. It is at a very early stage of development and at the moment implements only the core non-Gaussian AR1 model described in Blaauw and Christen (2011). There is as yet no ability 
to convert between ^14^C ages and calibrated calendar ages.

There are currently just three exported functions: 

* `make_stan_dat` prepares data, and parameter values to be passed as data, into the correct format for the Stan model 
* `stan_bacon` calls a pre-compiled Stan implementation of Bacon and estimates the age model
* `plot_stan_bacon` plots a sample from the posterior distribution of the estimated age model

**baconr** also includes the example data, core MSB2K, included in the existing C++ implementation of [Bacon v2.2](http://www.chrono.qub.ac.uk/blaauw/bacon.html)

The motivation for creating this package is to make use of the Bacon age modelling routine more flexible, and to make further development 
of the age model itself easier, by coding it in a widely used higher-level probabilistic programming language.



*  Blaauw, Maarten, and J. Andrés Christen. 2011. Flexible Paleoclimate Age-Depth Models Using an Autoregressive Gamma Process. Bayesian Analysis 6 (3): 457-74. doi:10.1214/ba/1339616472.

*  Stan Development Team. 2016. Stan Modeling Language Users Guide and Reference Manual, Version 2.14.0.   http://mc-stan.org



## Installation

**baconr** can be installed directly from github


```r
if (!require("devtools")) {
  install.packages("devtools")
}

devtools::install_github("andrewdolman/baconr")
```


## Development

We aim to develop this package following the [Guidelines for R packages providing interfaces to Stan](https://cran.r-project.org/web/packages/rstantools/vignettes/developer-guidelines.html)


## Example



```r
library(baconr)
library(knitr)
library(ggplot2)
library(tidyr)
library(dplyr)

opts_chunk$set(echo=TRUE, message = FALSE, warning = FALSE, cache = TRUE,
               fig.width = 6, fig.pos = "H", dpi = 150, autodep = TRUE)
```


```r
# Get number of sections K, so that they will be ~ 5cm
K_for_5cm <- round(diff(range(MSB2K$depth)) / 5)

make_stan_dat(depth = MSB2K$depth, 
  obs_age = MSB2K$age, 
  obs_err = MSB2K$error,
  K = K_for_5cm, nu = 6,
  acc_mean = 20, acc_var = "default",
  mem_mean = 0.7, mem_strength = 4)
```

```
## $depth
##  [1]  1.5  4.5  8.5 12.5 14.5 14.5 14.5 17.5 20.5 21.5 21.5 22.5 28.5 31.5
## [15] 32.5 33.5 34.5 37.5 38.5 41.5 43.5 46.5 47.5 48.5 49.5 50.5 52.5 53.5
## [29] 54.5 55.5 58.5 59.5 64.5 70.5 71.5 73.5 75.5 77.5 79.5 99.5
## 
## $obs_age
##  [1] 4128 4106 4046 4184 4076 4107 4097 4177 4220 4281 4374 4493 4452 4616
## [15] 4662 4743 4638 4810 4757 4839 4913 4880 4989 5070 4993 5115 5026 5242
## [29] 5159 5130 5238 5293 5293 5368 5498 5588 5514 5535 5644 5885
## 
## $obs_err
##  [1] 65 60 59 58 62 61 58 53 59 64 64 62 52 64 64 67 67 67 82 59 65 57 70
## [24] 66 67 79 51 64 50 66 65 38 54 51 69 55 57 52 77 45
## 
## $K
## [1] 20
## 
## $nu
## [1] 6
## 
## $acc_mean
## [1] 20
## 
## $acc_var
## [1] 266.6667
## 
## $mem_mean
## [1] 0.7
## 
## $mem_strength
## [1] 4
## 
## $N
## [1] 40
## 
## $mem_alpha
## [1] 2.8
## 
## $mem_beta
## [1] 1.2
## 
## $c
##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
## 
## $delta_c
## [1] 5.22
## 
## $c_depth_bottom
##  [1]   5.22  10.44  15.66  20.88  26.10  31.32  36.54  41.76  46.98  52.20
## [11]  57.42  62.64  67.86  73.08  78.30  83.52  88.74  93.96  99.18 104.40
## 
## $c_depth_top
##  [1]  0.00  5.22 10.44 15.66 20.88 26.10 31.32 36.54 41.76 46.98 52.20
## [12] 57.42 62.64 67.86 73.08 78.30 83.52 88.74 93.96 99.18
## 
## $which_c
##  [1]  1  1  2  3  3  3  3  4  4  5  5  5  6  7  7  7  7  8  8  8  9  9 10
## [24] 10 10 10 11 11 11 11 12 12 13 14 14 15 15 15 16 20
```


```r
fit <- stan_bacon(
  depth = MSB2K$depth, 
  obs_age = MSB2K$age, 
  obs_err = MSB2K$error,
  K = K_for_5cm, nu = 6,
  acc_mean = 20, acc_var = "default",
  mem_mean = 0.7, mem_strength = 4,
  iter = 2000, chains = 4)
```

```
## 
## SAMPLING FOR MODEL 'bacon' NOW (CHAIN 1).
## 
## Chain 1, Iteration:    1 / 2000 [  0%]  (Warmup)
## Chain 1, Iteration:  200 / 2000 [ 10%]  (Warmup)
## Chain 1, Iteration:  400 / 2000 [ 20%]  (Warmup)
## Chain 1, Iteration:  600 / 2000 [ 30%]  (Warmup)
## Chain 1, Iteration:  800 / 2000 [ 40%]  (Warmup)
## Chain 1, Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Chain 1, Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Chain 1, Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Chain 1, Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Chain 1, Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Chain 1, Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Chain 1, Iteration: 2000 / 2000 [100%]  (Sampling)
##  Elapsed Time: 1.415 seconds (Warm-up)
##                0.786 seconds (Sampling)
##                2.201 seconds (Total)
## 
## 
## SAMPLING FOR MODEL 'bacon' NOW (CHAIN 2).
## 
## Chain 2, Iteration:    1 / 2000 [  0%]  (Warmup)
## Chain 2, Iteration:  200 / 2000 [ 10%]  (Warmup)
## Chain 2, Iteration:  400 / 2000 [ 20%]  (Warmup)
## Chain 2, Iteration:  600 / 2000 [ 30%]  (Warmup)
## Chain 2, Iteration:  800 / 2000 [ 40%]  (Warmup)
## Chain 2, Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Chain 2, Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Chain 2, Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Chain 2, Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Chain 2, Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Chain 2, Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Chain 2, Iteration: 2000 / 2000 [100%]  (Sampling)
##  Elapsed Time: 1.519 seconds (Warm-up)
##                1.037 seconds (Sampling)
##                2.556 seconds (Total)
## 
## 
## SAMPLING FOR MODEL 'bacon' NOW (CHAIN 3).
## 
## Chain 3, Iteration:    1 / 2000 [  0%]  (Warmup)
## Chain 3, Iteration:  200 / 2000 [ 10%]  (Warmup)
## Chain 3, Iteration:  400 / 2000 [ 20%]  (Warmup)
## Chain 3, Iteration:  600 / 2000 [ 30%]  (Warmup)
## Chain 3, Iteration:  800 / 2000 [ 40%]  (Warmup)
## Chain 3, Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Chain 3, Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Chain 3, Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Chain 3, Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Chain 3, Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Chain 3, Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Chain 3, Iteration: 2000 / 2000 [100%]  (Sampling)
##  Elapsed Time: 1.524 seconds (Warm-up)
##                0.789 seconds (Sampling)
##                2.313 seconds (Total)
## 
## 
## SAMPLING FOR MODEL 'bacon' NOW (CHAIN 4).
## 
## Chain 4, Iteration:    1 / 2000 [  0%]  (Warmup)
## Chain 4, Iteration:  200 / 2000 [ 10%]  (Warmup)
## Chain 4, Iteration:  400 / 2000 [ 20%]  (Warmup)
## Chain 4, Iteration:  600 / 2000 [ 30%]  (Warmup)
## Chain 4, Iteration:  800 / 2000 [ 40%]  (Warmup)
## Chain 4, Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Chain 4, Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Chain 4, Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Chain 4, Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Chain 4, Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Chain 4, Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Chain 4, Iteration: 2000 / 2000 [100%]  (Sampling)
##  Elapsed Time: 1.442 seconds (Warm-up)
##                0.717 seconds (Sampling)
##                2.159 seconds (Total)
```




```r
print(fit$fit, par = c("c_ages"))
```

```
## Inference for Stan model: bacon.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##               mean se_mean    sd    2.5%     25%     50%     75%   97.5%
## c_ages[1]  4039.88    0.68 38.00 3957.44 4016.82 4042.67 4065.56 4106.66
## c_ages[2]  4079.28    0.48 30.30 4012.99 4060.52 4080.06 4099.54 4134.67
## c_ages[3]  4110.93    0.44 27.71 4055.91 4093.26 4111.04 4128.93 4164.56
## c_ages[4]  4149.50    0.46 28.81 4097.25 4129.38 4147.66 4168.41 4210.27
## c_ages[5]  4294.94    0.62 38.97 4219.86 4268.50 4294.87 4320.38 4372.01
## c_ages[6]  4434.93    1.35 70.84 4310.50 4382.69 4432.61 4482.90 4578.78
## c_ages[7]  4608.68    0.65 41.35 4529.86 4581.21 4608.70 4636.49 4687.50
## c_ages[8]  4750.93    0.89 44.59 4662.08 4721.01 4752.01 4781.50 4835.49
## c_ages[9]  4853.97    0.63 39.61 4776.34 4827.93 4853.50 4880.85 4932.24
## c_ages[10] 4961.58    0.61 38.70 4888.53 4934.46 4960.88 4987.93 5038.01
## c_ages[11] 5093.46    0.75 40.09 5019.35 5065.35 5092.30 5120.86 5174.14
## c_ages[12] 5227.27    0.54 34.07 5157.39 5203.94 5228.20 5251.04 5291.26
## c_ages[13] 5299.24    0.63 39.64 5225.53 5272.55 5297.66 5324.95 5382.50
## c_ages[14] 5370.13    0.91 48.06 5281.47 5336.57 5368.11 5401.75 5471.60
## c_ages[15] 5498.39    0.58 36.83 5421.66 5474.38 5499.72 5523.55 5566.43
## c_ages[16] 5569.93    0.57 36.27 5502.53 5544.79 5568.63 5594.33 5644.10
## c_ages[17] 5654.47    0.96 60.99 5550.54 5610.23 5649.78 5693.51 5782.91
## c_ages[18] 5733.10    1.09 68.92 5603.52 5683.47 5732.68 5780.12 5868.20
## c_ages[19] 5810.04    1.04 66.09 5673.27 5767.58 5812.43 5855.87 5931.03
## c_ages[20] 5887.31    0.79 50.27 5793.45 5856.24 5884.41 5916.76 5993.22
## c_ages[21] 5990.81    1.48 93.39 5851.08 5925.01 5976.40 6038.06 6210.41
##            n_eff Rhat
## c_ages[1]   3086    1
## c_ages[2]   4000    1
## c_ages[3]   4000    1
## c_ages[4]   4000    1
## c_ages[5]   4000    1
## c_ages[6]   2747    1
## c_ages[7]   4000    1
## c_ages[8]   2509    1
## c_ages[9]   4000    1
## c_ages[10]  4000    1
## c_ages[11]  2835    1
## c_ages[12]  4000    1
## c_ages[13]  4000    1
## c_ages[14]  2765    1
## c_ages[15]  4000    1
## c_ages[16]  4000    1
## c_ages[17]  4000    1
## c_ages[18]  4000    1
## c_ages[19]  4000    1
## c_ages[20]  4000    1
## c_ages[21]  4000    1
## 
## Samples were drawn using NUTS(diag_e) at Thu Apr 06 22:50:23 2017.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```



```r
set.seed(20170406)
plot_stan_bacon(fit, 100)
```

![](readme_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

