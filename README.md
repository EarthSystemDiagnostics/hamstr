---
title: "baconr: An rstan implementation of the Bayesian Age-Depth model *Bacon* (Blaauw and Christen, 2011)."
output: 
  html_document: 
    keep_md: yes
---

-------------------------------

**baconr** implements the Bayesian Age~Depth model, **Bacon**, in the **Stan** probabilistic programming language. It is at a very early stage of development and at the moment implements only the core non-Gaussian AR1 model described in Blaauw and Christen (2011). Functions from the R package **Bchron** (Parnell 2016) can be used to calibrate ^14^C ages to calendar ages.


There are currently just three exported functions: 

* `make_stan_dat` prepares data and parameter values into the correct format to be passed to the Stan model 
* `stan_bacon` calls a pre-compiled Stan implementation of Bacon and estimates the age model
* `plot_stan_bacon` plots a sample from the posterior distribution of the estimated age model

**baconr** also includes the example data, core MSB2K, included in the existing C++ implementation of [Bacon v2.2](http://www.chrono.qub.ac.uk/blaauw/bacon.html)

The motivation for creating this package is to make use of the Bacon age modelling routine more flexible, and to make further development 
of the age model itself easier, by coding it in a widely used higher-level probabilistic programming language.


*  Parnell, Andrew. 2016. Bchron: Radiocarbon Dating, Age-Depth Modelling, Relative Sea Level Rate Estimation, and Non-Parametric Phase Modelling. R package version 4.2.6. https://CRAN.R-project.org/package=Bchron

*  Blaauw, Maarten, and J. Andrés Christen. 2011. Flexible Paleoclimate Age-Depth Models Using an Autoregressive Gamma Process. Bayesian Analysis 6 (3): 457-74. doi:10.1214/ba/1339616472.

*  Stan Development Team. 2016. Stan Modeling Language Users Guide and Reference Manual, Version 2.14.0.   http://mc-stan.org




## Installation

**baconr** can be installed directly from Github


```r
if (!require("devtools")) {
  install.packages("devtools")
}

devtools::install_github("andrewdolman/baconr", args = "--preclean", build_vignettes = FALSE)
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
library(Bchron)

opts_chunk$set(echo=TRUE, message = FALSE, warning = FALSE, cache = TRUE,
               fig.width = 6, fig.pos = "H", dpi = 300, autodep = TRUE)
```

### Data and parameters

First convert 14^C ages to calendar ages. Functions from the R package **Bchron** can be used for this.


```r
cal.ages <- Bchron::BchronCalibrate(ages=MSB2K$age,
                            ageSds=MSB2K$error,
                            calCurves=rep("intcal13", nrow(MSB2K)),
                            ids=paste0("Date-", 1:nrow(MSB2K)))

SimplifyCalAge <- function(ages, densities){
  
  wts <- densities / sum(densities)
  
  M <- sum(wts > 0)

  w.mean <- sum(ages * wts)
  w.std <- sqrt(sum(wts * (ages-w.mean)^2) / ((M-1)/M * sum(wts)))
  
  return(c("mean" = w.mean, "std" = w.std))
  
}


MSB2K$age.cal <- sapply(cal.ages, function(x) SimplifyCalAge(x$ageGrid, x$densities)["mean"])
MSB2K$age.cal.sd <- sapply(cal.ages, function(x) SimplifyCalAge(x$ageGrid, x$densities)["std"])
```


the output from `make_stan_dat` is shown here to illustrate the required data and parameter format, but a separate call to `make_stan_dat` is not normally necessary as it is used internally by `stan_bacon`


```r
# Get number of sections K, so that they will be ~ 5cm
K_for_5cm <- round(diff(range(MSB2K$depth)) / 5)

stan_dat <- make_stan_dat(depth = MSB2K$depth, 
  obs_age = MSB2K$age.cal, 
  obs_err = MSB2K$age.cal.sd,
  K = K_for_5cm, nu = 6,
  acc_mean = 20, acc_alpha = 1.5,
  mem_mean = 0.7, mem_strength = 4)

stan_dat
```

```
## $depth
##  [1]  1.5  4.5  8.5 12.5 14.5 14.5 14.5 17.5 20.5 21.5 21.5 22.5 28.5 31.5
## [15] 32.5 33.5 34.5 37.5 38.5 41.5 43.5 46.5 47.5 48.5 49.5 50.5 52.5 53.5
## [29] 54.5 55.5 58.5 59.5 64.5 70.5 71.5 73.5 75.5 77.5 79.5 99.5
## 
## $obs_age
##  [1] 4663.006 4643.940 4559.684 4707.835 4606.703 4644.705 4634.061
##  [8] 4705.752 4737.626 4843.905 4989.204 5140.893 5101.359 5330.612
## [15] 5404.113 5467.301 5368.963 5523.004 5474.133 5560.100 5661.345
## [22] 5619.289 5743.136 5808.139 5746.404 5848.643 5780.562 6033.162
## [29] 5903.081 5861.310 6028.906 6078.348 6078.738 6151.047 6296.178
## [36] 6372.028 6318.868 6338.488 6437.480 6706.671
## 
## $obs_err
##  [1] 100.20567 104.14839 110.60003  83.46879 112.46756 104.35427 105.81326
##  [8]  80.32471  85.63910 105.20423 113.63604 108.04005 115.90679 127.29059
## [15]  95.85990  84.52935 116.69397  83.28079  94.58266  74.86960  77.45291
## [22]  66.45828  85.77970  74.27373  84.25610  97.50202  73.56497  88.55855
## [29]  75.77078  86.64969  90.17830  64.11373  76.04572  81.85147  75.21436
## [36]  49.59601  56.20994  45.74972  85.31270  50.88403
## 
## $hiatus_depth
## NULL
## 
## $hiatus_length
## NULL
## 
## $hiatus_shape
## [1] 1
## 
## $hiatus_interval
## [1] 0.1
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
## $acc_alpha
## [1] 1.5
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
## $acc_beta
## [1] 0.075
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




### Fit the bacon model with `stan_bacon`


```r
fit <- stan_bacon(
  depth = MSB2K$depth, 
  obs_age = MSB2K$age.cal, 
  obs_err = MSB2K$age.cal.sd,
  K = K_for_5cm, nu = 6,
  acc_mean = 20, acc_alpha = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  iter = 2000, chains = 3)
```

```
## 
## SAMPLING FOR MODEL 'bacon' NOW (CHAIN 1).
## 
## Gradient evaluation took 5.1e-05 seconds
## 1000 transitions using 10 leapfrog steps per transition would take 0.51 seconds.
## Adjust your expectations accordingly!
## 
## 
## Iteration:    1 / 2000 [  0%]  (Warmup)
## Iteration:  200 / 2000 [ 10%]  (Warmup)
## Iteration:  400 / 2000 [ 20%]  (Warmup)
## Iteration:  600 / 2000 [ 30%]  (Warmup)
## Iteration:  800 / 2000 [ 40%]  (Warmup)
## Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Iteration: 2000 / 2000 [100%]  (Sampling)
## 
##  Elapsed Time: 2.18345 seconds (Warm-up)
##                1.53961 seconds (Sampling)
##                3.72306 seconds (Total)
## 
## 
## SAMPLING FOR MODEL 'bacon' NOW (CHAIN 2).
## 
## Gradient evaluation took 1.8e-05 seconds
## 1000 transitions using 10 leapfrog steps per transition would take 0.18 seconds.
## Adjust your expectations accordingly!
## 
## 
## Iteration:    1 / 2000 [  0%]  (Warmup)
## Iteration:  200 / 2000 [ 10%]  (Warmup)
## Iteration:  400 / 2000 [ 20%]  (Warmup)
## Iteration:  600 / 2000 [ 30%]  (Warmup)
## Iteration:  800 / 2000 [ 40%]  (Warmup)
## Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Iteration: 2000 / 2000 [100%]  (Sampling)
## 
##  Elapsed Time: 2.05763 seconds (Warm-up)
##                1.5288 seconds (Sampling)
##                3.58643 seconds (Total)
## 
## 
## SAMPLING FOR MODEL 'bacon' NOW (CHAIN 3).
## 
## Gradient evaluation took 1.7e-05 seconds
## 1000 transitions using 10 leapfrog steps per transition would take 0.17 seconds.
## Adjust your expectations accordingly!
## 
## 
## Iteration:    1 / 2000 [  0%]  (Warmup)
## Iteration:  200 / 2000 [ 10%]  (Warmup)
## Iteration:  400 / 2000 [ 20%]  (Warmup)
## Iteration:  600 / 2000 [ 30%]  (Warmup)
## Iteration:  800 / 2000 [ 40%]  (Warmup)
## Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Iteration: 2000 / 2000 [100%]  (Sampling)
## 
##  Elapsed Time: 2.22362 seconds (Warm-up)
##                1.48273 seconds (Sampling)
##                3.70634 seconds (Total)
```



```r
options(width = 85)
print(fit$fit, par = c("R", "w", "c_ages"))
```

```
## Inference for Stan model: bacon.
## 3 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=3000.
## 
##               mean se_mean    sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
## R             0.73    0.00  0.18    0.29    0.63    0.78    0.87    0.94  1799    1
## w             0.30    0.01  0.22    0.00    0.09    0.28    0.47    0.72  1725    1
## c_ages[1]  4571.78    0.93 50.73 4463.70 4540.89 4575.39 4605.93 4662.86  3000    1
## c_ages[2]  4611.96    0.80 43.61 4519.83 4585.77 4613.60 4641.45 4691.18  3000    1
## c_ages[3]  4657.49    0.72 39.53 4577.14 4631.48 4657.85 4683.84 4733.11  3000    1
## c_ages[4]  4718.55    0.75 40.88 4640.47 4691.31 4718.34 4744.71 4801.53  3000    1
## c_ages[5]  4864.85    1.02 55.71 4755.74 4826.32 4863.32 4902.85 4970.28  3000    1
## c_ages[6]  5094.61    1.74 85.60 4920.92 5041.33 5096.38 5150.37 5261.55  2430    1
## c_ages[7]  5317.74    0.96 52.69 5212.87 5283.63 5317.63 5353.55 5418.70  3000    1
## c_ages[8]  5468.79    0.78 42.79 5387.30 5440.42 5468.53 5496.87 5551.85  3000    1
## c_ages[9]  5588.78    0.73 39.84 5505.82 5563.06 5590.21 5615.51 5666.17  3000    1
## c_ages[10] 5711.90    0.67 36.89 5635.88 5688.33 5712.99 5736.86 5781.00  3000    1
## c_ages[11] 5848.83    0.71 39.01 5769.01 5823.03 5850.31 5874.09 5924.35  3000    1
## c_ages[12] 5980.77    0.73 39.73 5900.93 5955.43 5980.16 6006.59 6059.57  3000    1
## c_ages[13] 6085.76    0.82 44.96 5999.67 6055.20 6087.12 6115.30 6173.87  3000    1
## c_ages[14] 6184.12    0.86 46.95 6082.58 6156.28 6186.72 6216.31 6269.04  3000    1
## c_ages[15] 6296.79    0.57 31.34 6234.16 6276.20 6297.75 6317.95 6355.54  3000    1
## c_ages[16] 6376.09    0.58 31.81 6315.44 6354.82 6375.18 6396.30 6442.39  3000    1
## c_ages[17] 6460.69    0.91 50.00 6367.76 6426.74 6457.64 6490.57 6568.26  3000    1
## c_ages[18] 6542.83    1.07 58.57 6427.63 6503.54 6543.48 6581.17 6659.80  3000    1
## c_ages[19] 6624.42    1.06 58.29 6505.01 6587.96 6626.66 6662.92 6738.17  3000    1
## c_ages[20] 6709.57    0.96 52.80 6608.30 6677.26 6707.02 6739.43 6822.53  3000    1
## c_ages[21] 6809.23    1.51 82.70 6667.54 6754.56 6800.79 6853.60 7000.61  3000    1
## 
## Samples were drawn using NUTS(diag_e) at Mon Sep  3 16:46:25 2018.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```


### Plot the estimated age model


```r
set.seed(20170406)
plot_stan_bacon(fit, 1000)
```

![](readme_files/figure-html/bacon_defaults-1.png)<!-- -->


### Fit again with stronger prior on accumulation rates



```r
fit2 <- stan_bacon(
  depth = MSB2K$depth, 
  obs_age = MSB2K$age.cal, 
  obs_err = MSB2K$age.cal.sd,
  K = K_for_5cm, nu = 6,
  acc_mean = 20, acc_alpha = 25,
  mem_mean = 0.7, mem_strength = 4,
  iter = 2000, chains = 3)
```



```r
plot_stan_bacon(fit2, 1000)
```

![](readme_files/figure-html/stronger_prior-1.png)<!-- -->


### Fit again with strong prior for lower memory



```r
fit3 <- stan_bacon(
  depth = MSB2K$depth, 
  obs_age = MSB2K$age.cal, 
  obs_err = MSB2K$age.cal.sd,
  K = K_for_5cm, nu = 6,
  acc_mean = 20, acc_alpha = 1.5,
  mem_mean = 0.1, mem_strength = 4,
  iter = 2000, chains = 3)
```



```r
plot_stan_bacon(fit3, 1000)
```

![](readme_files/figure-html/weaker_prior-1.png)<!-- -->


### Add tephras at 25 and 90 cm



```r
teph <- data.frame(depth = c(25, 87.5), age.cal = c(5000, 6500), age.cal.sd = c(3, 5))
MSB2K.2 <- bind_rows(MSB2K, teph) %>% 
  arrange(age.cal)
  

fit4 <- stan_bacon(
  depth = MSB2K.2$depth, 
  obs_age = MSB2K.2$age.cal, 
  obs_err = MSB2K.2$age.cal.sd,
  K = K_for_5cm, nu = 6,
  acc_mean = 20, acc_alpha = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  iter = 2000, chains = 3)
```



```r
plot_stan_bacon(fit4, 1000) +
  geom_pointrange(data = teph,
                  aes(x = depth, y = age.cal, ymax = age.cal + age.cal.sd, ymin = age.cal - age.cal.sd),
                  group = NA, colour = "Blue", alpha = 0.5) 
```

![](readme_files/figure-html/tephras-1.png)<!-- -->


### Add hiatus


```r
# Add 1000 years to all depths below 40 cm
hi_ages <- MSB2K$age
hi_ages[MSB2K$depth > 40] <- hi_ages[MSB2K$depth > 40] + 1000

fit_hi <- stan_bacon(
  depth = MSB2K$depth, 
  obs_age = hi_ages, 
  obs_err = MSB2K$age.cal.sd,
  hiatus_depth = c(40),
  hiatus_length = c(1000),
  K = 20, nu = 6,
  acc_mean = 20, acc_alpha = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  iter = 2000, chains = 3)
```


```r
plot_stan_bacon(fit_hi)
```

![](readme_files/figure-html/hiatus-1.png)<!-- -->



### Distribution of sediment accumulation rates


```r
age.fit <- fit
age.mod <- rstan::extract(age.fit$fit)

bp.dat <- age.mod$x[1:3000,] %>% 
  as_tibble() %>% 
  tibble::rownames_to_column("Rep") %>% 
  gather(Depth, value, -Rep) %>% 
  mutate(Depth = as.numeric(gsub("V", "", Depth)),
         Depth = Depth * age.fit$data$delta_c,
         Rep = as.numeric(Rep))

# bp.dat.bw <- bp.dat %>% 
#   group_by(Depth) %>% 
#   summarise(mean.val = mean(value),
#             upr = quantile(value, 0.9),
#             lwr = quantile(value, 0.1))
# 
# bp.dat.bw %>% 
#   ggplot(aes(x = Depth)) + 
#   geom_pointrange(aes(y = mean.val, ymax = upr, ymin = lwr))
```


```r
bp.dat %>% 
  ggplot(aes(x = Depth, y = 1/value, group = Depth)) + 
  geom_violin() +
  #geom_boxplot() +
  #geom_point(alpha = 0.15)
  scale_x_continuous("Depth [cm]")+
  scale_y_continuous("Sediment accumulation rate [cm/yr]",
                     trans = "log10", breaks = c(0.01, 0.05, 0.1, 0.5)) +
  annotation_logticks(sides = "l") +
  expand_limits(y = 0.01) + 
  theme_bw()
```

![](readme_files/figure-html/acc_rates-1.png)<!-- -->


```r
bp.dat %>% 
  ggplot(aes(x = Depth, y = value, group = Depth)) + 
  geom_violin() +
  #geom_boxplot() +
  #geom_point(alpha = 0.15)
  scale_x_continuous("Depth [cm]")+
  scale_y_continuous("Sedimentation interval [yr/cm]",
                     trans = "log10", breaks = c(1, 5, 10, 20, 50)) +
  annotation_logticks(sides = "l") +
#  expand_limits(y = 0.01) + 
  theme_bw()
```

![](readme_files/figure-html/acc_intervals-1.png)<!-- -->






