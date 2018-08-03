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
## Gradient evaluation took 0 seconds
## 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
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
##  Elapsed Time: 2.045 seconds (Warm-up)
##                1.357 seconds (Sampling)
##                3.402 seconds (Total)
## 
## 
## SAMPLING FOR MODEL 'bacon' NOW (CHAIN 2).
## 
## Gradient evaluation took 0 seconds
## 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
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
##  Elapsed Time: 1.937 seconds (Warm-up)
##                1.355 seconds (Sampling)
##                3.292 seconds (Total)
## 
## 
## SAMPLING FOR MODEL 'bacon' NOW (CHAIN 3).
## 
## Gradient evaluation took 0 seconds
## 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
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
##  Elapsed Time: 1.943 seconds (Warm-up)
##                1.334 seconds (Sampling)
##                3.277 seconds (Total)
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
## R             0.73    0.00  0.18    0.29    0.63    0.78    0.86    0.94  2070    1
## w             0.29    0.01  0.22    0.00    0.09    0.27    0.47    0.71  1839    1
## c_ages[1]  4573.50    0.88 48.00 4473.82 4543.32 4575.54 4607.04 4659.60  3000    1
## c_ages[2]  4612.87    0.77 41.99 4528.62 4584.65 4614.57 4641.71 4689.95  3000    1
## c_ages[3]  4656.90    0.72 39.34 4577.84 4630.75 4657.62 4684.27 4730.64  3000    1
## c_ages[4]  4717.58    0.76 41.53 4637.12 4690.90 4716.21 4744.41 4800.38  3000    1
## c_ages[5]  4862.61    1.00 54.75 4760.58 4824.77 4862.12 4899.54 4971.74  3000    1
## c_ages[6]  5096.96    1.71 84.40 4924.89 5040.93 5100.29 5151.61 5262.94  2435    1
## c_ages[7]  5319.55    0.95 52.22 5216.42 5285.60 5319.26 5354.62 5422.30  3000    1
## c_ages[8]  5470.25    0.76 41.72 5386.68 5444.00 5470.40 5496.72 5550.31  3000    1
## c_ages[9]  5589.82    0.69 38.00 5510.56 5565.14 5591.18 5615.51 5661.10  3000    1
## c_ages[10] 5712.59    0.66 36.32 5639.86 5688.01 5714.03 5737.88 5779.46  3000    1
## c_ages[11] 5848.11    0.69 37.93 5771.57 5823.01 5849.39 5873.08 5921.27  3000    1
## c_ages[12] 5981.24    0.72 39.42 5903.19 5954.51 5982.41 6006.95 6059.26  3000    1
## c_ages[13] 6087.21    0.80 43.74 6002.39 6058.34 6087.48 6116.77 6173.36  3000    1
## c_ages[14] 6185.54    0.87 47.85 6084.94 6155.00 6189.17 6218.59 6271.70  3000    1
## c_ages[15] 6297.44    0.58 31.81 6230.28 6276.77 6298.33 6319.18 6358.12  3000    1
## c_ages[16] 6375.13    0.58 31.94 6313.38 6353.34 6374.42 6396.71 6440.18  3000    1
## c_ages[17] 6459.12    0.91 49.94 6367.28 6425.51 6456.82 6490.32 6562.01  3000    1
## c_ages[18] 6540.40    1.08 59.24 6421.17 6501.19 6541.78 6577.11 6662.44  3000    1
## c_ages[19] 6621.52    1.09 59.52 6499.68 6583.91 6622.73 6659.22 6736.90  3000    1
## c_ages[20] 6707.87    0.98 53.87 6605.34 6673.48 6706.02 6740.03 6819.63  3000    1
## c_ages[21] 6805.96    1.50 82.06 6666.92 6748.84 6798.72 6854.39 6986.32  3000    1
## 
## Samples were drawn using NUTS(diag_e) at Fri Aug 03 15:27:28 2018.
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






