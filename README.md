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
##  [1]  1.5  4.5  8.5 12.5 14.5 14.5 14.5 17.5 20.5 21.5 21.5 22.5 28.5 31.5 32.5
## [16] 33.5 34.5 37.5 38.5 41.5 43.5 46.5 47.5 48.5 49.5 50.5 52.5 53.5 54.5 55.5
## [31] 58.5 59.5 64.5 70.5 71.5 73.5 75.5 77.5 79.5 99.5
## 
## $obs_age
##  [1] 4663.006 4643.940 4559.684 4707.835 4606.703 4644.705 4634.061 4705.752
##  [9] 4737.626 4843.905 4989.204 5140.893 5101.359 5330.612 5404.113 5467.301
## [17] 5368.963 5523.004 5474.133 5560.100 5661.345 5619.289 5743.136 5808.139
## [25] 5746.404 5848.643 5780.562 6033.162 5903.081 5861.310 6028.906 6078.348
## [33] 6078.738 6151.047 6296.178 6372.028 6318.868 6338.488 6437.480 6706.671
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
##  [1]  0.00  5.22 10.44 15.66 20.88 26.10 31.32 36.54 41.76 46.98 52.20 57.42
## [13] 62.64 67.86 73.08 78.30 83.52 88.74 93.96 99.18
## 
## $which_c
##  [1]  1  1  2  3  3  3  3  4  4  5  5  5  6  7  7  7  7  8  8  8  9  9 10 10 10
## [26] 10 11 11 11 11 12 12 13 14 14 15 15 15 16 20
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
## Chain 1: 
## Chain 1: Gradient evaluation took 5.6e-05 seconds
## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.56 seconds.
## Chain 1: Adjust your expectations accordingly!
## Chain 1: 
## Chain 1: 
## Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
## Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
## Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
## Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
## Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
## Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
## Chain 1: 
## Chain 1:  Elapsed Time: 2.30366 seconds (Warm-up)
## Chain 1:                1.48962 seconds (Sampling)
## Chain 1:                3.79328 seconds (Total)
## Chain 1: 
## 
## SAMPLING FOR MODEL 'bacon' NOW (CHAIN 2).
## Chain 2: 
## Chain 2: Gradient evaluation took 2.1e-05 seconds
## Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.21 seconds.
## Chain 2: Adjust your expectations accordingly!
## Chain 2: 
## Chain 2: 
## Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
## Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
## Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
## Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
## Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
## Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
## Chain 2: 
## Chain 2:  Elapsed Time: 2.0355 seconds (Warm-up)
## Chain 2:                1.39621 seconds (Sampling)
## Chain 2:                3.43171 seconds (Total)
## Chain 2: 
## 
## SAMPLING FOR MODEL 'bacon' NOW (CHAIN 3).
## Chain 3: 
## Chain 3: Gradient evaluation took 1.9e-05 seconds
## Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.19 seconds.
## Chain 3: Adjust your expectations accordingly!
## Chain 3: 
## Chain 3: 
## Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
## Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
## Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
## Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
## Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
## Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
## Chain 3: 
## Chain 3:  Elapsed Time: 2.13185 seconds (Warm-up)
## Chain 3:                1.57365 seconds (Sampling)
## Chain 3:                3.70551 seconds (Total)
## Chain 3:
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
## R             0.72    0.00  0.18    0.27    0.61    0.77    0.86    0.94  2230    1
## w             0.29    0.00  0.22    0.00    0.08    0.25    0.47    0.71  2113    1
## c_ages[1]  4573.37    0.90 49.81 4473.83 4541.43 4576.07 4606.96 4667.31  3081    1
## c_ages[2]  4613.41    0.68 43.11 4528.33 4584.73 4614.60 4642.34 4696.39  4067    1
## c_ages[3]  4658.10    0.67 40.90 4577.42 4630.62 4658.66 4686.04 4738.92  3708    1
## c_ages[4]  4718.25    0.75 40.96 4641.13 4690.69 4717.79 4746.04 4800.16  2944    1
## c_ages[5]  4861.94    1.10 56.56 4755.90 4822.12 4861.27 4899.77 4978.31  2640    1
## c_ages[6]  5097.06    1.63 86.24 4920.02 5044.34 5098.30 5150.41 5267.04  2789    1
## c_ages[7]  5318.39    0.85 53.10 5215.31 5283.12 5318.30 5354.00 5421.20  3885    1
## c_ages[8]  5470.69    0.79 42.35 5386.77 5442.80 5470.64 5499.19 5552.13  2892    1
## c_ages[9]  5589.73    0.67 39.06 5507.39 5564.69 5591.47 5614.61 5666.04  3427    1
## c_ages[10] 5712.15    0.65 37.17 5633.94 5687.87 5713.69 5737.98 5779.94  3224    1
## c_ages[11] 5847.78    0.72 38.47 5770.39 5823.53 5848.40 5873.44 5921.15  2816    1
## c_ages[12] 5980.56    0.66 38.27 5903.00 5955.96 5980.12 6005.99 6055.01  3333    1
## c_ages[13] 6086.04    0.75 43.58 5999.57 6057.65 6085.84 6114.66 6172.38  3391    1
## c_ages[14] 6183.92    0.85 46.82 6081.15 6155.00 6186.88 6215.17 6266.01  3058    1
## c_ages[15] 6297.86    0.54 30.80 6234.84 6278.00 6298.45 6318.69 6356.24  3222    1
## c_ages[16] 6375.34    0.55 31.93 6315.57 6353.19 6374.60 6396.09 6442.08  3428    1
## c_ages[17] 6460.90    0.87 50.46 6371.78 6426.00 6458.80 6491.04 6572.70  3387    1
## c_ages[18] 6541.88    0.99 58.64 6428.86 6502.35 6543.08 6579.05 6663.10  3491    1
## c_ages[19] 6624.08    1.03 59.00 6502.84 6587.20 6624.11 6661.81 6737.75  3307    1
## c_ages[20] 6707.15    0.90 51.83 6606.28 6674.77 6705.87 6738.50 6814.59  3284    1
## c_ages[21] 6805.24    1.30 85.12 6667.69 6748.63 6797.16 6849.93 7005.70  4274    1
## 
## Samples were drawn using NUTS(diag_e) at Mon Feb 10 17:03:53 2020.
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






