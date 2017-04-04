# baconr: An rstan implementation of the *Bacon* Age-Depth model of Blaauw and Christen (2011).
Andrew M. Dolman  
`r Sys.Date()`  


```r
library(tidyverse)
```

```
## Loading tidyverse: ggplot2
## Loading tidyverse: tibble
## Loading tidyverse: tidyr
## Loading tidyverse: readr
## Loading tidyverse: purrr
## Loading tidyverse: dplyr
```

```
## Conflicts with tidy packages ----------------------------------------------
```

```
## filter(): dplyr, stats
## lag():    dplyr, stats
```

```r
library(baconr)
```

```
## Loading required package: Rcpp
```

```r
library(knitr)

opts_chunk$set(echo=FALSE, message = FALSE, warning = FALSE, cache = TRUE,
               fig.width = 7, fig.pos = "H", dpi = 150, autodep = TRUE)
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
##  Elapsed Time: 1.74342 seconds (Warm-up)
##                0.995789 seconds (Sampling)
##                2.73921 seconds (Total)
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
##  Elapsed Time: 1.67576 seconds (Warm-up)
##                1.13479 seconds (Sampling)
##                2.81055 seconds (Total)
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
##  Elapsed Time: 1.78697 seconds (Warm-up)
##                1.00955 seconds (Sampling)
##                2.79652 seconds (Total)
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
##  Elapsed Time: 1.95335 seconds (Warm-up)
##                1.04124 seconds (Sampling)
##                2.99459 seconds (Total)
```


```
## Inference for Stan model: bacon.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##               mean se_mean    sd    2.5%     25%     50%     75%   97.5%
## c_ages[1]     0.00    0.00  0.00    0.00    0.00    0.00    0.00    0.00
## c_ages[2]    39.51    0.55 30.40    3.11   16.99   32.43   53.97  117.33
## c_ages[3]    74.17    0.68 37.77   18.38   46.00   67.79   95.79  161.63
## c_ages[4]   114.69    0.81 44.59   41.70   82.53  110.51  141.76  214.58
## c_ages[5]   263.27    0.88 55.64  158.00  225.87  262.35  298.69  379.28
## c_ages[6]   390.79    1.57 79.80  245.14  334.32  387.54  444.66  553.96
## c_ages[7]   550.23    1.02 57.72  438.28  510.83  550.41  588.34  669.66
## c_ages[8]   693.19    1.04 59.06  579.15  652.79  693.20  733.44  810.11
## c_ages[9]   791.51    0.96 55.77  688.14  753.06  789.82  827.81  906.49
## c_ages[10]  882.22    1.07 56.36  776.35  842.42  882.05  919.23  995.35
## c_ages[11] 1018.74    0.98 55.04  916.66  980.62 1017.36 1054.99 1131.50
## c_ages[12] 1128.23    0.89 52.39 1028.76 1092.41 1127.37 1162.50 1233.61
## c_ages[13] 1237.98    0.84 53.12 1136.87 1202.52 1236.65 1272.64 1347.34
## c_ages[14] 1295.10    0.96 56.18 1192.94 1255.99 1292.74 1330.96 1412.09
## c_ages[15] 1380.38    1.00 57.08 1273.44 1341.15 1378.93 1417.98 1496.49
## c_ages[16] 1485.31    0.93 53.18 1386.17 1449.81 1483.46 1518.81 1593.96
## c_ages[17] 1559.91    0.92 58.11 1453.43 1519.61 1557.27 1597.36 1679.43
## c_ages[18] 1631.48    1.23 72.90 1501.26 1579.96 1627.20 1679.30 1787.52
## c_ages[19] 1704.33    1.33 76.74 1556.44 1651.91 1703.20 1756.41 1857.89
## c_ages[20] 1776.97    1.24 73.97 1628.28 1728.82 1777.19 1826.39 1917.93
## c_ages[21] 1851.12    0.96 60.49 1736.13 1810.63 1848.90 1888.59 1974.95
##            n_eff Rhat
## c_ages[1]   4000  NaN
## c_ages[2]   3042    1
## c_ages[3]   3106    1
## c_ages[4]   3059    1
## c_ages[5]   4000    1
## c_ages[6]   2594    1
## c_ages[7]   3217    1
## c_ages[8]   3245    1
## c_ages[9]   3375    1
## c_ages[10]  2759    1
## c_ages[11]  3182    1
## c_ages[12]  3463    1
## c_ages[13]  4000    1
## c_ages[14]  3412    1
## c_ages[15]  3239    1
## c_ages[16]  3257    1
## c_ages[17]  4000    1
## c_ages[18]  3502    1
## c_ages[19]  3343    1
## c_ages[20]  3549    1
## c_ages[21]  4000    1
## 
## Samples were drawn using NUTS(diag_e) at Tue Apr  4 17:06:12 2017.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

