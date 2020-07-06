// Bacon as a Stan model
// 20.06.2020 Andrew Dolman
// Add addtional layer to alphas.

// ideas for parameter renaming
// acc_mean_record = record_acc_mean
// acc_mean_sections = section_acc_mean

// hyper-hyper-parameters
// acc_mean_hhp = record_prior_acc_mean_mean (mean of the prior on acc_mean)
// acc_mean_shape_hhp = record_prior_acc_mean_shape (shape of the prior on acc_mean)

// acc_shape_mean_hhp = record_prior_acc_shape_mean (mean of the prior on acc_mean)
// acc_shape_shape_hhp = record_prior_acc_shape_shape (shape of the prior on acc_mean)

// hyper-parameters
// acc_shape_sections_hp

// mem_mean_hp
// mem_strength_hp




data {
  int<lower=0, upper=1> inflate_errors;
  int<lower=0> N;
  int<lower=0> K;  // no of fine sections
  int<lower=0> K1; // no of coarse sections
  int<lower=0> nu; // degrees of freedom of t error distribution
  vector[N] depth;
  vector[N] obs_age;
  vector[N] obs_err;
  vector[K] c_depth_bottom;
  vector[K] c_depth_top;
  int whichK1[K]; // index fine sections to their parent coarse sections
  int which_c[N]; // index observations to their fine section
  real<lower = 0> delta_c; // width of each fine section
  
  // hyperparameters
  
  // parameters for the prior distribution of the overall mean acc rate
  real<lower = 0> record_prior_acc_mean_mean;
  real<lower = 0> record_prior_acc_mean_shape;
  
  // parameters for the prior distribution of the shape of the distribution of section means
  real<lower = 0> record_prior_acc_shape_mean;
  real<lower = 0> record_prior_acc_shape_shape;
  
  real<lower = 0> section_acc_shape; // shape of the within section innovations
  
  real<lower = 0> mem_mean;
  real<lower = 0> mem_strength;
  
}
transformed data{
  
  // transform shape and mean to shape and beta (rate) of gamma dist
  real<lower=0> record_prior_acc_mean_beta = record_prior_acc_mean_shape / record_prior_acc_mean_mean;
  real<lower=0> record_prior_acc_shape_beta  = record_prior_acc_shape_shape / record_prior_acc_shape_mean;
  
  
  // transform mean and strength of memory beta distribution to alpha and beta 
  real<lower=0> mem_alpha = mem_strength * mem_mean;
  real<lower=0> mem_beta = mem_strength * (1-mem_mean);
  
  
}
parameters {
  real<lower = 0, upper = 1> R;
  vector<lower = 0>[K] alpha;
  real age0;
  real<lower = 0> record_acc_mean;
  real<lower = 0> record_acc_shape;
  vector<lower=0>[K1] section_acc_mean;
  
  // measurement error inflation factor
  // these have length 0 if inflate_errors == 0
  // parameters are in scope so model runs, but zero length so nothing sampled
  real<lower = 0> infl_mean[inflate_errors];
  real<lower = 0> infl_shape[inflate_errors];
  vector<lower = 0>[inflate_errors ? N : 0] infl;
  
}
transformed parameters{
  
  real<lower = 0> record_acc_mean_beta;
  
  real<lower = 0, upper = 1> w;
  
  vector[K] x;
  vector[K+1] c_ages;
  vector[N] Mod_age;
  vector[K1] section_acc_mean_beta;
  
  
  // the inflated observation errors
  vector[N] obs_err_infl;
  
  if (inflate_errors == 1){
    obs_err_infl = obs_err + infl .* obs_err;
  } else {
    obs_err_infl = obs_err;
  }
  
  
  record_acc_mean_beta = record_acc_shape / record_acc_mean;
  
  section_acc_mean_beta = section_acc_shape ./ section_acc_mean;
  
  w = R^(delta_c);
  
  x[1] = alpha[1];
  
  // call to sampling function is already vectorised below so there is no
  // advantage to vectorising here
  // x[2:K] = w * x[1:(K-1)] + (1-w) * alpha[2:K];
  
  for(i in 2:K){
    x[i] = w*x[i-1] + (1-w)*alpha[i];
  }
  
  // Get the Mod_ages
  c_ages[1] = age0;
  c_ages[2:(K+1)] = age0 + cumulative_sum(x * delta_c);
  Mod_age = c_ages[which_c] + x[which_c] .* (depth - c_depth_top[which_c]);
}

model {
  
  // parameters for the hierarchical prior on the section means
  record_acc_mean ~ gamma(record_prior_acc_mean_shape, record_prior_acc_mean_beta);
  record_acc_shape ~ gamma(record_prior_acc_shape_shape, record_prior_acc_shape_beta);
  
  // the K1 section means, parametrised by the record mean and a user specified shape (alpha)
  section_acc_mean ~ gamma(record_prior_acc_mean_shape, record_acc_mean_beta);
  
  // the Gamma distributed innovations
  // parametrised from the section mean acc rates and a user supplied shape (alpha)
  alpha ~ gamma(section_acc_shape, section_acc_mean_beta[whichK1]);
  
  // loosely model the first age as being anywhere between 0 and the youngest data point
  //age0 ~ normal(min(obs_age), 100);
  
  // the memory parameters
  R ~ beta(mem_alpha, mem_beta);
  
  // the observation error inflation model
  if (inflate_errors == 1){
    infl_mean ~ gamma(1.5, 1);
    infl_shape ~ gamma(1.5, 1);
    infl ~ gamma(infl_shape[1], infl_shape[1] / infl_mean[1]);
  } 
  
  obs_age ~ student_t(nu, Mod_age, obs_err_infl);
  
  
}
