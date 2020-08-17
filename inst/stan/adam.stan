// Bacon as a Stan model
// 03.08.2020 Andrew Dolman
// Add flexible addtional layer to alphas.

data {
  int<lower=0, upper=1> inflate_errors;
  int<lower=0> N;
  int<lower=0> K_fine;  // no of fine sections
  int<lower=0> K_tot;  // total no of gamma parameters
  int<lower=0> nu; // degrees of freedom of t error distribution
  vector[N] depth;
  vector[N] obs_age;
  vector[N] obs_err;
  vector[K_fine] c_depth_bottom;
  vector[K_fine] c_depth_top;
  int parent[K_tot]; // index fine sections to their parent coarse sections
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

  int<lower = 1> first_K_fine = K_tot - K_fine+1;
  
  real<lower=0> mean_obs_err = mean(obs_err);

}
parameters {
  real<lower = 0, upper = 1> R;
  vector<lower = 0>[K_tot] alpha;
  real age0;
  //real<lower = 0> record_acc_mean;
  real<lower = 0> shape;
  //vector<lower = 0>[K_tot] shape;
  //vector<lower=0>[K1] section_acc_mean;

  // measurement error inflation factor
  // these have length 0 if inflate_errors == 0
  // parameters are in scope so model runs, but zero length so nothing sampled
  real<lower = 0> infl_mean[inflate_errors];
  real<lower = 0> infl_shape[inflate_errors];
  vector<lower = 0>[inflate_errors ? N : 0] infl;
  //vector<lower = 0>[inflate_errors ? 1 : 0] infl_sigma;
  real<lower = 0> infl_sigma[inflate_errors];
 

}
transformed parameters{

  //real<lower = 0> record_acc_mean_beta;

  real<lower = 0, upper = 1> w;

  vector[K_fine] x;
  vector[K_fine+1] c_ages;
  vector[N] Mod_age;
  vector[K_tot] beta;


  // the inflated observation errors
  vector[N] obs_err_infl;
 
  if (inflate_errors == 1){
    
   for (n in 1:N)
   
   obs_err_infl[n] = obs_err[n] + infl_sigma[1] * infl[n];
   //obs_err_infl[n] = obs_err[n] + obs_err[n] * infl[n];
   
  } else {
    obs_err_infl = obs_err;
  }

  beta[1] = shape / alpha[1];
  beta[2:K_tot] = section_acc_shape ./ alpha[parent[2:K_tot]];

  //section_acc_mean_beta = section_acc_shape ./ section_acc_mean;

  w = R^(delta_c);

  // only the "fine" alphas

  x[1] = alpha[first_K_fine];

  // call to sampling function is already vectorised below so there is no
  // advantage to vectorising here
  // x[2:K_fine] = w * x[1:(K_fine-1)] + (1-w) * alpha[2:K_fine];

  for(i in 2:K_fine){
    x[i] = w*x[i-1] + (1-w)*alpha[i + first_K_fine -1];
  }

  // Get the Mod_ages
  c_ages[1] = age0;
  c_ages[2:(K_fine+1)] = age0 + cumulative_sum(x * delta_c);
  Mod_age = c_ages[which_c] + x[which_c] .* (depth - c_depth_top[which_c]);
}

model {

  // parameters for the hierarchical prior on the section means
  //alpha[1] ~ gamma(record_prior_acc_mean_shape, record_prior_acc_mean_beta);
  
  alpha[1] ~ normal(0, 10*record_prior_acc_mean_mean);
  
  shape ~ gamma(record_prior_acc_shape_shape, record_prior_acc_shape_beta);


  // the Gamma distributed innovations
  // betas have already been indexed to correct parent
  alpha[2:K_tot] ~ gamma(section_acc_shape, beta[2:K_tot]);

  // maybe need to model the shapes also - at the moment only the first shape is modelled
  // the rest are fixed

  // the memory parameters
  R ~ beta(mem_alpha, mem_beta);

  // the observation error inflation model
  if (inflate_errors == 1){
    infl_mean ~ gamma(1, 1);
    infl_shape ~ gamma(1, 1);
    infl ~ gamma(infl_shape[1], infl_shape[1] / infl_mean[1]);
    
    infl_sigma ~ normal(0, mean_obs_err);
  }

  obs_age ~ student_t(nu, Mod_age, obs_err_infl);

}
