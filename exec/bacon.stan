// Bacon as a Stan model
// 25.03.2017 Andrew Dolman
data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> nu;
  vector[N] depth;
  vector[N] obs_age;
  vector[N] obs_err;
  vector[K] c_depth_bottom;
  vector[K] c_depth_top;
  int which_c[N];
  real<lower = 0> delta_c;

  // hyperparameters
  real<lower = 0> acc_mean;
  real<lower = 0> acc_var;

  real<lower = 0> mem_mean;
  real<lower = 0> mem_strength;
}
transformed data{
// transformed hyperparameters
  real<lower=0> acc_alpha;
  real<lower=0> acc_beta;

  real<lower=0> mem_alpha;
  real<lower=0> mem_beta;

// depths corresponding to c_ages
//  vector[K+1] c_depths;

  acc_alpha = acc_mean^2 / acc_var;
  acc_beta = acc_mean / acc_var;

  mem_alpha = mem_strength * mem_mean;
  mem_beta = mem_strength * (1-mem_mean);
  
//  c_depths = append_row(c_depth_top[1], c_depth_bottom);

}
parameters {
  real<lower = 0, upper = 1> R;
  vector<lower = 0>[K] alpha;
  real<lower = 0> age0;
}
transformed parameters{

  real<lower = 0, upper = 1> w;

  vector[K] x;
  vector[K+1] c_ages;
  vector[N] Mod_age;

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
  age0 ~ normal(obs_age[1], 100);
  R ~ beta(mem_alpha, mem_beta);
  alpha ~ gamma(acc_alpha, acc_beta);
  //obs_age ~ normal(Mod_age, obs_err);
  obs_age ~ student_t(nu, Mod_age, obs_err);
  }
