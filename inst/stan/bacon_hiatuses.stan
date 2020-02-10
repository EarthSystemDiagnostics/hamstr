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
  vector[K] delta_c;
  
  // hiatus data
  int<lower=0> N_hiatus;
  vector[N_hiatus] hiatus_depth;
  vector[N_hiatus] hiatus_length;
  vector[K] hiatus_factor;
  real<lower = 0> hiatus_shape;
  vector[N_hiatus] hiatus_beta;
  int which_c_hiatus[N_hiatus];


  // hyperparameters
  real<lower = 0> acc_mean;

  real<lower = 0> acc_alpha;
  real<lower = 0> acc_beta;

  real<lower = 0> mem_mean;
  real<lower = 0> mem_strength;
}
transformed data{
// transformed hyperparameters
  real<lower=0> mem_alpha;
  real<lower=0> mem_beta;

// depths corresponding to c_ages
//  vector[K+1] c_depths;

  mem_alpha = mem_strength * mem_mean;
  mem_beta = mem_strength * (1-mem_mean);
  
//  c_depths = append_row(c_depth_top[1], c_depth_bottom);

}
parameters {
  real<lower = 0, upper = 1> R;
  vector<lower = 0>[K] alpha;
  real<lower = 0> age0;
  vector<lower = 0>[N_hiatus] hiatus;
  
}
transformed parameters{

  vector<lower = 0, upper = 1>[K] w;

  vector[K] x;
  vector[K] hiatus_vec = rep_vector(0, K);
  vector[K+1] c_ages;
  vector[N] Mod_age;

 for(i in 1:K){
  w[i] = R^delta_c[i];
 }
  // pretty sure now this is the correct scaling 
  // if the correlation between 2 1cm layers is 0.5
  // the correlation between 2 0.1cm layers would be 0.5^0.1 = 0.933
  // after 10 0.1 cm layers the correlation would be 0.933^10 = 0.5

  x[1] = alpha[1];

  // multiply w by "hiatus_factor" 1 when no hiatus, 0 when hiatus, to break the 
  // autocorrelation between layers when there was an hiatus
  for(i in 2:K){
    x[i] = (w[i]*hiatus_factor[i])*x[i-1] + (1-(w[i]*hiatus_factor[1]))*alpha[i];
    }
  
  // call to sampling function is already vectorised below so there is only a small 
  // advantage to vectorising here
  //x[2:K] = w * x[1:(K-1)] + (1-w) * alpha[2:K];

  // Get the Mod_ages
  c_ages[1] = age0;
  hiatus_vec[which_c_hiatus] = hiatus;
  c_ages[2:(K+1)] = age0 + cumulative_sum(x .* delta_c[1:K]) + cumulative_sum(hiatus_vec);
  Mod_age = c_ages[which_c] + x[which_c] .* (depth - c_depth_top[which_c]);
}
model {
  age0 ~ normal(obs_age[1], 100);
  R ~ beta(mem_alpha, mem_beta);
  alpha ~ gamma(acc_alpha, acc_beta);
  hiatus ~ gamma(hiatus_shape, hiatus_beta);
  //obs_age ~ normal(Mod_age, obs_err);
  obs_age ~ student_t(nu, Mod_age, obs_err);
  }
  
 
  
  
