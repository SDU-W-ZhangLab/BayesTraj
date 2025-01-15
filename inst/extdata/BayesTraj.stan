data {
  
  int<lower = 1> N; 
  int<lower = 1> n_branch; 
  int<lower = 1> G;
  
  int b1_G_s;
  int b1_G_t;
  
  int b2_G_s;
  int b2_G_t;
  
  int b3_G_s;
  int b3_G_t;
  
  int b4_G_s;
  int b4_G_t;
  
  vector[G] lst1;  
  vector[G] lst2;
  vector[G] lst3;
  vector[G] lst4;
 
  matrix[N, G] y;
  
  real b1_mu0_s_mu[b1_G_s];
  real b1_mu0_s_sd[b1_G_s];
  real b2_mu0_s_mu[b2_G_s];
  real b2_mu0_s_sd[b2_G_s];
  real b3_mu0_s_mu[b3_G_s];
  real b3_mu0_s_sd[b3_G_s];
  real b4_mu0_s_mu[b4_G_s];
  real b4_mu0_s_sd[b4_G_s];
  
  real b1_mu0_t_mu[b1_G_t];
  real b1_mu0_t_sd[b1_G_t];
  real b2_mu0_t_mu[b2_G_t];
  real b2_mu0_t_sd[b2_G_t];
  real b3_mu0_t_mu[b3_G_t];
  real b3_mu0_t_sd[b3_G_t];
  real b4_mu0_t_mu[b4_G_t];
  real b4_mu0_t_sd[b4_G_t];
  
  real b1_k_s_mu[b1_G_s];
  real b1_k_s_sd[b1_G_s];
  real b1_k_t_mu[b1_G_t];
  real b1_k_t_sd[b1_G_t];
  
  real b2_k_s_mu[b2_G_s];
  real b2_k_s_sd[b2_G_s];
  real b2_k_t_mu[b2_G_t];
  real b2_k_t_sd[b2_G_t];
  
  real b3_k_s_mu[b3_G_s];
  real b3_k_s_sd[b3_G_s];
  real b3_k_t_mu[b3_G_t];
  real b3_k_t_sd[b3_G_t];
  
  real b4_k_s_mu[b4_G_s];
  real b4_k_s_sd[b4_G_s];
  real b4_k_t_mu[b4_G_t];
  real b4_k_t_sd[b4_G_t];
  
  real lower_mu0_s;
  real lower_mu0_t;
  real upper_mu0_s;
  real upper_mu0_t;
  
  real lower_k_s;
  real lower_k_t;
  real lower_t0;
  real upper_t0;
  real upper_b1_t0;
  real lower_b2_t0;
  
  real alpha;
  real eplison;
}
parameters {
	
	vector<lower = lower_mu0_s,upper = upper_mu0_s>[b1_G_s] b1_mu0_s;
	vector<lower = lower_mu0_s,upper = upper_mu0_s>[b2_G_s] b2_mu0_s;
	vector<lower = lower_mu0_s,upper = upper_mu0_s>[b3_G_s] b3_mu0_s;
	vector<lower = lower_mu0_s,upper = upper_mu0_s>[b4_G_s] b4_mu0_s;
	
	vector<lower = lower_mu0_t,upper = upper_mu0_t>[b1_G_t] b1_mu0_t;
	vector<lower = lower_mu0_t,upper = upper_mu0_t>[b2_G_t] b2_mu0_t;
	vector<lower = lower_mu0_t,upper = upper_mu0_t>[b3_G_t] b3_mu0_t;
	vector<lower = lower_mu0_t,upper = upper_mu0_t>[b4_G_t] b4_mu0_t;
	
	vector<lower = lower_k_s>[b1_G_s] b1_k_s;
	vector<lower = lower_k_t>[b1_G_t] b1_k_t;
	vector<lower = lower_k_s>[b2_G_s] b2_k_s;
	vector<lower = lower_k_t>[b2_G_t] b2_k_t;
	vector<lower = lower_k_s>[b3_G_s] b3_k_s;
	vector<lower = lower_k_t>[b3_G_t] b3_k_t;
	vector<lower = lower_k_s>[b4_G_s] b4_k_s;
	vector<lower = lower_k_t>[b4_G_t] b4_k_t;
	
	vector<lower = lower_t0, upper = upper_b1_t0>[b1_G_s] b1_t0;
	vector<lower = lower_b2_t0, upper = upper_t0>[b1_G_t] b2_t0;
	
	vector<lower = lower_t0, upper = upper_t0>[G] new_t01;
	vector<lower = lower_t0, upper = upper_t0>[G] new_t02;
	vector<lower = lower_t0, upper = upper_t0>[G] new_t03;
	vector<lower = lower_t0, upper = upper_t0>[G] new_t04;

	real<lower = 0> mu_hyper;
	real<lower = 0, upper = 1> t[N]; 
	real<lower = 0> phi[1]; 
	vector[2] beta;

	simplex[n_branch] theta;

}
transformed parameters {

	matrix[n_branch, G] all_mu0;
	matrix[n_branch, G] all_k;
	matrix[n_branch, G] all_t0;
	
	matrix[N, G] muu;
	matrix<lower = 0>[N, G] sdd;

	vector<lower = lower_t0, upper = upper_t0>[G] b_t0;
	
	for(i in 1:G){
		int S = 1;
		int T = 1;
		if(lst1[i] == 0){
			b_t0[i] = b1_t0[S];
			S += 1;
		}else{
			b_t0[i] = b2_t0[T];
			T += 1;
		}
	}

	array[N] vector[n_branch] branch;
	
	for(i in 1:N) {
	  for(j in 1:n_branch){
		  branch[i][j] =0;
	  }
	}

	for(l in 1:n_branch){
	
		int S = 1;
		int T = 1;
		vector[G] mu0;
		vector[G] k;
		vector[G] t0;
		vector[G] type;

		if(l == 1){
			type = lst1;
			t0 = new_t01;
			
			for(i in 1:G){
				if(lst1[i] == 0){
					mu0[i] = b1_mu0_s[S];
					k[i] = b1_k_s[S];
					S += 1;
				}else{
					mu0[i] = b1_mu0_t[T];
					k[i] = b1_k_t[T];
					T += 1;
				}
			}
			
			
		}else if(l == 2){
			type = lst2;
			t0 = new_t02;

			for(i in 1:G){
				if(lst2[i] == 0){
					mu0[i] = b2_mu0_s[S];
					k[i] = b2_k_s[S];
					S += 1;
				}else{
					mu0[i] = b2_mu0_t[T];
					k[i] = b2_k_t[T];
					T += 1;
				}
			}

		}else if(l == 3){
			type = lst3;
			t0 = new_t03;

			for(i in 1:G){
				if(lst3[i] == 0){
					mu0[i] = b3_mu0_s[S];
					k[i] = b3_k_s[S];
					S += 1;
				}else{
					mu0[i] = b3_mu0_t[T];
					k[i] = b3_k_t[T];
					T += 1;
				}
			}
		}else if(l == 4){
			type = lst4;
			t0 = new_t04;
			
			for(i in 1:G){
				if(lst4[i] == 0){
					mu0[i] = b4_mu0_s[S];
					k[i] = b4_k_s[S];
					S += 1;
				}else{
					mu0[i] = b4_mu0_t[T];
					k[i] = b4_k_t[T];
					T += 1;
				}
			}
		}
	
	for(i in 1:N) {
		for(g in 1:G){
			if(type[g] == 0){
				muu[i,g] = 2 * mu0[g] / (1 + exp(- k[g] * (t[i] - t0[g])));
				sdd[i,g] = sqrt( alpha * (1 + phi[1]) * muu[i,g] + eplison);
			}else{
				muu[i,g] = 2 * mu0[g] * exp(- 10 * k[g] * (t[i] - t0[g])^2);
				sdd[i,g] = sqrt( alpha * (1 + phi[1]) * muu[i,g] + eplison);
			}
		}
	}
		
	all_mu0[l] = to_row_vector(mu0);
	all_k[l] = to_row_vector(k);
	all_t0[l] = to_row_vector(t0);
	
//likelihood

	for (i in 1:N) {
		
		for(g in 1:G){
			if(y[i,g] == 0) {
				branch[i][l] += log_sum_exp(bernoulli_logit_lpmf(1 | beta[1] + beta[2] * muu[i,g]),
									  bernoulli_logit_lpmf(0 | beta[1] + beta[2] * muu[i,g]) + 
									  normal_lpdf(y[i,g] | muu[i,g], sdd[i,g])); 
			} else {
				branch[i][l] += bernoulli_logit_lpmf(0 | beta[1] + beta[2] * muu[i,g]) + 
								normal_lpdf(y[i,g] | muu[i,g], sdd[i,g]);
			}
		}
	
	}
}


	array[N] vector[n_branch] log_all;
		
	vector[n_branch] lps;
	for (i in 1:N) {
		lps = log(theta);
		for (l in 1:n_branch){
			lps[l] += branch[i][l];
		}
		log_all[i] = lps;
	}
}
model {
	t ~ uniform(0,1);
	
	b1_mu0_s ~ normal(b1_mu0_s_mu,b1_mu0_s_sd);
	b2_mu0_s ~ normal(b2_mu0_s_mu,b2_mu0_s_sd);
	b3_mu0_s ~ normal(b3_mu0_s_mu,b3_mu0_s_sd);
	b4_mu0_s ~ normal(b4_mu0_s_mu,b4_mu0_s_sd);
	
	b1_mu0_t ~ normal(b1_mu0_t_mu,b1_mu0_t_sd);
	b2_mu0_t ~ normal(b2_mu0_t_mu,b2_mu0_t_sd);
	b3_mu0_t ~ normal(b3_mu0_t_mu,b3_mu0_t_sd);
	b4_mu0_t ~ normal(b4_mu0_t_mu,b4_mu0_t_sd);
	
	
	b1_k_s ~ normal(b1_k_s_mu, b1_k_s_sd);
	b1_k_t ~ normal(b1_k_t_mu, b1_k_t_sd);
	b2_k_s ~ normal(b2_k_s_mu, b2_k_s_sd);
	b2_k_t ~ normal(b2_k_t_mu, b2_k_t_sd);
	b3_k_s ~ normal(b3_k_s_mu, b3_k_s_sd);
	b3_k_t ~ normal(b3_k_t_mu, b3_k_t_sd);
	b4_k_s ~ normal(b4_k_s_mu, b4_k_s_sd);
	b4_k_t ~ normal(b4_k_t_mu, b4_k_t_sd);
	
	
	b1_t0 ~ beta(5,5);	
	b2_t0 ~ beta(5,5);	

	phi ~ gamma(12, 4);
	beta ~ normal(0, 0.1);
	theta ~ dirichlet(rep_vector(1000.0,n_branch));

    for(i in 1:G){
		new_t01[i] ~ normal(b_t0[i],0.01);
		new_t02[i] ~ normal(b_t0[i],0.01);
		new_t03[i] ~ normal(b_t0[i],0.01);
		new_t04[i] ~ normal(b_t0[i],0.01);
	}

	for (n in 1:N) {
		target += log_sum_exp(log_all[n]);
	}

}
generated quantities {
	array[N] vector[n_branch] category_probs;
	for (n in 1:N) {
		vector[n_branch] lps1 = log_all[n];
		lps1 = lps1 - max(lps1); 
		category_probs[n] = exp(lps1) / sum(exp(lps1)); 
	}
}

