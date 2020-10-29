// saved as multi-curves-general.stan
// fit a linear mixed effects model with subject specific curves y = X*beta + Z*u + Xs*b1 + ... + Xs*bJ + Zs*v1 + ... + Zs*vJ
data {
  int<lower = 1> J;		/* total number of subjects */
  int<lower = 1> N;		/* total number of samples */
  int<lower = 1> ncol_X;	/* number of column of X design matrix */
  int<lower = 1> ncol_Z;	/* number of column of Z design matrix */
  int<lower = 1> ncol_Xs;	/* number of column of Xs design matrix */
  int<lower = 1> ncol_Zs;	/* number of column of Zs design matrix */

  real y[N];			    /* response */
  int<lower = 1, upper = J> sub[N]; /* subject index */
  matrix[N, ncol_X] X;		    /* design mat X */
  matrix[N, ncol_Z] Z;		    /* design mat Z */
  matrix[N, ncol_Xs] Xs;		    /* subject-specific design mat Xs */
  matrix[N, ncol_Zs] Zs;		    /* subject-specific design mat Zs */
}

parameters {
  vector[ncol_X] beta;		/* population X coefs */
  vector[ncol_Z] u;		/* population Z coefs */
  vector[ncol_Xs] b[J];		/* individual X coefs */
  vector[ncol_Zs] v[J];		/* individual Z coefs */

  /* /\* vanilla correlation matrix *\/ */
  /* corr_matrix[ncol_Xs] cor_b;	/\* individual X corr *\/ */

  /* scaled cholesky correlation matrix */
  /* cholesky_factor_corr[ncol_Xs] L_cor_b;	/\* individual X cholesky corr *\/ */
  /* vector<lower = 0>[ncol_Xs] scl_b;	/\* individual X scale of corr *\/ */

  /* vanilla covariance matrix */
  cov_matrix[ncol_Xs] cov_b;	/* individual X cov */
  
  /* /\* scaled covariance matrix *\/ */
  /* cov_matrix[ncol_Xs] cov_b;	   /\* individual X cov *\/ */
  /* vector<lower = 0>[ncol_Xs] scl_b; /\* individual X scale of cov *\/ */
  
  /* real<lower = 0> sig_u;		    /\* population Z sd *\/ */
  /* real<lower = 0> sig_v;		    /\* individual Z sd *\/ */
  /* real<lower = 0> sig_eps;		    /\* Gaussian error sd *\/ */
  real<lower = 0> sig2_u;		    /* population Z var */
  real<lower = 0> sig2_v;		    /* individual Z var */
  real<lower = 0> sig2_eps;		    /* Gaussian error var */
}

transformed parameters {
/*   /\* real<lower = 0> sig2_eps = exp(log_sig2_eps); /\\* error variance *\\/ *\/ */
  real<lower = 0> sig_u = sqrt(sig2_u);	    /* population Z sd */
  real<lower = 0> sig_v = sqrt(sig2_v);	    /* individual Z sd */
  real<lower = 0> sig_eps = sqrt(sig2_eps); /* error sd */
}

model {
  sig2_u ~ inv_gamma(0.001, 0.001);
  sig2_v ~ inv_gamma(0.001, 0.001);
  sig2_eps ~ inv_gamma(0.001, 0.001);

  /* prior for cov_b */
  cov_b ~ inv_wishart(ncol_Xs + 1, diag_matrix(rep_vector(1, ncol_Xs)));

  /* /\* prior for L_cor_b *\/ */
  /* L_cor_b ~ lkj_corr_cholesky(1); */

  /* prior for scl_b */
  /* uniform */
  
  /* prior for beta */
  /* uniform */

  /* prior for u */
  u ~ normal(0, sig_u);
  
  /* prior for b */
  /* b ~ multi_normal_cholesky(rep_vector(0, ncol_X), diag_pre_multiply(scl_b, L_cor_b)); */
  b ~ multi_normal(rep_vector(0, ncol_Xs), cov_b);

  /* prior for v */
  for (j in 1:J) {
    v[j] ~ normal(0, sig_v);
  }

  /* likelihood */
  for (n in 1:N) {
    y[n] ~ normal(X[n] * beta + Xs[n] * b[sub[n]] +
  		  Z[n] * u + Zs[n] * v[sub[n]], sig_eps);
  }
}


