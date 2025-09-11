data {
  // data dimensions
  int<lower=1> N_in; // obs
  int<lower=1> N_out; // obs
  int<lower=1> K; // choices

  // data
  int<lower=1,upper=K> y_in[N_in];  // observed choice

  // fixed effects
  int<lower=1> Xw;  // num fixed effects
  int<lower=1> Xn_in;  // num non-zero pairs
  int<lower=1,upper=N_in> X1_in[Xn_in];  // row of non-zero
  int<lower=1,upper=Xw> X2_in[Xn_in]; // col of non-zero
  int<lower=1> Xn_out;
  int<lower=1,upper=N_out> X1_out[Xn_out];
  int<lower=1,upper=Xw> X2_out[Xn_out];

  // year
  vector[N_in] year_in;
  vector[N_out] year_out;

  // random effect dimensions
  int<lower=1> S; // num states
  int<lower=1> R; // num races
  int<lower=1> P; // num places

  // random effect data
  int<lower=1,upper=S> state_in[N_in];  // states
  int<lower=1,upper=S> state_out[N_out];  // states
  int<lower=1,upper=R> race_in[N_in];   // races
  int<lower=1,upper=R> race_out[N_out];   // races
  int<lower=1,upper=P> place_in[N_in];  // places
  int<lower=1,upper=P> place_out[N_out];  // places
}

transformed data {
  int Km1;
  Km1 <- K - 1;
}

parameters {
  // intercepts
  row_vector[Km1] alpha;

  // fixed effects
  vector[Km1] beta0[Xw];  // causes
  vector[Km1] beta1;      // year

  // random effects
  real<lower=1e-5> sigma_S[Km1];
  vector[S] pi_S[Km1];  // state
  real<lower=1e-5> sigma_R[Km1];
  vector[R] pi_R[Km1];  // race
  real<lower=1e-5> sigma_P[Km1];
  vector[P] pi_P[Km1];  // place
}

model {
  // predicted probabilities of each
  vector[K] mu_in[N_in];

  // priors
  alpha ~ normal(0.0, 1.0E3);
  for (k in 1:Km1) {
    beta0[k] ~ normal(0.0, 1.0);
    pi_S[k] ~ normal(0.0, sigma_S[k]);
    pi_R[k] ~ normal(0.0, sigma_R[k]);
    pi_P[k] ~ normal(0.0, sigma_P[k]);
  }
  beta1 ~ normal(0.0, 1.0);

  // predictions
  {
    matrix[N_in,Km1] bX;
    int s;
    int r;
    int p;
    int x1;
    int x2;
    vector[4] summand;
    // fill in fixed effects
      bX <- rep_matrix(alpha, N_in) + (year_in * beta1');
      for (x in 1:Xn_in) {
        x1 <- X1_in[x];
        x2 <- X2_in[x];
        for (k in 1:Km1) {
          bX[x1,k] <- bX[x1,k] + beta0[x2,k];
        }
      }
    // fill in predictions
      for (n in 1:N_in) {
        s <- state_in[n];
        r <- race_in[n];
        p <- place_in[n];
        mu_in[n,1] <- 0;
        for (k in 1:Km1) {
          summand[1] <- bX[n,k];
          summand[2] <- pi_S[k,s];
          summand[3] <- pi_R[k,r];
          summand[4] <- pi_P[k,p];
          mu_in[n,k+1] <- sum(summand);
        }
      }
  }

  // data likelihood
  for (n in 1:N_in) {
    y_in[n] ~ categorical_logit(mu_in[n]);
    #y[n] ~ categorical(softmax(mu[n]));
  }
}

generated quantities {
  vector[K] mu_out[N_out];
  {
    matrix[N_out,Km1] bX;
    int s;
    int r;
    int p;
    int x1;
    int x2;
    vector[4] summand;
    vector[K] tmp;
    // fill in fixed effects
      bX <- rep_matrix(alpha, N_out) + (year_out * beta1');
      for (x in 1:Xn_out) {
        x1 <- X1_out[x];
        x2 <- X2_out[x];
        for (k in 1:Km1) {
          bX[x1,k] <- bX[x1,k] + beta0[x2,k];
        }
      }
    // fill in predictions
      for (n in 1:N_out) {
        s <- state_out[n];
        r <- race_out[n];
        p <- place_out[n];
        tmp[1] <- 0;
        for (k in 1:Km1) {
          summand[1] <- bX[n,k];
          summand[2] <- pi_S[k,s];
          summand[3] <- pi_R[k,r];
          summand[4] <- pi_P[k,p];
          tmp[k+1] <- sum(summand);
        }
        mu_out[n] <- softmax(tmp);
      }
  }
}
