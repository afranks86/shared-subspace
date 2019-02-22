// -*- mode: C -*-
data {

  int<lower=1> M;  // Length Y - 1 (e.g. number of angles)
  real<lower=0> lam[M+1];  // lambda's are quadratic coefficients
  real gamma[M+1];  // linear coefficients 
  
}
parameters {
  real<lower=0, upper=pi()> phi[M-1];
  real<lower=0, upper=2*pi()> phiLast;
}
transformed parameters {

  real<lower=-1, upper=1> Y[M+1];
  real ym;
  real ylast;

  for(i in 1:(M-1) ){
    real yi;
    yi <- 1;
    for( j in 1:(i-1) ) {
      yi <- yi*sin(phi[j]);
    }
    Y[i]  <- yi*cos(phi[i]);
  }

  ym <- 1;
  for( j in 1:(M-1) ) {
    ym <- ym*sin(phi[j]);
  }
  Y[M]  <- ym*cos(phiLast);
  Y[M+1] <- ym*sin(phiLast);
  
}
model {
  for(i in 1:(M-1)) {
    increment_log_prob(lam[i]*Y[i]^2+gamma[i]*Y[i]+(M-i)*log(sin(phi[i])));
  }
  increment_log_prob(lam[M]*Y[M]^2+gamma[M]*Y[M]);
  increment_log_prob(lam[M+1]*Y[M+1]^2+gamma[M+1]*Y[M+1]);
}
