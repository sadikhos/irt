##########################################################################################
#
#   Title:  Longitudinal analysis of the simulated IRT data using total scores
#           and IRT models
#   Author: Shamil Sadikhov
#   Date:   Dec 2014 
#   Comments:
#
##########################################################################################

source("irt_func.R")

set.seed = "12345"

require(nlme)
require(rstan)
require(reshape2)
require(abind)
require(mvtnorm)
# require(varComp)
require(latticeExtra)
require(coda)
# require(glmer2stan)

stdize <- function(x, narm = TRUE) (x - mean(x, na.rm = narm))/sd(x, na.rm = narm)

##########################################################################################
# TOTAL SCORES
# AR(1) covariance, population-averaged model
# Delta - treatment effect at each time point
# Data are in wide format, 2D
##########################################################################################

irt_pa_totscr_delta1 <- '
data {
     int<lower=0> T;                // number of categorical time points
     int<lower=0> N;                // number of subjects
     vector[T] Y[N];                  // data of total scores
     int trt[N];                    // treatment arm
     int n_trt;                     // number of arms
}

parameters {
     real<lower=0.0, upper=1.0> rho;        // Correlation parameter
     vector[T] mu_y[n_trt];

}

transformed parameters {
}

model{
   real sigsq_y;
   matrix[T, T] var_y;
   matrix[T, T] L_y;
   real m_mu_y;
   real var_mu_y;


// hyperparameters for mean of total score at each time point

    m_mu_y <- 0;
    var_mu_y <- 1;
    for (t in 1:T)
       for (j in 1:n_trt)
          mu_y[j, t] ~ normal(m_mu_y, var_mu_y);

// AR(1) structure for the covariance of thetas
   sigsq_y <- 1.0;
   var_y[1,1] <- sigsq_y;

   for (t in 2:T) {
      var_y[t,t] <- sigsq_y;
      for (j in 1:(t-1)) {
        var_y[t,j] <- sigsq_y*pow(rho, t-j);
        var_y[j,t] <- var_y[t,j];
      }
   }

   L_y <- cholesky_decompose(var_y);

   for (i in 1:N)
      Y[i] ~ multi_normal_cholesky(mu_y[trt[i]], L_y);
     //  print(Y[1]);

}

generated quantities {
    vector[T-1] deltaout[n_trt-1];             // Treatment effect
// Code the treatment effect - difference in differences
   for (t in 1:(T-1))
      for (l in 1:(n_trt - 1))
         deltaout[l, t] <-  mu_y[l+1, t+1] - mu_y[1, t+1] - (mu_y[l+1, 1] - mu_y[1, 1]);              
}
'

#compile the model 
c_irt_pa_totscr_delta1 <- stan_model(model_code = 'irt_pa_totscr_delta1')


##########################################################################################
# TOTAL SCORES
# Random coefficients model
# Delta - treatment effect at each time point
# Data are in long format
##########################################################################################

irt_rc_totscr_delta1 <- '
data{
    int N;
    real Y[N];
    real trt[N];
    real time[N];
    int Subject[N];
    real trt_X_time[N];
    int N_Subject;
}

transformed data{
    vector[2] zeros_Subject;
    for ( i in 1:2 ) zeros_Subject[i] <- 0;
}

parameters{
    real Intercept;
    real beta_trt;
    real beta_time;
    real beta_trt_X_time;
    real<lower=0> sigma;
    vector[2] vary_Subject[N_Subject];
    corr_matrix[2] Omega;    // correlation mat for random int an slope
    vector<lower=0>[2] Subject_sd;  // subject sds
}

transformed parameters {
    matrix[2, 2] D;
    D <- diag_matrix(Subject_sd);
}

model{
    real vary[N];
    real glm[N];
    matrix[2,2] L;
    matrix[2,2] DL;
    matrix[2,2] Sigma_Subject; 

    // Priors
    Intercept ~ normal(0, 100);
    beta_trt ~ normal(0, 100);
    beta_time ~ normal(0, 100);
    beta_trt_X_time ~ normal(0, 100);
    sigma ~ uniform(0, 100);
    Omega ~ lkj_corr(2.0);
    Subject_sd ~ cauchy(0,2);

    L <- cholesky_decompose(Omega);
    DL <- D * L;


    // Random effects
    for ( j in 1:N_Subject ) vary_Subject[j] ~ multi_normal_cholesky(zeros_Subject, DL);
    // HLM
    for ( i in 1:N ) {
        vary[i] <- vary_Subject[Subject[i],1]
                    + vary_Subject[Subject[i],2] * time[i];
        glm[i] <- vary[i] + Intercept
                   + beta_trt * trt[i]
                   + beta_time * time[i]
                   + beta_trt_X_time * trt_X_time[i];
    }
    Y ~ normal( glm , sigma );
}

generated quantities{
    real deltaout;
    real dev;
    real vary[N];
    real glm[N];
    cov_matrix[2] Sigma_Subject;
    dev <- 0;
    for ( i in 1:N ) {
        vary[i] <- vary_Subject[Subject[i],1]
                    + vary_Subject[Subject[i],2] * time[i];
        glm[i] <- vary[i] + Intercept
                   + beta_trt * trt[i]
                   + beta_time * time[i]
                   + beta_trt_X_time * trt_X_time[i];
        dev <- dev + (-2) * normal_log( Y[i] , glm[i] , sigma );
    }
    deltaout <- beta_trt_X_time*4 - beta_trt;

  Sigma_Subject <- D * Omega * D;
}
'

#compile the model 
c_irt_rc_totscr_delta1 <- stan_model(model_code = 'irt_rc_totscr_delta1')



##########################################################################################
# Multiple Group MLIRT model, AR(1) covariance, population-averaged model
# Delta - treatment effect at each time point
##########################################################################################

mlirt_pa_bin_delta1 <- '
data {
  int n_trt;                     // number of arms
  real m_mu_theta;               // prior for mean of theta
  real var_mu_theta;             // prior for variance of mean of theta
  int<lower=0> T;                // number of categorical time points
  int<lower=0> N;                // number of subjects
  int<lower=0> J;                // number of items
  int Y[N, J, T];                // data, 3D array of integers
  int trt[N];                    // treatment arm
}

parameters {
    vector[J] kappa;                        // intercept 
    vector[T] theta[N];                     // ability : theta[N,T] - array of vectors
    real<lower=0.0, upper=10.0> alpha[J];   // discrimination
    real<lower=0.0, upper=1.0> rho;         // Correlation parameter
    vector[T] mu_theta_raw[n_trt];
    real <lower=0.0, upper = 5.0> mu_alpha; // Prior for item discrim.
    real mu_kappa;                          // prior for item diff.
    real var_alpha;
    real var_kappa;
}

transformed parameters {
   vector[T] mu_theta[n_trt];
// prior for mu_theta
   mu_theta[1, 1] <- 0.0;  // fix to 0 for identifiability

// set the values of mu_theta before sampling - stan complains otherwise
   for (j in 2:n_trt)
      mu_theta[j, 1] <- mu_theta_raw[j, 1];
   for (t in 2:T)
      for (j in 1:n_trt)
         mu_theta[j, t] <- mu_theta_raw[j, t];

}

model{
   real sigsq_theta;
   matrix[T, T] var_theta;
   matrix[T, T] L_theta;
   real p[N, J, T, 2];
   real Q[N, J, T, 1];
   real m_mu_alpha;
   real m_mu_kappa;
   real var_mu_alpha;
   real var_mu_kappa;

// hyperpriors for item parameters
   m_mu_alpha <- 1;
   m_mu_kappa <- 0;
   var_mu_alpha <- 10000;
   var_mu_kappa <- 10000;

// hyperparameters for item parameters;
   mu_alpha ~ normal(m_mu_alpha, sqrt(var_mu_alpha));
   mu_kappa ~ normal(m_mu_kappa, sqrt(var_mu_kappa));
   var_alpha ~ inv_gamma(1, 1);    // inv gamma hyperprior
   var_kappa ~ inv_gamma(1, 1);

// prior for discrimination, truncated below 0
   alpha ~ normal(mu_alpha, sqrt(var_alpha));

// prior for item difficulty
   kappa ~ normal(mu_kappa, sqrt(var_kappa));

 // hyperparameters for ability parameters
    for (t in 1:T)
       for (j in 1:n_trt)
          mu_theta_raw[j, t] ~ normal(m_mu_theta, var_mu_theta);

// AR(1) structure for the covariance of thetas
   sigsq_theta <- 1.0;
   var_theta[1,1] <- sigsq_theta;
   for (t in 2:T) {
      var_theta[t,t] <- sigsq_theta;
      for (j in 1:(t-1)) {
        var_theta[t,j] <- sigsq_theta*pow(rho, t-j);
        var_theta[j,t] <- var_theta[t,j];
      }
   }

   L_theta <- cholesky_decompose(var_theta);
   for (i in 1:N)
      theta[i] ~ multi_normal_cholesky(mu_theta[trt[i]], L_theta);


  for (t in 1:T){
    for (i in 1:N) {
      for (j in 1:J) {
          Y[i, j, t] ~  bernoulli_logit(alpha[j]*(theta[i, t] - kappa[j]));
                     }
                   }
                }
}

generated quantities {
    vector[T-1] deltaout[n_trt-1];             // Treatment effect
// Code the treatment effect - difference in differences
   for (t in 1:(T-1))
      for (l in 1:(n_trt - 1))
         deltaout[l, t] <-  mu_theta[l+1, t+1] - mu_theta[1, t+1] - (mu_theta[l+1, 1] - mu_theta[1, 1]);              
}
'

#compile the model 
c_mlirt_pa_bin_delta1 <- stan_model(model_code = 'mlirt_pa_bin_delta1')



##########################################################################################
#
#      Multiple group MLIRT random coefficients model, binomial data
#
##########################################################################################

mlirt_rc_bin_delta1 <- '
data {
  int n_trt;                     // number of arms
  int<lower=0> T;                // number of categorical time points
  int<lower=0> N;                // number of subjects
  int<lower=0> J;                // number of items
  int Y[N, J, T];                // data, 3D array of integers
  int trt[N];                    // treatment arm coded as integer, 1 = control arm
}

parameters {
    vector[J] kappa;                        // intercept 
    real<lower=0.0, upper=5.0> alpha[J];    // discrimination parameter
 //   real<lower= -1.0, upper= 1.0> rho;    // Correlation parameter between random int and slope
    vector[N] gamma0;                       // random intercept
    vector[N] gamma1;                       // random slope
    vector[n_trt] mu_gamma0_raw;            // intercept prior mean
    vector[n_trt] mu_gamma1;                // slope prior mean
    real <lower=0.0, upper=5.0> var_gamma1; // prior on variance of slope
    real m_mu_gamma0;                       // hyperparameters for random slope mean
    real m_mu_gamma1;
    real <lower=0.0, upper = 5.0> mu_alpha; // Prior for item discrim.
    real mu_kappa;                          // prior for item diff.
    real var_alpha;
    real var_kappa;
}

transformed parameters {
   vector[n_trt] mu_gamma0;
   mu_gamma0[1] <- 0.0;        // fix to 0 for identifiability, intercept, control arm mean

// set the values of mu_gamma0 before sampling - stan complains otherwise
   for (j in 2:n_trt)
      mu_gamma0[j] <- mu_gamma0_raw[j];
}

model{
   real var_mu_gamma1;
   vector[T] theta[N];                     // ability : theta[N,T] - array of vectors
   real m_mu_alpha;
   real m_mu_kappa;
   real var_mu_alpha;
   real var_mu_kappa;


// hyperpriors for item parameters
   m_mu_alpha <- 1;
   m_mu_kappa <- 0;
   var_mu_alpha <- 10000;
   var_mu_kappa <- 10000;

// hyperparameters for item parameters;
   mu_alpha ~ normal(m_mu_alpha, sqrt(var_mu_alpha));
   mu_kappa ~ normal(m_mu_kappa, sqrt(var_mu_kappa));
   var_alpha ~ inv_gamma(1, 1);    // inv gamma hyperprior
   var_kappa ~ inv_gamma(1, 1);

// prior for discrimination, truncated below 0
   alpha ~ normal(mu_alpha, sqrt(var_alpha));

// prior for item difficulty
   kappa ~ normal(mu_kappa, sqrt(var_kappa));


// hyperpriors for random intercept and slope for the ability parameters
    //  var_mu_gamma0 ~ ;
    var_mu_gamma1 <- 10000.0;          
    m_mu_gamma0 ~ normal(0, 100.0);  
    m_mu_gamma1 ~ normal(0, 100.0);

 

// priors for random intercept and slope for ability parameters by group
   for (j in 1:n_trt) {
       mu_gamma0_raw[j] ~ normal(m_mu_gamma0, 1.0);
       mu_gamma1[j] ~ normal(m_mu_gamma1, sqrt(var_mu_gamma1));
   }

   // var_gamma1 ~ inv_gamma(0.01, 0.01);    // inv gamma hyperprior
   var_gamma1 ~ gamma(1, 1);    // gamma hyperprior

   for (i in 1:N){
      gamma0[i] ~ normal(mu_gamma0[trt[i]], 1.0);
      gamma1[i] ~ normal(mu_gamma1[trt[i]], sqrt(var_gamma1));
   }

   for (i in 1:N){
      for (t in 1:T){
         theta[i, t] <- gamma0[i] + gamma1[i]*(t - 1);
      }
   }

   
  for (t in 1:T){
    for (i in 1:N) {
      for (j in 1:J) {
          Y[i, j, t] ~  bernoulli_logit(alpha[j]*(theta[i, t] - kappa[j]));

      }
    }
  }
}

generated quantities {
    vector[n_trt - 1] deltaout;             // Treatment effect
// Code the treatment effect - difference in differences
   for (l in 1:(n_trt - 1))
      deltaout[l] <-  mu_gamma1[l+1] - mu_gamma1[1];              
}
'

#compile the model and sample in two steps
c_mlirt_rc_bin_delta1 <- stan_model(model_code = 'mlirt_rc_bin_delta1')



##########################################################################################
###                    Simulating power, type I error etc.
##########################################################################################

sim.stan.fun <- function(n.s = n.subj
                       , amat = irt.pars1$a
                       , bmat = irt.pars1$d
                       , cmat = irt.pars1$c
                       , mtheta = mu.vec
                       , m.true = mu.true             # true effect on total score scale
                       , vcov = V.theta
                       , trend = "linear"
                       , type = "3pl"
                       , stancmod1 = c_irt_pa_totscr_delta1   # compiled code for the PA model
                       , stancmod2 = c_irt_rc_totscr_delta1   # compiled code for RC model
                       , stancmod3 = c_mlirt_pa_bin_delta1    # compiled code for MLIRT PA model
                       , stancmod4 = c_mlirt_rc_bin_delta1    # compiled code for MLIRT RC model
                       , niter = 100                # number of MCMC samples
                       , rngseed = 12345
                         ) {
    sim.data.theta1 <- lirtgen(n.s = n.subj, amat = irt.pars1$a, bmat = irt.pars1$d, cmat = irt.pars1$c, mtheta = mu.vec, vcov = V.theta, trend = "linear", type = "3pl")
    sim.data1 <- sim.data.theta1[[1]]
    n.subj <- 1000
    n.time <- 5

    sim.data1$totscrsum <- apply(as.matrix(sim.data1), 1, sum)
    sim.data1$totscr <- stdize(sim.data1$totscrsum)
    sim.data1$usubjid <- rep(1:n.subj, n.time)
    sim.data1$trt <- rep(rep(1:2, each=floor(n.subj/2)), n.time)
    sim.data1$time <- rep(1:n.time, each = n.subj)
    sim.data1 <- sim.data1[order(sim.data1$usubjid, sim.data1$time), ]

    ######################################################################
                                       #
                                        #                       TOTAL SCORES
                                        #
    ######################################################################

                                        # random coef model
                                        # Prepare data
    sim.data1.melt <- melt(sim.data1, id.vars = c("usubjid", "trt", "time"), measure.var = "totscr")

    Y2 <- sim.data1.melt$value
    Subject <- sim.data1.melt$usubjid
    N_Subject <- length(unique(Subject))
    N <- dim(sim.data1.melt)[1]
    time <- sim.data1.melt$time - 1
    trt <- sim.data1.melt$trt
    n.trt <- length(levels(as.factor(trt)))    
    trt_X_time <- trt*time

    stan_dat2 <- list(N = N
                     ,Y = Y2
                     ,time = time
                     ,trt = trt
                     ,trt_X_time = trt_X_time
                     ,n_trt = n.trt
                     ,Subject = Subject
                     ,N_Subject = N_Subject
                      )

    mod.rc <- sampling(stancmod2, data = stan_dat2, iter = niter, chains = 1, seed = rngseed)
    mod.rc.res <- extract(mod.rc, inc_warmup = FALSE)
    p.deltaout.rc <- mod.rc.res$deltaout
    mod.rc.ci <- coda::HPDinterval(coda::as.mcmc(as.vector(p.deltaout.rc)))     

    
                                        # Categorical time model, AR(1)
                                        #Prepare data
    sim.data1.melt <- melt(sim.data1, id.vars = c("usubjid", "trt", "time"), measure.var = "totscr")
    sim.data1.long <- dcast(sim.data1.melt, formula = usubjid + trt  ~ time + variable)

    Y <- sim.data1.long[, -(1:2)]   # response data
    trt <- sim.data1.long[, 2]
    T <- dim(Y)[2]       # time
    N <- dim(Y)[1]       # sample size
    n.trt <- length(levels(as.factor(trt)))   # number of arms
                                        # data for stan()
    stan_dat1 <- list(T = T
                     ,N = N
                     ,Y = Y
                     ,trt = trt
                     ,n_trt = n.trt)

    mod.pa <- sampling(stancmod1, data = stan_dat1, iter = niter, chains = 1, seed = rngseed)
    mod.pa.res <- extract(mod.pa, inc_warmup = FALSE)
    p.deltaout.pa <- mod.pa.res$deltaout
    mod.pa.ci <- coda::HPDinterval(coda::as.mcmc(as.vector(p.deltaout.pa)))

    ######################################################################
                                        #
                                        #                      MLIRT MODELS
                                        #
    ######################################################################
    sim.data1 <- sim.data.theta1[[1]]
    
    ######################
                                        #    RC MLIRT model
    ######################
                                        # prepare data
    sim.data1$usubjid <- rep(1:n.subj, n.time)
    sim.data1$group <- rep(rep(1:2, each=floor(n.subj/2)), n.time)
    sim.data1$time <- rep(1:n.time, each = n.subj)

    sim.data1.melt <- melt(sim.data1, id.vars = c("usubjid", "group", "time"))
    sim.data1.long <- dcast(sim.data1.melt, formula = usubjid + group + time ~ variable)

    sim.data1.long <- sim.data1.long[order(sim.data1.long$time, sim.data1.long$group,  sim.data1.long$usubjid), ]
    Y <- abind(split(sim.data1.long[-(1:3)], sim.data1.long$time), along = 3)
    T <- dim(Y)[3]    # time
    J <- dim(Y)[2]    # n.variables 
    N <- dim(Y)[1]    # sample size
    trt <- rep(1:2, each = n.subj/2)
    n.trt <- nlevels(as.factor(trt))    # number of arms

                                        # data for stan()
    mlirt_dat2 <- list(T = T
                      ,N = N
                      ,J = J 
                      ,Y = Y
                      ,trt = trt
                      ,n_trt = n.trt)

    initf <- function() {
        list(gamma0 = rnorm(N), gamma1 = rnorm(N), alpha = runif(J), kappa = rnorm(J),
             mu_gamma0_raw= as.array(c(0.0, 0.0)), mu_gamma1= .2, var_gamma1 = .2, m_mu_gamma0=.2, m_mu_gamma1=.2,
             mu_alpha = 1, mu_kappa = 1, var_alpha = 1, var_kappa = 1)
    }
    
    mod.mlirt.rc <- sampling(stancmod4, data = mlirt_dat2, iter = niter, chains = 1) #, init = initf)
    mlirt.rc.res <- extract(mod.mlirt.rc, inc_warmup = FALSE)
    p.mlirt.deltaout.rc <- mlirt.rc.res$deltaout
    mod.mlirt.rc.ci <- coda::HPDinterval(coda::as.mcmc(as.vector(p.mlirt.deltaout.rc)))

    ######################
                                        #    PA MLIRT model
    ######################


                                        # priors
    m.mu.theta <- 0.0
    var.mu.theta <- 1.0

                                        # data for stan()
    mlirt_dat1 <- list(T = T
                      ,N = N
                      ,J = J 
                      ,Y = Y
                      ,m_mu_theta = m.mu.theta
                      ,var_mu_theta = var.mu.theta
                      ,trt = trt
                      ,n_trt = n.trt)

    mod.mlirt.pa <- sampling(stancmod3, data = mlirt_dat1, iter = niter, chains = 1)    
    mlirt.pa.res <- extract(mod.mlirt.pa, inc_warmup = FALSE)
    p.mlirt.deltaout.pa <- mlirt.pa.res$deltaout
    mod.mlirt.pa.ci <- coda::HPDinterval(coda::as.mcmc(as.vector(p.mlirt.deltaout.pa)))
    
    return(c(reject.totscr.rc = prod(mod.rc.ci) > 0,
             reject.totscr.pa = prod(mod.pa.ci) > 0,
             reject.mlirt.rc = prod(mod.mlirt.rc.ci) > 0,
             reject.mlirt.pa = prod(mod.mlirt.pa.ci) > 0,
             bias.totscr.rc = m.true - mean(p.deltaout.rc),
             bias.totscr.pa = m.true - mean(p.deltaout.pa),
             bias.mlirt.rc = m.true - mean(p.mlirt.deltaout.rc),
             bias.mlirt.pa = m.true - mean(p.mlirt.deltaout.pa)))
}


#### Simulate Item parameters

numberItems <- 11
numberItemLevels <- 2
n.time <- 5  # number of time points
times <- 1:n.time

a1 <- matrix(rlnorm(numberItems, meanlog=.2, sdlog=.2))
a.mat <- cbind(a1,a1,a1,a1,a1)

# for binary or ordinal data
d1 <- matrix(rnorm(numberItems*(numberItemLevels - 1)), numberItems)
# d1 <- apply(d, 1, sort, decreasing=TRUE)

d.mat <- cbind(d1, d1, d1, d1, d1)


# guessing parameter, set to 0 to get 2PL model
c.mat <- matrix(rep(0, numberItems*n.time), ncol=n.time)

irt.pars1 <- list(a = a.mat, d=d.mat, c=c.mat)

##################################################
# Generate a covariance matrix for theta over time
##################################################
# AR(1)
rho <- 0.5 
sigma <- 1 
times <- 1:n.time
H <- abs(outer(times, times, "-")) # matrix of powers
V.theta <- sigma * rho^H


## Trt effect
n.subj <- 1000
mu.vec <- c(0, 0.2)
mu.true <- mu.vec[2] - mu.vec[1]

# test
sim.res1a <- replicate(1, expr = try(sim.stan.fun(mtheta = mu.vec, niter = 1000)), simplify = FALSE)


##########################################################################################
#
###                                  Parallelize!
#
##########################################################################################


parReplicate <- function(cl, n, expr, simplify=TRUE, USE.NAMES=TRUE) {
    parSapply(cl, integer(n), function(i, ex) eval(ex, envir=.GlobalEnv),
              substitute(expr), simplify=simplify, USE.NAMES=USE.NAMES)
}

### The arguments simplify and USE.NAMES are compatible with sapply rather than replicate, but they make it a better wrapper around parSapply.

require(parallel)

n.clus <- 8
n.sim <- 8
rngseed <- 12345

cl <- makePSOCKcluster(n.clus)
clusterCall(cl, function() {library(reshape2); library(rstan); library(abind); library(mvtnorm)})
clusterExport(cl, c("parReplicate", "stdize", "sim.stan.fun", "irf", "lirtgen", "V.theta", "irt.pars1", "n.subj", "mu.vec",
                    "mu.true", "c_irt_rc_totscr_delta1", "c_irt_pa_totscr_delta1", "c_mlirt_pa_bin_delta1", "c_mlirt_rc_bin_delta1", "rngseed"))

sim.res1 <- parReplicate(cl, n.sim, expr = try(sim.stan.fun(mtheta = mu.vec, niter=1000)), simplify = FALSE)

sim.res <- apply(do.call("rbind", sim.res1), 2, function(x) mean(as.numeric(x), na.rm = TRUE))
# sim.res <- apply( sim.res1, 1, function(x) mean(as.numeric(x), na.rm = TRUE))

stopCluster(cl)

simres <- list(sim.res, sim.res1, mu.vec, mu.true)
save(simres, file = "simres_all.RData")
