# ------------------------------------------------------------------------------
# Simulation Study for 
#   "Bridging Biomarker Measurements to Identify Cost-Effective Biomarkers"
#
# Logistic Regression
#
# Date Created: March 4, 2025
#
# Author: Trevor Thomson
#         tthomson@fredhutch.org
# ------------------------------------------------------------------------------

rm(list =ls())
set.seed(5347547)

library(MASS)
library(pbapply)
library(nleqslv)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(stringr)
library(numDeriv) # for "grad"
library(cubature) # for "cuhre"
library(survey) # for "svyglm"

N1_Setting = 1
MeasurementErrorVariance_Setting = 1
Setting_theta = 2

# Specify the initial sizes of the datasets
N0 = 250 # <- We will consider increments of 50


# Specify the cohort sample size for population 1
if (N1_Setting == 1) N1 = 5000 # <- initial cohort sample size for population 1
if (N1_Setting == 2) N1 = 10000 # <- initial cohort sample size for population 1
if (N1_Setting == 3) N1 = 15000 # <- initial cohort sample size for population 1


# Specify the size of the first phase sampling under the case-cohort study design
n1 = 0.04 * N1

# Specify various sample sizes for the paired sample
n0_vec = seq(100, N0, by = 50)

# Specify beta
Beta = c(-3.6, -1)

# Specify the mean and variance of U for samples 0 and 1
mu_U1_vec = c(0,0)
sigma2_U1_vec = c(1,1)

# Specify the value of theta0 and theta1
if (Setting_theta == 1) theta = c(0, 1)
if (Setting_theta == 2) theta = c(0, 0.6)

# Based on the earlier specifications, define other parameters
mu_U2_vec = theta[1] + theta[2] * mu_U1_vec
sigma2_U2_vec = theta[2]^2 * sigma2_U1_vec
Cov_U1U2_vec = theta[2] * sigma2_U1_vec

# Specify the variance of measurement errors across cohorts
if (MeasurementErrorVariance_Setting == 1){
  sigma2_eps1_vec = sigma2_U1_vec * (1.05-1) # Ratio of sigma2_X / sigma2_U = 1.05
  sigma2_eps2_vec = sigma2_U1_vec * (1.05-1) # Ratio of sigma2_X / sigma2_U = 1.05
}
if (MeasurementErrorVariance_Setting == 2){
  sigma2_eps1_vec = sigma2_U1_vec * (1.25-1) # Ratio of sigma2_X / sigma2_U = 1.25
  sigma2_eps2_vec = sigma2_U1_vec * (1.25-1) # Ratio of sigma2_X / sigma2_U = 1.25
}
if (MeasurementErrorVariance_Setting == 3){
  sigma2_eps1_vec = sigma2_U1_vec * (1.50-1) # Ratio of sigma2_X / sigma2_U = 1.50
  sigma2_eps2_vec = sigma2_U1_vec * (1.50-1) # Ratio of sigma2_X / sigma2_U = 1.50
}

# Define the variance of X1 and X2
Var_X1_vec = sigma2_U1_vec + sigma2_eps1_vec
Var_X2_vec = sigma2_U2_vec + sigma2_eps2_vec
Cov_X1X2_vec = Cov_U1U2_vec


# Number of simulations
M = 1000

# ---------------------------------------------------------------------------
#           TABLE OF CONTENTS
#
# (0) Define Functions
#
#   (0.1) P(Y = 1 | X2) Function
#
#   (0.2) Log Likelihood Function Paired Sample Data
#
#   (0.3) Estimate Mean and Variance from Normal Distribution: 
#         Log-Likelihood Function
#
#   (0.4) Estimate Mean and Variance from Normal Distribution: 
#         Score Function
#
#   (0.5) P(Y = 1 | X2) Function
#
# (1) Initialize vectors and matrices
#
# ------------------------------------------------------------------------
#                     RUN SIMULATION
# ------------------------------------------------------------------------
#
# (2) Generate Covariates and the Response for each Cohort
#
# (3) Estimate Parameters from Paired Sample
#
# (4) Estimate Parameters with Population 1
#
# (5) Estimate the Variance of Thetahat
#
#
# ------------------------------------------------------------------------
#                     END OF SIMULATION
# ------------------------------------------------------------------------
#
# ---------------------------------------------------------------------------

# ... Source C++ functions
sourceCpp('Cpp_Functions_Logistic_Github.cpp')
f = function(x){
  val = 1 + exp(-x)
  val2 = 1 / val
  return(val2)
}


# ---------------------------------------------------------------------------
# (0) Define Functions
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# (0.1) P(Y = 1 | X2) Function
# ---------------------------------------------------------------------------

# Assuming that the response generating model is a logistic model,
# We have P(Y = 1 | X2) = \int 1/(exp(-beta0 - beta1 u_1)) dF(u_1|X2),
#   where f(u_1|X2) is a normal pdf
# We numerically approximate this integral to obtain a value for P(Y = 1 | X)


# --- Input:
# Beta: Vector of regression parameters (beta0, beta1)
# X2: Vector for covariate
# mu_U1: Scalar for mean of U_1
# mu_U2: Scalar for mean of U_2
# sigma2_U1: Scalar for Var(U_1)
# Cov_U1U2: Scalar for Cov(U_1, U_2)
# Var_X2: Scalar for Var(X_2)
#
# --- Output:
# A vector containing P(Y = 1 | X)

Prob_Y1_GivenX2_Logistic_func = function(Beta,
                                         X2,
                                         mu_U1,
                                         mu_U2,
                                         sigma2_U1,
                                         Cov_U1U2,
                                         Var_X2){
  
  beta0 = Beta[1]
  beta1 = Beta[2]
  
  # compute \int \pi(u_1; \beta) dF(u_1|X2) 
  #   pi(u; \beta) = 1/[1+exp(-beta0 - beta1*u)]
  #   U_1|X2 ~ N(mu_U1 + Cov_U1U2 / Var_X2 * (X_j - mu_Uj),  sigma2_U1 - Cov_U1Uj^2 / Var_Xj) 
  Integral_val = sapply(1:length(X2),
                        function(i){
                          x = X2[i]
                          # Define the integrand $\pi(u; \beta) f(u|x)$, 
                          #    where f(u|x) is normal
                          Integrand = function(u){
                            lin.pred = beta0 + beta1*u
                            prob_model = 1/(1+exp(-lin.pred))
                            mu = mu_U1 + Cov_U1U2/Var_X2 * (x - mu_U2)
                            sigma2 = sigma2_U1 - Cov_U1U2^2/Var_X2
                            sigma = ifelse(sigma2 >= 0,
                                           sqrt(sigma2),
                                           0)
                            dens.val = dnorm(u, mean = mu, sd = sigma)
                            val = prob_model * dens.val
                            return(val)
                          }
                          # Approximate the integral "Integrand" with "cuhre"
                          #    where we integral is wrt to "u"
                          val = cuhre(f = Integrand, 
                                      nComp = 1,
                                      lowerLimit = -Inf, 
                                      upperLimit = Inf)$integral
                          
                          # return "val"
                          return(val)
                        })
  
  if (sum(Integral_val) == 0) Integral_val = rep(NA, length(X2))
  
  # Return "Integral_val"
  return(Integral_val)
}


# ---------------------------------------------------------------------------
#   (0.2) Log Likelihood Function Paired Sample Data
# ---------------------------------------------------------------------------

# Specify the log-likelihood function for paired sample with lower limits of detection

# --- Input
#
# mu_1: Scalar for the mean of X_1
# sigma2_1: Scalar for the variance of X_1
# theta: Vector of length 2 for parameters theta0 and theta1
# sigma2_eps1: Scalar for the measurement error variance of X_1
# sigma2_eps2: Scalar for the measurement error variance of X_2
# X1_vec: Vector of X_1 data
# X2_vec: Vector of X_2 data
#
# --- Output:
# logL: log-likelihood value

LogLikFcn_PairedSample = function(mu_1,
                                  sigma2_1,
                                  theta,
                                  sigma2_eps1 = sigma2_eps1_vec[1],
                                  sigma2_eps2 = sigma2_eps2_vec[1],
                                  X1_vec,
                                  X2_vec){
  
  # Define parameters
  theta0 = theta[1]
  theta1 = theta[2]
  mu_2 = theta0 + theta1 * mu_1
  sigma2_U1 = sigma2_1 - sigma2_eps1
  Cov_X1X2 = theta1*sigma2_U1
  sigma2_U2 = theta1^2*sigma2_U1
  sigma2_2 = sigma2_U2 + sigma2_eps2
  
  # Define mu
  mu = c(mu_1, mu_2)
  
  # Define Sigma
  Sigma = matrix(c(sigma2_1, Cov_X1X2, Cov_X1X2, sigma2_2),
                 nrow = 2, ncol = 2)
  
  val1 = mvtnorm::dmvnorm(cbind(X1_vec,
                                X2_vec), 
                          mean = mu, 
                          sigma = Sigma,
                          log = T)
  
  val = sum(val1)
  return(val)
}


# ---------------------------------------------------------------------------
#   (0.3) Estimate Mean and Variance from Normal Distribution: 
#         Log-Likelihood Function
# ---------------------------------------------------------------------------

# Let X ~ N(mu, sigma^2). Suppose we observe W = max(X, a), where a is known.
#   Given a random sample of {W_1, ..., W_n}, the ith unit's contribution to the
#   likelihood function is 
#   phi(W_i; mu, sigma^2)
#   where phi(x;mu,sigma^2) is the PDF evaluated at x
#   of the normal distribution with mean mu and variance sigma^2.

# --- INPUT:
# W: Data of left-censored observations
# mu: mean parameter
# sigma2: variance parameter
# IPW: Inverse probability weights (Vector) (Default is 1)

# --- OUTPUT:
# U(delta): value of estimating function evaluated at delta 
#       (Scalar if "Sum == T"; Vector of length n otherwise)

LogLik_LeftCensored_Normal = function(W,
                                      mu, 
                                      IPW = 1,
                                      sigma2,
                                      LLOD = -Inf){
  
  # Compute log(dnorm(W[i], mu, sigma2)) for each i; we will update this value for the units with a value of the LLOD
  valA = dnorm(W, 
               mean = mu, 
               sd = sqrt(sigma2), 
               log = T) 
  
  # Multiply each observation by its weight
  val_final = valA * IPW
  
  # return "sum(val_final)"
  val = sum(val_final)
  return(val)
}


# ---------------------------------------------------------------------------
#   (0.4) Estimate Mean and Variance from Normal Distribution: 
#         Score Function
# ---------------------------------------------------------------------------

# Let X ~ N(mu, sigma^2). Suppose we observe W = max(X, a), where a is known.
#   The function below corresponds to the score function

# --- INPUT:
# W: Data of left-censored observations
# mu: mean parameter
# sigma2: variance parameter
# IPW: Inverse probability weights (Vector) (Default is 1)
# LLOD: Lower limit of detection (Default is "-Inf")
# Sum: Do we want to take the sum of each component of the estimating function?
#       (Logical; default is "T")

# --- OUTPUT:
# U(delta): value of estimating function evaluated at delta 
#       (Scalar if "Sum == T"; Vector of length n otherwise)


ScoreFcn_LeftCensored_Normal = function(W,
                                        mu, 
                                        sigma2,
                                        IPW = 1,
                                        LLOD = -Inf,
                                        Sum = T){
  
  # For units with "W[i]  > LLOD", their contribution to the score function is...
  valA = (W - mu) / sigma2
  valB = 0.5/sigma2 * ( (W-mu)^2 / sigma2 - 1)
  
  # Weight each contribution to the score function
  valA = IPW * valA
  valB = IPW * valB
  
  # Combine "valA" and "valB" together
  val = cbind(valA, valB)
  if (Sum == T) val = as.numeric(colSums(val))
  return(val)
}



# ---------------------------------------------------------------------------
# (1) Initialize vectors and matrices
# ---------------------------------------------------------------------------

# Note: beta is true regression coefficient
#       gamma is coefficient using noisy data

# Mean of Y
Mean_Y_D1_vec = rep(NA, M)

# Sample Size for population 1
SampleSize_D1_vec = rep(NA, M)

# Estimated effects for beta = (beta0, beta1)'
Betahat_Logistic_True_D1_mat = matrix(NA, nrow = M, ncol = 2)
Betahat_Logistic_MLE_D1_mat = matrix(NA, nrow = M, ncol = 2)
Var_Betahat_Logistic_True_D1_mat = matrix(NA, nrow = M, ncol = 2)
Var_Betahat_Logistic_MLE_D1_mat = matrix(NA, nrow = M, ncol = 2)

# Estimated effects for theta = (theta0, theta1)' and variance estimator
MeanVarX1Paired_MultiSampleSize_list = lapply(1:length(n0_vec),
                                              function(i) matrix(NA, nrow = M, ncol = 2))
MeanVarX2Paired_MultiSampleSize_list = lapply(1:length(n0_vec),
                                              function(i) matrix(NA, nrow = M, ncol = 2))
Thetahat_MultiSampleSize_list = lapply(1:length(n0_vec),
                                       function(i) matrix(NA, nrow = M, ncol = 2))
Var_Thetahat_MultiSampleSize_list = Thetahat_MultiSampleSize_list

# Mean and variance parameters for each sample
Estimated_Mean_D1_vec = rep(NA, M)
Estimated_Variance_D1_vec = rep(NA, M)

# Estimated probabilities P(Y=1|X) = \int pi(u) dF(u|X) over a grid of X values
x_grid = qnorm(c(0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.5, 0.8, 0.9, 0.925, 0.95, 0.975, 0.99),
               mean = 0,
               sd = sqrt(mu_U1_vec[2] + sigma2_eps1_vec[2]))
Prob_Y1_GivenX2_Estimate_list = lapply(1:length(n0_vec), function(i) matrix(NA, nrow = M, ncol = length(x_grid)))
Prob_Y1_GivenX2_Variance_list = lapply(1:length(n0_vec), function(i) matrix(NA, nrow = M, ncol = length(x_grid)))

# Get the probability based on the logistic regression with effect Gamma2
Prob_Gammahat_Logistic_X2Direct_D1_mat = matrix(NA, nrow = M, ncol = length(x_grid))
Var_Prob_Gammahat_Logistic_X2Direct_D1_mat = matrix(NA, nrow = M, ncol = length(x_grid))

# ---------------------------------------------------------------------------
#                     RUN SIMULATION
# ---------------------------------------------------------------------------


for (m in 1:M){
  
  seed_initial = 5347547 + m
  set.seed(seed_initial)
  
  print(paste("Working on Simulation", m, "out of", M))
  
  # ---------------------------------------------------------------------------
  # (2) Generate Covariates and the Response for each Cohort
  # ---------------------------------------------------------------------------
  
  # -----------------------------------------
  #       Data "D0"
  # -----------------------------------------
  
  # Generate true biomarker values
  U1_D0 = rnorm(N0, 
                mean = mu_U1_vec[1],
                sd = sqrt(sigma2_U1_vec[1]))
  
  U2_D0 = theta[1] + theta[2] * U1_D0
  
  # Generate measurement error
  MeasurementError_D0 = MASS::mvrnorm(n = N0, 
                                      mu = c(0, 0),
                                      Sigma = matrix(c(sigma2_eps1_vec[1], 0,
                                                       0, sigma2_eps2_vec[1]),
                                                     nrow = 2, ncol = 2) )
  
  # Obtain observed biomarker values
  X1_D0 = U1_D0 + MeasurementError_D0[,1]
  X2_D0 = U2_D0 + MeasurementError_D0[,2]
  
  
  # Have the data in lists
  U1_D0_list = lapply(1:length(n0_vec), 
                      function(i) U1_D0[1:n0_vec[i]])
  U2_D0_list = lapply(1:length(n0_vec), 
                      function(i) U2_D0[1:n0_vec[i]])
  X1_D0_list = lapply(1:length(n0_vec), 
                      function(i) X1_D0[1:n0_vec[i]])
  X2_D0_list = lapply(1:length(n0_vec), 
                      function(i) X2_D0[1:n0_vec[i]])
  
  
  # Generate "D0_list"
  D0_list = lapply(1:length(n0_vec),
                   function(i){
                     data.frame(X1 = X1_D0_list[[i]],
                                X2 = X2_D0_list[[i]],
                                Prob = 1)
                   })
  
  
  # -----------------------------------------
  #       Data "D1"
  # -----------------------------------------
  
  # Generate true biomarker values
  U1_D1 = rnorm(N1, 
                mean = mu_U1_vec[2],
                sd = sqrt(sigma2_U1_vec[2]))
  U2_D1 = theta[1] + theta[2] * U1_D1
  
  # Generate measurement error
  MeasurementError_D1 = MASS::mvrnorm(n = N1, 
                                      mu = c(0, 0),
                                      Sigma = matrix(c(sigma2_eps1_vec[2], 0,
                                                       0, sigma2_eps2_vec[2]),
                                                     nrow = 2, ncol = 2) )
  
  # Obtain observed biomarker values
  X1_D1 = U1_D1 + MeasurementError_D1[,1]
  X2_D1 = U2_D1 + MeasurementError_D1[,2]
  
  # Generate Y
  Y_D1 = sapply(1:N1,
                function(i) rbinom(1, size = 1, prob = f(Beta[1] + Beta[2] * U1_D1[i])))
  
  # Proceed to sample units into the subcohort
  index_SubCohort_D1 = sample(1:N1,
                              size = n1)
  
  # Include all cases as well
  index_SubCohort_D1 = c(index_SubCohort_D1,
                         which(Y_D1 == 1))
  index_SubCohort_D1 = unique(index_SubCohort_D1)
  
  # Generate "D1"
  D1 = data.frame(Y = Y_D1,
                  X1 = X1_D1,
                  X2 = X2_D1,
                  U1 = U1_D1,
                  Strata = 1,
                  Prob = n1 / N1)
  D1 = D1[index_SubCohort_D1,]
  
  # Update the probability if Y == 1
  D1[which(D1$Y == 1),]$Prob = 1
  D1$Weight = 1/D1$Prob
  
  
  Mean_Y_D1_vec[m] = mean(D1$Y)
  SampleSize_D1_vec[m] = dim(D1)[1]
  
  
  # ---------------------------------------------------------------------------
  # (3) Estimate Parameters from Paired Sample
  # ---------------------------------------------------------------------------
  
  A_theta_list = list() # Initialize
  
  for (i in 1:length(n0_vec)){
    
    # Extract the paired sample data
    D0_data = D0_list[[i]]
    
    # Define the negative log-likelihood function for paired sample
    OF_PairedSample = function(par){
      -LogLikFcn_PairedSample(mu_1 = par[1],
                              sigma2_1 = par[2],
                              theta = par[3:4],
                              X1_vec = D0_list[[i]]$X1,
                              X2_vec = D0_list[[i]]$X2)   
    }
    
    # Set the starting point
    start = c(mu_U1_vec[1],
              sigma2_U1_vec[1] + sigma2_eps1_vec[1],
              theta)
    
    # Minimize the objective function
    MLEs_Paired = suppressWarnings(nlm(f = OF_PairedSample,
                                       p = start,
                                       hessian = T))
    
    # Store the hessian
    A_theta_list[[i]] = MLEs_Paired$hessian
    
    # Obtain the gradient
    ScoreFunction_Paired_use = function(par,
                                        Sum = F){
      val_mat = sapply(1:length(D0_list[[i]]$X1),
                       function(j){
                         f = function(par){
                           -LogLikFcn_PairedSample(mu_1 = par[1],
                                                   sigma2_1 = par[2],
                                                   theta = par[3:4],
                                                   X1_vec = D0_list[[i]]$X1[j],
                                                   X2_vec = D0_list[[i]]$X2[j])
                         }
                         val = grad(f, par)
                         return(val)
                       })
      val_mat = t(val_mat)
      val = val_mat # Initalize
      if (Sum == T) val = apply(val, 2, sum)
      return(val)
    }
    
    # Obtain an estimate of the variance
    A_alpha_D0 = MLEs_Paired$hessian
    B_alpha_D0 = crossprod(ScoreFunction_Paired_use(MLEs_Paired$estimate))
    Var_alpha_D0 = MASS::ginv(A_alpha_D0) %*% B_alpha_D0 %*% t(MASS::ginv(A_alpha_D0))
    Var_Thetahat_MultiSampleSize_list[[i]][m,] = diag(Var_alpha_D0)[3:4]
    
    # Store the estimate for mu, sigma2, theta0 and theta1
    MeanVarX1Paired_MultiSampleSize_list[[i]][m,] = MLEs_Paired$estimate[1:2]
    Thetahat_MultiSampleSize_list[[i]][m,] = MLEs_Paired$estimate[3:4]
    
    # Estimate the variance of X2 with paired sample
    EstFcn_alpha_D0 = function(par){
      ScoreFcn_LeftCensored_Normal(W = D0_list[[i]]$X2,
                                   mu = par[1],
                                   sigma2 = par[2],
                                   IPW = 1/D0_list[[i]]$Prob) 
    }
    
    # Set the starting value
    start = c(mu_U2_vec[1], Var_X2_vec[1])
    
    # Solve the estimating equation
    mod_alpha_D0 = nleqslv(x = start,
                           fn = EstFcn_alpha_D0,
                           jacobian = T,
                           control = list(trace = 0,
                                          maxit = 500))
    
    MeanVarX2Paired_MultiSampleSize_list[[i]][m,] = mod_alpha_D0$x
  }
  
  
  # ---------------------------------------------------------------------------
  # (4) Estimate Parameters with Population 1
  # ---------------------------------------------------------------------------
  
  
  # -----------------------------------------
  # Obtain Direct Estimates
  # -----------------------------------------
  
  svyDesign = svydesign(id = ~1,
                        weights = ~Weight,
                        data = D1)
  
  # Estimate the coefficients with true covariate U1, and store results
  mod_U1_logit_D1 = suppressWarnings( svyglm(Y ~ U1,
                                             data = D1,
                                             design = svyDesign,
                                             family=binomial(link='logit')) )
  Betahat_Logistic_True_D1_mat[m,] = as.numeric(mod_U1_logit_D1$coefficients)
  Var_Betahat_Logistic_True_D1_mat[m,] = as.numeric(diag(summary(mod_U1_logit_D1)$cov.unscaled))
  
  
  # Pretend we have X2 in our dataset, proceed to estimate the regression parameters
  mod_X2_logit_D1 = suppressWarnings( svyglm(Y ~ X2,
                                             data = D1,
                                             design = svyDesign,
                                             family=binomial(link='logit')) )
  
  # Estimate P(Y = 1 | x) here, obtain Monte Carlo variance
  BB = 1000
  set.seed(seed_initial)
  lin.pred_X2 = mod_X2_logit_D1$coefficients[1] + mod_X2_logit_D1$coefficients[2] * x_grid
  Prob_Gammahat_Logistic_X2Direct_D1_mat[m,] = 1/(1+exp(-lin.pred_X2))
  
  Gammahat_Logistic_B = MASS::mvrnorm(BB, 
                                      mu = mod_X2_logit_D1$coefficients, 
                                      Sigma = summary(mod_X2_logit_D1)$cov.unscaled)
  
  Predicted_Probs_Gammahat_B = sapply(1:BB,
                                      function(b){
                                        gamma.b = Gammahat_Logistic_B[b,]
                                        val = 1/(1+exp(-gamma.b[1] - gamma.b[2] * x_grid))
                                        return(val)
                                      })
  Predicted_Probs_Gammahat_B = t(Predicted_Probs_Gammahat_B)
  Var_Prob_Gammahat_Logistic_X2Direct_D1_mat[m,] = diag(var(Predicted_Probs_Gammahat_B))
  
  # -----------------------------
  # Estimate \hat{\alpha}_1
  # -----------------------------
  
  # Estimate mean and variance for X1
  EstFcn_alpha_D1 = function(par){
    ScoreFcn_LeftCensored_Normal(W = D1$X1,
                                 mu = par[1],
                                 sigma2 = par[2],
                                 IPW = 1/D1$Prob) 
  }
  
  # Test the estimating function
  start = c(mu_U1_vec[2], Var_X2_vec[2])
  EstFcn_alpha_D1(start)
  
  # Solve the estimating equation
  mod_alpha_D1 = nleqslv(x = start,
                         fn = EstFcn_alpha_D1,
                         jacobian = T,
                         control = list(trace = 0,
                                        maxit = 500))
  
  mod_alpha_D1$x; c(mean(D1$X1), var(D1$X1)); start
  
  # Store the results
  Estimated_Mean_D1_vec[m] = mod_alpha_D1$x[1]
  Estimated_Variance_D1_vec[m] = mod_alpha_D1$x[2]
  
  # -----------------------------------------
  # Estimate \hat{\beta}
  # -----------------------------------------
  
  
  # Specify a function for the estimating function (Logistic Regression)
  EstFcn_beta_Logistic_D1 = function(par){
    set.seed(seed_initial)
    val_EF = ScoreFunction_cpp(beta = par,
                               Y_vec = D1$Y,
                               Xstar_vec = D1$X1,
                               Omega_vec = 1 / D1$Prob,
                               LLOD = -Inf,
                               mu_U = mod_alpha_D1$x[1],
                               mu_X = mod_alpha_D1$x[1],
                               sigma2_U = mod_alpha_D1$x[2] - sigma2_eps1_vec[2],
                               sigma2_X = mod_alpha_D1$x[2],
                               B = 500)
    val = apply(val_EF, 2, sum)
    return(val)
  }
  
  # Test the estimating function
  EstFcn_beta_Logistic_D1(Beta)
  
  # Solve the estimating equation
  mod_beta_Logistic_D1 = nleqslv(x = Beta,
                                 fn = EstFcn_beta_Logistic_D1,
                                 jacobian = T,
                                 control = list(trace = 0,
                                                maxit = 500,
                                                ftol = 1e-8))
  # mod_beta_Logistic_D1$x; as.numeric(mod_U1_logit_D1$coefficients); Beta
  
  # Store the results
  Betahat_Logistic_MLE_D1_mat[m,] = mod_beta_Logistic_D1$x
  
  
  # ---------------------------------------------------------------------------
  # (6) Estimate the Variance of Thetahat
  # ---------------------------------------------------------------------------
  
  
  
  # Run a for loop for each paired data sample size
  for (i in 1:length(n0_vec)){
    
    A_logistic_func = function(par){
      set.seed(seed_initial)
      alpha1 = par[1:2]
      alpha2 = par[3:4]
      beta = par[5:6]
      omega = par[7]
      
      # Compute the estimating equation for the nuisance parameter
      val1 = ScoreFcn_LeftCensored_Normal(W = D1$X1,
                                          mu = alpha1[1],
                                          sigma2 = alpha1[2],
                                          IPW = 1/D1$Prob,
                                          Sum = T)
      val2 = ScoreFcn_LeftCensored_Normal(W = D0_list[[i]]$X2,
                                          mu = alpha2[1],
                                          sigma2 = alpha2[2],
                                          IPW = 1,
                                          Sum = T)
      
      # Compute the estimating equation for beta
      val3 = ScoreFunction_cpp(beta = beta,
                               Y_vec = D1$Y,
                               Xstar_vec = D1$X1,
                               Omega_vec = 1 / D1$Prob,
                               LLOD = -Inf,
                               mu_U = alpha1[1],
                               mu_X = alpha1[1],
                               sigma2_U = alpha1[2] - sigma2_eps1_vec[2],
                               sigma2_X = alpha1[2],
                               B = 500)
      val3_sum = apply(val3, 2, sum)
      
      val4 = rep(omega - dim(D1)[1] / N1, dim(D1)[1])
      val4_sum = sum(val4)
      
      val = c(val1,val2,val3_sum, val4_sum)
      return(val)
    }
    
    A_2_logit = jacobian(A_logistic_func, 
                         x = c(mod_alpha_D1$x,
                               MeanVarX2Paired_MultiSampleSize_list[[i]][m,],
                               mod_beta_Logistic_D1$x,
                               dim(D1)[1]/N1))
    
    # Obtain the score function for theta
    ScoreFunction_Paired_use = function(param,
                                        Sum = F){
      val_mat = sapply(1:length(D0_list[[i]]$X1),
                       function(j){
                         f = function(param){
                           -LogLikFcn_PairedSample(mu_1 = param[1],
                                                   sigma2_1 = param[2],
                                                   theta = param[3:4],
                                                   X1_vec = D0_list[[i]]$X1[j],
                                                   X2_vec = D0_list[[i]]$X2[j])
                         }
                         val = grad(f, param)
                         return(val)
                       })
      val_mat = t(val_mat)
      return(val_mat)
    }
    
    # Create a block diagonal matrix 
    A_logit_mat = Matrix::bdiag(A_theta_list[[i]], A_2_logit)
    A_logit_mat = as.matrix(A_logit_mat)
    
    
    B_logistic_func = function(par){
      set.seed(seed_initial)
      alpha0 = par[1:2]
      theta = par[3:4]
      alpha1 = par[5:6]
      alpha2 = par[7:8]
      beta = par[9:10]
      omega = par[11]
      
      # Compute the estimating equation for beta
      val1 = ScoreFunction_Paired_use(param = c(alpha0,
                                                theta))
      
      # Compute the estimating equation for the nuisance parameter
      val2 = ScoreFcn_LeftCensored_Normal(W = D1$X1,
                                          mu = alpha1[1],
                                          sigma2 = alpha1[2],
                                          IPW = 1/D1$Prob,
                                          Sum = F)
      val3 = ScoreFcn_LeftCensored_Normal(W = D0_list[[i]]$X2,
                                          mu = alpha2[1],
                                          sigma2 = alpha2[2],
                                          IPW = 1,
                                          Sum = F)
      
      # Compute the estimating equation for beta
      val4 = ScoreFunction_cpp(beta = beta,
                               Y_vec = D1$Y,
                               Xstar_vec = D1$X1,
                               Omega_vec = 1 / D1$Prob,
                               LLOD = -Inf,
                               mu_U = alpha1[1],
                               mu_X = alpha1[1],
                               sigma2_U = alpha1[2] - sigma2_eps1_vec[2],
                               sigma2_X = alpha1[2],
                               B = 500)
      
      # Initialize a matrix to compile results together
      m0 = n0_vec[i]
      m1 = dim(D1)[1]
      mat = matrix(0, 
                   nrow = m0 + m1,
                   ncol = length(par))
      mat[1:m0,1:4] = val1
      mat[(m0+1):(m0+m1),5:6] = val2
      mat[1:m0,7:8] = val3
      mat[(m0+1):(m0+m1),9:10] = val4
      mat[,11] = 0 # for omega
      return(mat)
    }
    
    B_logit_mat = crossprod(B_logistic_func(c(MeanVarX1Paired_MultiSampleSize_list[[i]][m,],
                                              Thetahat_MultiSampleSize_list[[i]][m,],
                                              mod_alpha_D1$x,
                                              MeanVarX2Paired_MultiSampleSize_list[[i]][m,],
                                              mod_beta_Logistic_D1$x,
                                              dim(D1)[1]/N1)))
    Var_Theta_logit_mat = Var_beta_logit_D1 = MASS::ginv(A_logit_mat) %*% B_logit_mat %*% t(MASS::ginv(A_logit_mat))
    
    Var_Theta_logit_vec = diag(Var_Theta_logit_mat)
    
    # Store the estimated variance of beta
    Var_Betahat_Logistic_MLE_D1_mat[m,] = Var_Theta_logit_vec[9:10]
    
    # Run the Monte Carlo algorithm
    # Generate 1000 realizations from the sampling distribution of Thetahat
    BB = 1500
    set.seed(seed_initial)
    Thetahat_Logistic = c(MeanVarX1Paired_MultiSampleSize_list[[i]][m,],
                          Thetahat_MultiSampleSize_list[[i]][m,],
                          mod_alpha_D1$x,
                          MeanVarX2Paired_MultiSampleSize_list[[i]][m,],
                          mod_beta_Logistic_D1$x,
                          dim(D1)[1]/N1)
    Thetahat_Logistic_B = MASS::mvrnorm(BB, 
                                        mu = Thetahat_Logistic, 
                                        Sigma = Var_Theta_logit_mat)
    
    # Obtain our point estimate for P(Y = 1 | X)
    Prob_Y1_GivenX2_Estimate_list[[i]][m,] = Prob_Y1_GivenX2_Logistic_func(Beta = Thetahat_Logistic[9:10],
                                                                           X2 = x_grid,
                                                                           mu_U1 = Thetahat_Logistic[5],
                                                                           mu_U2 = Thetahat_Logistic[7],
                                                                           sigma2_U1 = Thetahat_Logistic[6] - sigma2_eps1_vec[2],
                                                                           Cov_U1U2 = (Thetahat_Logistic[6] - sigma2_eps1_vec[2])*Thetahat_Logistic[4],
                                                                           Var_X2 = Thetahat_Logistic[8])
    
    # Run an sapply statement to obtain the estimate for P(Y=1|X)
    Est_Probabilities_Logistic = sapply(1:BB,
                                        function(b){
                                          Thetahat.b = Thetahat_Logistic_B[b,]
                                          thetahat.b = Thetahat.b[3:4]
                                          Beta.b = Thetahat.b[9:10]
                                          mu_1.b = Thetahat.b[5]
                                          sigma2_1.b = Thetahat.b[6]
                                          mu_2.b = Thetahat.b[7]
                                          sigma2_2.b = Thetahat.b[8]
                                          val = Prob_Y1_GivenX2_Logistic_func(Beta = Beta.b,
                                                                              X2 = x_grid,
                                                                              mu_U1 = mu_1.b,
                                                                              mu_U2 = mu_2.b,
                                                                              sigma2_U1 = sigma2_1.b - sigma2_eps1_vec[2],
                                                                              Cov_U1U2 = (sigma2_1.b - sigma2_eps1_vec[2])*thetahat.b[2],
                                                                              Var_X2 = sigma2_2.b)
                                          return(val)
                                        })
    Est_Probabilities_Logistic = t(Est_Probabilities_Logistic)
    Prob_Y1_GivenX2_Variance_list[[i]][m,] = diag(var(Est_Probabilities_Logistic, na.rm = T))
    
  }  # <- End of for loop for paired sample data size
  
  # ------------------------------------------------------------------------
  #                     END OF SIMULATION
  # ------------------------------------------------------------------------
  
}