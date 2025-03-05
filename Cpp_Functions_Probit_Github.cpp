//#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//#include <truncnorm.h>

using namespace Rcpp;



// ---------------------------------------------------------------------------------------------------------------------------------------------

// Write a function that computes pi(U; beta)

// Input: beta = (beta0, beta1): Vector of parameters
//        U_i: Covariate (Scalar)

// [[Rcpp::export]] 
double pi_U_beta_probit_cpp(NumericVector beta,
                            double U_i){
  
  // Extract beta0 and beta1 from beta
  double beta0_pi = beta(0);
  double beta1_pi = beta(1);
  // Compute the linear predictor
  double x = beta0_pi + beta1_pi * U_i;
  // return "Phi(beta0 + beta1 * U_i)"
  double val = 0.5 * (1 + std::erf(x / std::sqrt(2.0)));
  return(val);
}


// ---------------------------------------------------------------------------------------------------------------------------------------------

// Write a function that generates variables a normal distribution with mean mu and variance sigma2
// Input: Xstar: Scalar covariate
//        mu_U: Mean of U (Scalar)
//        mu_X: Mean of X (Scalar)
//        sigma2_U: Variance of U (Scalar)
//        sigma2_X: Variance of X (Scalar)
//        B: Number of iterations for the Monte Carlo approximation (Scalar; default is 10000)

// [[Rcpp::export]]
NumericVector GenerateRVs_func(double Xstar,
                               double mu_U,
                               double mu_X,
                               double sigma2_U,
                               double sigma2_X,
                               int B = 10000){
  
  NumericVector val(B); // Initialize
  
  // Generate random variables from [U|X]: normal with mean and variance defined below
  double mu_U_X = mu_U + sigma2_U / sigma2_X * (Xstar - mu_X);
  double sigma2_U_X = sigma2_U - (sigma2_U * sigma2_U / sigma2_X);
  
  // Generate B observations from [U|X]
  val = rnorm(B, mu_U_X, sqrt(sigma2_U_X));
  
  
  return(val);
}


// ---------------------------------------------------------------------------------------------------------------------------------------------

// Write a function for the score function corresponding to the observed 
//    data likelihood function

// Input: beta: Vector of parameters
//        Y_vec: Response (Vector)
//        Xstar_vec: Noisy Covariate (Vector)
//        Omega_vec: IPWs (Vector)
//        LLOD: Lower limit of detection (Scalar)
//        mu_U: Mean of U (Scalar)
//        mu_X: Mean of X (Scalar)
//        sigma2_U: Variance of U (Scalar)
//        sigma2_X: Variance of X (Scalar)
//        B: Number of iterations for the Monte Carlo approximation (Scalar; default is 10000)
//        MeasurementErrorLogi: Do we want to account for measurement error (Logical)


// [[Rcpp::export]]
NumericMatrix ScoreFunction_cpp(NumericVector beta,
                                NumericVector Y_vec,
                                NumericVector Xstar_vec,
                                NumericVector Omega_vec,
                                double LLOD,
                                double mu_U,
                                double mu_X,
                                double sigma2_U,
                                double sigma2_X,
                                int B = 10000,
                                bool MeasurementErrorLogi = true){
  
  // Get the number of units
  int n = Y_vec.size();
  
  // Run a for loop over each unit, 
  //    Generate observations from [U | X; alpha]
  //    Take a Monte Carlo expectation
  NumericMatrix val(n,2); // Initialize
  
  for (int i = 0 ; i < n ; ++i){ 
    
    // Extract elements from "Y", "Xstar", and "Omega_vec"
    double Y_i = Y_vec(i);
    double Xstar_i = Xstar_vec(i);
    double IPW_i = Omega_vec(i);
    
    // Break into cases depending if we want to account for measurement error or not
    if (MeasurementErrorLogi == true){
      
      // Proceed to generate B random variables from [U | Xstar]
      NumericVector U_i_B = GenerateRVs_func(Xstar_i,
                                             mu_U,
                                             mu_X,
                                             sigma2_U,
                                             sigma2_X,
                                             B); // generated from normal distribution [U|X]
      
      // Run a loop to compute [Y_i | U_{ib}]
      double val_Z = 0; // initialize
      
      for (int b = 0 ; b < B ; ++b){ 
        
        // Compute pi(beta; U_i[b])
        double U_ib = U_i_B(b);
        double prob = 0; // Initialize
        prob = pi_U_beta_probit_cpp(beta, U_ib);
        
        
        // Compute [Y_i | U_{ib}] and update "val_Z"
        val_Z += ( (Y_i * prob) + ( (1-Y_i) * (1-prob) ) );
      }
      
      // Update "val_Z"
      val_Z = val_Z / B;
      
      // Run another for loop to get components we want
      double val1 = 0; // initialize
      double val2 = 0; // initialize
      
      for (int b = 0 ; b < B ; ++b){ 
        
        // Compute pi(beta; U_i[b])
        double U_ib = U_i_B(b);
        double prob = 0; // Initialize
        prob = pi_U_beta_probit_cpp(beta, U_ib);
        
        
        // Get the two components of the estimating function
        val1 += (2*Y_i-1) * prob * (1 - prob) / val_Z;
        val2 += (2*Y_i-1) * prob * (1 - prob) * U_ib / val_Z;
      }
      
      // Store "val1 * IPW_i" and "val2 * IPW_i"
      val(i,0) = IPW_i * val1 / B;
      val(i,1) = IPW_i * val2 / B;
      
    } else{ // "MeasurementErrorLogi == T"
      
      double val1 = 0; // initialize
      double val2 = 0; // initialize
      
      // No need to do any integration, proceed to compute contribution to the score function
      double prob = 0; // Initialize
      prob = pi_U_beta_probit_cpp(beta, Xstar_i);
      
      double val_Z = (Y_i * prob) + ( (1-Y_i) * (1-prob) );
      
      // Get the two components of the estimating function
      val1 = (2*Y_i-1) * prob * (1 - prob) / val_Z;
      val2 = (2*Y_i-1) * prob * (1 - prob) * Xstar_i / val_Z;
      // val1 = (Y_i - prob);
      // val2 = (Y_i - prob) * Xstar_i;
      
      // Store val1 and val2
      val(i,0) = IPW_i * val1; 
      val(i,1) = IPW_i * val2;
      
    }
  } // end of for loop (1:n)
  return(val);
}