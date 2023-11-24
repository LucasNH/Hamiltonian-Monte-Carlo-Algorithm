grad_U = function(q, type, alpha = 50, beta = 50){
  # U represents a distribution, and grad_U is the gradient of U.
  # q: represents a sample from the distribution of U.
  # type: the type of distribution the user is sampling from.
  # alpha: the alpha parameter if the user samples from a beta dist.
  # beta: the beta parameter if the user samples from a beta dist.
  if(type == "beta"){
    return(-(alpha/q - beta/(1-q)))
  } else if (type == "norm"){
    return(q)
  }
}

HMC = function(n_samples, U_type, epsilon, L, current_q,
               alpha = 50, beta = 50, mu = 0, sigma = 1){
  # n_samples: the number of samples we want to obtain.
  # U_type: the type of distribution we are sampling from. (For now, beta/normal)
  # epsilon: the stepsize.
  # L: number of steps per stepsize.
  # current_q: the start value for the representative sample of our target dist.
  # alpha, beta, mu, sigma: parameters of the U_type.
  
  q = current_q # q represents a sample from our target distribution.
  D = numeric(length = n_samples) # D contains our n_samples
  
  for(j in 1:n_samples){
    # Here, our kinetic energy is a standard normal distribution.
    p = rnorm(n = 1, mean = 0, sd = 1) 
    current_p = p
    ######################################################
    # Leapfrog Algorithm for solving the dynamical system
    # Make a half step for momentum at the beginning
    p = p - epsilon * grad_U(q, U_type, alpha, beta) / 2
    # Alternate full steps for position and momentum
    for (i in 1:L){
      # Make a full step for the position
      q = q + epsilon*p
      # Make a full step for the momentum, except at end of trajectory
      if (i!=L){ p = p - epsilon * grad_U(q, U_type, alpha, beta)}
    }
    # Make a half step for momentum at the end.
    p = p - epsilon * grad_U(q, U_type, alpha, beta) / 2
    ######################################################
    # Evaluate potential and kinetic energies at start and end of trajectory
    if(U_type == "beta"){
      current_U = -dbeta(current_q, alpha, beta, log = TRUE)
      proposed_U = -dbeta(q, alpha, beta, log = TRUE)
    } else if (U_type == "norm"){
      current_U = -dnorm(current_q, mu, sigma, log = TRUE)
      proposed_U = -dnorm(q, mu, sigma, log = TRUE)
    }
    current_K = sum(current_p^2)/2
    proposed_K = sum(p^2)/2
    # This part corresponds with the acceptance-rejection method.
    if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)){
      D[j] = q
      current_q = q # accept, need to change q value
    } 
    else {
      D[j] = current_q
    } 
  }
  return(D)
}


# below are some one-off simulations that work with our code.
par(mfrow=c(1,2))
set.seed(012920)
HMC_beta_samples = HMC(n_samples = 2000, U_type = "beta", epsilon = 0.01, L = 2, 
                  current_q = 0.3, alpha = 50, beta = 50)
hist(HMC_beta_samples, breaks = 'scott', freq = FALSE,
     border = "#ffffff",
     col = "#81ddff",
     xlab = "Values of HMC Samples",
     main = "Histogram of HMC Samples (Beta Dist)")
curve(dbeta(x, 50, 50), col = "#FF6666", add = TRUE, lwd = 2)
set.seed(012920)
HMC_norm_samples = HMC(n_samples = 2000, U_type = "norm", epsilon = 0.2, L = 2, 
                       current_q = 0.3, mu = 0, sigma = 1)
hist(HMC_norm_samples, breaks = 'scott', freq = FALSE,
     border = "#ffffff",
     col = "#f79ee1",
     xlab = "Values of HMC Samples",
     main = "Histogram of HMC Samples (Normal Dist)")
curve(dnorm(x,0, 1), col = "#6699FF", add = TRUE, lwd = 2)



