grad_U <- function(q, type, alpha = 50, beta = 50) {
  #' U represents a distribution, and grad_U is the gradient of U.
  #'
  #' @param q represents a sample from the distribution of U.
  #' @param type the type of distribution the user is sampling from.
  #' @param alpha the alpha parameter if the user samples from a beta dist.
  #' @param beta the beta parameter if the user samples from a beta dist.

  if ("beta" == type) {
    return(-(alpha / q - beta / (1 - q)))
  } else if ("norm" == type) {
    return(q)
  }
}

HMC <- function(n_samples, U_type, epsilon, L, current_q, alpha = 50, beta = 50,
                mu = 0, sigma = 1) {
  #' @param n_samples the number of samples we want to obtain.
  #' @param U_type the type of distribution we are sampling from. (For now,
  #' beta/normal)
  #' @param epsilon the stepsize.
  #' @param L number of steps per stepsize.
  #' @param current_q the start value for the representative sample of our
  #' target dist.
  #' @param alpha parameter of the U_type,
  #' @param beta parameter of the U_type,
  #' @param mu parameter of the U_type,
  #' @param sigma parameter of the U_type.

  q <- current_q # q represents a sample from our target distribution.
  D <- numeric(length = n_samples) # D contains our n_samples

  for (j in 1:n_samples){
    # Here, our kinetic energy is a standard normal distribution.
    p <- rnorm(n = 1, mean = 0, sd = 1)
    current_p <- p

    ######################################################
    # Leapfrog Algorithm for solving the dynamical system
    # Make a half step for momentum at the beginning
    p <- p - epsilon * grad_U(q, U_type, alpha, beta) / 2
    # Alternate full steps for position and momentum
    for (i in 1:L) {
      # Make a full step for the position
      q <- q + epsilon * p
      # Make a full step for the momentum, except at end of trajectory
      if (i != L) {
        p <- p - epsilon * grad_U(q, U_type, alpha, beta)
      }
    }
    # Make a half step for momentum at the end.
    p <- p - epsilon * grad_U(q, U_type, alpha, beta) / 2
    ######################################################

    # Evaluate potential and kinetic energies at start and end of trajectory
    if ("beta" == U_type) {
      current_U <- -dbeta(current_q, alpha, beta, log = TRUE)
      proposed_U <- -dbeta(q, alpha, beta, log = TRUE)
    } else if ("norm" == U_type) {
      current_U <- -dnorm(current_q, mu, sigma, log = TRUE)
      proposed_U <- -dnorm(q, mu, sigma, log = TRUE)
    }
    current_K <- sum(current_p^2) / 2
    proposed_K <- sum(p^2) / 2
    # This part corresponds with the acceptance-rejection method.
    if (runif(1) < exp(current_U - proposed_U + current_K - proposed_K)) {
      D[j] <- q
      current_q <- q # accept, need to change q value
    } else {
      D[j] <- current_q
    }
  }
  return(D)
}