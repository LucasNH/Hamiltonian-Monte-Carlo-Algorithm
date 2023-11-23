# This file implements the Hamiltonian Monte Carlo Algorithm

source("leapfrog.R")

leapfrog <- leapfrog

hamiltonian_monte_carlo <- function(n_samples, potential, initial_position,
  initial_potential = NULL, initial_potential_grad = NULL, path_len = 1,
  step_size = 0.1, integrator = leapfrog, max_energy_change = 1000.0,
  do_reject = TRUE
) {
  #' Run Hamiltonian Monte Carlo sampling.
  #'
  #' @param n_samples Number of samples to return
  #' @param negative_log_prob The negative log probability from which to sample
  #' @param initial_position A place from which to start sampling.
  #' @param path_len How long each integration path is. Smaller is faster and
  #' more correlated (L).
  #' @param step_size How long each integration step is. Smaller is slower and
  #' more accurate (epsilon).
  #' @param integrator Integrator to use, from `leapfrog.R`
  #' @param do_reject Turn off metropolis correction. Not valid MCMC if False!
  initial_position <- array(initial_position)
  negative_log_prob <- function(q) potential(q)[1]  # NOQA
  dVdq <- function(q) potential(q)[2]  # NOQA

  # collect all our samples in a list
  samples <- list(initial_position)
  sample_positions <- list()
  sample_momentums <- list()
  accepted <- list()
  p_accepts <- list()

  # If initial_position is a 10d vector and n_samples is 100, we want 100 x 10
  # momentum draws. We can do this in one call to rnorm, and iterate over rows
  size <- c(n_samples, dim(initial_position)[1])

  # Keep a single object for momentum resampling
  momentum <- rnorm(size, 0, 1)
  for (p0 in momentum) {
    # Integrate over our path to get a new position and momentum
    result <- integrator(samples[length(samples)], p0, dVdq,
      path_len = path_len, step_size = step_size
    )
    q_new <- result$q_new
    p_new <- result$p_new
    positions <- result$positions
    momentums <- result$momentums

    sample_positions <- c(sample_positions, list(positions))
    sample_momentums <- c(sample_momentums, list(momentums))

    # Check Metropolis acceptance criterion
    start_log_p <-
      negative_log_prob(samples[length(samples)]) - sum(dnorm(p0, 0, 1))
    new_log_p <- negative_log_prob(q_new) - sum(dnorm(p_new, 0, 1))
    energy_change <- start_log_p - new_log_p
    p_accept <- exp(energy_change)

    if (runif(1) < p_accept) {
      samples <- c(samples, list(q_new))
      accepted <- c(accepted, TRUE)
    } else {
      if (do_reject) {
        samples <- c(samples, list(samples[length(samples)]))
      } else {
        samples <- c(samples, list(q_new))
      }
      accepted <- c(accepted, FALSE)
    }
    p_accepts <- c(p_accepts, p_accept)
  }
  return(list(
    samples = array(samples[-1]),
    sample_positions = array(sample_positions),
    sample_momentums = array(sample_momentums),
    accepted = array(accepted),
    p_accepts = array(p_accepts)
  ))
}