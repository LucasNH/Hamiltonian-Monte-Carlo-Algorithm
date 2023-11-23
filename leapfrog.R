# This file implements the Leapfrog method for numerical integration

leapfrog <- function(q, p, dVdq, path_len, step_size) {
  #' Leapfrog integrator for Hamiltonian Monte Carlo.
  #'
  #' Returns a new position and momentum
  #'
  #' @param q Initial position
  #' @param p Initial momentum
  #' @param dVdq Gradient of the velocity
  #' @param path_len How long to integrate for (L)
  #' @param step_size How long each integration step should be (epsilon)

  # Making copy to avoid mutation
  positions <- list(q)
  momentums <- list(p)
  stages <- list(list(q, p))

  velocity <- dVdq(q)
  for (i in seq_len(round(path_len / step_size))) {
    p <- p - step_size * velocity / 2  # half step
    stages <- c(stages, list(list(q, p)))
    q <- q + step_size * p  # whole step
    stages <- c(stages, list(list(q, p)))
    positions <- c(positions, list(q))
    velocity <- dVdq(q)
    p <- p - step_size * velocity / 2  # half step
    stages <- c(stages, list(list(q, p)))
    momentums <- c(momentums, list(p))
  }

  # momentum flip at end of steps
  return(list(q, -p, array(positions), array(momentums), array(stages)))
}