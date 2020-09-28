FHT_parametric_survival <- function(time, mu, y0) {
  return(pnorm((mu*time + y0)/sqrt(time)) - exp(-2*y0*mu)*pnorm((mu*time - y0)/sqrt(time)))
}

brier_score <- function(ts, deltas, censoring_estimates, t_star, estimated_survival_t_star) {
  ### ts, deltas, censoring_estimated are vectors where index i of each vector correspond to observation i
  before_or_at <- ts <= t_star
  after <- ts > t_star
  before_or_at_and_observed <- before_or_at * deltas
  censoring_estimate_t_star <- censoring_estimates[i]
  left_term <- ((estimated_survival_t_star^2)*before_or_at_and_observed)/censoring_estimates
  right_term <- ((1-estimated_survival_t_star)^2*after)/censoring_estimate_t_star
  return(mean(left_term + right_term))
}

calculate_at_risk <- function(times, survival_times) {
  at_risk <- sapply(times, function(ti) { sum(survival_times>=ti) })
  return(at_risk)
}

kaplan_meier_estimate_of_censoring <- function(time, times, delta) {
  indexes_to_look_at <- which(times <= time)
  if (sum(indexes_to_look_at) == 0) {
    return(1)
  }
  times_to_look_at <- times[indexes_to_look_at]
  deltas_to_look_at <- delta[indexes_to_look_at]
  at_risk <- calculate_at_risk(times_to_look_at, times)
  terms <- 1 - (1-deltas_to_look_at)/at_risk
  return(prod(terms))
}


brier_score_many_observations_fht <- function(times, delta, y0s, mus) {
  # Assumes these are vectors, where index i correspond to observation i for each vector
  # Assumes times are sorted by time ascending
  censoring_estimates <- sapply(times, function(time) {
    kaplan_meier_estimate_of_censoring(time, times, delta)
  })
  n <- length(times)
  
  # Drop the last observation
  times <- times[-n]
  delta <- delta[-n]
  y0s <- y0s[-n]
  mus <- mus[-n]
  censoring_estimates <- censoring_estimates[-n]
  N <- length(times)

  # Generate matrix of estimated survival probabilities, where
  # element [i, j] corresponds to individual i's probability at time j
  estimated_survival_probabilities <- matrix(nrow=N, ncol=N)
  for (j in 1:N) {
    estimated_survival_probabilities[, j] <- FHT_parametric_survival(time=times[j], mu=mus, y0=y0s)
  }
  
  brier_scores <- rep(NA, N)
  for (i in 1:N) {
    brier_scores[i] <- brier_score(
      ts=times,
      deltas=delta,
      censoring_estimates=censoring_estimates,
      t_star=current_time,
      estimated_survival_t_star=estimated_survival_probabilities[, i]
    )
  }
  return(
    data.frame(times=times, brier_scores=brier_scores)
  )
}

brier_score_many_observations <- function(times, delta, estimated_survival_probabilities) {
  # Assumes times and delta are vectors, where index i correspond to observation i for each vector
  # Assumes times are sorted by time ascending
  # Assumes estimated_survival_probabilities is a matrix of estimated survival probabilities,
  # where element [i, j] corresponds to individual i's probability at time j
  censoring_estimates <- sapply(times, function(time) {
    kaplan_meier_estimate_of_censoring(time, times, delta)
  })
  n <- length(times)
  
  # Drop the last observation
  times <- times[-n]
  delta <- delta[-n]
  censoring_estimates <- censoring_estimates[-n]
  estimated_survival_probabilities <- estimated_survival_probabilities[-n, -n]
  N <- length(times)
  
  brier_scores <- rep(NA, N)
  for (i in 1:N) {
    brier_scores[i] <- brier_score(
      ts=times,
      deltas=delta,
      censoring_estimates=censoring_estimates,
      t_star=current_time,
      estimated_survival_t_star=estimated_survival_probabilities[, i]
    )
  }
  return(
    data.frame(times=times, brier_scores=brier_scores)
  )
}

