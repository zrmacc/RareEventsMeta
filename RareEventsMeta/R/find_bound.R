# Updated: 2021-03-16

# -----------------------------------------------------------------------------

#' Initialize Mu
#' 
#' Initializes the search for the boundary mu. If `mu0` is provided, directly
#' assesses the p-value. If not, a moment estimate of mu is obtained, and the
#' search is initialized from the boundary of the confidence interval for the
#' moment estimator.
#' 
#' @param data Study sizes and event counts.
#' @param mu_to_pval_map Function that accepts a value of mu and
#'   returns a p-value.
#' @param lower_bound Logical, searches for the lower bound? 
#' @param mu0 Optional initial value for mu.
#' @param step_size Distance between successive estimates of mu.
#' @param t1e Type I error level.
#' @return List containing the current 'bound', the 'mc' results, 
#'   and the 'pval' associated with the current bound.

InitMu <- function(
  data,
  mu_to_pval_map,
  lower_bound = FALSE,
  mu0 = NULL, 
  step_size = 0.0002,
  t1e = 0.05
) {
  
  # Try initial value of mu.
  if (!is.null(mu0)) {
    current_bound <- mu0
    current_mc <- mu_to_pval_map(mu = current_bound)
    current_p <- current_mc$p
    if (current_p <= t1e) {
      stop("Search was initiated outside the boundary.\n")
    }
  } else {
    
    # Observed moment estimators.
    obs_est <- MomentEst(
      data$size_1, 
      data$events_1, 
      data$size_2, 
      data$events_2
    )
    
    # Moment estimator CI.
    crit <- qnorm(p = 1 - t1e / 2)
    if (lower_bound) {
      mom_lower <- obs_est$mu - crit * obs_est$mu_se2
      current_bound <- max(step_size, mom_lower)
    } else {
      mom_upper <- obs_est$mu + crit * obs_est$mu_se2
      current_bound <- min(mom_upper, 1 - step_size)
    }
    current_mc <- mu_to_pval_map(mu = current_bound)
    current_p <- current_mc$p
    
    if (current_p <= t1e) {
      current_bound <- obs_est$mu
      current_p <- 1
    }
  }
  
  # Output.
  out <- list(
    bound = current_bound,
    mc = current_mc,
    pval = current_p
  )
  return(out)
}

# -----------------------------------------------------------------------------

#' Extend Boundary
#' 
#' Extends the retention region by moving mu along the boundary 
#' by a single step.
#' 
#' @param mu_to_pval_map Function that accepts a value of mu and
#'   returns a p-value.
#' @param search_results Current search results.
#' @param keep_history Keep search history?
#' @param lower_bound Logical, searches for the lower bound? 
#' @param step_size Distance between successive estimates of mu.
#' @return List containing the current 'bound', the 'mc' results, 
#'   and the 'pval' associated with the current bound.

ExtendBoundary <- function(
  mu_to_pval_map,
  search_results,
  keep_history = FALSE,
  lower_bound = FALSE,
  step_size = 0.0002
) {
  
  # Propose new bound.
  if (lower_bound) {
    new_bound <- max(search_results$bound - step_size, 0)
  } else {
    new_bound <- min(search_results$bound + step_size, 1)
  }
  
  # Generate p-value.
  new_mc <- mu_to_pval_map(new_bound)
  new_p <- new_mc$p
  if (keep_history) {
    new_mc <- rbind(search_results$mc, new_mc)
  }
  
  # Output.
  out <- list(
    bound = new_bound,
    mc = new_mc,
    pval = new_p
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Check Boundary Region.
#' 
#' Check a grid of mu and nu values in the region of the final value of mu
#' from the search along the boundary.
#' 
#' @param ab_to_pval_map Function that accepts a vector contaiing alpha and 
#'   beta then returns a p-value.
#' @param search_results Current search results.
#' @param lower_bound Logical, searches for the lower bound? 
#' @param mu_num_extra_steps Number of steps for mu to check beyond boundary.
#' @param nu_num_extra_steps Number of steps for nu to check beyond boundary.
#' @param step_size Distance between successive estimates of mu.
#' @param t1e Type I error level.
#' @return Data.frame containing (alpha, beta) values in the check region
#'   where the p-value exceeded the type I error.

CheckBoundaryRegion <- function(
  ab_to_pval_map,
  search_results,
  lower_bound = FALSE,
  mu_num_extra_steps = 10,
  nu_num_extra_steps = 10,
  step_size = 0.0002,
  t1e = 0.05
) {
  
  # Return an empty data.frame if no extra search steps were requested.
  if ((mu_num_extra_steps == 0) | (nu_num_extra_steps == 0)) {
    return(data.frame())
  }
  
  current_bound <- search_results$bound
  # Extra values of mu.
  if (lower_bound) {
    extra_mu_values <- seq(
      from = current_bound - step_size * mu_num_extra_steps,
      to = current_bound - step_size, 
      by = step_size)
  } else {
    extra_mu_values <- seq(
      from = current_bound + step_size,
      to = current_bound + step_size * mu_num_extra_steps, 
      by = step_size)
  }
  
  # Extra values of nu.
  search_grid <- lapply(X = extra_mu_values, function(x) {
    return(NuSeq(num_nu_vals = nu_num_extra_steps, mu = x))
  })
  search_grid <- do.call(rbind, search_grid)
  
  # Run MC for each (alpha, beta) pair.
  mc_results <- apply(
    X = search_grid, 
    MARGIN = 1, 
    FUN = ab_to_pval_map
  )
  mc_results <- do.call(rbind, mc_results)
  rownames(mc_results) <- NULL
  
  # Filter to results with p-values > the type I error.
  mc_results <- mc_results[mc_results$p >= t1e, ]
  return(mc_results)
}


# -----------------------------------------------------------------------------

#' Bounds for Mu.
#' 
#' Finds the lower/upper boundary for mu via grid search.
#' 
#' @param size_1 Size of study 1.
#' @param events_1 Events in study 1.
#' @param size_2 Size of study 2.
#' @param events_2 Events in study 2.
#' @param lower_bound Logical, searches for the lower bound? 
#' @param reps MC replicates for p-value estimation.
#' @param t1e Type I error level.
#' @param mu0 Optional initial value for mu.
#' @param step_size Distance between successive estimates of mu.
#' @param maxit Maximum number of iterations to perform.
#' @param keep_history Keep search history?
#' @param mu_num_extra_steps Number of steps for mu to check beyond boundary.
#' @param nu_num_extra_steps Number of steps for nu to check beyond boundary.
#' @importFrom stats qnorm
#' @export
#' @return List containing the `search_results` and the output of the
#'   `boundary_region_check`.
#' @examples 
#' set.seed(2013)
#' data <- GenData(
#'   total_studies = 10,
#'   n1 = 100,
#'   n2 = 100,
#'   alpha2 = 10,
#'   beta2 = 10
#' )
#' # Note: use more `reps` and smaller `step_size` for more accurate results.
#' 
#' # Find lower bound.
#' FindBound(
#'   size_1 = data$size_1,
#'   events_1 = data$events_1,
#'   size_2 = data$size_2,
#'   events_2 = data$events_2,
#'   reps = 50,
#'   lower_bound = TRUE,  # Set to FALSE for upper bound.
#'   step_size = 0.02
#' )

FindBound <- function(
  size_1,
  events_1,
  size_2,
  events_2,
  lower_bound = TRUE, 
  reps = 500,
  t1e = 0.05,
  mu0 = NULL,
  step_size = 0.0002,
  maxit = 100,
  keep_history = FALSE,
  mu_num_extra_steps = 10,
  nu_num_extra_steps = 10
) {
  
  # Bundle data.
  data <- data.frame(size_1, events_1, size_2, events_2)
  
  # Function that maps directly from mu to output of RunMC.
  mu_to_pval_map <- function(mu) {
    out <- RunMC(
          size_1 = size_1,
          events_1 = events_1,
          size_2 = size_2,
          events_2 = events_2,
          reps = reps,
          mu = mu
        )
    return(out)
  }
  
  # Function that maps from (alpha, beta) to output of RunMC.
  ab_to_pval_map <- function(ab) {
      out <- RunMC(
        size_1 = size_1,
        events_1 = events_1,
        size_2 = size_2,
        events_2 = events_2,
        reps = reps,
        alpha = ab[1],
        beta = ab[2]
      )
    return(out)
  }
  
  # -------------------------------------------------------
  
  # Initialize search.
  search_results <- InitMu(
    data = data,
    mu_to_pval_map = mu_to_pval_map,
    lower_bound = lower_bound,
    mu0 = mu0,
    step_size = step_size,
    t1e = t1e
  )
  
  # Extend boundary.
  for(i in 1:maxit) {
    search_results <- ExtendBoundary(
      mu_to_pval_map = mu_to_pval_map,
      search_results = search_results,
      keep_history = keep_history,
      lower_bound = lower_bound,
      step_size = step_size
    )
    if (search_results$pval <= t1e) {break}
  }

  # Send warning if max iteration was reached.
  if (i == maxit & search_results$pval > t1e) {
    warning("Maximum iterations performed without reaching lower bound.
            Consider increasing the step_size or initializing the search.\n")
  }
  
  # Check boundary region.
  boundary_check <- CheckBoundaryRegion(
    ab_to_pval_map = ab_to_pval_map,
    search_results,
    lower_bound = lower_bound,
    mu_num_extra_steps = mu_num_extra_steps,
    nu_num_extra_steps = nu_num_extra_steps,
    step_size = step_size,
    t1e = t1e
  )
  
  # Output.
  out <- list(
    search_results = search_results,
    boundary_region_check = boundary_check
  )
  
  return(out)
}
