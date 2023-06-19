# Purpose: Moment estimator. 

# -----------------------------------------------------------------------------
# Data preparation.
# -----------------------------------------------------------------------------

#' Prepare Data for Simulation
#' 
#' @param size_1 Size of study 1.
#' @param events_1 Events in study 1.
#' @param size_2 Size of study 2.
#' @param events_2 Events in study 2.
#' @param study Optional study identifier.
#' 
#' @importFrom stats dhyper
#' @return Data.frame containing
#' \itemize{
#'   \item `study` identifier.
#'   \item `events` down-sampled number of events.
#'   \item `total_events` total number of events.
#'   \item `weight` based on the hypergeometric distribution.
#'   \item `n1` and `n2`, original study sizes.
#' }

PrepData <- function(
    size_1,
    events_1,
    size_2,
    events_2,
    study
) {
  
  # Loop over studies.
  loop <- function(i) {
    
    ## More patients in group 1.
    if (size_1[i] >= size_2[i]) {
      if (events_1[i] == 0 & events_2[i] != 0) {
        
        # Only 1 way to sample zero events from group 1.
        weight <- 1
        events <- events_1[i]
        total_events <- events_1[i] + events_2[i]
        
      } else {
        
        # Resample n2 patients and record weight associated with each value of
        # number of events observed according to the hypergeometric distribution.
        weight <- dhyper(
          0:events_1[i], events_1[i],
          size_1[i] - events_1[i], size_2[i]
        )
        events <- 0:events_1[i]
        total_events <- events_2[i] + events
        
      }
      
      ## More patients in group 2.
    } else {
      if (events_1[i] != 0 & events_2[i] == 0) {
        
        # Only one way to sample zero events in group 2.
        weight <- 1
        events <- events_1[i]
        total_events <- events_1[i] + events_2[i]
        
      } else {
        
        # Resample n1 patients and record weight associated with each value of
        # number of events observed according to the hypergeometric distribution
        weight <- dhyper(
          0:events_2[i], events_2[i],
          size_2[i] - events_2[i], size_1[i]
        )
        events <- events_1[i]
        total_events <- 0:events_2[i] + events
      }
      
    }
    
    out <- data.frame(
      cbind(
        "study" = study[i],
        "events" = events, 
        "total_events" = total_events, 
        "weight" = weight, 
        "n1" = size_1[i], 
        "n2" = size_2[i]
      )
    )
    
    ## If we have zero in the smaller group we employ the following
    ## rescaling of the probabilities as a continuity correction
    is_zero <- (out$total_events == 0)
    if (any(is_zero)) {
      out <- out[!is_zero, ]
      out$weight <- out$weight / sum(out$weight)
    }
    
    # Remove studies with weight zero.
    out <- out[out$weight > 0, ]
    
    # Return.
    return(out)
  }
  
  out <- lapply(seq_len(length(size_1)), loop)
  out <- do.call(rbind, out)
  
  return(out)
}


# -----------------------------------------------------------------------------
# Moment estimation.
# -----------------------------------------------------------------------------

#' Moment Estimation
#' 
#' @param size_1 Size of study 1.
#' @param events_1 Events in study 1.
#' @param size_2 Size of study 2.
#' @param events_2 Events in study 2.
#' @param study Optional study identifier.
#'  
#' @export
#' @return List containing:
#' \itemize{
#'   \item `mu`, estimated first moment.
#'   \item `mu_se2`, square standard error of estimated first moment.
#'   \item `mu_cc`, estimated first moment with continuity correction.
#'   \item `mu2`, estimated second moment.
#'   \item `nu`, estimate of nu.
#' }

MomentEst <- function(
    size_1,
    events_1,
    size_2,
    events_2,
    study = NULL,
    corrected = TRUE,
    weighted = TRUE
) {
  
  # Create study identifier if not provided.
  if (is.null(study)) {
    study <- seq_len(length(size_1))
  }
  studies <- length(study)
  
  
  # Prepare data.
  if(weighted){
    
  data <- PrepData(
    size_1 = size_1,
    events_1 = events_1,
    size_2 = size_2,
    events_2 = events_2,
    study = study
  )
  
  }else{
    
    data <- data.frame(cbind(study = 1:studies,
                  events = events_1,
                  total_events = events_1 + events_2,
                  weight = rep(1, studies),
                  n1 = size_1,
                  n2 = size_2
                  ))
  }
  
  # Row first moment.
  mu <- sum(data$weight * data$events / data$total_events) / studies
  
  # Continuity correction.
  if(corrected){
    
    #all(data$events == 0) | all(data$total_events == 1)
    data$total_events <- data$total_events + 1/(data$n1 + data$n2)
    data$events <- data$events + 1/data$n1
    
  }else{
    
    data$total_events <- data$total_events  
    data$events <- data$events  
    
  }

  
  # Continuity corrected first moment.
  mu_cc <- sum(data$weight * data$events / data$total_events) / studies
  
  # Estimate nu using first and second moments.
  num <- (sum(data$weight * (data$events / data$total_events)^2) / studies -
            sum(data$weight * mu_cc / data$total_events) / studies)
  
  denum <- sum(data$weight * (1 - 1 / (data$total_events))) / studies
  
  
  mu2 <- num / denum
  nu <- max(0, mu2 - mu_cc^2)

  
  # Standard error of first moment estimator.
  se2 <- sum(
    data$weight * (mu_cc * (1 - mu_cc) / data$total_events) + 
      data$weight * (1 - 1 / data$total_events) * nu
  ) / studies^2
  
 
  # Output.
  out <- list(
    "mu" = mu,
    "mu_se2" = se2,
    "mu_cc" = mu_cc,
    "mu2" = mu2,
    "nu" = nu
  )
  return(out)
}
