#' Comparison methods.
CompMethods <- function(data){

  # ------------------------------------------------ #
  # Comparison to existing fixed effects approaches. #
  # ------------------------------------------------ #

  # MH odds ratio with continuity correction (include DZ studies).
  or <- tryCatch(metabin(data[,"events_1"], data[, "size_1"],
                         data[,"events_2"], data[, "size_2"],
                         sm = "OR",
                         allstudies = TRUE,
                         control=list(stepadj=0.25, maxiter=1000)), # improve convergence of Fisher scoring
                 error = function(e){
                   return(rep(NA, 2))
                 })
  if(is.na(or[1])){
    or_MH_cc <- or
  }else{
    or_MH_cc <- c(or$lower.fixed, or$upper.fixed)
  }


  # MH odds ratio without continuity correction.
  or <- tryCatch(metabin(data[,"events_1"], data[, "size_1"],
                         data[,"events_2"], data[, "size_2"],
                         sm = "OR",
                         MH.exact = TRUE,
                         control=list(stepadj=0.25, maxiter=1000)),
                 error = function(e){
                   return(rep(NA, 2))
                 })

  if(is.na(or[1])){
    or_MH <- or
  }else{
    or_MH <- c(or$lower.fixed, or$upper.fixed)
  }

  # Peto method for odds ratio, fixed effects.
  or <- tryCatch(metabin(data[,"events_1"], data[, "size_1"],
                         data[,"events_2"], data[, "size_2"],
                         sm= "OR",
                         method = "Peto",
                         control=list(stepadj=0.25, maxiter=1000)),
                 error = function(e){
                   return(rep(NA, 2))
                 })
  if(is.na(or[1])){
    or_peto_fixed <- or
  }else{
    or_peto_fixed <- c(or$lower.fixed, or$upper.fixed)
  }


  # ------------------------------------------------- #
  # Comparison to existing random effects approaches. #
  # ------------------------------------------------- #

  # DL method for odds ratio with continuity correction.
  or <- tryCatch(metabin(data[,"events_1"], data[, "size_1"],
                         data[,"events_2"], data[, "size_2"],
                         sm = "OR",
                         allstudies = TRUE,
                         random = TRUE,
                         control=list(stepadj=0.25, maxiter=1000)),
                 error = function(e){
                   return(rep(NA, 2))
                 })

  if(is.na(or[1])){
    or_dl <- or
  }else{
    or_dl <- c(or$lower.random, or$upper.random)
  }

  # Peto method for odds ratio, random effects.
  or <- tryCatch(metabin(data[,"events_1"], data[, "size_1"],
                         data[,"events_2"], data[, "size_2"],
                         sm= "OR",
                         method = "Peto",
                         control=list(stepadj=0.25, maxiter=1000)),
                 error = function(e){
                   return(rep(NA, 2))
                 })
  if(is.na(or[1])){
    or_peto_rand <- or
  }else{
    or_peto_rand <- c(or$lower.random, or$upper.random)
  }

  all_CIs <- rbind(or_MH_cc,
                   or_MH,
                   or_peto_fixed,
                   or_peto_rand,
                   or_dl
  )

  all_CIs_e <- cbind(all_CIs,
                     sapply(1:nrow(all_CIs), function(xx)
                       IncludeNull(all_CIs[xx, ])))

  return(all_CIs_e)
}
