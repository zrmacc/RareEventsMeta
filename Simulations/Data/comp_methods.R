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
    or_MH_cc <-  c(exp(or$lower.fixed), exp(or$upper.fixed),
                   or$pval.fixed,
                   exp(or$upper.fixed) - exp(or$lower.fixed))
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
    or_MH <-  c(exp(or$lower.fixed), exp(or$upper.fixed),
                or$pval.fixed,
                exp(or$upper.fixed) - exp(or$lower.fixed))
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
    or_peto_fixed <-  c(exp(or$lower.fixed), exp(or$upper.fixed),
                        or$pval.fixed,
                        exp(or$upper.fixed) - exp(or$lower.fixed))
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
    or_dl <- c(exp(or$lower.random), exp(or$upper.random),
               or$pval.random,
               exp(or$upper.random) - exp(or$lower.random))
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
    or_peto_rand <- c(exp(or$lower.random), exp(or$upper.random),
                      or$pval.random,
                      exp(or$upper.random) - exp(or$lower.random))
  }

  all_CIs <- rbind(or_MH_cc,
                   or_MH,
                   or_peto_fixed,
                   or_peto_rand,
                   or_dl
  )

  all_CIs_e <- cbind(all_CIs,
                     sapply(1:nrow(all_CIs), function(xx)
                       IncludeNull(all_CIs[xx, ],
                                   null_val = 1)))

  return(all_CIs_e)
}


#' Comparison methods.
IncludeNull <- function(CI, null_val = log(1)){

  lower_less <- I(CI[1] <= null_val) * 1
  upper_more <- I(CI[2] >= null_val) * 1

  return(lower_less * upper_more)
}

