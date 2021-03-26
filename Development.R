# Note: reinstall package and run devtools::check() after modifying
# scripts in RareEventsMeta/R. 
library(RareEventsMeta)

# Settings.
t1e <- 0.05
step_size <- 0.005
maxit <- 10
reps <- 50
lower_bound <- TRUE
keep_history <- TRUE
mu_extra_steps <- 3
nu_extra_steps <- 3

# Generate data.
set.seed(2013)
data <- GenData(
  total_studies = 10,
  n1 = 100,
  n2 = 100,
  alpha2 = 10,
  beta2 = 10
)

study <- data$study
size_1 <- data$size_1
events_1 <- data$events_1
size_2 <- data$size_2
events_2 <- data$events_2

# Moment estimator
moments <- MomentEst(size_1, events_1, size_2, events_2)

# Lower bound
lower <- FindBound(
  size_1 = size_1,
  events_1 = events_1,
  size_2 = size_2,
  events_2 = events_2,
  lower_bound = TRUE,
  keep_history = FALSE,
  reps = reps,
  t1e = t1e,
  step_size = step_size,
  maxit = 50,
  mu_extra_steps = 5,
  nu_extra_steps = 5
)

# Upper bound
upper <- FindBound(
  size_1 = size_1,
  events_1 = events_1,
  size_2 = size_2,
  events_2 = events_2,
  lower_bound = FALSE,
  keep_history = FALSE,
  reps = reps,
  t1e = t1e,
  step_size = step_size,
  maxit = 50,
  mu_extra_steps = 5,
  nu_extra_steps = 5
)

# Confidence interval.
ci <- ExactConfInt(
  size_1 = size_1,
  events_1 = events_1,
  size_2 = size_2,
  events_2 = events_2,
  reps = 50,
  step_size = 0.02,
  mu_extra_steps = 5,
  nu_extra_steps = 5
)
show(ci)
