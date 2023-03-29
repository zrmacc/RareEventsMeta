# Rosiglitazone analysis

setwd('~/Documents/GitHub/RareEventsMeta/Simulations/Data')
my_data <- read.csv('rosiglitazone.csv')


n1 <-  my_data[, 'n.t']
e1_mi <- my_data[, 'e.mi.t']
e1_cvd <- my_data[, 'e.death.t']

n2 <-  my_data[, 'n.c']
e2_mi <- my_data[, 'e.mi.c']
e2_cvd <- my_data[, 'e.death.c']

mi_data <- cbind(n1, e1_mi, n2, e2_mi)
cvd_data <- cbind(n1, e1_cvd, n2, e2_cvd)
colnames(mi_data) <- colnames(cvd_data) <- c("size_1", "events_1",
                                             "size_2", "events_2")
# Library.
library(RareEventsMeta)

# Library for comparison methods.
library(meta)
source("comp_methods.R")

##### ##### ##### #####
##### MI Analysis #####
##### ##### ##### #####

# Run comparison methods.
mi_comp <- CompMethods(mi_data)
exp(mi_comp)

mi_data_dzr <- mi_data[ ! ((mi_data[, 'events_1'] == 0) &
                             (mi_data[, 'events_2'] == 0)),]
dim(mi_data_dzr)

# Run XRRmeta - change defaults for this and adjust nu stepsize, add pvalue option.
t0 <- proc.time()
mi_xrrmeta <- ExactConfInt(
  events_1 = mi_data_dzr[, 'events_1'],
  size_1 = mi_data_dzr[, 'size_1'],
  events_2 = mi_data_dzr[, 'events_2'],
  size_2 = mi_data_dzr[, 'size_2'],
  reps = 2000,
  step_size = 0.001,
  maxit = 500,
  mu_extra_steps = 10,
  nu_extra_steps = 10
)
t1 <- proc.time()
elapsed <- t1-t0
mi_xrrmeta
elapsed
#####  ##### #####  #####
#####  CVD Analysis #####
#####  ##### #####  #####

# Run comparison methods.
cvd_comp <- CompMethods(cvd_data)
exp(cvd_comp)

cvd_data_dzr <- cvd_data[ ! ((cvd_data[, 'events_1'] == 0) &
                               (cvd_data[, 'events_2'] == 0)),]
dim(cvd_data_dzr)
# Run XRRmeta.
t0 <- proc.time()
cvd_xrrmeta <- ExactConfInt(
  events_1 = cvd_data_dzr[, 'events_1'],
  size_1 = cvd_data_dzr[, 'size_1'],
  events_2 = cvd_data_dzr[, 'events_2'],
  size_2 =cvd_data_dzr[, 'size_2'],
  reps = 2000,
  step_size = 0.001,
  maxit = 500,
  mu_extra_steps = 10,
  nu_extra_steps = 10
)
t1 <- proc.time()
elapsed <- t1-t0
cvd_xrrmeta
