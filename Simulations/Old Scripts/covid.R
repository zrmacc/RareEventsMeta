# Covid analysis

setwd('~/Documents/GitHub/RareEventsMeta/Simulations/Data')


my_data<- read.csv('covid_death.csv')
n1 <-  my_data[, 'size_m']
e1_death <- my_data[, 'deaths_m']
n2 <-  my_data[, 'size_f']
e2_death <- my_data[, 'deaths_f']
death_data <- cbind(n1, e1_death, n2, e2_death)

my_data<- read.csv('covid_itu.csv')
n1 <-  my_data[, 'size_m']
e1_itu <- my_data[, 'itu_m']
n2 <-  my_data[, 'size_f']
e2_itu <- my_data[, 'itu_f']
itu_data <- cbind(n1, e1_itu, n2, e2_itu)

colnames(death_data) <- colnames(itu_data) <- c("size_1", "events_1",
                                             "size_2", "events_2")

dim(death_data)
dim(itu_data)

# Library for comparison methods.
library(meta)
source("comp_methods.R")

##### ##### ##### ########
##### Death Analysis #####
##### ##### ##### ########

# Run comparison methods.
death_comp <- CompMethods(death_data)
death_comp

death_data_dzr <- death_data[ ! ((death_data[, 'events_1'] == 0) &
                                   (death_data[, 'events_2'] == 0)),]
dim(death_data_dzr)

# Run XRRmeta - change defaults for this and adjust nu stepsize, add pvalue option.
t0 <- proc.time()
death_xrrmeta <- ExactConfInt(
  events_1 = death_data_dzr[, 'events_1'],
  size_1 = death_data_dzr[, 'size_1'],
  events_2 = death_data_dzr[, 'events_2'],
  size_2 = death_data_dzr[, 'size_2'],
  reps = 2000,
  step_size = 0.001,
  maxit = 500,
  mu_extra_steps = 10,
  nu_extra_steps = 10
)
t1 <- proc.time()
elapsed <- t1-t0
death_xrrmeta
elapsed


#####  ##### #####  #####
#####  itu Analysis #####
#####  ##### #####  #####

# Run comparison methods.
itu_comp <- CompMethods(itu_data)
itu_comp

itu_data_dzr <- itu_data[ ! ((itu_data[, 'events_1'] == 0) &
                               (itu_data[, 'events_2'] == 0)),]
dim(itu_data_dzr)
# Run XRRmeta.
t0 <- proc.time()
itu_xrrmeta <- ExactConfInt(
  events_1 = itu_data_dzr[, 'events_1'],
  size_1 = itu_data_dzr[, 'size_1'],
  events_2 = itu_data_dzr[, 'events_2'],
  size_2 =itu_data_dzr[, 'size_2'],
  reps = 2000,
  step_size = 0.001,
  maxit = 500,
  mu_extra_steps = 10,
  nu_extra_steps = 10
)
t1 <- proc.time()
elapsed <- t1-t0
itu_xrrmeta
