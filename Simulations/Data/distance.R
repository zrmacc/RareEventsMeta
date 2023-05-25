# Covid analysis

setwd('~/Documents/GitHub/RareEventsMeta/Simulations/Data')


my_data<- read.csv('distance.csv')
n1 <-  my_data[, 'Further_Total']
e1 <- my_data[, 'Further_Events']
n2 <-  my_data[, 'Shorter_Total']
e2 <- my_data[, 'Shorter_Events']

distance_data <- cbind(n1, e1, n2, e2)
dim(distance_data)
colnames(distance_data) <- c("size_1", "events_1",
                                "size_2", "events_2")


# Library for comparison methods.
library(meta)
source("comp_methods.R")

##### ##### ##### ########
##### distance Analysis #####
##### ##### ##### ########

# Run comparison methods.
distance_comp <- CompMethods(distance_data)
distance_comp

distance_data_dzr <- distance_data[ ! ((distance_data[, 'events_1'] == 0) &
                                   (distance_data[, 'events_2'] == 0)),]
dim(distance_data_dzr)

# Run XRRmeta - change defaults for this and adjust nu stepsize, add pvalue option.

# Library.
library(RareEventsMeta)
MomentEst(events_1 = distance_data_dzr[, 'events_1'],
          size_1 = distance_data_dzr[, 'size_1'],
          events_2 = distance_data_dzr[, 'events_2'],
          size_2 = distance_data_dzr[, 'size_2'])

t0 <- proc.time()
distance_xrrmeta <- ExactConfInt(
  events_1 = distance_data_dzr[, 'events_1'],
  size_1 = distance_data_dzr[, 'size_1'],
  events_2 = distance_data_dzr[, 'events_2'],
  size_2 = distance_data_dzr[, 'size_2'],
  reps = 2000,
  step_size = 0.001,
  maxit = 500,
  mu_extra_steps = 10,
  nu_extra_steps = 10
)
t1 <- proc.time()
elapsed <- t1-t0
distance_xrrmeta
elapsed

# https://onlinelibrary.wiley.com/doi/full/10.1002/jrsm.1521
# https://www.thelancet.com/article/S0140-67362031142-9/fulltext