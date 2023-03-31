my_studies <- 48

# 0.003 - 24, 48, 96
# 0.01, 0.03 - 12, 24, 48

setwd("/Users/jgrons/Documents/GitHub/RareEventsMeta/Simulations/Rscripts")

# Setting 1A
my_alpha <- 1.45
my_multiplier <- 1
my_rate <- 0.003
source('comp_only_jesse.R')

# Setting 1B
my_alpha <- 1.45
my_multiplier <- 1
my_rate <- 0.01
source('comp_only_jesse.R')

# Setting 1C
my_alpha <- 1.45
my_multiplier <- 1
my_rate <- 0.03
source('comp_only_jesse.R')

# Setting 1C
my_alpha <- 1.1
my_multiplier <- 1.5
my_rate <- 0.003
source('comp_only_jesse.R')

# Setting 1D
my_alpha <- 1.1
my_multiplier <- 1.5
my_rate <- 0.01
source('comp_only_jesse.R')

# Setting 1E
my_alpha <- 1.1
my_multiplier <- 1.5
my_rate <- 0.03
source('comp_only_jesse.R')

# Setting 2A
my_alpha <- 145
my_multiplier <- 1
my_rate <- 0.003
source('comp_only_jesse.R')

# Setting 2B
my_alpha <- 145
my_multiplier <- 1
my_rate <- 0.01
source('comp_only_jesse.R')

# Setting 2C
my_alpha <- 145
my_multiplier <- 1
my_rate <- 0.03
source('comp_only_jesse.R')

# Setting 2D
my_alpha <- 110
my_multiplier <- 1.5
my_rate <- 0.003
source('comp_only_jesse.R')

# Setting 2E
my_alpha <- 110
my_multiplier <- 1.5
my_rate <- 0.01
source('comp_only_jesse.R')

# Setting 2R
my_alpha <- 110
my_multiplier <- 1.5
my_rate <- 0.03
source('comp_only_jesse.R')


# alpha <- 145
# beta <- 145
# calc_params <- function(alpha, beta){
#
#   mu <- alpha/(alpha + beta)
#   het <- mu * (1-mu) * 1/(alpha + beta + 1)
#   return(c(mu, het))
# }
# calc_params(alpha, beta)


