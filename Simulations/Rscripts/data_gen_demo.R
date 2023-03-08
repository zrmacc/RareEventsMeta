# OLD VERSION OF THE CODE

total_studies <- 500
rate1 <- 0.015
alpha1 <- 1.44

alpha2 <-3
beta2 <- 3

n1 <- rep(500, total_studies)
n2 <- rep(1000, total_studies)
# Baseline event rate in group 1.
base_rate <- rgamma(total_studies, alpha1, alpha1 / rate1)

# Events in group 1.
y2 <- rpois(total_studies, n2 * base_rate)

# Events in group 2.
rr <- rbeta(total_studies, alpha2, beta2)
# problem: rr/1-rr is not gamma - its beta prime.
y1 <- rpois(total_studies, n1 * base_rate * rr / (1 - rr))


mean(base_rate * rr / (1 - rr))
mean(base_rate)
mean(rr)


# PROPOSED NEW VERSION OF THE CODE

total_studies <- 500
rate1 <- 0.015
alpha1 <- 1.44

alpha2 <-3
beta2 <- 3

n1 <- rep(500, total_studies)
n2 <- rep(1000, total_studies)
# Baseline event rate in group 1.
base_rate_2 <- rgamma(total_studies, alpha2, alpha2 / rate1)
y2 <- rpois(total_studies, n2 * base_rate_2)

# Events in group 2.
base_rate_1 <- rgamma(total_studies, beta2, alpha2 / rate1)
y1 <- rpois(total_studies, n1 * base_rate_1)


mean(base_rate_1)
mean(base_rate_2)
mean(base_rate_1/(base_rate_1 + base_rate_2))

