# For SLURM
args <- commandArgs(trailingOnly = TRUE)
K_org <- as.numeric(args[[1]])
alpha <- as.numeric(args[[2]])
beta <- as.numeric(args[[3]])
base_rate <- as.numeric(args[[4]])
jj <- as.numeric(args[[5]])
# AI <- Sys.getenv('SLURM_ARRAY_TASK_ID')
# kk <- as.numeric(AI);


# To run on desktop:
# setwd('~/Desktop/cp_code')
K_org <- 24
alpha <- 6.5
beta <- 6.5
base_rate <- 0.006
jj <- 2
kk <- 1
i <- 1


# if(base_rate == 0.01){base_rate = 0.05}

# import package
# source('meta_fun_052818.R')

## if(K_org == 8){K_org = 12}else if(K_org == 12){K_org = 16}else if(K_org == 16){K_org = 24}else if(K_org == 24){K_org = 48}
set.seed(92047)

## study specific sizes
if (jj == 2) {
  N2 <- c(
    357, 391, 774, 213, 232, 43, 121, 110, 382, 284, 294, 563, 278, 418,
    395, 203, 104, 212, 138, 196, 122, 175, 56, 39, 561, 116, 148, 231, 89,
    168, 116, 1172, 706, 204, 288, 254, 314, 162, 442, 394, 2635, 1456, 101, 232,
    70, 25, 196, 676
  )
  N1 <- N2
  N1 <- N1[1:K_org]
  N2 <- N2[1:K_org]
} else {
  N2 <- c(
    176, 207, 185, 109, 116, 47, 124, 114, 384, 135, 302, 142, 279, 212, 198,
    106, 99, 107, 139, 96, 120, 173, 58, 38, 276, 111, 143, 242, 88, 172,
    61, 377, 325, 185, 280, 272, 154, 160, 112, 124, 2634, 2895, 51, 115, 75,
    24, 195, 225
  )
  N1 <- c(
    357, 391, 774, 213, 232, 43, 121, 110, 382, 284, 294, 563, 278, 418,
    395, 203, 104, 212, 138, 196, 122, 175, 56, 39, 561, 116, 148, 231, 89,
    168, 116, 1172, 706, 204, 288, 254, 314, 162, 442, 394, 2635, 1456, 101, 232,
    70, 25, 196, 676
  )
  N1 <- N1[1:K_org]
  N2 <- N2[1:K_org]
}


res.track <- NULL # store order of the p values
ests.all <- NULL # store mom ests
K.all <- NULL # store overall K due to double zero d
inds.all <- NULL # store dropped studies due to data generation
res.all <- NULL # results of reject/accept
bb <- 2000 # number of MC replications
numb <- 50 # number of runs
max.p <- NULL # max pvalue

for (i in ((kk - 1) * numb + 1):((kk - 1) * numb + numb)) {
  set.seed(i)

  # generate the data
  data.obj <- gen.data(K_org, N1, N2, alpha, beta, base_rate = base_rate)

  # data with DZ studies removed
  data <- data.obj$data
  K <- data.obj$K
  summary(data)

  # skip an iteration if number of events too large in either group
  inds <- c(which(data[, 5] > data[, 4]), which(data[, 3] > data[, 2]))
  inds
  inds.all <- c(inds.all, I(length(inds) > 0))


  if (length(inds) >= 1) {
    next
  }



  if (!is.null(nrow(data))) {
    ptm <- proc.time()
    ests <- MOM_est(data, K)
    proc.time() - ptm
  } else {
    ests <- rep(NA, 6)
  }

  if (sum(N1 == N2) == length(N1)) {
    mu.type <- "one2one"
  } else {
    mu.type <- "unbal"
  }

  mom.o <- ests
  m <- alpha / (alpha + beta)

  # get alpha and beta along the boundary
  v1 <- m^2 * (1 - m) / (1 + m)
  v2 <- m * (1 - m)^2 / (2 - m)
  v.t0 <- min(v1, v2)
  my.var <- seq(1e-6, v.t0, length.out = 15)
  # my.var = v.t0
  l <- length(my.var)
  res.t.all <- NULL
  for (i in 1:l) {
    v.t <- my.var[i]
    a.t <- m * ((m * (1 - m) - v.t) / v.t) * (1 + 1e-6)
    b.t <- (1 - m) * ((m * (1 - m) - v.t) / v.t) * (1 + 1e-6)
    mu.t <- a.t / (a.t + b.t)

    # run the Monte Carlo
    res.t <- run_MC(bb, mom.o, data[, "Grp1_Ev"], data[, "Grp1_Ev"] + data[, "Grp2_Ev"],
      data[, "Grp1_Sz"], data[, "Grp2_Sz"], a.t, b.t,
      mu.type = mu.type
    )
    res.t.all <- rbind(res.t.all, res.t)
  }
  res.0 <- I(colMeans(I(cbind(res.t.all[, 5], res.t.all[, 6]) >= 0.05)) > 0)
  res.all <- rbind(res.all, res.0)

  res.track <- rbind(res.track, c(
    max(which(res.t.all[, 5] == max(res.t.all[, 5]))),
    max(which(res.t.all[, 6] == max(res.t.all[, 6])))
  ))
  ests.all <- rbind(ests.all, unlist(mom.o))
  K.all <- c(K.all, K)
  res.t.all
  max.p <- rbind(max.p, c(max(res.t.all[, 5]), max(res.t.all[, 6])))
}

# library(matrixStats)
# mean(K.all)
# colMeans(ests.all)
# colVars(ests.all)
# alpha/(alpha+beta)
# alpha/(alpha+beta)*(1-alpha/(alpha+beta))/(alpha+beta+1)
# colMeans(res.all)
# mean(inds.all)


# setwd('~/Desktop')
write(t(res.all),
  file =
    paste0(
      kk, "balan-", jj, "CP-", K_org, "-alpha-", alpha, "-beta-", beta,
      "-baseR-", base_rate, ".txt"
    ),
  append = T, ncol = ncol(res.all)
)



write(t(res.track),
  file =
    paste0(
      kk, "balan-", jj, "resTRACK-", K_org, "-alpha-", alpha, "-beta-", beta,
      "-baseR-", base_rate, ".txt"
    ),
  append = T, ncol = ncol(res.track)
)


write(K.all,
  file =
    paste0(
      kk, "balan-", jj, "kALL-", K_org, "-alpha-", alpha, "-beta-", beta,
      "-baseR-", base_rate, ".txt"
    ),
  append = T, ncol = 1
)


write(t(ests.all),
  file =
    paste0(
      kk, "balan-", jj, "Ests-", K_org, "-alpha-", alpha, "-beta-", beta,
      "-baseR-", base_rate, ".txt"
    ),
  append = T, ncol = ncol(ests.all)
)

write(t(max.p),
  file =
    paste0(
      kk, "balan-", jj, "MAXp-", K_org, "-alpha-", alpha, "-beta-", beta,
      "-baseR-", base_rate, ".txt"
    ),
  append = T, ncol = 2
)
