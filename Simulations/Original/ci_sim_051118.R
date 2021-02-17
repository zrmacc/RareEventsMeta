# # For SLURM
# args = commandArgs(trailingOnly = TRUE)
# K_org = as.numeric(args[[1]]);
# alpha = as.numeric(args[[2]]);
# beta = as.numeric(args[[3]]);
# base_rate = as.numeric(args[[4]]);
# jj = as.numeric(args[[5]]);
#
#
# AI <- Sys.getenv('SLURM_ARRAY_TASK_ID')
# kk <- as.numeric(AI);


# To run on desktop:
# setwd('~/Desktop/ci_code')
K_org <- 12
alpha <- 6.5
beta <- 6.5
base_rate <- 0.006
jj <- 3
kk <- 1
i <- 1


# read in functions
# source('meta_fun_052818.R')
library(parallel)
devtools::install(pkg = "exactinf4meta", reload = TRUE)
library("exactinf4meta")


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



CIs.all <- NULL
ests.all <- NULL
K.all <- NULL
inds.all <- NULL
tot <- NULL
CI.all <- NULL
time.all <- NULL
bb <- 500
numb <- 5
for (i in ((kk - 1) * numb + 1):((kk - 1) * numb + numb)) {
  set.seed(i)
  # generate the data
  data.obj <- gen.data(K_org, N1, N2, alpha, beta, base_rate = base_rate)

  # data with DZ studies removed
  data <- data.obj$data
  K <- data.obj$K
  summary(data)


  mean(data[, "Grp1_Ev"] + data[, "Grp2_Ev"])
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

  mom.o <- ests
  # Get the normal-based CI
  CI.dsl <- c(mom.o$muhat - 1.96 * sqrt(mom.o$v_muhat), mom.o$muhat + 1.96 * sqrt(mom.o$v_muhat))
  CI.dsl
  # Calculate the number of cores
  no_cores <- 2

  # Initiate cluster
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, library("exactinf4meta"))

  ptm <- proc.time()
  # run the perturbation
  bs.res <- parSapply(
    cl, 1:2,
    function(kk) {
      if (kk == 1) {
        res <- outer_bound(bb, mom.o, data[, "Grp1_Ev"], data[, "Grp1_Ev"] + data[, "Grp2_Ev"],
          data[, "Grp1_Sz"], data[, "Grp2_Sz"], CI.dsl,
          mesh.size = 0.01
        )
      } else {
        res <- lower_bound(bb, mom.o, data[, "Grp1_Ev"], data[, "Grp1_Ev"] + data[, "Grp2_Ev"],
          data[, "Grp1_Sz"], data[, "Grp2_Sz"], CI.dsl,
          mesh.size = 0.01
        )
      }
    }
  )
  ptm2 <- proc.time() - ptm
  stopCluster(cl)


  res_CI <- do.call(rbind, bs.res)
  CI_norm_res <- res_CI[which(res_CI[, "pval.norm"] >= 0.05), ]
  CI_norm <- c(min(CI_norm_res[, "mu"]), max(CI_norm_res[, "mu"]))
  CI_norm




  CIs.all <- rbind(CIs.all, CI_norm)
  ests.all <- rbind(ests.all, unlist(mom.o))
  K.all <- c(K.all, K)
  tot <- rbind(tot, c(max(res_CI[, c("epsN_upper")]), max(res_CI[, c("epsN_lower")])))
  time.all <- c(time.all, ptm2[3])
}



# colMeans(ests.all)
# alpha/(alpha+beta)
# alpha/(alpha+beta)*(1-alpha/(alpha+beta))*(1/(alpha+beta+1))
# alpha/(alpha+beta)*(1-alpha/(alpha+beta))*(1/(alpha+beta+1)) + (alpha/(alpha+beta))^2
# colSds(ests.all)
# mean(inds.all)

write(t(CIs.all),
  file =
    paste0(
      kk, "balan-", jj, "CI-", K_org, "-alpha-", alpha, "-beta-", beta,
      "-baseR-", base_rate, ".txt"
    ),
  append = T, ncol = ncol(CIs.all)
)



write(t(tot),
  file =
    paste0(
      kk, "balan-", jj, "TOT-", K_org, "-alpha-", alpha, "-beta-", beta,
      "-baseR-", base_rate, ".txt"
    ),
  append = T, ncol = ncol(tot)
)



write((time.all),
  file =
    paste0(
      kk, "balan-", jj, "TIME-", K_org, "-alpha-", alpha, "-beta-", beta,
      "-baseR-", base_rate, ".txt"
    ),
  append = T, ncol = 1
)




#
#
# res2 = lower_bound(bb, mom.o, data[, 'Grp1_Ev'], data[, 'Grp1_Ev'] + data[, 'Grp2_Ev'],
#                   data[, 'Grp1_Sz'], data[, 'Grp2_Sz'], CI.dsl, mesh.size = 0.01)
#
#
