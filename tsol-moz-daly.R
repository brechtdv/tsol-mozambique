### TSOL MOZAMBIQUE
### -- disability adjusted life years

## attach packages
library(mc2d)
library(bd)

## import helper functions
source("helpers.R")


## -------------------------------------------------------------------------
## SETTINGS
##

## set seed to allow for reproducibility
set.seed(123)

## number of iterations
n <- 1e5


## -------------------------------------------------------------------------
## PARAMETER VALUES 
##

## POPULATION

## define population matrix - Agonia district
pop_mx <-
  matrix(c(30777, 45925, 59408, 15156, 5065,
           31990, 45317, 71449, 18588, 6653),
         ncol = 2)

## EPILEPSY

## active epilepsy prevalence, Beta distribution
ae_prev_alpha <-   22
ae_prev_beta  <- 1701

## proportion treated, Beta distribution
ep_trt_alpha <-  22
ep_trt_beta  <- 129

## epilepsy duration, Fixed
ep_dur <-
  matrix(c(1.4, 2.0, 3.6, 2.8, 1.6,
           1.6, 3.1, 5.9, 6.0, 2.8),
         ncol = 2)

## epilepsy onset, Fixed
ep_ons <- c(2.50, 9.95, 26.99, 51.94, 73.60)

## epilepsy treated disability weight, Uniform distribution 
dsw_ep_tr_min <- 0.047
dsw_ep_tr_max <- 0.106

## epilepsy non-treated disability weight, Uniform distribution 
dsw_ep_nt_min <- 0.279
dsw_ep_nt_max <- 0.572

## epilepsy mortality rate, Gamma distribution
ep_mrt_n   <-     1800  # nationwide epilepsy deaths
ep_mrt_pop <- 22383000  # nationwide pop

## age at death, Fixed
ep_aad <- c(2.50, 10.0, 30.0, 53.50, 77.5)

## NEUROCYSTICERCOSIS

## proportion PWE with NCC, Uniform parameters 
prop_ncc_min <- 0.427
prop_ncc_max <- 0.592

## HEADACHE

## proportion migraine headache, tension-type headache, no headache
## .. among people with active epilepsy prev, Dirichlet dist
ha_prev <- c(37, 57, 57)

## headache durations, fixed
ha_dur <-
  matrix(c(4.28, 4.28, 4.28, 2.90, 4.75,
           4.28, 4.28, 4.28, 2.90, 4.75),
         ncol = 2)

## headache onset, Fixed
ha_ons <- c(2.50, 9.95, 26.99, 51.94, 73.60)

## GBD 2013 disability weights, Uniform distribution
dsw_mgr_min <- 0.294
dsw_mgr_max <- 0.588

tt_dsw_min <- 0.022
tt_dsw_max <- 0.057


## -------------------------------------------------------------------------
## SIMULATIONS
##

## active epilepsy prevalence (population wise)
prev_ep_act <- rbeta(n, ae_prev_alpha, ae_prev_beta)
mean_ci(prev_ep_act)

## proportion PWE with NCC, P(NCC|EP)
prop_ncc <- runif(n, prop_ncc_min, prop_ncc_max)
mean_ci(prop_ncc)

## proportion treated, P(TR|EP)=P(TR|EP,NCC)
prop_ep_trt <- rbeta(n, ep_trt_alpha, ep_trt_beta)
mean_ci(prop_ep_trt)

## epilepsy disability weights
dsw_ep_tr <- runif(n, dsw_ep_tr_min, dsw_ep_tr_max)
dsw_ep_nt <- runif(n, dsw_ep_nt_min, dsw_ep_nt_max)
mean_ci(dsw_ep_tr); mean_ci(dsw_ep_nt)

## epilepsy mortality rate
mrt_ep <- rgamma(n, ep_mrt_n, ep_mrt_pop)
mean_ci(mrt_ep)

## proportions headache among people with epilepsy, P(HA|EP)
## .. 1/migraine, 2/tension, 3/no
prop_ha <- rdirichlet(n, ha_prev)
prop_mgr <- prop_ha[, 1]
prop_tns <- prop_ha[, 2]
mean_ci(prop_mgr); mean_ci(prop_tns)

## migraine disability weights
dsw_mgr <- runif(n, dsw_mgr_min, dsw_mgr_max)
mean_ci(dsw_mgr)

## tension-type headache disability weights
dsw_tns <- runif(n, tt_dsw_min, tt_dsw_max)
mean_ci(dsw_tns)


## -------------------------------------------------------------------------
## CALCULATIONS
##

## prevalence NCC-associated active epilepsy
prev_ncc_ep_act <- prev_ep_act * prop_ncc
mean_ci(prev_ncc_ep_act)

## prevalence of NCC-migraine in NCC-epilepsy patients
prev_ncc_ep_mgr <- prev_ncc_ep_act * prop_mgr
mean_ci(prev_ncc_ep_mgr)

## prevalence of NCC-tension type headache in NCC-epilepsy patients
prev_ncc_ep_tns <- prev_ncc_ep_act * prop_tns
mean_ci(prev_ncc_ep_tns)


## -------------------------------------------------------------------------
## EPILEPSY CASES
##

## NCC-epilepsy incidence rate, age-sex specific [5*2 matrix]
## .. prevalence / duration
inc_ncc_ep_act <- apply(ep_dur, 1:2, function(x) prev_ncc_ep_act / x)
str(inc_ncc_ep_act)

## NCC-epilepsy incident cases, age-sex specific [10 rows]
N_ncc_ep_act <- apply(inc_ncc_ep_act, 1, function(x) x * pop_mx)
t(apply(N_ncc_ep_act, 1, mean_ci))

## NCC-epilepsy incident cases - treated, age-sex specific [10 cols]
N_ncc_ep_act_tr <- t(apply(N_ncc_ep_act, 1, function(x) x * prop_ep_trt))
t(apply(N_ncc_ep_act_tr, 1, mean_ci))

## NCC-epilepsy incident cases - non-treated, age-sex specific [10 cols]
N_ncc_ep_act_nt <- t(apply(N_ncc_ep_act, 1, function(x) x * (1-prop_ep_trt)))
t(apply(N_ncc_ep_act_nt, 1, mean_ci))


## -------------------------------------------------------------------------
## EPILEPSY DEATHS
##

## Epilepsy deaths, age-sex specific [10 rows]
M_ep <- t(apply(t(c(pop_mx)), 2, function(x) x * mrt_ep))
apply(M_ep, 1, mean_ci)

## NCC-epilepsy deaths, age-sex specific [10 rows]
M_ncc_ep <- t(apply(M_ep, 1, function(x) x * prop_ncc))
apply(M_ncc_ep, 1, mean_ci)


## -------------------------------------------------------------------------
## HEADACHE CASES
##

## NCC-epilepsy-migraine incidence rate, age-sex specific [5*2 matrix]
inc_ncc_ep_mgr <- apply(ha_dur, 1:2, function(x) prev_ncc_ep_mgr / x)
str(inc_ncc_ep_mgr)

## NCC-epilepsy-ttha incidence rate, age-sex specific [5*2 matrix]
inc_ncc_ep_tns <- apply(ha_dur, 1:2, function(x) prev_ncc_ep_tns / x)
str(inc_ncc_ep_tns)

## NCC-epilepsy-migraine incident cases, age-sex specific [10 rows]
N_ncc_ep_mgr <- apply(inc_ncc_ep_mgr, 1, function(x) x * pop_mx)
t(apply(N_ncc_ep_mgr, 1, mean_ci))

## NCC-epilepsy-ttha incident cases, age-sex specific [10 rows]
N_ncc_ep_tns <- apply(inc_ncc_ep_tns, 1, function(x) x * pop_mx)
t(apply(N_ncc_ep_tns, 1, mean_ci))


## -------------------------------------------------------------------------
## DISABILITY ADJUSTED LIFE YEARS
##

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## SOCIAL WEIGHTING - select one of these
K <- 0; r <- 0     # DALY[0;0]
#K <- 0; r <- 0.03  # DALY[0;0.03]
#K <- 1; r <- 0     # DALY[1;0]
#K <- 1; r <- 0.03  # DALY[1;0.03]

## RESIDUAL LIFE EXPECTANCY - select 'who' or 'gbd'
ep_rle <- rle(ep_aad, "who")

##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

## YEARS LIVED WITH DISABILITY -- epilepsy treated cases
yld_ep_tr <- array(dim = c(5, 2, n))

yld_ep_tr[1, 1, ] <-
  burden(N = N_ncc_ep_act_tr[1, ], DW = dsw_ep_tr,
         A = ep_ons[1], L = ep_dur[1, 1], K = K, r = r, a = ep_ons[1])
yld_ep_tr[2, 1, ] <-
  burden(N = N_ncc_ep_act_tr[2, ], DW = dsw_ep_tr,
         A = ep_ons[2], L = ep_dur[2, 1], K = K, r = r, a = ep_ons[2])
yld_ep_tr[3, 1, ] <-
  burden(N = N_ncc_ep_act_tr[3, ], DW = dsw_ep_tr,
         A = ep_ons[3], L = ep_dur[3, 1], K = K, r = r, a = ep_ons[3])
yld_ep_tr[4, 1, ] <-
  burden(N = N_ncc_ep_act_tr[4, ], DW = dsw_ep_tr,
         A = ep_ons[4], L = ep_dur[4, 1], K = K, r = r, a = ep_ons[4])
yld_ep_tr[5, 1, ] <-
  burden(N = N_ncc_ep_act_tr[5, ], DW = dsw_ep_tr,
         A = ep_ons[5], L = ep_dur[5, 1], K = K, r = r, a = ep_ons[5])

yld_ep_tr[1, 2, ] <-
  burden(N = N_ncc_ep_act_tr[6, ], DW = dsw_ep_tr,
         A = ep_ons[1], L = ep_dur[1, 2], K = K, r = r, a = ep_ons[1])
yld_ep_tr[2, 2, ] <-
  burden(N = N_ncc_ep_act_tr[7, ], DW = dsw_ep_tr,
         A = ep_ons[2], L = ep_dur[2, 2], K = K, r = r, a = ep_ons[2])
yld_ep_tr[3, 2, ] <-
  burden(N = N_ncc_ep_act_tr[8, ], DW = dsw_ep_tr,
         A = ep_ons[3], L = ep_dur[3, 2], K = K, r = r, a = ep_ons[3])
yld_ep_tr[4, 2, ] <-
  burden(N = N_ncc_ep_act_tr[9, ], DW = dsw_ep_tr,
         A = ep_ons[4], L = ep_dur[4, 2], K = K, r = r, a = ep_ons[4])
yld_ep_tr[5, 2, ] <-
  burden(N = N_ncc_ep_act_tr[10, ], DW = dsw_ep_tr,
         A = ep_ons[5], L = ep_dur[5, 2], K = K, r = r, a = ep_ons[5])

yld_ep_tr_all <- apply(yld_ep_tr, 3, sum)
mean_ci(yld_ep_tr_all)


## YEARS LIVED WITH DISABILITY -- epilepsy not-treated cases
yld_ep_nt <- array(dim = c(5, 2, n))

yld_ep_nt[1, 1, ] <-
  burden(N = N_ncc_ep_act_nt[1, ], DW = dsw_ep_nt,
         A = ep_ons[1], L = ep_dur[1, 1], K = K, r = r, a = ep_ons[1])
yld_ep_nt[2, 1, ] <-
  burden(N = N_ncc_ep_act_nt[2, ], DW = dsw_ep_nt,
         A = ep_ons[2], L = ep_dur[2, 1], K = K, r = r, a = ep_ons[2])
yld_ep_nt[3, 1, ] <-
  burden(N = N_ncc_ep_act_nt[3, ], DW = dsw_ep_nt,
         A = ep_ons[3], L = ep_dur[3, 1], K = K, r = r, a = ep_ons[3])
yld_ep_nt[4, 1, ] <-
  burden(N = N_ncc_ep_act_nt[4, ], DW = dsw_ep_nt,
         A = ep_ons[4], L = ep_dur[4, 1], K = K, r = r, a = ep_ons[4])
yld_ep_nt[5, 1, ] <-
  burden(N = N_ncc_ep_act_nt[5, ], DW = dsw_ep_nt,
         A = ep_ons[5], L = ep_dur[5, 1], K = K, r = r, a = ep_ons[5])

yld_ep_nt[1, 2, ] <-
  burden(N = N_ncc_ep_act_nt[6, ], DW = dsw_ep_nt,
         A = ep_ons[1], L = ep_dur[1, 2], K = K, r = r, a = ep_ons[1])
yld_ep_nt[2, 2, ] <-
  burden(N = N_ncc_ep_act_nt[7, ], DW = dsw_ep_nt,
         A = ep_ons[2], L = ep_dur[2, 2], K = K, r = r, a = ep_ons[2])
yld_ep_nt[3, 2, ] <-
  burden(N = N_ncc_ep_act_nt[8, ], DW = dsw_ep_nt,
         A = ep_ons[3], L = ep_dur[3, 2], K = K, r = r, a = ep_ons[3])
yld_ep_nt[4, 2, ] <-
  burden(N = N_ncc_ep_act_nt[9, ], DW = dsw_ep_nt,
         A = ep_ons[4], L = ep_dur[4, 2], K = K, r = r, a = ep_ons[4])
yld_ep_nt[5, 2, ] <-
  burden(N = N_ncc_ep_act_nt[10, ], DW = dsw_ep_nt,
         A = ep_ons[5], L = ep_dur[5, 2], K = K, r = r, a = ep_ons[5])

yld_ep_nt_all <- apply(yld_ep_nt, 3, sum)
mean_ci(yld_ep_nt_all)


## YEARS LIVED WITH DISABILITY -- epilepsy
yld_ep <- yld_ep_tr + yld_ep_nt
yld_ep_all <- apply(yld_ep, 3, sum)
mean_ci(yld_ep_all)


## YEARS OF LIFE LOST -- epilepsy
yll_ep <- array(dim = c(5, 2, n))

yll_ep[1, 1, ] <-
  burden(N = M_ncc_ep[1, ], DW = 1,
         A = ep_aad[1], L = ep_rle[1], K = K, r = r, a = ep_aad[1])
yll_ep[2, 1, ] <-
  burden(N = M_ncc_ep[2, ], DW = 1,
         A = ep_aad[2], L = ep_rle[2], K = K, r = r, a = ep_aad[2])
yll_ep[3, 1, ] <-
  burden(N = M_ncc_ep[3, ], DW = 1,
         A = ep_aad[3], L = ep_rle[3], K = K, r = r, a = ep_aad[3])
yll_ep[4, 1, ] <-
  burden(N = M_ncc_ep[4, ], DW = 1,
         A = ep_aad[4], L = ep_rle[4], K = K, r = r, a = ep_aad[4])
yll_ep[5, 1, ] <-
  burden(N = M_ncc_ep[5, ], DW = 1,
         A = ep_aad[5], L = ep_rle[5], K = K, r = r, a = ep_aad[5])

yll_ep[1, 2, ] <-
  burden(N = M_ncc_ep[6, ], DW = 1,
         A = ep_aad[1], L = ep_rle[1], K = K, r = r, a = ep_aad[1])
yll_ep[2, 2, ] <-
  burden(N = M_ncc_ep[7, ], DW = 1,
         A = ep_aad[2], L = ep_rle[2], K = K, r = r, a = ep_aad[2])
yll_ep[3, 2, ] <-
  burden(N = M_ncc_ep[8, ], DW = 1,
         A = ep_aad[3], L = ep_rle[3], K = K, r = r, a = ep_aad[3])
yll_ep[4, 2, ] <-
  burden(N = M_ncc_ep[9, ], DW = 1,
         A = ep_aad[4], L = ep_rle[4], K = K, r = r, a = ep_aad[4])
yll_ep[5, 2, ] <-
  burden(N = M_ncc_ep[10, ], DW = 1,
         A = ep_aad[5], L = ep_rle[5], K = K, r = r, a = ep_aad[5])

yll_ep_all <- apply(yll_ep, 3, sum)
mean_ci(yll_ep_all)


## DISABILITY-ADJUSTED LIFE YEARS -- epilepsy
daly_ep <- yld_ep + yll_ep
daly_ep_all <- apply(daly_ep, 3, sum)
mean_ci(daly_ep_all)


## YEARS LIVED WITH DISABILITY -- migraine
yld_ha_mgr <- array(dim = c(5, 2, n))

yld_ha_mgr[1, 1, ] <-
  burden(N = N_ncc_ep_mgr[1, ], DW = dsw_mgr,
         A = ha_ons[1], L = ha_dur[1, 1], K = K, r = r, a = ha_ons[1])
yld_ha_mgr[2, 1, ] <-
  burden(N = N_ncc_ep_mgr[2, ], DW = dsw_mgr,
         A = ha_ons[2], L = ha_dur[2, 1], K = K, r = r, a = ha_ons[2])
yld_ha_mgr[3, 1, ] <-
  burden(N = N_ncc_ep_mgr[3, ], DW = dsw_mgr,
         A = ha_ons[3], L = ha_dur[3, 1], K = K, r = r, a = ha_ons[3])
yld_ha_mgr[4, 1, ] <-
  burden(N = N_ncc_ep_mgr[4, ], DW = dsw_mgr,
         A = ha_ons[4], L = ha_dur[4, 1], K = K, r = r, a = ha_ons[4])
yld_ha_mgr[5, 1, ] <-
  burden(N = N_ncc_ep_mgr[5, ], DW = dsw_mgr,
         A = ha_ons[5], L = ha_dur[5, 1], K = K, r = r, a = ha_ons[5])

yld_ha_mgr[1, 2, ] <-
  burden(N = N_ncc_ep_mgr[6, ], DW = dsw_mgr,
         A = ha_ons[1], L = ha_dur[1, 2], K = K, r = r, a = ha_ons[1])
yld_ha_mgr[2, 2, ] <-
  burden(N = N_ncc_ep_mgr[7, ], DW = dsw_mgr,
         A = ha_ons[2], L = ha_dur[2, 2], K = K, r = r, a = ha_ons[2])
yld_ha_mgr[3, 2, ] <-
  burden(N = N_ncc_ep_mgr[8, ], DW = dsw_mgr,
         A = ha_ons[3], L = ha_dur[3, 2], K = K, r = r, a = ha_ons[3])
yld_ha_mgr[4, 2, ] <-
  burden(N = N_ncc_ep_mgr[9, ], DW = dsw_mgr,
         A = ha_ons[4], L = ha_dur[4, 2], K = K, r = r, a = ha_ons[4])
yld_ha_mgr[5, 2, ] <-
  burden(N = N_ncc_ep_mgr[10, ], DW = dsw_mgr,
         A = ha_ons[5], L = ha_dur[5, 2], K = K, r = r, a = ha_ons[5])

yld_ha_mgr_all <- apply(yld_ha_mgr, 3, sum)
mean_ci(yld_ha_mgr_all)


## YEARS LIVED WITH DISABILITY -- tension-type headached
yld_ha_tns <- array(dim = c(5, 2, n))

yld_ha_tns[1, 1, ] <-
  burden(N = N_ncc_ep_tns[1, ], DW = dsw_tns,
         A = ha_ons[1], L = ha_dur[1, 1], K = K, r = r, a = ha_ons[1])
yld_ha_tns[2, 1, ] <-
  burden(N = N_ncc_ep_tns[2, ], DW = dsw_tns,
         A = ha_ons[2], L = ha_dur[2, 1], K = K, r = r, a = ha_ons[2])
yld_ha_tns[3, 1, ] <-
  burden(N = N_ncc_ep_tns[3, ], DW = dsw_tns,
         A = ha_ons[3], L = ha_dur[3, 1], K = K, r = r, a = ha_ons[3])
yld_ha_tns[4, 1, ] <-
  burden(N = N_ncc_ep_tns[4, ], DW = dsw_tns,
         A = ha_ons[4], L = ha_dur[4, 1], K = K, r = r, a = ha_ons[4])
yld_ha_tns[5, 1, ] <-
  burden(N = N_ncc_ep_tns[5, ], DW = dsw_tns,
         A = ha_ons[5], L = ha_dur[5, 1], K = K, r = r, a = ha_ons[5])

yld_ha_tns[1, 2, ] <-
  burden(N = N_ncc_ep_tns[6, ], DW = dsw_tns,
         A = ha_ons[1], L = ha_dur[1, 2], K = K, r = r, a = ha_ons[1])
yld_ha_tns[2, 2, ] <-
  burden(N = N_ncc_ep_tns[7, ], DW = dsw_tns,
         A = ha_ons[2], L = ha_dur[2, 2], K = K, r = r, a = ha_ons[2])
yld_ha_tns[3, 2, ] <-
  burden(N = N_ncc_ep_tns[8, ], DW = dsw_tns,
         A = ha_ons[3], L = ha_dur[3, 2], K = K, r = r, a = ha_ons[3])
yld_ha_tns[4, 2, ] <-
  burden(N = N_ncc_ep_tns[9, ], DW = dsw_tns,
         A = ha_ons[4], L = ha_dur[4, 2], K = K, r = r, a = ha_ons[4])
yld_ha_tns[5, 2, ] <-
  burden(N = N_ncc_ep_tns[10, ], DW = dsw_tns,
         A = ha_ons[5], L = ha_dur[5, 2], K = K, r = r, a = ha_ons[5])

yld_ha_tns_all <- apply(yld_ha_tns, 3, sum)
mean_ci(yld_ha_tns_all)


## YEARS LIVED WITH DISABILITY -- headache
yld_ha <- yld_ha_mgr + yld_ha_tns
yld_ha_all <- apply(yld_ha, 3, sum)
mean_ci(yld_ha_all)


## YEARS LIVED WITH DISABILITY -- total
yld <- yld_ep + yld_ha
yld_all <- apply(yld, 3, sum)
mean_ci(yld_all)


## DISABILITY-ADJUSTED LIFE YEARS -- total
daly <- daly_ep + yld_ha
daly_all <- apply(daly, 3, sum)
mean_ci(daly_all)


###
### SUMMARIES
###

## cases, deaths
mean_ci(colSums(N_ncc_ep_act))
mean_ci(colSums(N_ncc_ep_mgr))
mean_ci(colSums(N_ncc_ep_tns))
mean_ci(colSums(M_ncc_ep))

## YLDs, YLLs, DALYs
mean_ci(yld_ep_all)
mean_ci(yld_ha_all)
mean_ci(yld_all)
mean_ci(yll_ep_all)
mean_ci(daly_all)
mean_ci(daly_ep_all)

## YLDs, YLLs, DALYs per 1000 population
mean_ci(1e3 * yld_ep_all / sum(pop_mx))
mean_ci(1e3 * yld_ha_all / sum(pop_mx))
mean_ci(1e3 * yld_all / sum(pop_mx))
mean_ci(1e3 * yll_ep_all / sum(pop_mx))
mean_ci(1e3 * daly_all / sum(pop_mx))
mean_ci(1e3 * daly_ep_all / sum(pop_mx))

## YLDs, YLLs, DALYs per incident epilepsy case
mean_ci(yld_ep_all / colSums(N_ncc_ep_act))
mean_ci(yld_ha_all / colSums(N_ncc_ep_act))
mean_ci(yld_all / colSums(N_ncc_ep_act))
mean_ci(yll_ep_all / colSums(N_ncc_ep_act))
mean_ci(daly_all / colSums(N_ncc_ep_act))
mean_ci(daly_ep_all / colSums(N_ncc_ep_act))

## YLD, YLL contribution
mean_ci(yld_ep_all / daly_ep_all)
mean_ci(yll_ep_all / daly_ep_all)

mean_ci(yld_ha_mgr_all / yld_ha_all)
mean_ci(yld_ha_tns_all / yld_ha_all)

mean_ci(yld_ep_all / daly_all)
mean_ci(yld_ha_all / daly_all)
mean_ci(yld_ha_mgr_all / daly_all)
mean_ci(yld_ha_tns_all / daly_all)
mean_ci(yld_all / daly_all)
mean_ci(yll_ep_all / daly_all)
