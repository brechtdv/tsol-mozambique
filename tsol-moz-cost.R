### TSOL MOZAMBIQUE
### -- economic impact

## attach packages
library(mc2d)
library(bd)

## import helper functions
source("helpers.R")

## run DALY calculation script
source("tsol-moz-daly.R")


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

## POPULATION SIZE

## human population - Angónia district
pop <- sum(pop_mx)

## human NCC cases - Angónia district
n_ncc <- pop * prev_ncc_ep_act
mean_ci(n_ncc)

## pig population - Angónia district
n_pigs_smallscale <- 20411


## HUMANS

## hospitalization probability, Beta distribution
p_hosp_alpha <- 38
p_hosp_beta <- 222

## hospitalization duration, Uniform distribution
stay_min <- 1
stay_max <- 21

## care seeking probabilities, Dirichlet distribution
p_care <- c(0.091, 0.227, 0.682)

## number of visits to healthcare provider, Uniform distribution
n_visit_med_min <- 1
n_visit_med_max <- 8

## number of visits to a traditional healer, Uniform distribution
n_year_heal_min <- 1
n_year_heal_max <- 10

## Phenobarbital treatment probability, Beta distribution
p_pheno_alpha <- 36
p_pheno_beta <- 46  

## loss of working time, Uniform distribution
loss_workingtime_min <- 1
loss_workingtime_max <- 30

## unemployment due to epilepsy, Beta distribution
unemployed_duetoepilepsy_alpha <-   6
unemployed_duetoepilepsy_beta  <- 261
  
## number of working days per year, Uniform distribution
working_days_min <- 220
working_days_max <- 312

## proportion of active population, Fixed
active <- 0.3969

## monthly salary, Gamma distribution
monthly_salary_shape <- 5.3
monthly_salary_rate  <- 0.06

## cost of visit to a medical doctor, Gamma distribution
price_med_shape <- 33.4 
price_med_rate  <- 11.1

## cost of hospitalization per day, Fixed
price_day_hosp <- 2.3

## cost of visiting traditional healer, Gamma distribution
price_heal_shape <- 8.5
price_heal_rate  <- 0.15

## price of Phenobarbital, Fixed
price_pheno <- 5

## number of working days per month, Fixed
n_working_day_by_month <- 26


## PIGS

## value of a pig, Gamma distribution
price_pigs_shape <- 16.5
price_pigs_rate  <-  0.32

## value reduction, Fixed
price_loss_pigs <- 0.5

## proportion of pigs sold per year, Fixed
pigs_sold <- 0.333
 
## porcine cysticercosis prevalence, Uniform distribution
prev_pigs_alpha <-  84
prev_pigs_beta  <- 577


## -------------------------------------------------------------------------
## SIMULATIONS / HUMANS
##

## probability hospitalization
p_hosp <- rbeta(n, p_hosp_alpha, p_hosp_beta)

## number of hospitalized NCC patients
n_hosp <- n_ncc * p_hosp

## duration of hospital stay
## .. per iteration, generate 'n_hosp' random durations, and sum them up
stay <- apply(t(n_hosp), 2, function(x) sum(runif(x, stay_min, stay_max)))

## patients not in hospital, seeking medical treatment
n_ncc2 <- n_ncc - n_hosp

## probility healthcare seeking
xyz <- rmultinom(n, n_ncc2, p_care)
n_heal <- xyz[1, ]
n_med <- xyz[2, ]
n_notreat <- xyz[3, ]

## number of visits to the medical doctor
n_visit_med <- runif(n, n_visit_med_min, n_visit_med_max)

## number of visits to the traditional healer 
n_visit_heal <- runif(n, n_year_heal_min, n_year_heal_max) 

## price of a traditional healer
price_heal <- rgamma(n, price_heal_shape, price_heal_rate)

## price of a visit to a medical doctor 
price_med <- rgamma(n, price_med_shape, price_med_rate)

## Phenobarbital use
p_pheno <- rbeta(n, p_pheno_alpha, p_pheno_beta)

## probability unemployed due to epilepsy
unemployed_duetoepilepsy <-
  rbeta(n, unemployed_duetoepilepsy_alpha, unemployed_duetoepilepsy_beta)

## working days per year
working_days <- runif(n, working_days_min, working_days_max)

## loss of work
loss_workingtime <- runif(n, loss_workingtime_min, loss_workingtime_max)

## monthly salary
monthly_salary <- rgamma(n, monthly_salary_shape, monthly_salary_rate)


## -------------------------------------------------------------------------
## SIMULATIONS / PIGS
##

## porcine cysticercosis prevalence
prev_pigs <- rbeta(n, prev_pigs_alpha, prev_pigs_beta)

## value of pig
price_pigs <- rgamma(n, price_pigs_shape, price_pigs_rate)


## -------------------------------------------------------------------------
## CALCULATIONS / HUMANS
##

## number of phenobarbital users
n_pheno <- (n_hosp + n_med) * p_pheno

## number of inactivity days
n_days1 <- n_ncc * active * unemployed_duetoepilepsy * working_days
n_days2 <- n_ncc * active * (1 - unemployed_duetoepilepsy) * loss_workingtime
n_days_inactivity <- n_days1 + n_days2


## -------------------------------------------------------------------------
## COSTS / HUMANS
##

## hospitalization costs
cost_hosp <- stay * price_day_hosp

## healthcare provider costs
cost_med <- n_visit_med * price_med 

## traditional healer costs
cost_heal <- n_visit_heal * price_heal

## medication costs
cost_medicine <- n_pheno * price_pheno

## productivity losses
cost_inactivity <-
  n_days_inactivity * monthly_salary / n_working_day_by_month


## -------------------------------------------------------------------------
## COSTS / PIGS
##

## number of infected pigs
n_pigs_infected <- n_pigs_smallscale * prev_pigs

## losses due to porcine cysticercosis
cost_pigs <-
  n_pigs_infected * pigs_sold * price_loss_pigs * price_pigs


## -------------------------------------------------------------------------
## COSTS / TOTAL
##

cost_total <-
  cost_hosp + cost_med + cost_heal + cost_medicine + cost_inactivity + 
  cost_pigs

cost_ncc <- cost_total - cost_pigs

cost_by_ncc <- cost_ncc / n_ncc


## -------------------------------------------------------------------------
## SUMMARIES
##

mean_ci(cost_total)
mean_ci(cost_ncc)
mean_ci(cost_pigs)

mean_ci(cost_by_ncc)

mean_ci(cost_hosp)
mean_ci(cost_med)
mean_ci(cost_heal)
mean_ci(cost_medicine)
mean_ci(cost_inactivity)

mean_ci(n_ncc)
mean_ci(n_pigs_infected)
mean_ci(n_pigs_infected / n_pigs_smallscale)

mean_ci(cost_hosp / cost_ncc)
mean_ci(cost_med / cost_ncc)
mean_ci(cost_heal / cost_ncc)
mean_ci(cost_medicine / cost_ncc)
mean_ci(cost_inactivity / cost_ncc)

mean_ci(cost_pigs / cost_total)
mean_ci(cost_inactivity / cost_total)
mean_ci((cost_hosp + cost_med + cost_heal + cost_medicine) /
          cost_total)
