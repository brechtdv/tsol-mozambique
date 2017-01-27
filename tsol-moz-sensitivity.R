### TSOL MOZAMBIQUE
### -- sensitivity analyses

## attach packages
library(ggplot2)

## import helper functions
source("helpers.R")


###
### DALYs
###

## load DALY script
source("tsol-moz-daly.R")

## total
sa_daly_pcc <-
sa_pcc(daly_all,
       cbind(prev_ep_act, prop_ncc, prop_ep_trt,
             prop_mgr, prop_tns, mrt_ep,
             dsw_ep_tr, dsw_ep_nt, dsw_mgr, dsw_tns))
sa_daly_pcc

## tornado graph
daly_items <-
c("E prevalence", "NCC prevalence in PWE", "E treatment proportion",
  "M prevalence in PWE", "TTHA prevalence in PWE", "E mortality rate",
  "E disability weight (treated)", "E disability weight (untreated)",
  "M disability weight", "TTHA disability weight")

tiff("moz-daly-tornado.tiff",
     6, 3, units = "in", res = 400, compress = "lzw")
tornado(sa_daly_pcc, daly_items)
graphics.off()


###
### COSTS
###

## load cost script
source("tsol-cost.R")

## pigs
sa_pcc(cost_pigs,
       cbind(prev_pigs, price_pigs))

## humans
sa_pcc(cost_ncc,
       cbind(prev_ep_act, prop_ncc,
             p_hosp, stay, n_heal, n_med, n_notreat,
             n_visit_med, n_visit_heal, price_heal, price_med, p_pheno,
             unemployed_duetoepilepsy, working_days, loss_workingtime,
             monthly_salary))

## total
sa_cost_pcc <-
sa_pcc(cost_total,
       cbind(prev_pigs, price_pigs,
             prev_ep_act, prop_ncc,
             p_hosp, stay, n_heal, n_med, n_notreat,
             n_visit_med, n_visit_heal, price_heal, price_med, p_pheno,
             unemployed_duetoepilepsy, working_days, loss_workingtime,
             monthly_salary))
sa_cost_pcc

## tornado graph
cost_items <-
c("Prevalence of PC", "Cost of a pig",
  "Prevalence of E", "NCC prevalence in PWE", 
  "Proportion of hospitalization", "Length of hospital stay",
  "visit healer", "visit medical doctor", "no treatment",
  "No of visits to a medical doctor", "No of visits to a healer",
  "Price traditional healer", "Price medical doctor",
  "Proportion phenobarbital use",
  "Proportion unemployment", "Number of working days", "Woking time lost", 
  "Monthly salary")

tiff("moz-cost-tornado.tiff",
     6, 3, units = "in", res = 400, compress = "lzw")
tornado(sa_cost_pcc, cost_items, p = 0.05)
graphics.off()
