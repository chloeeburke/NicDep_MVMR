# This script was created by Jasmine Khouja 27.04.22.
# Edited for this project by Chloe Burke January 2024.

# The script conducts univariable and multivariable MR exploring the effects of
# nicotine and non-nicotine constituents of tobacco smoke (measured by nicotine
# metabolite ratio [NMR] and cigarettes per day [CPD]) on lifetime depression.

# The SNPs used in this script were selected based on the findings from GSCAN
# (Liu et al., 2019) and Buchwald and colleagues.

################################################################################
##### Load packages #####
# Note: some packages already installed 
################################################################################
library(usethis)
library(remotes)
#install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
#devtools::install_github("mrcieu/ieugwasr") 
#install.packages("googleAuthR")
#ieugwasr::get_access_token() 
library(ieugwasr)
library(googleAuthR)
library(tidyverse)
library(stringr)
library(dplyr)
library(forestplot)
library(plyr)
library(gtable)
library(reshape)
library(gplots)
require(ggplot2)
library(ggplot2)
library(gridExtra)
library(grid)
library(extrafont)
library(plotly)
library(data.table)
library(curl)
#install.packages("MVMR")
#install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
library(MVMR)
library(MendelianRandomization)
library(simex)
library(writexl)

#########################
##### Set directory #####
#########################

setwd(" ")

################################################################################
##### Set SNP lists for NMR, CPD and both #####
################################################################################
NMR_SNPlist<- read.table("NMR_SNPlist.txt", header=FALSE)
NMR_SNPlist<-rename(NMR_SNPlist, c("SNP" = "V1"))
CPD_SNPlist<- read.table("CPD_SNPlist.txt", header=FALSE)
CPD_SNPlist<-rename(CPD_SNPlist, c("SNP" = "V1"))
SNPlist<- read.table("SNPlist.txt", header=FALSE)
SNPlist<-rename(SNPlist, c( "SNP"= "V1"))

################################################################################
##### Removing rs117090198  ##### 
# Note: Removed because the SNP has a very high p-value prior to the conditional independence analysis which is unusual and unexplained
################################################################################
# Authors ran the GWAS to identify conditionally independent SNPs, took top SNP and condition and then does it come past the threshold and have sig p-value
# Really high p-value (i.e., not looking like it is significant, that becomes significant after conditioned on)
# Authors couldn't explain why so excluded from instrument

NMR_SNPlist <- data.frame("SNP"= c(NMR_SNPlist[!NMR_SNPlist$SNP == "rs117090198", ])) 

################################################################################
#####Extract exposure data for MR of NMR and depression #####
################################################################################

nmr_dat <- read_exposure_data("META_nmr_Extended_r.ma",
                              clump = FALSE,
                              sep = "\t",
                              snp_col = "RSID",
                              beta_col = "BETA",
                              se_col = "SE",
                              eaf_col = "AF",
                              effect_allele_col = "ALT",
                              other_allele_col = "REF",
                              pval_col = "PVALUE",
                              samplesize_col = "N",
                              min_pval = 1e-200,
                              log_pval = FALSE
)

nmr_dat_mr <- format_data(
  nmr_dat,
  type = "exposure",
  snps = NMR_SNPlist$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)

################################################################################
##### Extract exposure data for MR of CPD and depression #####
################################################################################

cpd_dat <- read_exposure_data("CigarettesPerDay.WithoutUKB.txt",
                              clump = FALSE,
                              sep = "\t",
                              snp_col = "RSID",
                              beta_col = "BETA",
                              se_col = "SE",
                              eaf_col = "AF",
                              effect_allele_col = "ALT",
                              other_allele_col = "REF",
                              pval_col = "PVALUE",
                              samplesize_col = "N",
                              min_pval = 1e-200,
                              log_pval = FALSE
)

cpd_dat_mr <- format_data(
  cpd_dat,
  type = "exposure",
  snps = CPD_SNPlist$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)

################################################################################
##### Extract outcome data for MR #####
################################################################################

#########
## NMR ##
#########
#Never
mdd_dat_never_nmr <- read_outcome_data(
  "dep_never_imputed.txt",
  snps = nmr_dat_mr$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  min_pval = 1e-200,
  log_pval = FALSE
)
#Ever
mdd_dat_ever_nmr <- read_outcome_data(
  "dep_ever_imputed.txt",
  snps = nmr_dat_mr$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  min_pval = 1e-200,
  log_pval = FALSE
)
#Current
mdd_dat_current_nmr <- read_outcome_data(
  "dep_current_imputed.txt",
  snps = nmr_dat_mr$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  min_pval = 1e-200,
  log_pval = FALSE
)
#Former
mdd_dat_former_nmr <- read_outcome_data(
  "dep_former_imputed.txt",
  snps = nmr_dat_mr$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  min_pval = 1e-200,
  log_pval = FALSE
)

#########
## CPD ##
#########
#Never
mdd_dat_never_cpd <- read_outcome_data(
  "dep_never_imputed.txt",
  snps = cpd_dat_mr$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  min_pval = 1e-200,
  log_pval = FALSE
)
#Ever
mdd_dat_ever_cpd <- read_outcome_data(
  "dep_ever_imputed.txt",
  snps = cpd_dat_mr$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  min_pval = 1e-200,
  log_pval = FALSE
)
#Current
mdd_dat_current_cpd <- read_outcome_data(
  "dep_current_imputed.txt",
  snps = cpd_dat_mr$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  min_pval = 1e-200,
  log_pval = FALSE
)
#Former
mdd_dat_former_cpd <- read_outcome_data(
  "dep_former_imputed.txt",
  snps = cpd_dat_mr$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  min_pval = 1e-200,
  log_pval = FALSE
)

################################################################################
##### Convert odds ratios to log odds #####
# Note: the exposures are already on the log scale.
################################################################################
# To convert binary outcomes from BOLT-LMM files to log odds, the following calculation needs to be made: logOR<-beta(caseprevelance*(1-caseprevelance))

mdd_dat_never_nmr$beta.outcome <- as.numeric(as.character(
  mdd_dat_never_nmr$beta.outcome
))
mdd_dat_ever_nmr$beta.outcome <- as.numeric(as.character(
  mdd_dat_ever_nmr$beta.outcome
))
mdd_dat_current_nmr$beta.outcome <- as.numeric(as.character(
  mdd_dat_current_nmr$beta.outcome
))
mdd_dat_former_nmr$beta.outcome <- as.numeric(as.character(
  mdd_dat_former_nmr$beta.outcome
))
mdd_dat_never_cpd$beta.outcome <- as.numeric(as.character(
  mdd_dat_never_cpd$beta.outcome
))
mdd_dat_ever_cpd$beta.outcome <- as.numeric(as.character(
  mdd_dat_ever_cpd$beta.outcome
))
mdd_dat_current_cpd$beta.outcome <- as.numeric(as.character(
  mdd_dat_current_cpd$beta.outcome
))
mdd_dat_former_cpd$beta.outcome <- as.numeric(as.character(
  mdd_dat_former_cpd$beta.outcome
))

mdd_dat_never_nmr["beta.outcome"] <- (mdd_dat_never_nmr$beta.outcome/(0.20349*((1-0.20349))))
mdd_dat_never_nmr["se.outcome"] <- (mdd_dat_never_nmr$se.outcome/(0.20349*((1-0.20349))))

mdd_dat_ever_nmr["beta.outcome"] <- (mdd_dat_ever_nmr$beta.outcome/(0.25484*((1-0.25484))))
mdd_dat_ever_nmr["se.outcome"] <- (mdd_dat_ever_nmr$se.outcome/(0.25484*((1-0.25484))))

mdd_dat_current_nmr["beta.outcome"] <- (mdd_dat_current_nmr$beta.outcome/(0.30688*((1-0.30688))))
mdd_dat_current_nmr["se.outcome"] <- (mdd_dat_current_nmr$se.outcome/(0.30688*((1-0.30688))))

mdd_dat_former_nmr["beta.outcome"] <- (mdd_dat_former_nmr$beta.outcome/(0.23984*((1-0.23984))))
mdd_dat_former_nmr["se.outcome"] <- (mdd_dat_former_nmr$se.outcome/(0.23984*((1-0.23984))))

mdd_dat_never_cpd["beta.outcome"] <- (mdd_dat_never_cpd$beta.outcome/(0.20349*((1-0.20349))))
mdd_dat_never_cpd["se.outcome"] <- (mdd_dat_never_cpd$se.outcome/(0.20349*((0.20349))))

mdd_dat_ever_cpd["beta.outcome"] <- (mdd_dat_ever_cpd$beta.outcome/(0.25484*((1-0.25484))))
mdd_dat_ever_cpd["se.outcome"] <- (mdd_dat_ever_cpd$se.outcome/(0.25484*((1-0.25484))))

mdd_dat_current_cpd["beta.outcome"] <- (mdd_dat_current_cpd$beta.outcome/(0.30688*((1-0.30688))))
mdd_dat_current_cpd["se.outcome"] <- (mdd_dat_current_cpd$se.outcome/(0.30688*((1-0.30688))))

mdd_dat_former_cpd["beta.outcome"] <- (mdd_dat_former_cpd$beta.outcome/(0.23984*((1-0.23984))))
mdd_dat_former_cpd["se.outcome"] <- (mdd_dat_former_cpd$se.outcome/(0.23984*((1-0.23984))))

################################################################################
##### Harmonising #####
################################################################################
##NMR - none removed
outcome_mdd_dat_nmr_never <-harmonise_data(nmr_dat_mr, mdd_dat_never_nmr, action = 2)
outcome_mdd_dat_nmr_ever <-harmonise_data(nmr_dat_mr, mdd_dat_ever_nmr, action = 2)
outcome_mdd_dat_nmr_current <- harmonise_data(nmr_dat_mr, mdd_dat_current_nmr, action = 2)
outcome_mdd_dat_nmr_former <- harmonise_data(nmr_dat_mr, mdd_dat_former_nmr, action = 2)
##CPD - SNPs for being palindromic with intermediate allele frequencies: rs1737894, rs28438420
outcome_mdd_dat_cpd_never <-harmonise_data(cpd_dat_mr, mdd_dat_never_cpd, action = 2)
outcome_mdd_dat_cpd_ever <-harmonise_data(cpd_dat_mr, mdd_dat_ever_cpd, action = 2)
outcome_mdd_dat_cpd_current <-harmonise_data(cpd_dat_mr, mdd_dat_current_cpd, action = 2)
outcome_mdd_dat_cpd_former <-harmonise_data(cpd_dat_mr, mdd_dat_former_cpd, action = 2)

################################################################################
##### Find proxies to add to missing outcome SNPs #####
################################################################################
#Checking NMR - none identified
proxy_needed1 <- data.frame(setdiff(NMR_SNPlist$SNP, mdd_dat_never_nmr$SNP))
proxy_needed2 <- data.frame(setdiff(NMR_SNPlist$SNP, mdd_dat_ever_nmr$SNP))
proxy_needed3 <- data.frame(setdiff(NMR_SNPlist$SNP, mdd_dat_current_nmr$SNP))
proxy_needed4 <- data.frame(setdiff(NMR_SNPlist$SNP, mdd_dat_former_nmr$SNP))
#Checking CPD - none identified
proxy_needed5 <- data.frame(setdiff(CPD_SNPlist$SNP, mdd_dat_never_cpd$SNP))
proxy_needed6 <- data.frame(setdiff(CPD_SNPlist$SNP, mdd_dat_ever_cpd$SNP))
proxy_needed7 <- data.frame(setdiff(CPD_SNPlist$SNP, mdd_dat_current_cpd$SNP))
proxy_needed8 <- data.frame(setdiff(CPD_SNPlist$SNP, mdd_dat_former_cpd$SNP))
#No proxies needed, removing from environment
objects_to_remove <- ls(pattern = "^proxy_needed")
rm(list = objects_to_remove)
rm(objects_to_remove)

################################################################################
################################## MR ##########################################
################################################################################

################################################################################
##### Generate results inc. F and Q stats for heterogeneity for: #####
# NMR never
# NMR ever
# NMR current
# NMR former
# CPD never 
# CPD ever
# CPD current
# CPD former
################################################################################
# NMR never
result_nmr_never <- mr(
  outcome_mdd_dat_nmr_never,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
result_nmr_never <- generate_odds_ratios(result_nmr_never)
outcome_mdd_dat_nmr_never <- subset(outcome_mdd_dat_nmr_never, mr_keep)
mr_n_nmr <- TwoSampleMR::mr_ivw(
  outcome_mdd_dat_nmr_never$beta.exposure,
  outcome_mdd_dat_nmr_never$beta.outcome,
  outcome_mdd_dat_nmr_never$se.exposure,
  outcome_mdd_dat_nmr_never$se.outcome,
  parameters = default_parameters()
)
ptr_n_nmr <- data.frame(mr_n_nmr["Q"])
egger_n_nmr <- mr_egger_regression(
  outcome_mdd_dat_nmr_never$beta.exposure,
  outcome_mdd_dat_nmr_never$beta.outcome,
  outcome_mdd_dat_nmr_never$se.exposure,
  outcome_mdd_dat_nmr_never$se.outcome,
  parameters
)
ptr_n_nmr[2, 1] <- egger_n_nmr["Q"]
F <- abs(outcome_mdd_dat_nmr_never$beta.exposure^2 / outcome_mdd_dat_nmr_never$se.exposure^2)
mF <- mean(F)
ptr_n_nmr[3, 1] <- mF

# NMR ever
result_nmr_ever <- mr(
  outcome_mdd_dat_nmr_ever,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
result_nmr_ever <- generate_odds_ratios(result_nmr_ever)
outcome_mdd_dat_nmr_ever <- subset(outcome_mdd_dat_nmr_ever, mr_keep)
mr_e_nmr <- TwoSampleMR::mr_ivw(
  outcome_mdd_dat_nmr_ever$beta.exposure,
  outcome_mdd_dat_nmr_ever$beta.outcome,
  outcome_mdd_dat_nmr_ever$se.exposure,
  outcome_mdd_dat_nmr_ever$se.outcome,
  parameters = default_parameters()
)
ptr_e_nmr <- data.frame(mr_e_nmr["Q"])
egger_e_nmr <- mr_egger_regression(
  outcome_mdd_dat_nmr_ever$beta.exposure,
  outcome_mdd_dat_nmr_ever$beta.outcome,
  outcome_mdd_dat_nmr_ever$se.exposure,
  outcome_mdd_dat_nmr_ever$se.outcome,
  parameters
)
ptr_e_nmr[2, 1] <- egger_e_nmr["Q"]
F <- abs(outcome_mdd_dat_nmr_ever$beta.exposure^2 / outcome_mdd_dat_nmr_ever$se.exposure^2)
mF <- mean(F)
ptr_e_nmr[3, 1] <- mF

# NMR current
result_nmr_current <- mr(
  outcome_mdd_dat_nmr_current,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
result_nmr_current <- generate_odds_ratios(result_nmr_current)
outcome_mdd_dat_nmr_current <- subset(outcome_mdd_dat_nmr_current, mr_keep)
mr_c_nmr <- TwoSampleMR::mr_ivw(
  outcome_mdd_dat_nmr_current$beta.exposure,
  outcome_mdd_dat_nmr_current$beta.outcome,
  outcome_mdd_dat_nmr_current$se.exposure,
  outcome_mdd_dat_nmr_current$se.outcome,
  parameters = default_parameters()
)
ptr_c_nmr <- data.frame(mr_c_nmr["Q"])
egger_c_nmr <- mr_egger_regression(
  outcome_mdd_dat_nmr_current$beta.exposure,
  outcome_mdd_dat_nmr_current$beta.outcome,
  outcome_mdd_dat_nmr_current$se.exposure,
  outcome_mdd_dat_nmr_current$se.outcome,
  parameters
)
ptr_c_nmr[2, 1] <- egger_c_nmr["Q"]
F <- abs(outcome_mdd_dat_nmr_current$beta.exposure^2 / outcome_mdd_dat_nmr_current$se.exposure^2)
mF <- mean(F)
ptr_c_nmr[3, 1] <- mF

# NMR former
result_nmr_former <- mr(
  outcome_mdd_dat_nmr_former,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
result_nmr_former <- generate_odds_ratios(result_nmr_former)
outcome_mdd_dat_nmr_former <- subset(outcome_mdd_dat_nmr_former, mr_keep)
mr_f_nmr <- TwoSampleMR::mr_ivw(
  outcome_mdd_dat_nmr_former$beta.exposure,
  outcome_mdd_dat_nmr_former$beta.outcome,
  outcome_mdd_dat_nmr_former$se.exposure,
  outcome_mdd_dat_nmr_former$se.outcome,
  parameters = default_parameters()
)
ptr_f_nmr <- data.frame(mr_f_nmr["Q"])
egger_f_nmr <- mr_egger_regression(
  outcome_mdd_dat_nmr_former$beta.exposure,
  outcome_mdd_dat_nmr_former$beta.outcome,
  outcome_mdd_dat_nmr_former$se.exposure,
  outcome_mdd_dat_nmr_former$se.outcome,
  parameters
)
ptr_f_nmr[2, 1] <- egger_f_nmr["Q"]
F <- abs(outcome_mdd_dat_nmr_former$beta.exposure^2 / outcome_mdd_dat_nmr_former$se.exposure^2)
mF <- mean(F)
ptr_f_nmr[3, 1] <- mF

# CPD never
result_cpd_never <- mr(
  outcome_mdd_dat_cpd_never,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
result_cpd_never <- generate_odds_ratios(result_cpd_never)
outcome_mdd_dat_cpd_never <- subset(outcome_mdd_dat_cpd_never, mr_keep)
mr_n_cpd <- TwoSampleMR::mr_ivw(
  outcome_mdd_dat_cpd_never$beta.exposure,
  outcome_mdd_dat_cpd_never$beta.outcome,
  outcome_mdd_dat_cpd_never$se.exposure,
  outcome_mdd_dat_cpd_never$se.outcome,
  parameters = default_parameters()
)
ptr_n_cpd <- data.frame(mr_n_cpd["Q"])
egger_n_cpd <- mr_egger_regression(
  outcome_mdd_dat_cpd_never$beta.exposure,
  outcome_mdd_dat_cpd_never$beta.outcome,
  outcome_mdd_dat_cpd_never$se.exposure,
  outcome_mdd_dat_cpd_never$se.outcome,
  parameters
)
ptr_n_cpd[2, 1] <- egger_n_cpd["Q"]
F <- abs(outcome_mdd_dat_cpd_never$beta.exposure^2 / outcome_mdd_dat_cpd_never$se.exposure^2)
mF <- mean(F)
ptr_n_cpd[3, 1] <- mF

# CPD ever
result_cpd_ever <- mr(
  outcome_mdd_dat_cpd_ever,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
result_cpd_ever <- generate_odds_ratios(result_cpd_ever)
outcome_mdd_dat_cpd_ever <- subset(outcome_mdd_dat_cpd_ever, mr_keep)
mr_e_cpd <- TwoSampleMR::mr_ivw(
  outcome_mdd_dat_cpd_ever$beta.exposure,
  outcome_mdd_dat_cpd_ever$beta.outcome,
  outcome_mdd_dat_cpd_ever$se.exposure,
  outcome_mdd_dat_cpd_ever$se.outcome,
  parameters = default_parameters()
)
ptr_e_cpd <- data.frame(mr_e_cpd["Q"])
egger_e_cpd <- mr_egger_regression(
  outcome_mdd_dat_cpd_ever$beta.exposure,
  outcome_mdd_dat_cpd_ever$beta.outcome,
  outcome_mdd_dat_cpd_ever$se.exposure,
  outcome_mdd_dat_cpd_ever$se.outcome,
  parameters
)
ptr_e_cpd[2, 1] <- egger_e_cpd["Q"]
F <- abs(outcome_mdd_dat_cpd_ever$beta.exposure^2 / outcome_mdd_dat_cpd_ever$se.exposure^2)
mF <- mean(F)
ptr_e_cpd[3, 1] <- mF

# CPD current
result_cpd_current <- mr(
  outcome_mdd_dat_cpd_current,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
result_cpd_current <- generate_odds_ratios(result_cpd_current)
outcome_mdd_dat_cpd_current <- subset(outcome_mdd_dat_cpd_current, mr_keep)
mr_c_cpd <- TwoSampleMR::mr_ivw(
  outcome_mdd_dat_cpd_current$beta.exposure,
  outcome_mdd_dat_cpd_current$beta.outcome,
  outcome_mdd_dat_cpd_current$se.exposure,
  outcome_mdd_dat_cpd_current$se.outcome,
  parameters = default_parameters()
)
ptr_c_cpd <- data.frame(mr_c_cpd["Q"])
egger_c_cpd <- mr_egger_regression(
  outcome_mdd_dat_cpd_current$beta.exposure,
  outcome_mdd_dat_cpd_current$beta.outcome,
  outcome_mdd_dat_cpd_current$se.exposure,
  outcome_mdd_dat_cpd_current$se.outcome,
  parameters
)
ptr_c_cpd[2, 1] <- egger_c_cpd["Q"]
F <- abs(outcome_mdd_dat_cpd_current$beta.exposure^2 / outcome_mdd_dat_cpd_current$se.exposure^2)
mF <- mean(F)
ptr_c_cpd[3, 1] <- mF

# CPD former 
result_cpd_former <- mr(
  outcome_mdd_dat_cpd_former,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
result_cpd_former <- generate_odds_ratios(result_cpd_former)
outcome_mdd_dat_cpd_former <- subset(outcome_mdd_dat_cpd_former, mr_keep)
mr_f_cpd <- TwoSampleMR::mr_ivw(
  outcome_mdd_dat_cpd_former$beta.exposure,
  outcome_mdd_dat_cpd_former$beta.outcome,
  outcome_mdd_dat_cpd_former$se.exposure,
  outcome_mdd_dat_cpd_former$se.outcome,
  parameters = default_parameters()
)
ptr_f_cpd <- data.frame(mr_f_cpd["Q"])
egger_f_cpd <- mr_egger_regression(
  outcome_mdd_dat_cpd_former$beta.exposure,
  outcome_mdd_dat_cpd_former$beta.outcome,
  outcome_mdd_dat_cpd_former$se.exposure,
  outcome_mdd_dat_cpd_former$se.outcome,
  parameters
)
ptr_f_cpd[2, 1] <- egger_f_cpd["Q"]
F <- abs(outcome_mdd_dat_cpd_former$beta.exposure^2/outcome_mdd_dat_cpd_former$se.exposure^2)
mF <- mean(F)
ptr_f_cpd[3, 1] <- mF

################################################################################
##### Calculate I2GX #####
# Isq(y, s) where y = vector of effects and s = vector of standard errors
# Apply simulation extrapolation SIMEX corrections to MR-Egger analysis where
# I2GX estimates < 0.9
# This indicates the effect estimate is biased by 10% due to measurement error
################################################################################
##Manual ISQ function:
ISQ <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  ISQ        = (Q - (k-1))/Q
  ISQ        = max(0,ISQ)
  return(ISQ)
}
##NMR
#Never
BetaXG   = outcome_mdd_dat_nmr_never$beta.exposure
seBetaXG = outcome_mdd_dat_nmr_never$se.exposure 
seBetaYG <- outcome_mdd_dat_nmr_never$se.outcome
BXG  = abs(BetaXG)         # gene--exposure estimates are positive  
Isq_unweighted_nmr_never <- ISQ(BXG,seBetaXG) #unweighted = 0.99
Isq_weighted_nmr_never <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0.99
#Ever
BetaXG   = outcome_mdd_dat_nmr_ever$beta.exposure
seBetaXG = outcome_mdd_dat_nmr_ever$se.exposure 
seBetaYG <- outcome_mdd_dat_nmr_ever$se.outcome
BXG  = abs(BetaXG)         # gene--exposure estimates are positive  
Isq_unweighted_nmr_ever <- ISQ(BXG,seBetaXG) #unweighted = 0.99
Isq_weighted_nmr_ever <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0.99
#Current
BetaXG   = outcome_mdd_dat_nmr_current$beta.exposure
seBetaXG = outcome_mdd_dat_nmr_current$se.exposure 
seBetaYG <- outcome_mdd_dat_nmr_current$se.outcome
BXG  = abs(BetaXG)         # gene--exposure estimates are positive  
Isq_unweighted_nmr_current <- ISQ(BXG,seBetaXG) #unweighted = 0.99
Isq_weighted_nmr_current <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0.99
#Former
BetaXG   = outcome_mdd_dat_nmr_former$beta.exposure
seBetaXG = outcome_mdd_dat_nmr_former$se.exposure 
seBetaYG <- outcome_mdd_dat_nmr_former$se.outcome
BXG  = abs(BetaXG)         # gene--exposure estimates are positive  
Isq_unweighted_nmr_former <- ISQ(BXG,seBetaXG) #unweighted = 0.99
Isq_weighted_nmr_former <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0.99
##CPD
#Never
BetaXG   = outcome_mdd_dat_cpd_never$beta.exposure
seBetaXG = outcome_mdd_dat_cpd_never$se.exposure 
seBetaYG <- outcome_mdd_dat_cpd_never$se.outcome
BXG  = abs(BetaXG)         # gene--exposure estimates are positive  
Isq_unweighted_cpd_never <- ISQ(BXG,seBetaXG) #unweighted = 0.97
Isq_weighted_cpd_never <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0.97
#Ever
BetaXG   = outcome_mdd_dat_cpd_ever$beta.exposure
seBetaXG = outcome_mdd_dat_cpd_ever$se.exposure 
seBetaYG <- outcome_mdd_dat_cpd_ever$se.outcome
BXG  = abs(BetaXG)         # gene--exposure estimates are positive  
Isq_unweighted_cpd_ever <- ISQ(BXG,seBetaXG) #unweighted = 0.97
Isq_weighted_cpd_ever <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0.97
#Current
BetaXG   = outcome_mdd_dat_cpd_current$beta.exposure
seBetaXG = outcome_mdd_dat_cpd_current$se.exposure 
seBetaYG <- outcome_mdd_dat_cpd_current$se.outcome
BXG  = abs(BetaXG)         # gene--exposure estimates are positive  
Isq_unweighted_cpd_current <- ISQ(BXG,seBetaXG) #unweighted = 0.97
Isq_weighted_cpd_current <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0.97
#Former
BetaXG   = outcome_mdd_dat_cpd_former$beta.exposure
seBetaXG = outcome_mdd_dat_cpd_former$se.exposure 
seBetaYG <- outcome_mdd_dat_cpd_former$se.outcome
BXG  = abs(BetaXG)         # gene--exposure estimates are positive  
Isq_unweighted_cpd_former <- ISQ(BXG,seBetaXG) #unweighted = 0.97
Isq_weighted_cpd_former <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0.97

ISQ <- data.frame(c(1:8))
ISQ[1, 1] <- Isq_unweighted_nmr_never
ISQ[2, 1] <- Isq_unweighted_nmr_ever
ISQ[3, 1] <- Isq_unweighted_nmr_current
ISQ[4, 1] <- Isq_unweighted_nmr_former
ISQ[5, 1] <- Isq_unweighted_cpd_never
ISQ[6, 1] <- Isq_unweighted_cpd_ever
ISQ[7, 1] <- Isq_unweighted_cpd_current
ISQ[8, 1] <- Isq_unweighted_cpd_former

################################################################################
#########          Test for reverse causation                ############
################################################################################
### Steiger filtering ###

##Ever-smokers CPD - treating CPD as continuous
outcome_mdd_dat_cpd_ever$samplesize.outcome<- 160248
outcome_mdd_dat_cpd_ever$ncase.outcome<-40838
outcome_mdd_dat_cpd_ever$ncontrol.outcome<- 119410
outcome_mdd_dat_cpd_ever$units.outcome<-"log odds"
outcome_mdd_dat_cpd_ever$prevalence.outcome<-0.25
outcome_mdd_dat_cpd_ever$units.exposure<-"SD"
cpd_e_steiger <- subset(outcome_mdd_dat_cpd_ever, outcome_mdd_dat_cpd_ever$mr_keep==TRUE)
cpd_e_steiger <- steiger_filtering(cpd_e_steiger)

#Report the proportion of SNPs that are true
false <- length(cpd_e_steiger$steiger_dir[cpd_e_steiger$steiger_dir == FALSE])
true <- length(cpd_e_steiger$steiger_dir[cpd_e_steiger$steiger_dir == TRUE])
percent <- (true/(true+false))*100
print(percent) #100%

##Current-smokers CPD - treating CPD as continuous
outcome_mdd_dat_cpd_current$samplesize.outcome<- 35871
outcome_mdd_dat_cpd_current$ncase.outcome<-11008
outcome_mdd_dat_cpd_current$ncontrol.outcome<- 24863
outcome_mdd_dat_cpd_current$units.outcome<-"log odds"
outcome_mdd_dat_cpd_current$prevalence.outcome<-0.31
outcome_mdd_dat_cpd_current$units.exposure<-"SD"
cpd_c_steiger <- subset(outcome_mdd_dat_cpd_current, outcome_mdd_dat_cpd_current$mr_keep==TRUE)
cpd_c_steiger <- steiger_filtering(cpd_c_steiger)

#Report the proportion of SNPs that are true
false <- length(cpd_c_steiger$steiger_dir[cpd_c_steiger$steiger_dir == FALSE])
true <- length(cpd_c_steiger$steiger_dir[cpd_c_steiger$steiger_dir == TRUE])
percent <- (true/(true+false))*100
print(percent) #94.3% 
summary(cpd_c_steiger$steiger_dir)
# FALSE    TRUE 
#     3      50 

##Former-smokers CPD - treating CPD as continuous
outcome_mdd_dat_cpd_former$samplesize.outcome<- 124377
outcome_mdd_dat_cpd_former$ncase.outcome<-29831
outcome_mdd_dat_cpd_former$ncontrol.outcome<- 94546
outcome_mdd_dat_cpd_former$units.outcome<-"log odds"
outcome_mdd_dat_cpd_former$prevalence.outcome<-0.24
outcome_mdd_dat_cpd_former$units.exposure<-"SD"
cpd_f_steiger <- subset(outcome_mdd_dat_cpd_former, outcome_mdd_dat_cpd_former$mr_keep==TRUE)
cpd_f_steiger <- steiger_filtering(cpd_f_steiger)

#Report the proportion of SNPs that are true
false <- length(cpd_f_steiger$steiger_dir[cpd_f_steiger$steiger_dir == FALSE])
true <- length(cpd_f_steiger$steiger_dir[cpd_f_steiger$steiger_dir == TRUE])
percent <- (true/(true+false))*100
print(percent) #100%

##Never-smokers CPD - treating CPD as continuous
outcome_mdd_dat_cpd_never$samplesize.outcome<- 194881
outcome_mdd_dat_cpd_never$ncase.outcome<-39656
outcome_mdd_dat_cpd_never$ncontrol.outcome<- 155225
outcome_mdd_dat_cpd_never$units.outcome<-"log odds"
outcome_mdd_dat_cpd_never$prevalence.outcome<-0.20
outcome_mdd_dat_cpd_never$units.exposure<-"SD"
cpd_n_steiger <- subset(outcome_mdd_dat_cpd_never, outcome_mdd_dat_cpd_never$mr_keep==TRUE)
cpd_n_steiger <- steiger_filtering(cpd_n_steiger)

#Report the proportion of SNPs that are true
false <- length(cpd_n_steiger$steiger_dir[cpd_n_steiger$steiger_dir == FALSE])
true <- length(cpd_n_steiger$steiger_dir[cpd_n_steiger$steiger_dir == TRUE])
percent <- (true/(true+false))*100
print(percent) #100%

##Ever-smokers NMR - treating NMR as continuous
outcome_mdd_dat_nmr_ever$samplesize.outcome<- 160248
outcome_mdd_dat_nmr_ever$ncase.outcome<-40838
outcome_mdd_dat_nmr_ever$ncontrol.outcome<- 119410
outcome_mdd_dat_nmr_ever$units.outcome<-"log odds"
outcome_mdd_dat_nmr_ever$prevalence.outcome<-0.25
outcome_mdd_dat_nmr_ever$units.exposure<-"SD"
nmr_e_steiger <- subset(outcome_mdd_dat_nmr_ever, outcome_mdd_dat_nmr_ever$mr_keep==TRUE)
nmr_e_steiger <- steiger_filtering(nmr_e_steiger)

#Report the proportion of SNPs that are true
false <- length(nmr_e_steiger$steiger_dir[nmr_e_steiger$steiger_dir == FALSE])
true <- length(nmr_e_steiger$steiger_dir[nmr_e_steiger$steiger_dir == TRUE])
percent <- (true/(true+false))*100
print(percent) #100%

##Current-smokers NMR - treating NMR as continuous
outcome_mdd_dat_nmr_current$samplesize.outcome<- 35871
outcome_mdd_dat_nmr_current$ncase.outcome<-11008
outcome_mdd_dat_nmr_current$ncontrol.outcome<- 24863
outcome_mdd_dat_nmr_current$units.outcome<-"log odds"
outcome_mdd_dat_nmr_current$prevalence.outcome<-0.31
outcome_mdd_dat_nmr_current$units.exposure<-"SD"
nmr_c_steiger <- subset(outcome_mdd_dat_nmr_current, outcome_mdd_dat_nmr_current$mr_keep==TRUE)
nmr_c_steiger <- steiger_filtering(nmr_c_steiger)

#Report the proportion of SNPs that are true
false <- length(nmr_c_steiger$steiger_dir[nmr_c_steiger$steiger_dir == FALSE])
true <- length(nmr_c_steiger$steiger_dir[nmr_c_steiger$steiger_dir == TRUE])
percent <- (true/(true+false))*100
print(percent) #100%

##Former-smokers NMR - treating NMR as continuous
outcome_mdd_dat_nmr_former$samplesize.outcome<- 124377
outcome_mdd_dat_nmr_former$ncase.outcome<-29831
outcome_mdd_dat_nmr_former$ncontrol.outcome<- 94546
outcome_mdd_dat_nmr_former$units.outcome<-"log odds"
outcome_mdd_dat_nmr_former$prevalence.outcome<-0.24
outcome_mdd_dat_nmr_former$units.exposure<-"SD"
nmr_f_steiger <- subset(outcome_mdd_dat_nmr_former, outcome_mdd_dat_nmr_former$mr_keep==TRUE)
nmr_f_steiger <- steiger_filtering(nmr_f_steiger)

#Report the proportion of SNPs that are true
false <- length(nmr_f_steiger$steiger_dir[nmr_f_steiger$steiger_dir == FALSE])
true <- length(nmr_f_steiger$steiger_dir[nmr_f_steiger$steiger_dir == TRUE])
percent <- (true/(true+false))*100
print(percent) #100%

##Never-smokers NMR - treating NMR as continuous
outcome_mdd_dat_nmr_never$samplesize.outcome<- 194881
outcome_mdd_dat_nmr_never$ncase.outcome<-39656
outcome_mdd_dat_nmr_never$ncontrol.outcome<- 155225
outcome_mdd_dat_nmr_never$units.outcome<-"log odds"
outcome_mdd_dat_nmr_never$prevalence.outcome<-0.20
outcome_mdd_dat_nmr_never$units.exposure<-"SD"
nmr_n_steiger <- subset(outcome_mdd_dat_nmr_never, outcome_mdd_dat_nmr_never$mr_keep==TRUE)
nmr_n_steiger <- steiger_filtering(nmr_n_steiger)

#Report the proportion of SNPs that are true
false <- length(nmr_n_steiger$steiger_dir[nmr_n_steiger$steiger_dir == FALSE])
true <- length(nmr_n_steiger$steiger_dir[nmr_n_steiger$steiger_dir == TRUE])
percent <- (true/(true+false))*100
print(percent) #100%

### Re-running MR analyses which filtered instrument for CPD on MDD among current smokers
#Filtering
outcome_mdd_dat_cpd_current_filt <- subset(cpd_c_steiger, cpd_c_steiger$steiger_dir==TRUE)
#MR
#Results are very similar to analysis pre-filtered; not reported as current smokers are a supplementary analysis to ever and never groups
steiger_cpd_current <- mr(
  outcome_mdd_dat_cpd_current_filt,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
steiger_cpd_current <- generate_odds_ratios(steiger_cpd_current)
outcome_mdd_dat_cpd_current_filt <- subset(outcome_mdd_dat_cpd_current_filt, mr_keep)
mr_c_cpd_steiger <- TwoSampleMR::mr_ivw(
  outcome_mdd_dat_cpd_current_filt$beta.exposure,
  outcome_mdd_dat_cpd_current_filt$beta.outcome,
  outcome_mdd_dat_cpd_current_filt$se.exposure,
  outcome_mdd_dat_cpd_current_filt$se.outcome,
  parameters = default_parameters()
)
ptr_c_cpd_steiger <- data.frame(mr_c_cpd_steiger["Q"])
egger_c_cpd_steiger <- mr_egger_regression(
  outcome_mdd_dat_cpd_current_filt$beta.exposure,
  outcome_mdd_dat_cpd_current_filt$beta.outcome,
  outcome_mdd_dat_cpd_current_filt$se.exposure,
  outcome_mdd_dat_cpd_current_filt$se.outcome,
  parameters
)
ptr_c_cpd_steiger[2, 1] <- egger_c_cpd_steiger["Q"]
F <- outcome_mdd_dat_cpd_current_filt$beta.exposure^2 / outcome_mdd_dat_cpd_current_filt$se.exposure^2
mF <- mean(F)
ptr_c_cpd_steiger[3, 1] <- mF

################################################################################
################################   MVMR  #######################################
################################################################################

NMR_instr <- nmr_dat_mr[, c(
  "SNP",
  "effect_allele.exposure",
  "se.exposure",
  "pval.exposure",
  "beta.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure",
  "samplesize.exposure"
)]
NMR_instr <- NMR_instr[!NMR_instr$SNP == "rs117090198", ]

CPD_instr <- cpd_dat_mr[, c(
  "SNP",
  "effect_allele.exposure",
  "se.exposure",
  "pval.exposure",
  "beta.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure",
  "samplesize.exposure"
)]

INSTR <- do.call("rbind", list(NMR_instr, CPD_instr))

################################################################################
##### Check for overlapping SNPs between the exposure instruments #####
################################################################################

n_occur <- data.frame(table(INSTR$SNP))
n_occur[n_occur$Freq > 1, ]
INSTR[INSTR$SNP %in% n_occur$Var1[n_occur$Freq > 1], ] #N = 2, rs117824460 and rs56113850

################################################################################
##### Extract instruments from the NMR data set #####
################################################################################
nmr_dat_mvmr <- format_data(
  nmr_dat,
  type = "exposure",
  snps = INSTR$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)

##### Change name of GWAS and check n SNPs ####
nmr_dat_mvmr$id.exposure <- "1"
str(nmr_dat_mvmr) #SNPs = 58/61 (i.e., missing SNP = 1/59 due to two overlapping)

##### Check rsid of missing SNP ##### 
diff_values <- setdiff(n_occur$Var1, nmr_dat_mvmr$SNP) #SNP = rs4886550

################################################################################
##### Extract instruments from the CPD data set #####
################################################################################
cpd_dat_mvmr <- format_data(
  cpd_dat,
  type = "exposure",
  snps = INSTR$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)

##### Change name of GWAS and check n SNPs #####
cpd_dat_mvmr$id.exposure <- "2"
str(cpd_dat_mvmr) #SNPs = 59 (i.e., no missing SNPs, as two are overlapping)

##### Change name of GWAS and check n SNPs #####
diff_values <- setdiff(n_occur$Var1, cpd_dat_mvmr$SNP) #SNP = NA

################################################################################
##### Merge them to extract from the other exposure as 'outcome' #####
################################################################################

nmr_dat_mvmr_1 <- subset(nmr_dat_mvmr, select = c(
  "SNP", "effect_allele.exposure",
  "other_allele.exposure",
  "beta.exposure", "se.exposure",
  "pval.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure",
  "eaf.exposure",
  "samplesize.exposure"
))

cpd_dat_mvmr_1 <- subset(cpd_dat_mvmr, select = c(
  "SNP",
  "effect_allele.exposure",
  "other_allele.exposure",
  "beta.exposure",
  "se.exposure",
  "pval.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure",
  "eaf.exposure",
  "samplesize.exposure"
))

##### check structure is the same #####
str(nmr_dat_mvmr_1)
str(cpd_dat_mvmr_1)

################################################################################
##### Merge ####
################################################################################
exposures <- do.call("rbind", list(nmr_dat_mvmr_1, cpd_dat_mvmr_1))

################################################################################
##### Find proxies missing from either the NMR or CPD data set #####
# Note: rs4886550 is not available in the ref panel and is the only missing SNP
################################################################################
#LDProxy: "rs4886550 Variant is not in 1000G reference panel"
#https://ldlink.nih.gov/?tab=ldproxy
n_occur <- data.frame(table(exposures$SNP))
n_occur[n_occur$Freq < 2, ]
exposures[exposures$SNP %in% n_occur$Var1[n_occur$Freq < 2], ]

################################################################################
##### Clumping ####
################################################################################
##### Change all p-values for NMR to 1e-200 for clumping so that none are dropped #####
##### Save old p-values first #####
exposures$oldpvalues <- exposures$pval.exposure
exposures <- exposures %>%
  mutate(pval.exposure = if_else(
    exposures$SNP %in% NMR_instr$SNP,
    1e-201, pval.exposure
  ))

##### Clump the data #####
exposures$id.exposure[exposures$id.exposure == "2"] <- "1"
exposures <- clump_data(exposures, clump_kb = 500, clump_r2 = 0.1) #Different to pre-reg which is 0.001
str(exposures)

##### Add ID's back #####
exposures$id.exposure[exposures$samplesize.exposure < 6000] <- "1"
exposures$id.exposure[exposures$samplesize.exposure > 6000] <- "2"

##### Revert all p-values for NMR from 1e-200 #####
exposures$pval.exposure <- exposures$oldpvalues
exposures <- select(exposures, -c(oldpvalues))

##### Added rs56113850 back in ##### - slightly inconsistent but probably due to imbalanaced instruments
##### Removed due to LD with other NMR SNPs but conditionally independent #####
rs56113850 <- do.call("rbind", list(nmr_dat_mvmr_1, cpd_dat_mvmr_1))
rs56113850 <- rs56113850[rs56113850$SNP == "rs56113850", ]
exposures <- do.call("rbind", list(exposures, rs56113850))

##### Split again to harmonise based on exposure id #####
NMR <- split(exposures, exposures$id.exposure)[["1"]]
CPD <- split(exposures, exposures$id.exposure)[["2"]]

################################################################################
##### Harmonise NMR on CPD #####
################################################################################

names(NMR) <- gsub("exposure", "outcome", names(NMR))
NMR_CPD <- harmonise_data(CPD, NMR)
#Removing the following SNPs for being palindromic with intermediate allele frequencies: rs1737894

################################################################################
##### Keep only snps that are present across both exposures ####
# Note: they would have frequency 1 if only available in one dataset
################################################################################

n_occur <- data.frame(table(exposures$SNP))
n_occur[n_occur$Freq == 2, ]
exposures <- exposures[exposures$SNP %in% n_occur$Var1[n_occur$Freq == 2], ]
str(exposures)

################################################################################
##### Format exposures #####
################################################################################

##### Keep onlySNPs where MrKeep = TRUE #####
NMR_CPD <- NMR_CPD[NMR_CPD$mr_keep == TRUE, ]
str(NMR_CPD)

##### Split the tables - CPD #####
CPD_H <- subset(
  NMR_CPD,
  id.exposure == "2",
  select = c(
    SNP, exposure,
    id.exposure,
    effect_allele.exposure,
    other_allele.exposure,
    beta.exposure,
    se.exposure,
    pval.exposure,
    eaf.exposure
  )
)

##### Split the tables - NMR #####
NMR_H <- subset(NMR_CPD,
                id.outcome == "1",
                select = c(
                  SNP,
                  outcome,
                  id.outcome,
                  effect_allele.outcome,
                  other_allele.outcome,
                  beta.outcome,
                  se.outcome,
                  pval.outcome,
                  eaf.outcome
                )
)

##### Turn NMR from outcome to exposure to merge the data sets #####
names(NMR_H) <- gsub("outcome", "exposure", names(NMR_H))
Exposures_H <- merge(NMR_H, CPD_H, all = TRUE)
Exposures_H["Phenotype"] <- NA
Exposures_H$Phenotype[Exposures_H$id.exposure == 1] <- "NMR"
Exposures_H$Phenotype[Exposures_H$id.exposure == 2] <- "CPD"
str(Exposures_H)

################################################################################
##### Extract outcome data for MVMR #####
################################################################################

outcome_dat_never <- read_outcome_data(
  "dep_never_imputed.txt",
  snps = Exposures_H$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  min_pval = 1e-200,
  log_pval = FALSE
)


outcome_dat_ever <- read_outcome_data(
  "dep_ever_imputed.txt",
  snps = Exposures_H$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  min_pval = 1e-200,
  log_pval = FALSE
)


outcome_dat_current <- read_outcome_data(
  "dep_current_imputed.txt",
  snps = Exposures_H$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  min_pval = 1e-200,
  log_pval = FALSE
)


outcome_dat_former <- read_outcome_data(
  "dep_former_imputed.txt",
  snps = Exposures_H$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  min_pval = 1e-200,
  log_pval = FALSE
)

################################################################################
##### Convert odds ratios to log odds #####
################################################################################
outcome_dat_never$beta.outcome <- as.numeric(
  as.character(
    outcome_dat_never$beta.outcome
  )
)
outcome_dat_ever$beta.outcome <- as.numeric(
  as.character(
    outcome_dat_ever$beta.outcome
  )
)
outcome_dat_current$beta.outcome <- as.numeric(
  as.character(
    outcome_dat_current$beta.outcome
  )
)
outcome_dat_former$beta.outcome <- as.numeric(
  as.character(
    outcome_dat_former$beta.outcome
  )
)
outcome_dat_never["beta.outcome"] <- (outcome_dat_never$beta.outcome/(0.20349*((1-0.20349))))
outcome_dat_never["se.outcome"] <- (outcome_dat_never$se.outcome/(0.20349*((1-0.20349))))
outcome_dat_ever["beta.outcome"] <- (outcome_dat_ever$beta.outcome/(0.25484*((1-0.25484))))
outcome_dat_ever["se.outcome"] <- (outcome_dat_ever$se.outcome/(0.25484*((1-0.25484))))
outcome_dat_current["beta.outcome"] <- (outcome_dat_current$beta.outcome/(0.30688*((1-0.30688))))
outcome_dat_current["se.outcome"] <- (outcome_dat_current$se.outcome/(0.30688*((1-0.30688))))
outcome_dat_former["beta.outcome"] <- (outcome_dat_former$beta.outcome/(0.23984*((1-0.23984))))
outcome_dat_former["se.outcome"] <- (outcome_dat_former$se.outcome/(0.23984*((1-0.23984))))

################################################################################
##### Organise outcome #####
################################################################################
outcome_dat_never["Phenotype"] <- NA
outcome_dat_never$Phenotype <- "Depression"
outcome_dat_ever["Phenotype"] <- NA
outcome_dat_ever$Phenotype <- "Depression"
outcome_dat_current["Phenotype"] <- NA
outcome_dat_current$Phenotype <- "Depression"
outcome_dat_former["Phenotype"] <- NA
outcome_dat_former$Phenotype <- "Depression"

################################################################################
##### Harmonise with outcome #####
################################################################################
mvdat_never <- mv_harmonise_data(
  Exposures_H,
  outcome_dat_never,
  harmonise_strictness = 2
)
mvdat_ever <- mv_harmonise_data(
  Exposures_H,
  outcome_dat_ever,
  harmonise_strictness = 2
)
mvdat_current <- mv_harmonise_data(
  Exposures_H,
  outcome_dat_current,
  harmonise_strictness = 2
)
mvdat_former <- mv_harmonise_data(
  Exposures_H,
  outcome_dat_former,
  harmonise_strictness = 2
)

mvdat_never1 <- harmonise_data(Exposures_H, outcome_dat_never)
mvdat_never1 <- mvdat_never1[mvdat_never1$mr_keep == TRUE, ]
str(mvdat_never1)

mvdat_ever1 <- harmonise_data(Exposures_H, outcome_dat_ever)
mvdat_ever1 <- mvdat_ever1[mvdat_ever1$mr_keep == TRUE, ]
str(mvdat_ever1)

mvdat_current1 <- harmonise_data(Exposures_H, outcome_dat_current)
mvdat_current1 <- mvdat_current1[mvdat_current1$mr_keep == TRUE, ]
str(mvdat_current1)

mvdat_former1 <- harmonise_data(Exposures_H, outcome_dat_former)
mvdat_former1 <- mvdat_former1[mvdat_former1$mr_keep == TRUE, ]
str(mvdat_former1)

################################################################################
##### Find proxies to add to missing outcome SNPs #####
################################################################################
proxy_needed4 <- data.frame(setdiff(NMR_H$SNP, outcome_dat_never$SNP)) #0 missing SNPs
proxy_needed5 <- data.frame(setdiff(NMR_H$SNP, outcome_dat_ever$SNP)) #0 missing SNPs
proxy_needed6 <- data.frame(setdiff(NMR_H$SNP, outcome_dat_current$SNP)) #0 missing SNPs
proxy_needed7 <- data.frame(setdiff(NMR_H$SNP, outcome_dat_former$SNP)) #0 missing SNPs
rm(proxy_needed4, proxy_needed5, proxy_needed6, proxy_needed7)

###########################
##### Save dataframes #####
###########################

write.csv(
  mvdat_never1,
  "mvdat_never.csv",
  row.names = FALSE
)
write.csv(
  mvdat_ever1,
  "mvdat_ever.csv",
  row.names = FALSE
)

write.csv(
  mvdat_current1,
  "mvdat_current.csv",
  row.names = FALSE
)

write.csv(
  mvdat_former1,
  "mvdat_former.csv",
  row.names = FALSE
)

#Read in if running to test/reproduce plots
mvdat_never1 <- read.csv("mvdat_never.csv")
mvdat_ever1 <- read.csv("mvdat_ever.csv")
mvdat_current1 <- read.csv("mvdat_current.csv")
mvdat_former1 <- read.csv("mvdat_former.csv")

################################################################################
##### Run MVMR #####
################################################################################
##### IVW never #####
bX1 <- c(mvdat_never1$beta.exposure[mvdat_never1$id.exposure == 1])
bX2 <- c(mvdat_never1$beta.exposure[mvdat_never1$id.exposure == 2])
bY <- c(mvdat_never1$beta.outcome[mvdat_never1$id.exposure == 1])
bYse <- c(mvdat_never1$se.outcome[mvdat_never1$id.exposure == 1])

set.seed(1234)
mod.MVMR_never <- lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2)
se_theta1MI.random <- summary(lm(
  bY ~ bX1 + bX2 - 1,
  weights = bYse^-2
))$coef[1, 2] /
  min(summary(lm(bY ~ bX1 + bX2 - 1,
                 weights = bYse^-2
  ))$sigma, 1)

mod_n <- summary(mod.MVMR_never)

mod_n_or <- coef(summary(mod.MVMR_never))
colnames(mod_n_or) <- c("b", "se", "t", "p")
mod_n_or <- as.data.frame(mod_n_or)
mod_n_or <- generate_odds_ratios(mod_n_or)

##### Orientation NMR #####
##### As Egger analyses require the exposure betas to be positive,        #####
##### we first orient the betas to be positive for NMR, and then          #####
##### orient the betas to be positive for CPD. In the paper, we           #####
##### report the result for each exposure only with the right orientation #####

clist <- c("bX2", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX1>0,", var, ",", var, "*-1)")))
}
bX1 <- abs(bX1)

##### MVMR Egger #####
mod.MVMRME_never <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(
  bY ~ bX1 + bX2,
  weights = bYse^-2
))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_n <- summary(mod.MVMRME_never)
mod_ME_n_or <- data.frame(mod.MVMRME_never[["coefficients"]])
colnames(mod_ME_n_or) <- c("b", "se", "t", "p")
mod_ME_n_or <- as.data.frame(mod_ME_n_or)
mod_ME_n_or <- generate_odds_ratios(mod_ME_n_or)

##### Orientation CPD #####
clist <- c("bX1", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX2>0,", var, ",", var, "*-1)")))
}
bX2 <- abs(bX2)

##### MVMR Egger #####
mod.MVMRME_never_2 <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_n_2 <- summary(mod.MVMRME_never_2)

mod_ME_n_2_or <- data.frame(mod.MVMRME_never_2[["coefficients"]])
colnames(mod_ME_n_2_or) <- c("b", "se", "t", "p")
mod_ME_n_2_or <- as.data.frame(mod_ME_n_2_or)
mod_ME_n_2_or <- generate_odds_ratios(mod_ME_n_2_or)

##### IVW ever #####
bX1 <- c(mvdat_ever1$beta.exposure[mvdat_ever1$id.exposure == 1])
bX2 <- c(mvdat_ever1$beta.exposure[mvdat_ever1$id.exposure == 2])
bY <- c(mvdat_ever1$beta.outcome[mvdat_ever1$id.exposure == 1])
bYse <- c(mvdat_ever1$se.outcome[mvdat_ever1$id.exposure == 1])

set.seed(1234)
mod.MVMR_ever <- lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2)
se_theta1MI.random <- summary(lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2))$coef[1, 2] /
  min(summary(lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2))$sigma, 1)

mod_e <- summary(mod.MVMR_ever)

mod_e_or <- coef(summary(mod.MVMR_ever))
colnames(mod_e_or) <- c("b", "se", "t", "p")
mod_e_or <- as.data.frame(mod_e_or)
mod_e_or <- generate_odds_ratios(mod_e_or)

##### Orientation NMR #####
clist <- c("bX2", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX1>0,", var, ",", var, "*-1)")))
}
bX1 <- abs(bX1)

##### MVMR Egger #####
mod.MVMRME_ever <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_e <- summary(mod.MVMRME_ever)

mod_ME_e_or <- data.frame(mod.MVMRME_ever[["coefficients"]])
colnames(mod_ME_e_or) <- c("b", "se", "t", "p")
mod_ME_e_or <- as.data.frame(mod_ME_e_or)
mod_ME_e_or <- generate_odds_ratios(mod_ME_e_or)

##### Orientation CPD #####
clist <- c("bX1", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX2>0,", var, ",", var, "*-1)")))
}
bX2 <- abs(bX2)

##### MVMR Egger #####
mod.MVMRME_ever_2 <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm
                              (bY ~ bX1 + bX2, weights = bYse^-2))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_e_2 <- summary(mod.MVMRME_ever_2)


mod_ME_e_2_or <- data.frame(mod.MVMRME_ever_2[["coefficients"]])
colnames(mod_ME_e_2_or) <- c("b", "se", "t", "p")
mod_ME_e_2_or <- as.data.frame(mod_ME_e_2_or)
mod_ME_e_2_or <- generate_odds_ratios(mod_ME_e_2_or)


##### IVW current #####
bX1 <- c(mvdat_current1$beta.exposure[mvdat_current1$id.exposure == 1])
bX2 <- c(mvdat_current1$beta.exposure[mvdat_current1$id.exposure == 2])
bY <- c(mvdat_current1$beta.outcome[mvdat_current1$id.exposure == 1])
bYse <- c(mvdat_current1$se.outcome[mvdat_current1$id.exposure == 1])

set.seed(1234)
mod.MVMR_current <- lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2)
se_theta1MI.random <- summary(lm(
  bY ~ bX1 + bX2 - 1,
  weights = bYse^-2
))$coef[1, 2] /
  min(summary(lm(bY ~ bX1 + bX2 - 1,
                 weights = bYse^-2
  ))$sigma, 1)

mod_c <- summary(mod.MVMR_current)

mod_c_or <- coef(summary(mod.MVMR_current))
colnames(mod_c_or) <- c("b", "se", "t", "p")
mod_c_or <- as.data.frame(mod_c_or)
mod_c_or <- generate_odds_ratios(mod_c_or)

##### Orientation NMR #####
##### As Egger analyses require the exposure betas to be positive,        #####
##### we first orient the betas to be positive for NMR, and then          #####
##### orient the betas to be positive for CPD. In the paper, we           #####
##### report the result for each exposure only with the right orientation #####

clist <- c("bX2", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX1>0,", var, ",", var, "*-1)")))
}
bX1 <- abs(bX1)

##### MVMR Egger #####
mod.MVMRME_current <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(
  bY ~ bX1 + bX2,
  weights = bYse^-2
))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_c <- summary(mod.MVMRME_current)
mod_ME_c_or <- data.frame(mod.MVMRME_current[["coefficients"]])
colnames(mod_ME_c_or) <- c("b", "se", "t", "p")
mod_ME_c_or <- as.data.frame(mod_ME_c_or)
mod_ME_c_or <- generate_odds_ratios(mod_ME_c_or)

##### Orientation CPD #####
clist <- c("bX1", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX2>0,", var, ",", var, "*-1)")))
}
bX2 <- abs(bX2)

##### MVMR Egger #####
mod.MVMRME_current_2 <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_c_2 <- summary(mod.MVMRME_current_2)

mod_ME_c_2_or <- data.frame(mod.MVMRME_current_2[["coefficients"]])
colnames(mod_ME_c_2_or) <- c("b", "se", "t", "p")
mod_ME_c_2_or <- as.data.frame(mod_ME_c_2_or)
mod_ME_c_2_or <- generate_odds_ratios(mod_ME_c_2_or)

##### IVW former #####
bX1 <- c(mvdat_former1$beta.exposure[mvdat_former1$id.exposure == 1])
bX2 <- c(mvdat_former1$beta.exposure[mvdat_former1$id.exposure == 2])
bY <- c(mvdat_former1$beta.outcome[mvdat_former1$id.exposure == 1])
bYse <- c(mvdat_former1$se.outcome[mvdat_former1$id.exposure == 1])

set.seed(1234)
mod.MVMR_former <- lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2)
se_theta1MI.random <- summary(lm(
  bY ~ bX1 + bX2 - 1,
  weights = bYse^-2
))$coef[1, 2] /
  min(summary(lm(bY ~ bX1 + bX2 - 1,
                 weights = bYse^-2
  ))$sigma, 1)

mod_f <- summary(mod.MVMR_former)

mod_f_or <- coef(summary(mod.MVMR_former))
colnames(mod_f_or) <- c("b", "se", "t", "p")
mod_f_or <- as.data.frame(mod_f_or)
mod_f_or <- generate_odds_ratios(mod_f_or)

##### Orientation NMR #####
##### As Egger analyses require the exposure betas to be positive,        #####
##### we first orient the betas to be positive for NMR, and then          #####
##### orient the betas to be positive for CPD. In the paper, we           #####
##### report the result for each exposure only with the right orientation #####

clist <- c("bX2", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX1>0,", var, ",", var, "*-1)")))
}
bX1 <- abs(bX1)

##### MVMR Egger #####
mod.MVMRME_former <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(
  bY ~ bX1 + bX2,
  weights = bYse^-2
))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_f <- summary(mod.MVMRME_former)
mod_ME_f_or <- data.frame(mod.MVMRME_former[["coefficients"]])
colnames(mod_ME_f_or) <- c("b", "se", "t", "p")
mod_ME_f_or <- as.data.frame(mod_ME_f_or)
mod_ME_f_or <- generate_odds_ratios(mod_ME_f_or)

##### Orientation CPD #####
clist <- c("bX1", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX2>0,", var, ",", var, "*-1)")))
}
bX2 <- abs(bX2)

##### MVMR Egger #####
mod.MVMRME_former_2 <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_f_2 <- summary(mod.MVMRME_former_2)

mod_ME_f_2_or <- data.frame(mod.MVMRME_former_2[["coefficients"]])
colnames(mod_ME_f_2_or) <- c("b", "se", "t", "p")
mod_ME_f_2_or <- as.data.frame(mod_ME_f_2_or)
mod_ME_f_2_or <- generate_odds_ratios(mod_ME_f_2_or)


##### Format to analyse / cross check with MVMR package #####
bX1 <- c(mvdat_never1$beta.exposure[mvdat_never1$id.exposure == 1])
bX2 <- c(mvdat_never1$beta.exposure[mvdat_never1$id.exposure == 2])
bY <- c(mvdat_never1$beta.outcome[mvdat_never1$id.exposure == 1])
bYse <- c(mvdat_never1$se.outcome[mvdat_never1$id.exposure == 1])
bXse1 <- c(mvdat_never1$se.exposure[mvdat_never1$id.exposure == 1])
bXse2 <- c(mvdat_never1$se.exposure[mvdat_never1$id.exposure == 2])
df_never <- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)

bX1 <- c(mvdat_ever1$beta.exposure[mvdat_ever1$id.exposure == 1])
bX2 <- c(mvdat_ever1$beta.exposure[mvdat_ever1$id.exposure == 2])
bY <- c(mvdat_ever1$beta.outcome[mvdat_ever1$id.exposure == 1])
bYse <- c(mvdat_ever1$se.outcome[mvdat_ever1$id.exposure == 1])
bXse1 <- c(mvdat_ever1$se.exposure[mvdat_ever1$id.exposure == 1])
bXse2 <- c(mvdat_ever1$se.exposure[mvdat_ever1$id.exposure == 2])
df_ever <- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)

bX1 <- c(mvdat_current1$beta.exposure[mvdat_current1$id.exposure == 1])
bX2 <- c(mvdat_current1$beta.exposure[mvdat_current1$id.exposure == 2])
bY <- c(mvdat_current1$beta.outcome[mvdat_current1$id.exposure == 1])
bYse <- c(mvdat_current1$se.outcome[mvdat_current1$id.exposure == 1])
bXse1 <- c(mvdat_current1$se.exposure[mvdat_current1$id.exposure == 1])
bXse2 <- c(mvdat_current1$se.exposure[mvdat_current1$id.exposure == 2])
df_current <- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)

bX1 <- c(mvdat_former1$beta.exposure[mvdat_former1$id.exposure == 1])
bX2 <- c(mvdat_former1$beta.exposure[mvdat_former1$id.exposure == 2])
bY <- c(mvdat_former1$beta.outcome[mvdat_former1$id.exposure == 1])
bYse <- c(mvdat_former1$se.outcome[mvdat_former1$id.exposure == 1])
bXse1 <- c(mvdat_former1$se.exposure[mvdat_former1$id.exposure == 1])
bXse2 <- c(mvdat_former1$se.exposure[mvdat_former1$id.exposure == 2])
df_former <- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)

df_mvmr_never <- format_mvmr(
  df_never[, c(1, 3)],
  df_never[, 5],
  df_never[, c(2, 4)],
  df_never[, 6]
)
df_mvmr_ever <- format_mvmr(
  df_ever[, c(1, 3)],
  df_ever[, 5],
  df_ever[, c(2, 4)],
  df_ever[, 6]
)
df_mvmr_current <- format_mvmr(
  df_current[, c(1, 3)],
  df_current[, 5],
  df_current[, c(2, 4)],
  df_current[, 6]
)
df_mvmr_former <- format_mvmr(
  df_former[, c(1, 3)],
  df_former[, 5],
  df_former[, c(2, 4)],
  df_former[, 6]
)
##### Cross check IVW result with MVMR package #####
res_n <- ivw_mvmr(df_mvmr_never)
res_e <- ivw_mvmr(df_mvmr_ever)
res_c <- ivw_mvmr(df_mvmr_current)
res_f <- ivw_mvmr(df_mvmr_former)

################################################################################
##### Calculate F-statistic and covariance #####
# Note: >10 is strong
# correlation between NMR and CPD in Buchwald  = -0.019
# qhet_mvmr is used to adjust for covariance
# At present, CIs take substantial time to calculate - only generate if point estimates are substantially different or evidence MVMR is underpowered
# Compare effects with and without adjustment
# There is a small amount of sample overlap between exposure instruments (e.g., FinnTwin in both GSCAN and Buchwald et al., 2021)
# So setting covariance at 0 would be inappropriate
# https://wspiller.github.io/MVMR/articles/MVMR.html 
################################################################################

cov <- matrix(c(1, -0.019, -0.019, 1), nrow = 2, ncol = 2)

Xcovmat_n <- phenocov_mvmr(cov, df_mvmr_never[, c(6, 7)])
Fstat_n <- strength_mvmr(df_mvmr_never, gencov = Xcovmat_n)

Xcovmat_e <- phenocov_mvmr(cov, df_mvmr_ever[, c(6, 7)])
Fstat_e <- strength_mvmr(df_mvmr_ever, gencov = Xcovmat_e)

Xcovmat_c <- phenocov_mvmr(cov, df_mvmr_current[, c(6, 7)])
Fstat_c <- strength_mvmr(df_mvmr_current, gencov = Xcovmat_c)

Xcovmat_f <- phenocov_mvmr(cov, df_mvmr_former[, c(6, 7)])
Fstat_f <- strength_mvmr(df_mvmr_former, gencov = Xcovmat_f)

#Point estimates are all very similar
cov_adj_mvmr_n <- qhet_mvmr(df_mvmr_never, cov, CI = FALSE, iterations = 1000)
cov_adj_mvmr_e <- qhet_mvmr(df_mvmr_ever, cov, CI = FALSE, iterations = 1000)
cov_adj_mvmr_c <- qhet_mvmr(df_mvmr_current, cov, CI = FALSE, iterations = 1000)
cov_adj_mvmr_f <- qhet_mvmr(df_mvmr_former, cov, CI = FALSE, iterations = 1000)

################################################################################
##### Test for horizontal pleiotropy #####
# Note: Q should be *less* than the number of SNPs included
################################################################################

ptr_n <- pleiotropy_mvmr(df_mvmr_never, gencov = Xcovmat_n)
ptr_e <- pleiotropy_mvmr(df_mvmr_ever, gencov = Xcovmat_e)
ptr_c <- pleiotropy_mvmr(df_mvmr_current, gencov = Xcovmat_c)
ptr_f <- pleiotropy_mvmr(df_mvmr_former, gencov = Xcovmat_f)

################################################################################
##### Forest Plots and Tables #####
################################################################################

###################
#####  Tables #####
###################
# First generating table of main results (i.e., IVW and Egger for MVMR and supporting stats like I2GX)

Exposure <- c(
  "NMR", "NMR", "NMR", "NMR",
  "CPD", "CPD", "CPD", "CPD",
  "NMR", "NMR", "NMR", "NMR",
  "CPD", "CPD", "CPD", "CPD",
  "NMR", "NMR", "NMR", "NMR",
  "CPD", "CPD", "CPD", "CPD",
  "NMR", "NMR", "NMR", "NMR",
  "CPD", "CPD", "CPD", "CPD"
)
Smoking <- c(
  "Never", "Never", "Never", "Never", "Never", "Never", "Never", "Never",
  "Ever", "Ever", "Ever", "Ever", "Ever", "Ever", "Ever", "Ever",
  "Current", "Current", "Current", "Current", "Current", "Current", "Current", "Current",
  "Former", "Former", "Former", "Former", "Former", "Former", "Former", "Former"
)

# Converting to data frames
res_n <- data.frame(res_n)
res_e <- data.frame(res_e)
res_c <- data.frame(res_c)
res_f <- data.frame(res_f)

Method <- c(
  "MR-IVW", "MVMR-IVW", "MR-Egger", "MVMR-Egger",
  "MR-IVW", "MVMR-IVW", "MR-Egger", "MVMR-Egger",
  "MR-IVW", "MVMR-IVW", "MR-Egger", "MVMR-Egger",
  "MR-IVW", "MVMR-IVW", "MR-Egger", "MVMR-Egger",
  "MR-IVW", "MVMR-IVW", "MR-Egger", "MVMR-Egger",
  "MR-IVW", "MVMR-IVW", "MR-Egger", "MVMR-Egger",
  "MR-IVW", "MVMR-IVW", "MR-Egger", "MVMR-Egger",
  "MR-IVW", "MVMR-IVW", "MR-Egger", "MVMR-Egger"
)
OR <- c(
  result_nmr_never[3, "or"], mod_n_or[1, 7],
  result_nmr_never[1, "or"], mod_ME_n_or[2, 7],
  result_cpd_never[3, "or"], mod_n_or[2, 7],
  result_cpd_never[1, "or"], mod_ME_n_2_or[3, 7],
  result_nmr_ever[3, "or"], mod_e_or[1, 7],
  result_nmr_ever[1, "or"], mod_ME_e_or[2, 7],
  result_cpd_ever[3, "or"], mod_e_or[2, 7],
  result_cpd_ever[1, "or"], mod_ME_e_2_or[3, 7],
  result_nmr_current[3, "or"], mod_c_or[1, 7],
  result_nmr_current[1, "or"], mod_ME_c_or[2, 7],
  result_cpd_current[3, "or"], mod_c_or[2, 7],
  result_cpd_current[1, "or"], mod_ME_c_2_or[3, 7],
  result_nmr_former[3, "or"], mod_f_or[1, 7],
  result_nmr_former[1, "or"], mod_ME_f_or[2, 7],
  result_cpd_former[3, "or"], mod_f_or[2, 7],
  result_cpd_former[1, "or"], mod_ME_f_2_or[3, 7]
)

LCI <- c(
  result_nmr_never[3, "or_lci95"], mod_n_or[1, 8],
  result_nmr_never[1, "or_lci95"], mod_ME_n_or[2, 8],
  result_cpd_never[3, "or_lci95"], mod_n_or[2, 8],
  result_cpd_never[1, "or_lci95"], mod_ME_n_2_or[3, 8],
  result_nmr_ever[3, "or_lci95"], mod_e_or[1, 8],
  result_nmr_ever[1, "or_lci95"], mod_ME_e_or[2, 8],
  result_cpd_ever[3, "or_lci95"], mod_e_or[2, 8],
  result_cpd_ever[1, "or_lci95"], mod_ME_e_2_or[3, 8],
  result_nmr_current[3, "or_lci95"], mod_c_or[1, 8],
  result_nmr_current[1, "or_lci95"], mod_ME_c_or[2, 8],
  result_cpd_current[3, "or_lci95"], mod_c_or[2, 8],
  result_cpd_current[1, "or_lci95"], mod_ME_c_2_or[3, 8],
  result_nmr_former[3, "or_lci95"], mod_f_or[1, 8],
  result_nmr_former[1, "or_lci95"], mod_ME_f_or[2, 8],
  result_cpd_former[3, "or_lci95"], mod_f_or[2, 8],
  result_cpd_former[1, "or_lci95"], mod_ME_f_2_or[3, 8]
)

UCI <- c(
  result_nmr_never[3, "or_uci95"], mod_n_or[1, 9],
  result_nmr_never[1, "or_uci95"], mod_ME_n_or[2, 9],
  result_cpd_never[3, "or_uci95"], mod_n_or[2, 9],
  result_cpd_never[1, "or_uci95"], mod_ME_n_2_or[3, 9],
  result_nmr_ever[3, "or_uci95"], mod_e_or[1, 9],
  result_nmr_ever[1, "or_uci95"], mod_ME_e_or[2, 9],
  result_cpd_ever[3, "or_uci95"], mod_e_or[2, 9],
  result_cpd_ever[1, "or_uci95"], mod_ME_e_2_or[3, 9],
  result_nmr_current[3, "or_uci95"], mod_c_or[1, 9],
  result_nmr_current[1, "or_uci95"], mod_ME_c_or[2, 9],
  result_cpd_current[3, "or_uci95"], mod_c_or[2, 9],
  result_cpd_current[1, "or_uci95"], mod_ME_c_2_or[3, 9],
  result_nmr_former[3, "or_uci95"], mod_f_or[1, 9],
  result_nmr_former[1, "or_uci95"], mod_ME_f_or[2, 9],
  result_cpd_former[3, "or_uci95"], mod_f_or[2, 9],
  result_cpd_former[1, "or_uci95"], mod_ME_f_2_or[3, 9]
)

p <- c(
  result_nmr_never[3, "pval"], mod_n_or[1, 4],
  result_nmr_never[1, "pval"], mod_ME_n_or[2, 4],
  result_cpd_never[3, "pval"], mod_n_or[2, 4],
  result_cpd_never[1, "pval"], mod_ME_n_2_or[3, 4],
  result_nmr_ever[3, "pval"], mod_e_or[1, 4],
  result_nmr_ever[1, "pval"], mod_ME_e_or[2, 4],
  result_cpd_ever[3, "pval"], mod_e_or[2, 4],
  result_cpd_ever[1, "pval"], mod_ME_e_2_or[3, 4],
  result_nmr_current[3, "pval"], mod_c_or[1, 4],
  result_nmr_current[1, "pval"], mod_ME_c_or[2, 4],
  result_cpd_current[3, "pval"], mod_c_or[2, 4],
  result_cpd_current[1, "pval"], mod_ME_c_2_or[3, 4],
  result_nmr_former[3, "pval"], mod_f_or[1, 4],
  result_nmr_former[1, "pval"], mod_ME_f_or[2, 4],
  result_cpd_former[3, "pval"], mod_f_or[2, 4],
  result_cpd_former[1, "pval"], mod_ME_f_2_or[3, 4]
)
I2Gx <- c(
  ".", ".", ISQ[1, 1], ".",
  ".", ".", ISQ[5, 1], ".",
  ".", ".", ISQ[2, 1], ".",
  ".", ".", ISQ[6, 1], ".",
  ".", ".", ISQ[3, 1], ".",
  ".", ".", ISQ[7, 1], ".",
  ".", ".", ISQ[4, 1], ".",
  ".", ".", ISQ[8, 1], "."
)
Q <- c(
  ptr_n_nmr[1, 1], ptr_n[["Qstat"]], ptr_n_nmr[2, 1], ptr_n[["Qstat"]],
  ptr_n_cpd[1, 1], ptr_n[["Qstat"]], ptr_n_cpd[2, 1], ptr_n[["Qstat"]],
  ptr_e_nmr[1, 1], ptr_e[["Qstat"]], ptr_e_nmr[2, 1], ptr_e[["Qstat"]],
  ptr_e_cpd[1, 1], ptr_e[["Qstat"]], ptr_e_cpd[2, 1], ptr_e[["Qstat"]],
  ptr_c_nmr[1, 1], ptr_c[["Qstat"]], ptr_c_nmr[2, 1], ptr_c[["Qstat"]],
  ptr_c_cpd[1, 1], ptr_c[["Qstat"]], ptr_c_cpd[2, 1], ptr_c[["Qstat"]],
  ptr_f_nmr[1, 1], ptr_f[["Qstat"]], ptr_f_nmr[2, 1], ptr_f[["Qstat"]],
  ptr_f_cpd[1, 1], ptr_f[["Qstat"]], ptr_f_cpd[2, 1], ptr_f[["Qstat"]]
)
EggerI <- c(
  ".", ".", egger_n_nmr[["b_i"]], mod_ME_n_or[1, 7],
  ".", ".", egger_n_cpd[["b_i"]], mod_ME_n_2_or[1, 7],
  ".", ".", egger_e_nmr[["b_i"]], mod_ME_e_or[1, 7],
  ".", ".", egger_e_cpd[["b_i"]], mod_ME_e_2_or[1, 7],
  ".", ".", egger_c_nmr[["b_i"]], mod_ME_c_or[1, 7],
  ".", ".", egger_c_cpd[["b_i"]], mod_ME_c_2_or[1, 7],
  ".", ".", egger_f_nmr[["b_i"]], mod_ME_f_or[1, 7],
  ".", ".", egger_f_cpd[["b_i"]], mod_ME_f_2_or[1, 7]
)
EggerIp <- c(
  ".", ".", egger_n_nmr[["pval_i"]], mod_ME_n_or[1, 4],
  ".", ".", egger_n_cpd[["pval_i"]], mod_ME_n_2_or[1, 4],
  ".", ".", egger_e_nmr[["pval_i"]], mod_ME_e_or[1, 4],
  ".", ".", egger_e_cpd[["pval_i"]], mod_ME_e_2_or[1, 4],
  ".", ".", egger_c_nmr[["pval_i"]], mod_ME_c_or[1, 4],
  ".", ".", egger_c_cpd[["pval_i"]], mod_ME_c_2_or[1, 4],
  ".", ".", egger_f_nmr[["pval_i"]], mod_ME_f_or[1, 4],
  ".", ".", egger_f_cpd[["pval_i"]], mod_ME_f_2_or[1, 4]
)
F_stat <- c(
  ptr_n_nmr[3, 1], Fstat_n[1, 1], ptr_n_nmr[3, 1], Fstat_n[1, 1],
  ptr_n_cpd[3, 1], Fstat_n[1, 2], ptr_n_cpd[3, 1], Fstat_n[1, 2],
  ptr_e_nmr[3, 1], Fstat_e[1, 1], ptr_e_nmr[3, 1], Fstat_e[1, 1],
  ptr_e_cpd[3, 1], Fstat_e[1, 2], ptr_e_cpd[3, 1], Fstat_e[1, 2],
  ptr_c_nmr[3, 1], Fstat_c[1, 1], ptr_c_nmr[3, 1], Fstat_c[1, 1],
  ptr_c_cpd[3, 1], Fstat_c[1, 2], ptr_c_cpd[3, 1], Fstat_c[1, 2],
  ptr_f_nmr[3, 1], Fstat_f[1, 1], ptr_f_nmr[3, 1], Fstat_f[1, 1],
  ptr_f_cpd[3, 1], Fstat_f[1, 2], ptr_f_cpd[3, 1], Fstat_f[1, 2]
)

all_results <- data.frame(
  Exposure, Smoking, Method,
  OR, LCI, UCI, p,
  I2Gx, Q, EggerI, EggerIp, F_stat
)


write.csv(all_results,
          "MVMR_Results_NMR_CPD_MDD.csv",
          row.names = FALSE
)

# Then generating table of sensitivity analyses (i.e., weighted median and weighted mode)

Exposure <- c(
  "NMR", "NMR",
  "CPD", "CPD", 
  "NMR", "NMR", 
  "CPD", "CPD", 
  "NMR", "NMR", 
  "CPD", "CPD", 
  "NMR", "NMR", 
  "CPD", "CPD" 
)
Smoking <- c(
  "Never", "Never", "Never", "Never", 
  "Ever", "Ever", "Ever", "Ever", 
  "Current", "Current", "Current", "Current", 
  "Former", "Former", "Former", "Former" 
)

Method <- c(
  "Weighted median", "Weighted mode", 
  "Weighted median", "Weighted mode",
  "Weighted median", "Weighted mode",
  "Weighted median", "Weighted mode",
  "Weighted median", "Weighted mode",
  "Weighted median", "Weighted mode",
  "Weighted median", "Weighted mode",
  "Weighted median", "Weighted mode"
)

OR <- c(
  result_nmr_never[2, "or"], result_nmr_never[5, "or"],
  result_cpd_never[2, "or"], result_cpd_never[5, "or"],
  result_nmr_ever[2, "or"], result_nmr_ever[5, "or"],
  result_cpd_ever[2, "or"], result_cpd_ever[5, "or"],
  result_nmr_current[2, "or"],  result_nmr_current[5, "or"],
  result_cpd_current[2, "or"],   result_cpd_current[5, "or"],
  result_nmr_former[2, "or"], result_nmr_former[5, "or"],
  result_cpd_former[2, "or"], result_cpd_former[5, "or"]
)

LCI <- c(
  result_nmr_never[2, "or_lci95"], result_nmr_never[5, "or_lci95"],
  result_cpd_never[2, "or_lci95"], result_cpd_never[5, "or_lci95"],
  result_nmr_ever[2, "or_lci95"], result_nmr_ever[5, "or_lci95"],
  result_cpd_ever[2, "or_lci95"], result_cpd_ever[5, "or_lci95"],
  result_nmr_current[2, "or_lci95"],  result_nmr_current[5, "or_lci95"],
  result_cpd_current[2, "or_lci95"],   result_cpd_current[5, "or_lci95"],
  result_nmr_former[2, "or_lci95"], result_nmr_former[5, "or_lci95"],
  result_cpd_former[2, "or_lci95"], result_cpd_former[5, "or_lci95"]
)

UCI <- c(
  result_nmr_never[2, "or_uci95"], result_nmr_never[5, "or_uci95"],
  result_cpd_never[2, "or_uci95"], result_cpd_never[5, "or_uci95"],
  result_nmr_ever[2, "or_uci95"], result_nmr_ever[5, "or_uci95"],
  result_cpd_ever[2, "or_uci95"], result_cpd_ever[5, "or_uci95"],
  result_nmr_current[2, "or_uci95"],  result_nmr_current[5, "or_uci95"],
  result_cpd_current[2, "or_uci95"],   result_cpd_current[5, "or_uci95"],
  result_nmr_former[2, "or_uci95"], result_nmr_former[5, "or_uci95"],
  result_cpd_former[2, "or_uci95"], result_cpd_former[5, "or_uci95"]
)

p <- c(
  result_nmr_never[2, "pval"], result_nmr_never[5, "pval"],
  result_cpd_never[2, "pval"], result_cpd_never[5, "pval"],
  result_nmr_ever[2, "pval"], result_nmr_ever[5, "pval"],
  result_cpd_ever[2, "pval"], result_cpd_ever[5, "pval"],
  result_nmr_current[2, "pval"],  result_nmr_current[5, "pval"],
  result_cpd_current[2, "pval"],   result_cpd_current[5, "pval"],
  result_nmr_former[2, "pval"], result_nmr_former[5, "pval"],
  result_cpd_former[2, "pval"], result_cpd_former[5, "pval"]
)

median_mode_results <- data.frame(
  Exposure, Smoking, Method,
  OR, LCI, UCI, p
)


write.csv(median_mode_results,
          "MedianMode_Results_NMR_CPD_MDD.csv",
          row.names = FALSE
)

###################
#####  Plots  #####
###################

##Ever smokers
# Create dataframe of the effect estimates, lower and upper confidence intervals

cochrane_from_rmeta_ever <-
  structure(
    list(
      mean  = c(NA, all_results$OR[9:16]),
      lower = c(NA, all_results$LCI[9:16]),
      upper = c(NA, all_results$UCI[9:16])
    ),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -9L),
    class = "data.frame"
  )

OR_LCI_UCI <- str_c(
  round(all_results$OR, digits = 2),
  " (",
  round(all_results$LCI, digits = 2),
  ", ",
  round(all_results$UCI, digits = 2),
  ")"
)

tabletext_ever <- cbind(
  c("Exposure", "Nicotine Metabolite Ratio", NA, NA, NA, "Smoking Heaviness", NA, NA, NA),
  c(
    "Method", "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger",
    "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger"
  ),
  c("OR (95% CI)", OR_LCI_UCI[9:16]),
  c("P value", round(all_results$p[9:16], digits = 3))
)
tabletext_ever[, 4][tabletext_ever[, 4] == 0] <- "<0.001"

##Create local function for different colours in forest plot
fn <- local({
  i = 0
  no_lines <- sum(!is.na(cochrane_from_rmeta_ever$mean))
  b_clrs = c("red1", "red4", "red1", "red4", "red1", "red4", "red1", "red4")
  l_clrs = c("red1", "red4", "red1", "red4", "red1", "red4", "red1", "red4")
  
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
  }
})

##### Create graphs and write to file #####
pdf.options(reset = TRUE, onefile = FALSE)
pdf("MVMR_ever.pdf", width = 10, height = 10)

# Create the forest plot with smaller line width
forestplot(tabletext_ever,
           fn.ci_norm = fn,
           cochrane_from_rmeta_ever,
           new_page = TRUE,
           is.summary = c(TRUE, rep(FALSE, 9)),
           lineheight = unit(1, "cm"),
           graphwidth = unit(10, "cm"),
           boxsize = 0.09,
           clip = c(0.1, 10),
           xticks = c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0),
           xlog = FALSE,
           zero = 1,
           ci.vertices = TRUE,
           colgap = unit(4, "mm"),
           cex = 0.2,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9)))

dev.off()

##Current smokers
# Create dataframe of the effect estimates, lower and upper confidence intervals
cochrane_from_rmeta_current <-
  structure(
    list(
      mean  = c(NA, all_results$OR[17:24]),
      lower = c(NA, all_results$LCI[17:24]),
      upper = c(NA, all_results$UCI[17:24])
    ),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -9L),
    class = "data.frame"
  )

OR_LCI_UCI <- str_c(
  round(all_results$OR, digits = 2),
  " (",
  round(all_results$LCI, digits = 2),
  ", ",
  round(all_results$UCI, digits = 2),
  ")"
)

tabletext_current <- cbind(
  c("Exposure", "Nicotine Metabolite Ratio", NA, NA, NA, "Smoking Heaviness", NA, NA, NA),
  c(
    "Method", "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger",
    "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger"
  ),
  c("OR (95% CI)", OR_LCI_UCI[17:24]),
  c("P value", round(all_results$p[17:24], digits = 3))
)
tabletext_current[, 4][tabletext_current[, 4] == 0] <- "<0.001"

##Create local function for different colours in forest plot
fn <- local({
  i = 0
  no_lines <- sum(!is.na(cochrane_from_rmeta_current$mean))
  b_clrs = c("lightgreen", "darkgreen", "lightgreen", "darkgreen", "lightgreen", "darkgreen", "lightgreen", "darkgreen")
  l_clrs = c("lightgreen", "darkgreen", "lightgreen", "darkgreen", "lightgreen", "darkgreen", "lightgreen", "darkgreen")
  
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
  }
})

##### Create graphs and write to file #####
pdf.options(reset = TRUE, onefile = FALSE)
pdf("MVMR_current.pdf", width = 10, height = 10)

# Create the forest plot with smaller line width
forestplot(tabletext_current,
           fn.ci_norm = fn,
           cochrane_from_rmeta_current,
           new_page = TRUE,
           is.summary = c(TRUE, rep(FALSE, 9)),
           lineheight = unit(1, "cm"),
           graphwidth = unit(10, "cm"),
           boxsize = 0.09,
           clip = c(0.1, 10),
           xticks = c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0),
           xlog = FALSE,
           zero = 1,
           ci.vertices = TRUE,
           colgap = unit(4, "mm"),
           cex = 0.2,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9)))

dev.off()

##Former smokers
# Create dataframe of the effect estimates, lower and upper confidence intervals
cochrane_from_rmeta_former <-
  structure(
    list(
      mean  = c(NA, all_results$OR[25:32]),
      lower = c(NA, all_results$LCI[25:32]),
      upper = c(NA, all_results$UCI[25:32])
    ),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -9L),
    class = "data.frame"
  )

OR_LCI_UCI <- str_c(
  round(all_results$OR, digits = 2),
  " (",
  round(all_results$LCI, digits = 2),
  ", ",
  round(all_results$UCI, digits = 2),
  ")"
)

tabletext_former <- cbind(
  c("Exposure", "Nicotine Metabolite Ratio", NA, NA, NA, "Smoking Heaviness", NA, NA, NA),
  c(
    "Method", "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger",
    "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger"
  ),
  c("OR (95% CI)", OR_LCI_UCI[25:32]),
  c("P value", round(all_results$p[25:32], digits = 3))
)
tabletext_former[, 4][tabletext_former[, 4] == 0] <- "<0.001"

##Create local function for different colours in forest plot
fn <- local({
  i = 0
  no_lines <- sum(!is.na(cochrane_from_rmeta_former$mean))
  b_clrs = c("skyblue2", "darkblue", "skyblue2", "darkblue", "skyblue2", "darkblue", "skyblue2", "darkblue")
  l_clrs = c("skyblue2", "darkblue", "skyblue2", "darkblue", "skyblue2", "darkblue", "skyblue2", "darkblue")
  
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
  }
})

##### Create graphs and write to file #####
pdf.options(reset = TRUE, onefile = FALSE)
pdf("MVMR_former.pdf", width = 10, height = 10)

# Create the forest plot with smaller line width
forestplot(tabletext_former,
           fn.ci_norm = fn,
           cochrane_from_rmeta_former,
           new_page = TRUE,
           is.summary = c(TRUE, rep(FALSE, 9)),
           lineheight = unit(1, "cm"),
           graphwidth = unit(10, "cm"),
           boxsize = 0.09,
           clip = c(0.1, 10),
           xticks = c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0),
           xlog = FALSE,
           zero = 1,
           ci.vertices = TRUE,
           colgap = unit(4, "mm"),
           cex = 0.2,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9)))

dev.off()

##Never smokers
# Create dataframe of the effect estimates, lower and upper confidence intervals
cochrane_from_rmeta_never <-
  structure(
    list(
      mean  = c(NA, all_results$OR[1:8]),
      lower = c(NA, all_results$LCI[1:8]),
      upper = c(NA, all_results$UCI[1:8])
    ),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -9L),
    class = "data.frame"
  )

OR_LCI_UCI <- str_c(
  round(all_results$OR, digits = 2),
  " (",
  round(all_results$LCI, digits = 2),
  ", ",
  round(all_results$UCI, digits = 2),
  ")"
)

tabletext_never <- cbind(
  c("Exposure", "Nicotine Metabolite Ratio", NA, NA, NA, "Smoking Heaviness", NA, NA, NA),
  c(
    "Method", "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger",
    "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger"
  ),
  c("OR (95% CI)", OR_LCI_UCI[1:8]),
  c("P value", round(all_results$p[1:8], digits = 3))
)
tabletext_never[, 4][tabletext_never[, 4] == 0] <- "<0.001"

##Create local function for different colours in forest plot
fn <- local({
  i = 0
  no_lines <- sum(!is.na(cochrane_from_rmeta_never$mean))
  b_clrs = c("grey", "grey4", "grey", "grey4", "grey", "grey4", "grey", "grey4")
  l_clrs = c("grey", "grey4", "grey", "grey4", "grey", "grey4", "grey", "grey4")
  
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
  }
})

##### Create graphs and write to file #####
pdf.options(reset = TRUE, onefile = FALSE)
pdf("MVMR_never.pdf", width = 10, height = 10)

# Create the forest plot with smaller line width
forestplot(tabletext_never,
           fn.ci_norm = fn,
           cochrane_from_rmeta_never,
           new_page = TRUE,
           is.summary = c(TRUE, rep(FALSE, 9)),
           lineheight = unit(1, "cm"),
           graphwidth = unit(10, "cm"),
           boxsize = 0.09,
           clip = c(0.1, 10),
           xticks = c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0),
           xlog = FALSE,
           zero = 1,
           ci.vertices = TRUE,
           colgap = unit(4, "mm"),
           cex = 0.2,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9)))

dev.off()

################################################################################
##### Sensitivity Plots #####
################################################################################
#https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html#plots 

#Renaming exposure and outcome for plots
outcome_mdd_dat_nmr_ever$exposure <- gsub("exposure", "NMR", outcome_mdd_dat_nmr_ever$exposure)
outcome_mdd_dat_nmr_ever$outcome<- gsub("outcome", "MDD", outcome_mdd_dat_nmr_ever$outcome)
outcome_mdd_dat_cpd_ever$exposure <- gsub("exposure", "CPD", outcome_mdd_dat_cpd_ever$exposure)
outcome_mdd_dat_cpd_ever$outcome<- gsub("outcome", "MDD", outcome_mdd_dat_cpd_ever$outcome)

##NMR
# leave one out analyses
nmr_clmp_loo <- mr_leaveoneout(outcome_mdd_dat_nmr_ever)
nmr_clmp_loo <- mr_leaveoneout_plot(nmr_clmp_loo)
ggsave(nmr_clmp_loo[[1]], file="loo_nmr_ever.png", width=7, height=4)

#scatter plot
plot_result_nmr <- subset(result_nmr_ever, method != "Simple mode")
scatter_nmr <- mr_scatter_plot(plot_result_nmr, outcome_mdd_dat_nmr_ever)
ggsave(scatter_nmr[[1]], file="scatter_nmr_ever.png", width=7, height=7)

##CPD
# leave one out analyses
cpd_clmp_loo <- mr_leaveoneout(outcome_mdd_dat_cpd_ever)
cpd_clmp_loo <- mr_leaveoneout_plot(cpd_clmp_loo)
ggsave(cpd_clmp_loo[[1]], file="loo_cpd_ever.png", width=7, height=10)

#scatter plot
plot_result_cpd <- subset(result_cpd_ever, method != "Simple mode")
scatter_cpd<- mr_scatter_plot(plot_result_cpd, outcome_mdd_dat_cpd_ever)
ggsave(scatter_cpd[[1]], file="scatter_cpd_ever.png", width=7, height=7)