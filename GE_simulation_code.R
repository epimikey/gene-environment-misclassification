###*****************************************************************
###
### R code for analysis in:
###
### "Weisskopf M, Leung M. Misclassification Bias in the Assessment 
###  of Gene-By-Environment Interactions. Epidemiology. 2023."
###
###*****************************************************************

###*************
### N: Notes ###
###*************

#' @param iter number of iterations
#' @param n sample size of full underlying cohort from which case-control study is sampled
#' @param pE prevalence of the exposure
#' @param pG prevalence of the genetic factor
#' @param pD prevalence of the disease in doubly exposed (pE=0, pG=0)
#' @param or1 odds ratio for gene-environment dependence
#' @param or2 odds ratio for environmental exposure main effect
#' @param or3 odds ratio for gene-environment interaction
#' @param n_case number of cases to sample from underlying cohort
#' @param n_control number of controls to sample from underlying cohort
#' @param Sp specificity of exposure recall among cases

###*************************
### 1: Required Packages ###
###*************************

# install required packages
requiredPackages = c('Rlab', 'dplyr','svMisc')
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE))
    install.packages(p)
  library(p, character.only = TRUE)
}

###***************************
### 2: Simulation Function ###
###***************************

# set seed
set.seed(666)

### Set up data-generating function
geINT <- function(iter, n, pE, pG, pD, or1, or2, or3, n_case, n_control, Sp)
{
  # create results matrix (bias, coverage, pvalue)
  results <- data.frame(bias=rep(NA,iter), cover=rep(NA,iter), prev=rep(NA,iter))
  
  # generate data
  for (i in 1:iter){
    df <- data.frame("id" = 1:n) %>%
      mutate(G = rbern(n, plogis(log(pG/(1-pG)))),
             E = rbern(n, plogis(log(pE/(1-pE))+ log(or1)*G)),
             D = rbern(n, plogis(log(pD/(1-pD)) + log(or2)*E + log(or3)*E*G)))
    
    # sample cases
    cc_cases <- df %>% filter(D == 1) %>% sample_n(n_case, replace=T)
    cc_cases_exp <- cc_cases %>% filter(E == 1) 
    cc_cases_unexp <- cc_cases %>% filter(E == 0) %>% mutate(E = rbern(nrow(.), 1-Sp)) # induce recall bias
    cc_cases <- bind_rows(cc_cases_exp, cc_cases_unexp)
    
    # sample controls
    cc_control <- df %>% filter(D == 0) %>% sample_n(n_control, replace=T)
    cc <- bind_rows(cc_cases, cc_control)
    
    # fit a logistic model
    model1 <- glm(formula = D ~ G*E, family = binomial(link = "logit"), data = cc)
    
    # calculate bias and coverage for each parameter
    results[i,1] <- summary(model1)$coef[4,1] - log(or3)
    results[i,2] <- ((summary(model1)$coef[4,1]-1.96*summary(model1)$coef[4,2]) < log(or3) & (summary(model1)$coef[4,1]+1.96*summary(model1)$coef[4,2]) > log(or3))
    results[i,3] <- mean(df$E) # prevalence of exposure
    
    # monitor simulation progress
    progress(i,iter)
  }
  results_summary <- data.frame(bias = round(mean(results$bias),2), cover = round(mean(results$cover),2), prev = round(mean(results$prev),2),
                                or1 = or1, or2 = or2, or3 = or3, Sp = Sp)
  results_summary
}

###**********************
### 3: Example Inputs ###
###**********************

simResults <- geINT(iter=500, n=100000, pE=0.4, pG=0.4, pD=0.01, or1=1.5, or2=1.5, or3=1, n_case=1000, n_control=1000, Sp=0.8)
