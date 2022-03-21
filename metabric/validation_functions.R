require(fs)
require(dplyr)
require(stringr)
require(tidyverse)
require(survival)
require(survminer)
require(viridis)

test.sig <- function(p){
  for (i in 1:abs(length(p))){
    if (!is.na(p[i]) & p[i] >= 0.05){
      p[i] <- paste0(format(p[i], scientific = FALSE, digits = 2), ' (n.s)')
    } else {
      p[i] <- format(p[i], scientific = FALSE, digits = 5)
    }
  }
}  


get_landmark_analysis <- function(gene_var,fit_surv, dat, time_to_event, surv){
  landmark = time_to_event
  landmark_time = as.numeric(landmark)
  surv_sum = surv_summary(fit_surv, data = os_data)
  group_var = sym('group')
  n_surv = surv_sum %>%
    group_by(!!group_var, .drop = TRUE)%>%
    filter(time == min(time))%>%
    ungroup()%>%
    select(!!group_var, n.risk)%>%
    rename(all_risk = n.risk)
  surv_sum <- surv_sum%>%
    group_by(!!group_var, .drop = TRUE) %>%
    slice(which.min(abs(time - landmark_time)))%>%
    select(surv, !!group_var, n.risk)%>%
    left_join(., n_surv)%>%
    pivot_wider(names_from = !!group_var, values_from = c(surv, all_risk,
                                                          n.risk)) %>%
    mutate(surv_landmark_high_low = surv_high - surv_low) %>%
    mutate(time_point = landmark)
  return(surv_sum)
}

get_cox_df <- function(gene_var, os_data)
{ names_cox <- c('HazardRatio', 'WaldTest', 'LikelihoodRatio', 'LogRank')
fit_cox <- coxph(Surv(os_months, events) ~ group, data = os_data)
summary(fit_cox)
zph <- cox.zph(fit_cox)
cox_df = data.frame(gene = gene_var, zph$table)
col_names <- colnames(cox_df)
cox_df <- cox_df %>%
  mutate(param = rownames(.))%>%
  select(-gene)%>%
  pivot_wider(names_from = param, values_from = c(chisq, df, p))%>%
  mutate('gene' = gene_var)
cox_sum <- summary(fit_cox)
coeff <- cox_sum$coefficients
hr =  cox_sum$conf.int
colnames(hr) = c('HazardRatio_HighvsLow', 'HazardRatio_LowvsHigh', 
                 'conf_int_95_low', 'conf_int_95_high')
cox_test = c(cox_sum$waldtest[c(1,3)], cox_sum$logtest[c(1,3)], 
             cox_sum$sctest[c(1,3)])
names(cox_test) = c('wald_test_test_score', 'wald_test_p_value', 
                    'likelihood_ratio_test_test_score', 
                    'likelihood_ratio_test_pvalue', 'log_rank_test_score', 
                    'log_rank_pvalue')
cox_p = data.frame(t(cox_test))
cox_df <- cbind(cox_df, hr, cox_p)%>%
  select(gene, everything(.))
return(cox_df)
}
