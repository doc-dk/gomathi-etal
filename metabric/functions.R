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
  
  p
}

mytheme <- function(...){
  theme(panel.grid.minor=element_blank(),
        line = element_line(size = 1),
        plot.margin=unit(c(1,3,1,3), 'lines'),
        text = element_text(family = 'sans', size = 20, face = 'bold'),
        axis.title.x = element_text(size = 25, face = 'bold'),
        axis.title.y=element_text(size=25, face = 'bold', vjust=1, 
                                  angle=90, margin = margin(t = 0, r= 20, 
                                                            b = 0, l = 0)),
        axis.text = element_text(size=20, face = 'bold'), 
        axis.line = element_line(size=1.2),
        legend.direction = 'horizontal',
        legend.title=element_blank(),
        legend.text=element_blank(),
        legend.key=element_blank(),
        plot.title = element_text(size = 20, face = 'bold'),
        strip.text.y=element_text(size=20))}
bar_theme <- mytheme() + theme(legend.text=element_text(size = 16), 
                               legend.direction = 'vertical', 
                               legend.position = 'right', 
                               legend.title=element_text(size=18, face = 'plain'))


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

get_cox_info <- function(gene_var, gene_clin, type){
  cox_df <- as.data.frame(matrix(c(gene_var, rep('data_not_available', 29)),
                                 ncol = 30, nrow = 1))
  colnames(cox_df) <- c('gene_var', 'chisq', 'df', 'p')
  os_dat <- get_os_data_no_strata(gene_clin, type)
  os_data <- get_median_low_high_group(os_dat) %>%
    mutate(events = if_else(os_months > time_to_event, 0, events))
  cox_df = get_cox_df(gene_var, os_data)
   return (cox_df)
}

get_gene <- function(gene, metabric_expression){
  gene_exp <- metabric_expression %>%
    filter(Hugo_Symbol == gene) %>%
    dplyr::select(-Hugo_Symbol, -Entrez_Gene_Id) %>%
    t(.) %>%
    data.frame()
  if(ncol(gene_exp)>0){
    gene_exp <- gene_exp %>%
      rename(expression = '.')%>%
      mutate(sample_id = rownames(.))
  }
  gene_exp
}

get_gene_and_clin_dat <- function(gene_exp, met_clin_clean,filter_col_name = 'all', tbd_filter = FALSE){
  gene_exp_clinical <- gene_exp %>%
    mutate(expression = as.numeric(expression))%>%
    mutate(sample_id = gsub("\\.", "-", sample_id)) %>%
    mutate(patient_id = sample_id) %>%
    left_join(., met_clin_clean, by = 'patient_id') %>%
    mutate(expression = as.numeric(as.character(expression)))
  if(tbd_filter == TRUE){
    filter_col <- sym(filter_col_name)
    gene_exp_clinical <- filter(gene_exp_clinical, !!filter_col == TRUE)
  }
  return(gene_exp_clinical)
}

get_os_data_no_strata <- function(gene_clin, type){
  status_col <- case_when(type == 'disease_specific_survival' ~ 'vital_status',
                          type == 'disease_free_survival' ~ 'relapse_free_status',
                          type == 'all_cause_mortality'~ 'os_status')
  months_col <- if_else(grepl('^relapse', status_col),
                        'relapse_free_status_months',
                        'os_months')
  status <- sym(status_col)
  time <- sym(months_col)
  os_data <- gene_clin %>%
    mutate(os_months = as.numeric(as.character(!!time)))%>%
    mutate(status_events = str_to_lower(!!status))%>%
    mutate(events = if(grepl('^relapse', status_col)){
      if_else(grepl('^0', relapse_free_status), 0, 1)} else {
        if_else(status_events == 'died of disease'| status_events == 'deceased',
                1, 0)})%>%
    select(patient_id, expression, os_months, status_events, events)
  
  return (os_data)
}

get_median_low_high_group <- function(dat){
  cutpoint <- dat %>%
    select(expression)
  cut = summary(cutpoint)
  cut_med = str_split_fixed(cut[3], ":", 2)[2]
  dat <- dat %>%
    mutate(group = if_else(as.numeric(expression) <= as.numeric(cut_med), 'low', 'high'))
  
  return (dat)
}

get_landmark_analysis <- function(gene_var,fit_surv, dat, time_to_event, surv){
  landmark = time_to_event
  landmark_time = as.numeric(landmark)
  surv_sum = surv_summary(fit_surv, data = km_cat)
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

analysis.spearman <- function(dat_df, gene_var, factor_1, factor_2){
  # dat_df <- dat_df %>%
  #   mutate(factor_1 = ifelse(is.character(factor_1), as.numeric(factor_1), factor_1),
  #          factor_2 = ifelse(is.character(factor_2), as.numeric(factor_2), factor_2))
  colnames(dat_df) <- c('factor_1', 'factor_2', 'subtype')
  y_label <- max(dat_df$factor_2) + 0.15*max(dat_df$factor_2)
  x_label <- 0.7*(max(dat_df$factor_1))
  dat_pval = data.frame()
  for (type in unique(dat_df$subtype)){
    dat_test <- dat_df %>%
      filter(subtype == type)
    test <- cor.test(dat_test$factor_2, dat_test$factor_1, method = 'spearman',
                     exact = FALSE)
    pval_test <- test.sig(test$p.value)
    pval_test = ifelse(is.na(as.numeric(pval_test)),
                       pval_test,
                       as.numeric(pval_test))
    rho_val <- format(test$estimate, digits = 3)
    rho_val <- as.numeric(rho_val)
    n_dat <- nrow(dat_test)
    dat_type <- data.frame(subtype = type,
                           gene = gene_var,
                           N = n_dat,
                           spearman_rho = rho_val,
                           pval = pval_test)
    dat_pval = rbind(dat_pval, dat_type)
  }
  
  plot <- ggscatter(dat_df,
                    x = 'factor_1',
                    y = 'factor_2',
                    xlab = factor_1,
                    ylab = factor_2,
                    add = 'reg.line',
                    shape = 21,
                    size = 3,
                    title = 'Gene Expression Correlation',
                    fill = 'factor_2',
                    facet.by = 'subtype',
                    cor.coef = TRUE,
                    cor.coef.size = 7,
                    add.params = list(color='black'),
                    show.legend.text = TRUE,
                    show.legend.title = TRUE) +
    mytheme() + theme(axis.title.x = element_text(),
                      legend.position = 'right',
                      legend.direction = 'vertical',
                      legend.text = element_text(size = 16))
  dat <- plot + scale_fill_viridis(alpha = 0.5, option = 'inferno',
                                   guide = 'colorbar')
  
  # plot <- plot_annotate(plot, x_label, y_label, pval_test, rho_val, n_dat)
  # sum_dat <- data.frame(rho = rho_val)
  dat <- list(pval = pval_test, plot = plot, sum_dat = dat_pval)
  return(dat)
}

plot_annotate <- function(p, x_label, y_label,pval_test, rho_val, n_dat){
  p = p + annotate('text', x = x_label,y = y_label, size = 4, label =
                     paste('p = ', pval_test))+
    annotate('text', x = x_label,y = (y_label - 0.05*y_label), size = 4, label =
               paste('rho = ', rho_val))+
    annotate('text', x = x_label,y = (y_label - 0.1*y_label), size = 4, label =
               paste('n = ', n_dat))
  return(p)
}

gene_corr.2 <- function(gene_1, gene_2, expression_data, clinical_data){
  gene_1_expression <- get_gene(gene_1, expression_data)%>%
    rename(factor_1 = expression)
   clinical_data = met_idc_subtype%>%
    select(sample_id, subtype)
  dat_df <- get_gene (gene_2, expression_data)%>%
    right_join(gene_1_expression, ., by = 'sample_id')%>%
    mutate(sample_id = gsub('[[:punct:]]', '-', sample_id))%>%
    left_join(clinical_data, ., by = 'sample_id')%>%
    select(factor_1, expression, subtype)%>%
    rename(factor_2 = expression)%>%
    filter(complete.cases(.))
  dat <- analysis.spearman(dat_df, gene_var= gene, factor_1 = gene_1,
                           factor_2 = gene_2)
  
  return(dat)
}

clin_gene_analysis <- function(feature_col, gene_var, 
                               filter_col = 'all', 
                               met_idc_subtype){
  pval_test <- 'data_not_available'
  plot <- 'data_not_available'
  sum_dat <- 'data_not_available'
  sum_dat <- NA
  dat <- list(pval = pval_test, graph = plot, df_sum = sum_dat)
  gene_exp <- get_gene(gene_var, metabric_expression)
  gene_exp_clin <- get_gene_and_clin_dat(gene_exp, met_idc_subtype, 
                                         filter_col_name = 'all', 
                                         tbd_filter = FALSE)
  selected_dat <- get_selected_feature_pval(gene_exp_clin, feature_col)
  dat_df <- selected_dat$df %>%
    filter(!is.na(subtype))
  #print(colnames(gene_exp_clin))
  if (is.data.frame(dat_df)){
    n_dat <- nrow(dat_df)
    plot <- gene_var
    #    print(feature_col)
    feature <- sym(feature_col)
    ptest <- t.test(expression ~ clinical_feature, dat_df)
    pval_test <- ptest$p.value
    pval_test_type <- ptest$method
    n_y <- length(unique(dat_df$clinical_feature))
    label_text <- paste('p =', format(pval_test, digits = 4, scientific = TRUE),
                        '\nn = ', n_dat)
    y_label <- summary(dat_df$expression)[5]
    sum_dat <- dat_df %>%
      group_by(subtype, clinical_feature)%>%
      summarise(mean = mean(expression), median = median(expression), n_cases = n())%>%
      ungroup()%>%
      pivot_wider(names_from = clinical_feature,
                  values_from = c(mean, median, n_cases),
                  names_sep = '_')
    # title_text <- paste(str_to_upper(gsub("_", " ", feature_col)), ':',
    #                     gene_var, 'in\n', str_to_upper(gsub("_", " ",
    #                                                         filter_col)),
    #                     'subtype')
    x_text <- str_to_title(gsub("_", " ", feature_col))
    plot <- ggboxplot(dat_df,
                      x = 'clinical_feature',
                      y = 'expression',
                      add = c('mean', 'jitter'),
                      # data = feature_col,
                      ylab = 'gene expression',
                      xlab = x_text,
                      palette = viridis(n = n_y, alpha =0.5),
                      #fill = 'clinical_feature',
                      color = 'clinical_feature',
                      show.legend = FALSE,
                      title = gene_var,
                      facet.by = 'subtype',
                      legend = 'none',
                      order = sort(unique(dat$clinical_feature))
                      # ylim = c(0.5 * min(dat$expression), 3*(iqr))
    ) + mytheme() +
      theme(axis.title.x = element_text()) +
      stat_compare_means(method = 't.test', label.y = 8, size = 6)
      # annotate('text', x = (n_y), y = y_label, label = label_text, size = 4)
    #      scale_y_continuous(trans='log2')
    # sum_dat$pval <- pval_test
    # sum_dat$pval_test_type <- pval_test_type
    sum_dat <- left_join(sum_dat, selected_dat$pval_test, by = 'subtype')
    sum_dat$gene <- gene_var
    dat <- list(pval = pval_test, graph = plot, df_sum = sum_dat)
  }
}

get_selected_feature_pval <- function(gene_exp_clin, feature_col){
  feature <- sym(feature_col)
  dat_all <- gene_exp_clin %>%
    select(!!feature, expression, subtype)%>%
    rename(clinical_feature = !!feature)%>%
    filter(!is.na(clinical_feature))
  dat_count <- dat_all %>%
    count(clinical_feature)%>%
    filter(n>2)
  if (nrow(dat_count) == 2){
    df_pval = data.frame()
    for (type in unique(dat_all$subtype)){
      dat <- dat_all%>%
        filter(subtype == type)
      ptest <- t.test(expression ~ clinical_feature, dat)
      pval_test <- ptest$p.value
      pval_test_type <- ptest$method
      df_var = data.frame(pval = pval_test, 
                          pval_test = pval_test_type, 
                          subtype = type)
      df_pval = rbind(df_pval, df_var)
    }
  } else { df_pval <- NA
      #https://thestatsgeek.com/2013/09/28/the-t-test-and-robustness-to-non-normality/
    }
  selected_dat <- list(df = dat_all, pval_test = df_pval)
  return(selected_dat)
}
