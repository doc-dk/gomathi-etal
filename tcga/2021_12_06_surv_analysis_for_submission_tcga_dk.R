  # setwd('~/repos/TCGA_YAP/')
  # source('tcga_functions_survival.R')
  source('~/repos/TCGA_YAP/survival_anlaysis_var_tcga.R')
  # get tcga fpkm upper quartile normalized data from gdc using below function
  # and save as a '.RDS' object on file.
  # get_fpkm_data()
  ensembl_names <- get_ensembl_fpkm()
  # clinical data was also downloaded from 
  clin_file = 'brca_clinical_data_tcga_annotated_517_femalePrimaryIDC.xlsx'
  clin_sheet = 'Female_Primary_Final_IDC'
  clin_dat <- get_brca_tcga_clin_clean(clin_file, clin_sheet)
  # data is log2 normalized when fpkm data is retrieved
  filter_list <- c('all', 'ERPR', 'HER2', 'TNBC')
  survival_types = c('disease_specific_survival', 'overall_survival',
                     'disease_free_survival')
  gene_var<- gene_list[1]
  filter_col <- filter_list[1]
  type <- survival_types[1]
  time_to_event = 60
  fpkm_gene <- get_fpkm_gene(clin_dat, ensembl_names, fpkm, gene_var, filter_col)
  clin_gene = left_join(fpkm_gene, clin_dat, by = 'barcode')
  os_dat = tcga_get_os_data_no_strata(clin_gene, type)
  dat = get_median_low_high_group(os_dat) %>%
    mutate(events_os = if_else(os_months < time_to_event, events_os, 0))
  
    fit_surv <- do.call(survfit, list(formula = Surv(time = os_months,
                                                   event = events_os)~group,
                                    data = dat))
  # fit_diff <- do.call(survdiff, list(formula = Surv(time = os_months,
  #                                                   event = events_os)~group,
  #                                    data = os_data))
  title_text <- paste(gene_var, 'subset: ', str_to_upper(gsub("_"," ",filter_col)))
  y_label <- str_replace_all(type, "_", " ")%>%
    str_to_title()
  
  gene_plot = ggsurvplot(fit_surv,
                         data=dat,
                         xlab='Months',
                         break.x.by = 20,
                         ylab = y_label,
                         #title = 'Effect of EGFR deletion (GISTIC < 0) on 4 year survival',
                         conf.int=FALSE,
                         risk.table=FALSE,
                         risk.table.y.text=FALSE,
                         cumcensor = FALSE,
                         cumcensor.table.y.text = FALSE,
                         ggtheme=theme(line = element_line(size = 1, colour = 'black'),
                                       panel.grid.minor=element_blank(),
                                       panel.background = element_blank(),
                                       legend.background = element_blank(),
                                       axis.line.x = element_line(size = 0.5),
                                       axis.line.y.left = element_line(size = 0.5),
                                       text=element_text(size=16, 
                                                         face = 'bold', 
                                                         colour = 'black'),
                                       axis.text = element_text(size=16, 
                                                                face = 'bold', 
                                                                colour = 'black'),
                                       legend.key = element_blank()),
                         tables.height=0.2,
                         tables.theme=theme(panel.border = element_blank(),
                                            panel.grid.major=element_blank(),
                                            axis.title.x=element_blank(),
                                            axis.title.y=element_blank(),
                                            axis.text.x=element_blank(),
                                            axis.ticks.x=element_blank()),
                         pval = TRUE,
                         pval.method = TRUE,
                         pval.size = 5,
                         pval.color = 'black',
                         test.for.trend = FALSE,
                         legend.labs = c(paste0(gene_var, '-high'),
                                         paste0(gene_var, '-low')),
                         legend = c(0.8, 0.2),
                         censor.shape = 124,
                         censor.size = 3,
                         xlim = c(0,60),
                         fontsize = 3)

    # annotate(text, x = c(1.8, 0.8), y = c(0.2, 0.2), )
    gene_plot

  