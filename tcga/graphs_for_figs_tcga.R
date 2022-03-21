  setwd('~/repos/TCGA_YAP/')
  source('~/repos/TCGA_YAP/functions_survival.R')
  ensembl_names <- get_ensembl_fpkm()
  clin_file = '2021_07_16_brca_clinical_data_annotated_517_femalePrimaryIDC_dk_gv.xlsx'
  
  clin_sheet = 'Female_Primary_Final_IDC'
  
  clin_dat <- get_brca_tcga_clin_clean(clin_file, clin_sheet)
  # data is log2 normalized when fpkm data is retrieved
  # filters = filters[1:2]
  # all_genes = all_genes[1:5]
  survival_types = c('disease_specific_survival', 'overall_survival',
                     'disease_free_survival')
  gene_surv_list <- list(c('COL12A1', 'disease_free_survival', 'all', 'A'), 
                         c('ABCA8', 'overall_survival', 'ERPR','C'),
                         c('COL12A1', 'disease_free_survival', 'HER2', 'D'),
                         c('THBS1', 'overall_survival', 'TNBC', 'E'))
  dat_tcga <- data.frame()
    # events_col <- case_when(type == 'disease_specific_survival' ~ 'event_dss',
    #                         type == 'disease_free_survival' ~ 'event_dfs',
    #                         type == 'overall_survival'~ 'event_os')
    # days_col <- case_when(type == 'disease_specific_survival' ~ 'time_dss',
    #                       type == 'disease_free_survival' ~ 'time_dfs',
    #                       type == 'overall_survival'~ 'time_os')
    # events <- sym(events_col)
    # time <- sym(days_col)
    # # os_data = tcga_get_os_data_no_strata(gene_clin, type) %>%
    # #   rename(expression = fpkm, !!events, !!time)
    # os_data <- gene_clin %>%
    #   mutate(os_months = as.numeric(as.character(!!time))/30)%>%
    #   mutate(events_os = as.numeric(as.character(!!events)))%>%
    #   # mutate(events =  ifelse(os_months<time_to_event, events, 0))%>%
    #   select(barcode, events_os, os_months, subtype, fpkm)%>%
    #   rename(expression = fpkm)
  
    dat$type = type
    dat$data_set = 'TCGA'
    dat$subtype = filter_col
    dat$gene = gene_var
    dat_tcga <- rbind(dat_tcga, dat)
  gene_surv <- gene_surv_list[[4]]
  gene_var<- gene_surv[1]
  filter_col <- gene_surv[3]
  type <- gene_surv[2]
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
    saveRDS(gene_plot, 'h.RDS')
  fig_b <-  readRDS('b.RDS')
  fig_d <-  readRDS('d.RDS')
  fig_f = readRDS('f.RDS')
  fig_h =  readRDS('h.RDS')
  fig_a <- readRDS('a.RDS')
  fig_c <- readRDS('c.RDS')
  fig_e <- readRDS('e.RDS')
  fig_g <- readRDS('g.RDS')
  final_fig <- list()
  final_fig<- list( b = fig_b, 
                   d = fig_d, 
                   f = fig_f, 
                   h = fig_h)
  
  pdf('survival_tcga.pdf', width = 6, height = 3)
  final_fig
  dev.off()
  require(ggtext)
  fig <- ggarrange(plotlist = final_fig, ncol = 2)
  #   saveRDS(fig_g, 'g.RDS')
# a <- readRDS('a.RDS')  

   saveRDS(fig_h, 'h.RDS')
  