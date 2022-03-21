####Metabric Data Survival Analysis

source('functions.R')

# metabric source files downloaded from cbioportal.org on 25.07.2021 and extracted
#read metabric data from source files
metabric_expression <- read.table('data_expression_median.txt', header=TRUE, sep ='\t',
                                  stringsAsFactors = FALSE, fill=TRUE)
expression_data <- metabric_expression

# downloaded file brca_metabric_clinical_data.txt was manually read into excel 
# and saved as .xlsx due to an error when the downloaded txt file was read 
# directly by `read.table`

metabric_clinical <- readxl::read_xlsx('brca_metabric_clinical_data.xlsx', na = 'NA')

#modify column names and case, add subtype filters, add column for events code
colnames(metabric_clinical) <- str_to_lower(colnames(metabric_clinical))
colnames(metabric_clinical) <- str_replace_all(colnames(metabric_clinical), " ", "_")
colnames(metabric_clinical) <- gsub("\\(|\\)", "", colnames(metabric_clinical))
metabric_clinical <- data.frame(lapply(metabric_clinical, function(x)
  if(is.character(x)) {tolower(x)}else{x}))
metabric_clinical <- metabric_clinical %>%
mutate(grade_type = case_when(neoplasm_histologic_grade == '3' ~ 'high_grade',
                              neoplasm_histologic_grade == '1'| 
                                neoplasm_histologic_grade == '2' ~ 'low_grade'),
       age_class = case_when(age_at_diagnosis < 30 ~ '<30',
                             (age_at_diagnosis == 30| age_at_diagnosis > 30)
                             & age_at_diagnosis < 40 ~ '30-40',
                             (age_at_diagnosis == 40|age_at_diagnosis >40) &
                               age_at_diagnosis <50 ~ '40-50',
                             (age_at_diagnosis == 50|age_at_diagnosis >50)
                             & (age_at_diagnosis<60) ~ '50-60',
                             age_at_diagnosis > 60 ~ '>60'),
       age_young_old = case_when(age_at_diagnosis < 50|
                                   age_at_diagnosis == 50 ~ '<=50',
                                 age_at_diagnosis >50 ~ '>50'),
       pT = case_when(tumor_size <20| tumor_size == 20 ~ 'T1',
                                tumor_size >20 & (tumor_size < 50|
                                                    tumor_size == 50) ~ 'T2',
                                tumor_size >50 ~ 'T3'),
       pN = case_when(lymph_nodes_examined_positive >9 ~ 'N3',
                            lymph_nodes_examined_positive < 10 &
                              lymph_nodes_examined_positive >3   ~ 'N2',
                            lymph_nodes_examined_positive <4 &
                              lymph_nodes_examined_positive>0 ~ 'N1',
                            lymph_nodes_examined_positive == 0 ~ 'N0'),
       er_status = str_to_lower(er_status),
       pr_status = str_to_lower(pr_status),
       her2_status = str_to_lower(her2_status),
       tnbc = if_else(er_status == 'negative' & pr_status == 'negative'
                      & her2_status == 'negative', TRUE, FALSE),
       erpr_her2minus = if_else((er_status == 'positive'| 
                                   pr_status == 'positive') & 
                                  her2_status == 'negative',
                                TRUE, FALSE),
       all_her2_pos = if_else(her2_status == 'positive', TRUE, FALSE),
       non_tnbc = if_else(er_status == 'negative' & pr_status == 'negative'
                          & her2_status == 'negative', FALSE, TRUE),
       erpr_all = if_else(er_status == 'positive' | pr_status == 'positive',
                          TRUE, FALSE),
       only_her2_pos = if_else(er_status == 'negative' & pr_status == 'negative'
                               & her2_status == 'positive', TRUE, FALSE),
       all = TRUE,
       os_status_events = case_when(grepl("0", overall_survival_status) ~ 0,
                                    grepl("1", overall_survival_status) ~ 1),
       recurrence_status_events = case_when(grepl("0", relapse_free_status) ~ 0,
                                            grepl("1", relapse_free_status) ~ 1),
       ln_status = if_else(lymph_nodes_examined_positive >0,
                                  'positive', 'negative'),
       tumor_type = if_else(tumor_size > 50, 'large_tumor', 'small_tumor'),
       stage_type = case_when(tumor_stage < 3 & !is.na(tumor_stage) ~ 'early_stage',
                                tumor_stage == 3|tumor_stage>3 ~'late_stage'),
       stage_type = if_else(pT == 'T2' & pN == 'N1', 'late_stage', stage_type),
       stage = tumor_stage)
met_clin_clean <- metabric_clinical

## filter for only IDC cases ####
met_idc <- met_clin_clean %>%
  filter(tumor_other_histologic_subtype == 'ductal/nst')%>%
  filter(cancer_type == 'breast cancer')%>%
  rename(os_months = overall_survival_months) %>%
  mutate(os_status = str_split_fixed(overall_survival_status, ":", 2)[,2])%>%
  mutate(sample_id= str_to_upper(sample_id))%>%
  mutate(patient_id = str_to_upper(patient_id)) %>%
  mutate(vital_status = patient.s_vital_status)%>%
  mutate(subtype = 'IDC')
met_idc_subtype =met_idc%>%
  mutate(subtype = case_when(erpr_her2minus ~ 'ER/PR positive',
                             all_her2_pos ~ 'HER2 positive',
                             tnbc ~ 'TNBC'))%>%
  bind_rows(met_idc)



types_surv = c('disease_specific_survival', 'all_cause_mortality', 
               'disease_free_survival')
filter_list <- c('all', 'erpr_her2minus', 'all_her2_pos','only_her2_pos', 'tnbc','non_tnbc')

all_genes_x <- c('OLR1','ANKRD1','TNNT2','GJA5','MAMDC2','SRGN','PSG5', 'PSG4','WDR69','MFAP5','DAB2','PSG7','LHFP','FN1',
                 'FYN','LIFR','PLCB4','SULF1',' FBN1','ATP6V0A4','CAPN6','FSTL1','ADAMTS1','SLC44A2',
                 'PROS1','NAV3','PAPSS2','SLIT2',
                 'CCDC80','ABCA8','PAPPA','HEG1','CYP4F11','DNAH11','NEXN','COL12A1','CLIC5','TMEM2','TGFB2','TGM2','ITGB2',
                 'DAPK1','ADAMTS12','PLAC8','TRAM2','FBLN5','IGFBP3','CDH2','THBS1','NPPB','SPRR1B','FXYD3','SLITRK6','SULF2',
                 'S100A14','SERPINB13','FYB','FAT2','KRT15','S100A8','TP63','CLCA2')

gene_list <- str_remove(all_genes_x, " ")

# set filter type, survival type and gene from the types given above

filter_col <- filter_list[1]
gene_var <- gene_list[2]
type <- types_surv[1]

# run the below script to get km plot, logrank p value, number at risk and 
# probability of survival at 60 months in high and low gene expression groups and 
# cox hazard ratio

    gene_exp <- get_gene(metabric_expression = expression_data, gene = gene_var)
    if (ncol(gene_exp)>0){
      gene_clin <- get_gene_and_clin_dat (gene_exp, met_clin_clean = met_idc, 
                                          filter_col_name = filter_col, 
                                          tbd_filter = TRUE)
      
      time_to_event = 60
      gene_clin <- get_gene_and_clin_dat (gene_exp, met_clin_clean = met_idc, 
                                          filter_col_name = filter_col, 
                                          tbd_filter = TRUE)
      os_dat <- get_os_data_no_strata(gene_clin, type)
      os_data <- get_median_low_high_group(os_dat)
      dat <- os_data %>%
        mutate(events = if_else(os_months > time_to_event, 0, events))
      fit_surv <- do.call(survfit, list(formula = Surv(time = os_months,
                                                       event = events)~group,
                                        data = dat))
      fit_diff <- do.call(survdiff, list(formula = Surv(time = os_months,
                                                        event = events)~group,
                                         data = dat))
      title_text <- paste(gene_var, 'subset: ', str_to_upper(gsub("_"," ",filter_col)))
      y_label <- str_replace_all(type, "_", " ")%>%
        str_to_title()
      
      gene_plot = ggsurvplot(fit_surv,
                             data=dat,
                             xlab='Months',
                             break.x.by = 20,
                             ylab = y_label,
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
      km_pval <- surv_pvalue(fit_surv)%>%
        mutate(gene = gene_var)%>%
        mutate(km_stat_all = fit_diff$chisq)%>%
        select(gene, pval, km_stat_all, method)
      surv_sum = get_landmark_analysis(gene_var,fit_surv, dat, time_to_event, type)
      km_df <- bind_cols(km_pval, surv_sum)
      cox_df <- get_cox_info(gene_var, gene_clin, type)
      dat <- list(graph = gene_plot, km_pval = km_pval, cox = cox_df)
    }
    
    # correlation with YAP1 expression 
    gene_1 = 'YAP1'
    yap_corr = list()
    for (gene in gene_list){
      if (gene %in% expression_data$Hugo_Symbol){
        dat <- gene_corr.2(gene_1, gene, expression_data = expression_data, 
                           clinical_data = met_idc_subtype)
      } else { dat <- 'data_not_available'}
      yap_corr[[gene]] = dat
    }
    

## clinical data analysis
    features <- c('age_young_old', 'ln_status', 'tumor_type', 'stage_type')
    feature_col <- features[4]
    gene_var <- gene_list[2]
    dat_feature <- clin_gene_analysis(feature_col, gene_var, filter_col = 'all',met_idc_subtype)
    dat_feature$graph
    