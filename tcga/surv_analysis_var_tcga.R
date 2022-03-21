require(dplyr)
require(tidyverse)
require(survival)
require(survminer)
require(pander)
require(kableExtra)
require(viridis)
require(stringi)

get_fpkm_data <- function(){
  query.exp.hg38 <- GDCquery(project = "TCGA-BRCA", 
                             data.category = "Transcriptome Profiling", 
                             data.type = "Gene Expression Quantification", 
                             workflow.type = "HTSeq - FPKM-UQ")
  GDCdownload(query.exp.hg38)
  fpkm.uq.counts <- GDCprepare(query = query.exp.hg38, summarizedExperiment = FALSE)
  # Save an object to a file
  saveRDS(fpkm.uq.counts, file = "brca_fpkm_uq_counts.rds")
}
get_clin_dat <- function(clin_file, clin_sheet){
    clin_file = 'brca_clinical_data_tcga_annotated_517_femalePrimaryIDC.xlsx'
    clin_sheet = 'Female_Primary_Final_IDC'
    clinical_brca <<- readxl::read_excel(
    file.path(clin_file),
    sheet = clin_sheet, col_names = TRUE)
  clin_clean <- clinical_brca %>%
    select(bcr_patient_barcode...2, 'Final Pathology Diagnosis', er_status_by_ihc,
           pr_status_by_ihc, her2_status_by_ihc, subtype,
           PANCAN_CDR_residual_tumor, PANCAN_CDR_OS, PANCAN_CDR_OS.time,
           PANCAN_CDR_DSS, PANCAN_CDR_DSS.time, PANCAN_CDR_DFI,
           PANCAN_CDR_DFI.time,PANCAN_CDR_PFI, PANCAN_CDR_PFI.time) %>%
    rename(final_path = 'Final Pathology Diagnosis',
           barcode = bcr_patient_barcode...2,
           event_os = PANCAN_CDR_OS,
           time_os = PANCAN_CDR_OS.time,
           event_dss = PANCAN_CDR_DSS,
           time_dss = PANCAN_CDR_DSS.time,
           event_dfs = PANCAN_CDR_PFI,
           time_dfs = PANCAN_CDR_PFI.time)%>%
    filter(subtype != "NA")
  return(clin_clean)
}
get_ensembl_fpkm <- function(){
  fpkm <<- readRDS(file = "brca_fpkm_uq_counts.rds")
  ensemebl_file <- file.path('ensembl_gene_list','gene_data_12_01_2019_all.csv')
  ensembl_list <- read.csv(ensemebl_file, header= TRUE)
  ensembl <- fpkm %>%
    select(X1)%>%
    separate(X1, 1, into = c('ensembl_gene', 'version'), sep = "[.]",
             extra = 'merge', remove = FALSE)%>%
    left_join(., ensembl_list, by = 'ensembl_gene')
}

get_fpkm_gene <- function(dat, ensembl_names, fpkm, gene_var, subtype_type){
  if (filter_col=='all'){
    barcode_filtered = dat
  }else{
    barcode_filtered = dat %>% filter (subtype == subtype_type)
  }
  egene_val <- ensembl_names[ensembl_names$symbol == gene_var, 1]
  egene_val <- as.character(egene_val[complete.cases(egene_val),])
  gene.fpkm <- fpkm %>%
    filter(X1 == egene_val)
  barcode.fpkm = 'check_gene'
  if (nrow(gene.fpkm) > 0){
    gene.fpkm.sample <- gene.fpkm %>%
      t(.) %>%
      data.frame()%>%
      rename(fpkm = '.')%>%
      mutate(barcode = row.names(.)) %>%
      filter(barcode != 'X1')%>%
      rename(barcode_all = barcode)%>%
      mutate(fpkm = log2(as.numeric(as.character(fpkm))+1))
    barcode <- get_barcode(gene.fpkm.sample$barcode_all)
    barcode.fpkm <- right_join(barcode, gene.fpkm.sample, by = 'barcode_all')%>%
      filter(complete.cases(.))%>%
      select(-barcode_all)%>%
      average_remove_duplictes(., id_col_pos = 1, val_col_pos = 2)
    barcode.fpkm <- barcode.fpkm %>%
      inner_join(barcode_filtered)%>%
      average_remove_duplictes(., id_col_pos = 1, val_col_pos = 2)%>%
      select(barcode, fpkm)
    
  }
  return(barcode.fpkm)
}

get_barcode <- function(barcode_all){
  sample_barcode <- as.data.frame(barcode_all)%>%
    separate(1, into = c('1', '2', "3", "4", "5", "6", "7"), sep = "-",
             extra = 'merge', remove = FALSE) %>%
    unite(., "barcode",  c('1', '2', '3'), sep = '-')%>%
    filter(!grepl("11", `4`))%>%
    select(barcode, barcode_all)
  return(sample_barcode)
}

tcga_get_os_data_no_strata <- function(gene_clin, type){
  events_col <- case_when(type == 'disease_specific_survival' ~ 'event_dss',
                          type == 'disease_free_survival' ~ 'event_dfs',
                          type == 'overall_survival'~ 'event_os')
  days_col <- case_when(type == 'disease_specific_survival' ~ 'time_dss',
                        type == 'disease_free_survival' ~ 'time_dfs',
                        type == 'overall_survival'~ 'time_os')
  events <- sym(events_col)
  time <- sym(days_col)
  # os_data = tcga_get_os_data_no_strata(gene_clin, type) %>%
  #   rename(expression = fpkm, !!events, !!time)
  os_data <- gene_clin %>%
    mutate(os_months = as.numeric(as.character(!!time))/30)%>%
    mutate(events_os = as.numeric(as.character(!!events)))%>%
    # mutate(events =  ifelse(os_months<time_to_event, events, 0))%>%
    select(barcode, events_os, os_months, subtype, fpkm)%>%
    rename(expression = fpkm)
  return (os_data)
}

get_median_low_high_group <- function(dat){
  cutpoint <- dat %>%
    select(expression)
  cut = summary(cutpoint)
  cut_med = str_split_fixed(cut[3], ":", 2)[2]
  # dat <- dat %>%
  #   mutate(group = if_else(as.numeric(expression) <= as.numeric(cut_med), 'low', 'high'))
  dat <- dat %>%
    mutate(group = if_else(as.numeric(expression) <= as.numeric(cut_med), 'group_2', 'group_1'))
  return (dat)
}
