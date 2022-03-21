if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("cBioPortalData")

require('cBioPortalData')
study_id = 'brca_metabric'
cbio <- cBioPortal()
dat <- getStudies(cbio)
gene <- genePanels(cbio)
metabric <- cBioPortalData(cbio, studyId = 'brca_metabric', genePanelId = 'Affymetrix')
metabric <- cBioPortalData(cbio, studyId = study_id, genePanelId = 'Affymetrix', )