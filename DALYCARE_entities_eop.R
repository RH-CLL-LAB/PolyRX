## DALYCARE ENTITIES 
# Created:
## Sys.Date() # "2024-04-25"

## Loads RKKP dalycare_dx and divides patients into different LC entities 
source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_entities.R') #Also saved as #load_dalycare_icd10()
load_dataset(c('t_dalycare_diagnoses', 'patient'))

# first define
# 1. IPI please
# 2. all medicine, DX.FIRST.POLY and cohort (1)
# 3. ALL.DISEASE/PATIENT_OS2

# DX.FIRST.POLY = read_csv('/ngc/projects2/dalyca_r/chribr_r/PolyRx/data/DX_FIRST_POLY.csv') # Saved in CCI_polyRX script
# cohort = DX.FIRST.POLY$patientid # define cohort1 #added 28/10-24 
# 
# PATIENT_OS2 = read_csv2('/ngc/projects2/dalyca_r/chribr_r/PolyRx/data/PATIENT_OS2.csv') # from Table_3.R
# PATIENT_OS2 %>% n_patients()
# PATIENT_OS2$n.ATC %>% table

ENTITIES = t_dalycare_diagnoses %>% 
  filter_first_diagnosis(c(ICD10.HL, ICD10.DLBCL, ICD10.FL, ICD10.MZL, ICD10.MCL, ICD10.CLL, ICD10.LPL, ICD10.MM), str_contains = F) %>% 
  mutate(Disease = ifelse(diagnosis %in% ICD10.HL, 'cHL', NA),
         Disease = ifelse(diagnosis %in% ICD10.DLBCL, 'DLBCL', Disease),
         Disease = ifelse(diagnosis %in% ICD10.FL, 'FL', Disease),
         Disease = ifelse(diagnosis %in% ICD10.MZL, 'MZL', Disease),
         Disease = ifelse(diagnosis %in% ICD10.MCL, 'MCL', Disease),
         Disease = ifelse(diagnosis %in% ICD10.CLL, 'CLL', Disease),
         Disease = ifelse(diagnosis %in% ICD10.LPL, 'LPL', Disease),
         Disease = ifelse(diagnosis %in% ICD10.MM, 'MM', Disease))  %>%
  filter(time_dx_death >=0) %>% 
  left_join(IPI, c('patientid', 'Disease')) %>% 
  mutate(Disease = factor(Disease, levels = c('cHL', 'DLBCL',  'FL', 'MZL', 'MCL', 'CLL', 'LPL',  'MM'))) %>%   
  mutate(Age = diff_years(date_birth, date_diagnosis)) 

ENTITIES %>%  nrow_npatients()
ENTITIES$date_diagnosis %>% summary
## Go to Poly_CCI_eop.R