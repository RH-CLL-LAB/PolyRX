## Creates Poly Rx figures: tests Cox PolyRx ~ IPI
## Sys.Date() # "2024-04-26" 
## author: christian brieghel

source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R')

#### load data created PolyRx_CCI_eop.R ####
DX.FIRST.POLY = read_csv('/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/DX_FIRST_POLY.csv')
DX.FIRST.POLY_30 = read_csv('/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/DX_FIRST_POLY_30.csv')
DX.FIRST.ATC = read_csv('/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/DX.FIRST.ATC.csv')
DX.FIRST.ATC2 = read_csv('/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/DX_FIRST_ATC2.csv') 
DX.FIRST.CCI = read_csv('/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/DX_FIRST_CCI.csv') 
ENTITIES2 = read_csv('/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/ENTITIES2.csv') %>% 
  mutate(Disease = factor(Disease, c('cHL', 'DLBCL', 'FL', 'MCL', 'MZL', 'CLL', 'LPL', 'MM')))
IPI = read_csv2('/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/IPI_2025.csv') 

ENTITIES2 %>% n_patients()
cohort_1 = ENTITIES2$patientid

DX.FIRST.POLY %>% nrow_npatients()
cohort_1 %>% n_distinct()
DX.FIRST.CCI %>% nrow_npatients()
DX.FIRST.ATC %>% nrow_npatients() # all prescriptions
ENTITIES2 %>% nrow_npatients()

#### RKKP quality data ####
load_dataset(c('RKKP_CLL', 'RKKP_LYFO', 'RKKP_DaMyDa')) 

CLL_clean = RKKP_CLL %>% clean_RKKP_CLL()
LYFO_clean = RKKP_LYFO %>% clean_RKKP_LYFO()
MM_clean = RKKP_DaMyDa %>% clean_RKKP_DAMYDA()

HOSPITALS = bind_rows(CLL_clean %>% select(patientid, hospital_id),
                      LYFO_clean %>% select(patientid, hospital_id),
                      MM_clean %>% select(patientid, hospital_id)) %>% 
  distinct() %>% #44 with >1 hospital_id
  group_by(patientid) %>% 
  slice(1) %>% 
  ungroup()

#### OS data and cohort 2 ####
ALL.DISEASES = ENTITIES2 %>% 
  select(-sex.y) %>% 
  left_join(DX.FIRST.POLY %>% select(patientid, Polypharmacy, n.ATC)) %>% 
  left_join(DX.FIRST.CCI %>% select(patientid, CCI.score, CCI.2011.update)) %>% 
  mutate(IPI2 = factor(IPI, levels = c('Low', 'Intermediate', 'High', 'Very high')),
         IPS = factor(IPI, levels = c('Low', 'High')),
         IPI = recode_factor(IPI2,
                             Low  = 'Low', 
                             Intermediate = 'Intermediate',
                             High = 'High',
                             `Very high` = 'High'),
         Sex = factor(sex.x, levels = c('F', 'M')),
         PS = factor(PS, levels = c(0:4)),
         Polypharmacy = factor(Polypharmacy, levels = c('No', 'Yes')),
         ATC.all.cut = cut(n.ATC, c(-Inf, 3, 7, 11, Inf), labels = c('0-3', '4-7', '8-11', '>11')),
         CCI.f = cut(CCI.score, c(-Inf, 2, Inf), c('0-2', '>2')),
         CCI = CCI.2011.update,
         Time = time_dx_death) %>% 
  left_join(HOSPITALS, 'patientid') 

ALL.DISEASES.COX = ALL.DISEASES %>% 
  filter(!is.na(Age),
         !is.na(Sex),
         !is.na(IPI),
         !is.na(Polypharmacy),
         !is.na(CCI.2011.update),
         !is.na(Time)) 

#Create cohort2
cohort_2 = ALL.DISEASES.COX$patientid 

ALL.DISEASES %>% n_patients() == n_distinct(cohort_1)
ALL.DISEASES$Disease %>% table(exclude = NULL)

ALL.DISEASES.COX %>% nrow_npatients()
ALL.DISEASES.COX$Disease %>% table(exclude = NULL)

# write_csv(ALL.DISEASES.COX, '/ngc/projects2/dalyca_r/chribr_r/PolyRx/data/ALL_DISEASES_COX.csv') 
# write_csv(ALL.DISEASES, '/ngc/projects2/dalyca_r/chribr_r/PolyRx/data/ALL_DISEASES.csv') 

#### Hospitalization data ####
SP_DATASETS
load_dataset('SP_ADT_haendelser', cohort_1)
SP_ADT_haendelser_subset %>% head2
SP_ADT_haendelser_subset$event_type_name %>% table


FINAL_COHORT_ADT.DURATION = SP_ADT_haendelser_subset %>%
  distinct() %>%
  mutate(EVENT_ADT = recode_factor(event_type_name,
                                   INDLÆGGELSE = 'IND',
                                   UDSKRIVNING = 'UD')) %>%
  filter(EVENT_ADT  %in%  c('IND', 'UD')) %>%
  transmute(patientid, Date_ADT = effective_time, EVENT_ADT) %>%
  distinct()

FINAL_COHORT_ADT.DURATION2 = FINAL_COHORT_ADT.DURATION %>% 
  left_join(ALL.DISEASES %>%  select(patientid, ATC.all.cut, date_diagnosis, date_death_fu), by = 'patientid') %>%
  filter(Date_ADT > date_diagnosis,
         Date_ADT < date_death_fu) %>% 
  select(patientid,Date_ADT, EVENT_ADT) %>% 
  arrange(patientid, Date_ADT, EVENT_ADT) %>% 
  mutate(DATE = Date_ADT,
         n.F = as.numeric(factor(EVENT_ADT))) %>% 
  mutate(Same = ifelse(n.F == dplyr::lag(n.F), 'Yes', 'No')) %>% 
  filter(Same == 'No') %>% 
  group_by(patientid, Date_ADT, n.F) %>% 
  slice(1) %>% 
  ungroup() %>% 
  mutate(n = ceiling(row_number()/2)) %>% 
  group_by(patientid, n) %>% 
  mutate(TIME = diff_days(dplyr::lag(Date_ADT), Date_ADT)) %>% 
  ungroup() %>% 
  filter(!is.na(TIME)) %>% 
  group_by(patientid) %>% 
  mutate(N.Contacts = n(),
         SUM.in.hospital = sum(TIME)) %>% 
  ungroup()  %>% 
  filter(TIME>1.000) # > 24 h # time-lapse mins

FINAL_COHORT_ADT.DURATION3 = FINAL_COHORT_ADT.DURATION2 %>% 
  right_join(ALL.DISEASES %>% 
               select(patientid, Age, Sex, IPI, Disease, CCI.2011.update, n.ATC, ATC.all.cut, Polypharmacy, date_diagnosis, date_death_fu , time_dx_death, status, hospital_id), by = 'patientid') %>% 
  left_join(go_live(), 'hospital_id') %>%
  filter(hospital_id %in% c('ROS', 'HER', 'RH'),
         date_death_fu >= date_golive) %>% 
  filter(Date_ADT > date_golive | is.na(Date_ADT)) %>% 
  mutate(date_golive_diagnosis = if_else(date_diagnosis < date_golive, date_golive, date_diagnosis),
         date_ADT_death_fu = if_else(!is.na(Date_ADT), Date_ADT, date_death_fu),
         time_golive_death = diff_years(date_golive_diagnosis, date_death_fu),
         time_golive_ADT_death = diff_years(date_golive_diagnosis, date_ADT_death_fu),
         ADT = ifelse(is.na(Date_ADT), 0, 1), 
         TIME = ifelse(is.na(Date_ADT), 0, TIME),
         SUM.in.hospital = ifelse(is.na(Date_ADT), 0, SUM.in.hospital),
         N.Contacts = ifelse(is.na(Date_ADT), 0, N.Contacts),
         N.Contacts.per.year = N.Contacts/time_dx_death,
         N_row = 1) 

FINAL_COHORT_ADT.1.DURATION = FINAL_COHORT_ADT.DURATION3 %>% 
  group_by(patientid) %>% 
  arrange(time_golive_ADT_death) %>% 
  slice(1) %>% 
  ungroup()

utable(Disease ~ ATC.all.cut, FINAL_COHORT_ADT.1.DURATION)

# Create cohort_3
cohort_3 = FINAL_COHORT_ADT.1.DURATION$patientid 
cohort_3 %>% n_distinct()
cohort_1 %>% n_distinct() - cohort_3 %>% n_distinct()

FINAL_COHORT_ADT.1.DURATION$date_ADT_death_fu %>% summary
FINAL_COHORT_ADT.1.DURATION %>% nrow_npatients()
FINAL_COHORT_ADT.1.DURATION$date_golive_diagnosis %>% summary
FINAL_COHORT_ADT.1.DURATION$ADT %>% table

#### Severe infection data ####

SP_DATASETS
load_dataset('SP_AdministreretMedicin', cohort_1)
SP_AdministreretMedicin_subset %>% head2
SP.AB = SP_AdministreretMedicin_subset %>%
  ATC_AB(atc = atc)

SP.INF = SP.AB %>% 
  AE_infection() 

ALL.DISEASES_INF = ALL.DISEASES %>% 
  filter(patientid %in% cohort_3) %>% 
  left_join(SP.INF %>% select(patientid, date_infection, DIFF, n_inf, n_days_IVAB, date_received), 'patientid')
cohort_3 %>% n_distinct()
ALL.DISEASES_INF %>% nrow_npatients()
SP.INF$date_received %>% table() # == date_last_fu

ALL.DISEASES_INF$date_received %>% table 
ALL.DISEASES_INF.1 = ALL.DISEASES_INF %>% 
  left_join(go_live(), 'hospital_id') %>% 
  mutate(date_diagnosis_golive = if_else(date_diagnosis < date_golive, date_golive, date_diagnosis),
         date_infection = if_else(date_infection < date_diagnosis_golive, NA, date_infection),
         n_days_IVAB = ifelse(is.na(n_days_IVAB), 0, n_days_IVAB)) %>% 
  group_by(patientid) %>% 
  arrange(patientid, date_infection) %>%
  mutate(n_days_IVAB_sum = sum(n_days_IVAB)) %>% 
  slice(1) %>% 
  ungroup() %>%
  right_truncation(date_event = date_infection, 
                   date_start = date_diagnosis_golive,
                   date_truncation = '2023-03-31') %>% 
  mutate(n_inf_per_year = n_inf/time_dx_death,
         n_inf_per_year = ifelse(is.na(n_inf_per_year), 0, n_inf_per_year)) %>% 
  filter(time_death_fu_truncated >= 0) 

COX_ADT = ALL.DISEASES %>%
  right_join(FINAL_COHORT_ADT.1.DURATION %>% select(patientid, time_golive_ADT_death, ADT)) %>% 
  filter(!is.na(IPI))

# Create unique ATC 3rd level variable
DX.FIRST.ATC3 = DX.FIRST.ATC2 %>% 
  mutate(ATC3 = str_sub(atc, 1, 3)) %>% 
  select(patientid, ATC3) %>% 
  group_by(patientid, ATC3) %>% 
  mutate(n = n()) %>%
  slice(1) %>% 
  ungroup()  

# Create wide format
DX.FIRST.ATC3.OS = DX.FIRST.ATC3 %>% 
  filter(!is.na(ATC3)) %>% 
  group_by(ATC3) %>% 
  mutate(n = 1) %>% 
  ungroup() %>% 
  spread(ATC3, n) %>% 
  right_join(ALL.DISEASES %>% select(patientid, time_dx_death, status)) %>% 
  mutate(across(A01:V07, ~ factor(ifelse(is.na(.), 'No', 'Yes'), levels = c('No', 'Yes'))))
if(DX.FIRST.ATC3.OS %>% n_patients() == ALL.DISEASES %>% n_patients()){
  cat('Tranfomration successfull. Continue.')
}else{'Re-interate spread() and/or mutate(across())'}

# Create specific ATC3 groups
DX.FIRST.ATC3.SUMMARY = DX.FIRST.ATC3 %>% 
  group_by(ATC3) %>% 
  mutate(n.patients = n(),
         freq = n.patients/n_distinct(cohort_1)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(ATC3, n.patients, freq)

COX_ADT = ALL.DISEASES %>%
  right_join(FINAL_COHORT_ADT.1.DURATION %>% select(patientid, time_golive_ADT_death, ADT)) %>% 
  filter(!is.na(IPI))

diseases = ALL.DISEASES$Disease %>% levels()
ATC.CLASSES = names(DX.FIRST.ATC3.OS)[2:80] 

ALL.DISEASES_INF.1.ATC = ALL.DISEASES_INF.1 %>% 
  left_join(DX.FIRST.ATC3.OS %>% select(-time_dx_death, -status)) %>% 
  select(patientid, contains(ATC.CLASSES), time_death_fu_truncated, event, Age, sex =sex.x)
COX_INF = ALL.DISEASES %>%
  right_join(ALL.DISEASES_INF.1.ATC %>% select(patientid, time_death_fu_truncated, event)) %>% 
  filter(!is.na(IPI))

cohort_4 = bind_rows(ALL.DISEASES_INF.1 %>% select(patientid, Age, Sex, IPI, n.ATC, CCI.2011.update, Time = time_death_fu_truncated),
                     FINAL_COHORT_ADT.1.DURATION %>% select(patientid, Age, Sex, IPI, n.ATC, CCI.2011.update, Time = time_golive_ADT_death)) %>% 
  filter(!is.na(Age),
         !is.na(Sex),
         !is.na(IPI),
         !is.na(n.ATC),
         !is.na(CCI.2011.update),
         !is.na(Time)) %>% 
  pull(patientid) %>% 
  unique()

#### ORR data ####
CLL_clean = CLL_clean %>% 
  left_join(ALL.DISEASES, 'patientid')

LYFO_clean = LYFO_clean %>% 
  left_join(ALL.DISEASES %>% select(-date_death_fu), 'patientid') %>% 
  mutate(response = ifelse(response_1st_line %in% c('UNK', 'not_performed'), NA, response_1st_line),
         response = recode_factor(response,
                                  `CMR with residual mass` = 'CR',
                                  Cru = 'CR',
                                  `CMR and CR` = 'CR',
                                  PR = 'PR',
                                  SD = 'SD',
                                  PD = 'PD', 
                                  mors = 'death')) %>% 
  mutate(Aged70 = ifelse(Age > 70, '>70', '<70'),
         date_treatment_2nd_line_death = if_else(!is.na(date_treatment_2nd_line), date_treatment_2nd_line, date_death_fu),
         time_treatment_2nd_line = diff_years(date_chemo_end_1st_line, date_treatment_2nd_line_death))

MM_clean = MM_clean %>% 
  left_join(ALL.DISEASES %>% select(-date_diagnosis, -date_death_fu), 'patientid') %>% 
  mutate(date_treatment_2nd_line_death = if_else(!is.na(date_treatment_2nd_line_start), date_treatment_2nd_line_start, date_death_fu),
         time_treatment_2nd_line = diff_years(date_treatment_1st_line_start, date_treatment_2nd_line_death)) %>% 
  filter(treatment ==1,
         time_to_treatment <= 90)

if(cohort_1 %>% n_distinct() - DX.FIRST.ATC2 %>% n_patients() == DX.FIRST.POLY %>% filter(n.ATC==0) %>% n_patients()){
  cat('\nAll patients minus no. patients w/ prescriptions == no. w/o prescriptions\nYou`re good to go!')
}else{cat('\nSomething is wrong!\nRun line by line to detect error.')}

# save.image('/ngc/projects2/dalyca_r/chribr_r/PolyRx/data/Before_figures.RData')
# load('/ngc/projects2/dalyca_r/chribr_r/PolyRx/data/Before_coxmodels.RData')

# End of script
#### Go to: PolyRx_Figures_eop.R ####