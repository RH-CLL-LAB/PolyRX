## Calculates CCI score and polypharmacy
## Sys.Date() # "2024-04-25"
## author: christian brieghel

#### load data ####
load_dataset(c('SDS_ekokur', 'SDS_epikur', 'diagnoses_all'), ENTITIES$patientid)

Medicine_all = bind_rows(SDS_epikur_subset %>% transmute(patientid,
                                                      date = eksd,
                                                      atc,
                                                      source = 'LSR'),
                             SDS_ekokur_subset %>% transmute(patientid,
                                                      date = eksd,
                                                      atc,
                                                      source = 'LSR')) %>% 
  mutate(atc = gsub(':|_|\u008f\u008f|ÅÅ|XXX', '', atc),
         atc = gsub('Å', '', atc)) %>% 
  filter(! atc %in% c("N/A", "NUL", '', 'X')) # Time-lapse mins

date_max_atc = Medicine_all$date %>% max
ENTITIES2 = ENTITIES %>%  
  mutate(date_diagnosis < date_max_atc) # patients with prescription data
ENTITIES2 %>% nrow_npatients()
cohort_1 = ENTITIES2$patientid 

Medicine_all %>% nrow_npatients()
diagnoses_all_subset %>% nrow_npatients()
patient %>% nrow_npatients()

#### CCI ####
# Load "all_diagnoses" tables from DALY-CARE Core/public
diagnoses_all_subset %>% head2
ALL_ICD10_all = diagnoses_all_subset %>% 
  filter(patientid %in% cohort_1,
         date_diagnosis > as.Date('1970-01-01'),
         date_diagnosis <= date_max_atc) 
ALL_ICD10_all$date_diagnosis %>% summary
ALL_ICD10_all %>% nrow_npatients()

ALL_ICD10_all.1 = ALL_ICD10_all %>% 
  select(patientid, diagnosis) %>% 
  distinct() %>% #only distinct diagnoses 
  group_by(patientid) %>% 
  mutate(N = n()) %>% 
  slice(1) %>% 
  ungroup()

ALL_ICD10_all.1 %>% nrow_npatients()
ALL_ICD10_all.1$N %>%  summary
cat('Median no. of distinct diagnoses')
utable(~Q(N), ALL_ICD10_all.1) ## Suppl. result

DX.FIRST = ENTITIES2 %>% 
  group_by(patientid) %>% 
  arrange(date_diagnosis) %>% 
  slice(1) %>% 
  ungroup()
DX.FIRST %>% nrow_npatients()
ALL_ICD10_all %>% head2

DX.FIRST.CCI = DX.FIRST %>% 
  select(patientid, date_first_dx = date_diagnosis) %>% 
  left_join(ALL_ICD10_all %>% select(patientid, date_icd10 = date_diagnosis, icd10 = diagnosis)) %>% 
  filter(date_icd10 <= date_first_dx) %>% 
  CCI(patientid = patientid, icd10 = icd10, include_LC_score = T)

DX.FIRST.CCI$CCI.Cancer.and.Hem.score %>% unique() #Must be 2
DX.FIRST.CCI %>% names
utable(~ Q(CCI.2011.update), DX.FIRST.CCI)

#### PolyRX ####
Medicine_all.1 = Medicine_all %>%
  filter(patientid %in% cohort_1) %>% 
  distinct() %>%
  group_by(patientid) %>% 
  mutate(N = n()) %>% 
  slice(1) %>% 
  ungroup() %>% 
  right_join(patient, by = 'patientid') %>% 
  mutate(N = ifelse(is.na(N), 0, N))

## ANY ATC level
ggplot(Medicine_all) +
  geom_histogram(aes(date)) # atc coverage

DX.FIRST.ATC = DX.FIRST %>% 
  select(patientid, date_first_dx = date_diagnosis) %>% 
  filter(date_first_dx >= ymd('2002-01-01'), # t_dalycare_diagnosis first date
         date_first_dx <= as.Date(date_max_atc)) %>%  #LSR max date
  left_join(Medicine_all %>% distinct %>% select(patientid, date_atc = date, atc), 'patientid') %>% 
  distinct()

DX.FIRST.ATC %>% nrow # REPORT in Suppl Info #n.Rx
DX.FIRST.ATC$atc %>% n_distinct() #n.ATC among Rx
DX.FIRST.ATC %>% n_patients()  == n_distinct(cohort_1)
DX.FIRST.ATC %>% 
  left_join(ENTITIES2 %>% select(patientid, date_diagnosis), 'patientid') %>% 
  mutate(SAME = date_first_dx == date_diagnosis) %>% 
  pull(SAME) %>% 
  table # must be TRUE only

# LSR before first DX
DX.FIRST.ATC2 = DX.FIRST.ATC %>% 
  left_join(ENTITIES2 %>% select(patientid, date_diagnosis), 'patientid') %>% 
  mutate(Time = diff_days(date_first_dx, date_atc)) %>%
  filter(Time <= 0,
         Time > -365.25) %>% #0 + 365.25 
  left_join(DX.FIRST %>% select(patientid, date = date_diagnosis) %>% 
               filter(date <= as.Date(date_max_atc)), 'patientid') %>% 
  mutate(SAME = date_first_dx == date_diagnosis)

DX.FIRST.ATC2_30 = DX.FIRST.ATC %>% 
  left_join(ENTITIES2 %>% select(patientid, date_diagnosis), 'patientid') %>% 
  mutate(Time = diff_days(date_first_dx, date_atc)) %>% 
  filter(Time <= -30,
         Time > -365.25) %>% #as sensitivity analysis
  left_join(DX.FIRST %>% select(patientid, date = date_diagnosis) %>% 
              filter(date <= as.Date(date_max_atc)), 'patientid') %>% 
  mutate(SAME = date_first_dx == date_diagnosis)

DX.FIRST.ATC2$SAME %>% table # must be TRUE
DX.FIRST.ATC2 %>% nrow_npatients()

DX.FIRST.POLY = DX.FIRST.ATC2 %>% 
  group_by(patientid, atc) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(patientid) %>% 
  mutate(n.ATC = n()) %>% 
  slice(1) %>% 
  ungroup() %>% 
  right_join(DX.FIRST %>% select(patientid, date = date_diagnosis) %>% 
               filter(date <= as.Date(date_max_atc)), 'patientid') %>% 
  mutate(n.ATC = ifelse(is.na(n.ATC), 0, n.ATC), 
         Polypharmacy = ifelse(n.ATC < 5, 'No', 'Yes'))

# Sensititvity excluding 30 days up to diagnosis
DX.FIRST.POLY_30 = DX.FIRST.ATC2_30 %>% 
  group_by(patientid, atc) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(patientid) %>% 
  mutate(n.ATC = n()) %>% 
  slice(1) %>% 
  ungroup() %>% 
  right_join(DX.FIRST %>% select(patientid, date = date_diagnosis) %>% 
               filter(date <= as.Date(date_max_atc)), 'patientid') %>% 
  mutate(n.ATC = ifelse(is.na(n.ATC), 0, n.ATC), 
         Polypharmacy = ifelse(n.ATC < 5, 'No', 'Yes'))

DX.FIRST.POLY %>% nrow_npatients()
DX.FIRST.POLY$n.ATC %>% table(exclude = NULL)
utable(~ Q(n.ATC) + Polypharmacy, DX.FIRST.POLY) # report Suppl. Info
DX.FIRST.POLY %>% n_patients()
DX.FIRST.ATC2 %>% n_patients() #those w/o prescription not included
DX.FIRST.CCI %>% n_patients()

# write_csv(ENTITIES2, '/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/ENTITIES2.csv')
# write_csv(DX.FIRST.ATC, '/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/DX.FIRST.ATC.csv')
# write_csv(DX.FIRST.POLY, '/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/DX_FIRST_POLY.csv')
# write_csv(DX.FIRST.POLY_30, '/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/DX_FIRST_POLY_30.csv')
# write_csv(DX.FIRST.ATC2, '/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/DX_FIRST_ATC2.csv')
# write_csv(DX.FIRST.CCI, '/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/DX_FIRST_CCI.csv')


# SAVE = TRUE
if(SAVE){
  save.image(paste0(getwd(), 'After_all_data_is_loaded.RData')) 
}
# load(paste0(getwd(), 'After_all_data_is_loaded.RData')) #setwd()

clear_ram()

# End of script
#### Go to RolyRx_data_prep_eop.R ####