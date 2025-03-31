## This script plots figures and tables in PolyRx project 
## created Sys.Date() # "2024-09-05"
## author: christian brieghel 

##### FIGURE 1: CONSORT numbers #####
ENTITIES %>% nrow_npatients()
ENTITIES %>% n_patients() - ALL.DISEASES %>%  n_patients() # Step 1 excluded

ALL.DISEASES %>%  nrow_npatients()
ALL.DISEASES$Disease %>%  table()
ALL.DISEASES %>% n_patients() - ALL.DISEASES.COX %>%  n_patients() # Step 2 excluded

ALL.DISEASES.COX %>% n_patients()
cohort_2 %>% n_distinct()
ALL.DISEASES.COX$Disease %>% table()

FINAL_COHORT_ADT.1.DURATION %>% n_patients()
cohort_3 %>% n_distinct()
FINAL_COHORT_ADT.1.DURATION$Disease %>% table
ALL.DISEASES %>% n_patients() - FINAL_COHORT_ADT.1.DURATION %>%  n_patients() # Step 2 excluded

cohort_4 %>% n_distinct()
COX_ADT %>% n_patients()
COX_ADT$Disease %>% table
#Figure created with Lucid chart

#### FIGURE 2 ####
palette(BLOOD)
PLOT.2.A.OS = KM_plot(survfit(Surv(time_dx_death, status) ~ ATC.all.cut, data = ALL.DISEASES),
                      breaks = 2,
                      title = 'All',
                      labs = c('0-3', '4-7', '8-11', '>11'),
                      xlab = 'Time (years since diagnosis)',
                      xlim = c(0,10),
                      palette = c(1:4))

PLOT.2.B.ADMISSIONS = KM_plot(survfit(Surv(time_golive_ADT_death, ADT) ~ ATC.all.cut, FINAL_COHORT_ADT.1.DURATION),
                              fun = 'event',
                              title = 'Eastern Denmark',
                              labs = c('0-3', '4-7', '8-11', '>11'),
                              ylim = c(0,1),
                              ylab = '% admitted',
                              palette = BLOOD[1:4])

PLOT.2.C.INFECTIONS = KM_plot(survfit(Surv(time_death_fu_truncated, event) ~ ATC.all.cut, 
                                      ALL.DISEASES_INF.1),
                              fun = 'event',
                              title = 'Eastern Denmark',
                              labs = c('0-3', '4-7', '8-11', '>11'),
                              ylab = '% infected',
                              palette = BLOOD)



ggsave('/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/Figure_2ABC.png',
       arrange_ggsurvplots(list(PLOT.2.A.OS, PLOT.2.B.ADMISSIONS, PLOT.2.C.INFECTIONS), nrow = 1,  ncol = 3),
       height = 18,
       width = 42,
       units = 'cm',
       dpi = 300)

#### FIGURE 3 ####
COX_ADT = ALL.DISEASES %>%
  right_join(FINAL_COHORT_ADT.1.DURATION %>% select(patientid, time_golive_ADT_death, ADT)) %>% 
  filter(!is.na(IPI))

COX_INF = ALL.DISEASES %>%
  right_join(ALL.DISEASES_INF.1.ATC %>% select(patientid, time_death_fu_truncated, event)) %>% 
  filter(!is.na(IPI))

COX_ADT %>% n_patients()
COX_INF %>% n_patients()

PLOT.OS = ggforest(coxph(Surv(time_dx_death, status) ~ Age + Sex + IPI + ATC + CCI, ALL.DISEASES.COX %>% 
                           dplyr::rename(ATC = ATC.all.cut) %>% 
                           as.data.frame()), cpositions = c(0.02, 0.15, 0.3))

PLOT.ADT = ggforest(coxph(Surv(time_golive_ADT_death, ADT) ~ Age + Sex + IPI + ATC + CCI, data = COX_ADT %>% 
                            dplyr::rename(ATC = ATC.all.cut) %>% 
                            as.data.frame()), cpositions = c(0.02, 0.15, 0.3))


PLOT.INF = ggforest(coxph(Surv(time_death_fu_truncated, event) ~ Age + Sex + IPI + ATC + CCI, data = COX_INF %>% 
                            dplyr::rename(ATC = ATC.all.cut) %>% 
                            as.data.frame()), cpositions = c(0.02, 0.15, 0.3))

ggsave( '/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/Figure_3_COX_ABC.png',
        ggarrange(PLOT.OS, PLOT.ADT, PLOT.INF, nrow = 1, ncol = 3, labels = LETTERS[1:3]),
        height = 10,
        width = 16,
        dpi = 300)

#### FIGURE 4 ####
BOX.ADT.PLOT = ggplot(data = FINAL_COHORT_ADT.1.DURATION) + 
  geom_boxplot(aes(x = ATC.all.cut, y = N.Contacts.per.year, fill= ATC.all.cut)) +
  scale_fill_manual('No. of medications', values = BLOOD[1:4]) +
  theme_classic() +
  labs(x = 'No. of medications',
       y = 'Hospital admissions per patient-year') +
  ylim(0,10)

utable(ATC.all.cut ~ Q(SUM.in.hospital), FINAL_COHORT_ADT.1.DURATION.admitted %>% mutate(SUM.in.hospital =  round(SUM.in.hospital, 0))) ## REPORT!
BOX.ADT.PLOT2 = ggplot(data = FINAL_COHORT_ADT.1.DURATION, aes(x = ATC.all.cut, y = SUM.in.hospital)) + 
  geom_boxplot(aes(fill= ATC.all.cut)) +
  scale_fill_manual('No. of medications', values = BLOOD[1:4]) +
  theme_classic() +
  labs(x = 'No. of medications',
       y = 'Days in-hospital per patient') +
  ylim(0,100)


SP.AB2.1 = SP.INF %>% 
  group_by(patientid) %>% 
  summarize(n_infections = n(), 
            N_days_on_AB = sum(n_days_IVAB)) %>% 
  ungroup() %>% 
  right_join(FINAL_COHORT_ADT.1.DURATION %>% filter(patientid %in% cohort_3) %>% select(patientid) %>% distinct(), 'patientid') %>% 
  mutate(n_infections = ifelse(is.na(n_infections), 0, n_infections), 
         N_days_on_AB = ifelse(is.na(N_days_on_AB), 0, N_days_on_AB)) %>% 
  left_join(ALL.DISEASES, 'patientid')

utable(ATC.all.cut ~ Q(n_infections), SP.AB2.1 %>% filter(n_infections >0))
utable(ATC.all.cut ~ n_infections, SP.AB2.1 %>% filter(n_infections >0))

ALL.DISEASES_INF.1 %>% nrow_npatients()

BOX.INF.PLOT1 = ggplot(data = SP.AB2.1 %>% filter(n_infections >0), aes(x = ATC.all.cut, y = n_infections)) + 
  geom_boxplot(fill= BLOOD[1:4]) +
  theme_classic() +
  labs(x = 'No. of medications',
       y = 'Infection per patient (days)') +
  ylim(0,20)

utable(ATC.all.cut ~ N_days_on_AB, SP.AB2.1)
utable(ATC.all.cut ~ Q(N_days_on_AB), SP.AB2.1 %>% filter(N_days_on_AB >0))
utable(ATC.all.cut ~ Q(N_days_on_AB), SP.AB2.1 %>% filter(N_days_on_AB >0)) # report
BOX.INF.PLOT2 = ggplot(data = SP.AB2.1 %>% filter(N_days_on_AB >0), aes(x = ATC.all.cut, y = N_days_on_AB)) + 
  geom_boxplot(fill= BLOOD[1:4]) +
  theme_classic() +
  labs(x = 'No. of medications',
       y = 'Days on antimicrobial therapy per patient') +
  ylim(0,50)

BOX.INF.PLOT3 = ggplot(data = ALL.DISEASES_INF.1, aes(x = ATC.all.cut, y = n_inf_per_year)) + 
  geom_boxplot(fill= BLOOD[1:4]) +
  theme_classic() +
  labs(x = 'No. of medications',
       y = 'Infections per patient-year (n)') +
  ylim(0,5)

ggsave('/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/Figure_4_Infduration.png',
       ggarrange(BOX.ADT.PLOT, BOX.ADT.PLOT2, 
                 BOX.INF.PLOT3, BOX.INF.PLOT2, common.legend = T, ncol = 2, nrow = 2, labels = LETTERS[1:4]),
       height = 25,
       width = 22,
       units = 'cm',
       dpi = 300)

#### FIGURE 5 ####
PLOT.R.DLBCL = KM_plot(survfit(Surv(time_treatment_2nd_line, relapse_treatment) ~ ATC.all.cut, 
                               LYFO_clean %>% filter(Disease == 'DLBCL',
                                                     time_treatment_2nd_line >= 0,
                                                     !is.na(time_treatment_2nd_line))),
                       fun = 'event',
                       breaks = 2,
                       xlim = c(0,10),
                       pval = T,
                       xlab = 'Time (years from EoFLT)',
                       ylab = '% relapsed',
                       palette = BLOOD,
                       pval.coord = c(0.2, 0.95),
                       title = 'DLBCL',
                       labs = c('0-3', '4-7', '8-11', '>11')) 

PLOT.R.HL = KM_plot(survfit(Surv(time_treatment_2nd_line, relapse_treatment) ~ ATC.all.cut, 
                            LYFO_clean %>% filter(Disease == 'cHL',
                                                  time_treatment_2nd_line >= 0,
                                                  !is.na(time_treatment_2nd_line))),
                    fun = 'event',
                    breaks = 2,
                    xlim = c(0,10),
                    pval = T,
                    xlab = 'Time (years from EoFLT)',
                    ylab = '% relapsed',
                    palette = BLOOD,
                    pval.coord = c(0.2, 0.95),
                    title = 'cHL',
                    labs = c('0-3', '4-7', '8-11', '>11')) 

PLOT.R.MCL = KM_plot(survfit(Surv(time_treatment_2nd_line, relapse_treatment) ~ ATC.all.cut, 
                            LYFO_clean %>% filter(Disease == 'MCL',
                                                  time_treatment_2nd_line >= 0,
                                                  !is.na(time_treatment_2nd_line))),
                    fun = 'event',
                    breaks = 2,
                    xlim = c(0,10),
                    pval = T,
                    xlab = 'Time (years from EoFLT)',
                    ylab = '% relapsed',
                    palette = BLOOD,
                    pval.coord = c(0.2, 0.95),
                    title = 'MCL',
                    labs = c('0-3', '4-7', '8-11', '>11')) 

PLOT.R.MM = KM_plot(survfit(Surv(time_treatment_2nd_line, treatment_2nd_line) ~ ATC.all.cut, 
                             MM_clean %>% filter(time_treatment_2nd_line >= 0,
                                                 !is.na(time_treatment_2nd_line),
                                                 !is.na(ATC.all.cut))),
                     fun = 'event',
                     breaks = 2,
                     xlim = c(0,10),
                     pval = T,
                     xlab = 'Time (years from EoFLT)',
                     ylab = '% relapsed',
                     palette = BLOOD,
                     pval.coord = c(0.2, 0.95),
                     title = 'MM',
                     labs = c('0-3', '4-7', '8-11', '>11')) 

ggsave('/ngc/projects2/dalyca_r/chribr_r/PolyRx/eop_figures/Figure_5_TTNT.png',
       arrange_ggsurvplots(list(PLOT.R.HL, PLOT.R.MCL, PLOT.R.DLBCL, PLOT.R.MM), nrow = 2,  ncol = 2),
       height = 22,
       width = 17,
       dpi = 300)

#### Figure S-hospitalization-UNI-subtypes ####
PLOT.3.B.ADMISSIONS_SI = KM_plot(survfit(Surv(time_golive_ADT_death, ADT) ~ ATC.all.cut, FINAL_COHORT_ADT.1.DURATION),
                              fun = 'event',
                              title = 'All',
                              labs = c('0-3', '4-7', '8-11', '>11'),
                              ylim = c(0,1),
                              ylab = '% admitted',
                              palette = BLOOD[1:4])


list_KM_plots_hospitalization = list()
for(i in 1:length(diseases)){
  list_KM_plots_hospitalization[[i]] = KM_plot(survfit(Surv(time_golive_ADT_death, ADT) ~ ATC.all.cut, FINAL_COHORT_ADT.1.DURATION %>% 
                                                         filter(Disease ==diseases[i])),
                                                                     fun = 'event',
                                                                     title = diseases[i],
                                                                     labs = c('0-3', '4-7', '8-11', '>11'),
                                                                     ylab = '% admitted',
                                                                     palette = BLOOD[1:4])
}


ggsave(paste0(getwd(), '/Figure_S_hospitalization-UNI-subtypes.png'),
       arrange_ggsurvplots(list(PLOT.3.B.ADMISSIONS_SI, list_KM_plots_hospitalization[[3]], list_KM_plots_hospitalization[[6]],
                                list_KM_plots_hospitalization[[1]], list_KM_plots_hospitalization[[4]], list_KM_plots_hospitalization[[7]],
                                list_KM_plots_hospitalization[[2]], list_KM_plots_hospitalization[[5]], list_KM_plots_hospitalization[[8]]),
                           nrow  = 3,
                           ncol = 3),
       dpi = 300,
       height = 40,
       width = 40,
       units = 'cm')

utable(ATC.all.cut ~ Q(N.Contacts.per.year), FINAL_COHORT_ADT.1.DURATION)



##### Figure S-time to infection #####
list_KM_plots = list()
for (i in 1:length(diseases)) {print(diseases[i])
  list_KM_plots[[i]] = KM_plot(survfit(Surv(time_death_fu_truncated, event) ~ ATC.all.cut, 
                                       ALL.DISEASES_INF.1 %>% filter(Disease == diseases[i])),
                               fun = 'event',
                               title = diseases[i],
                               labs = c('0-3', '4-7', '8-11', '>11'),
                               ylab = '% infected',
                               palette = BLOOD)
}

PLOT.INF.ALL_SI = KM_plot(survfit(Surv(time_death_fu_truncated, event) ~ ATC.all.cut, 
                                  ALL.DISEASES_INF.1),
                          fun = 'event',
                          title = 'All',
                          labs = c('0-3', '4-7', '8-11', '>11'),
                          ylab = '% infected',
                          palette = BLOOD)

ggsave(paste0(getwd(), '/figure_S_time_to_infection.png'),
       arrange_ggsurvplots(list(PLOT.INF.ALL_SI, list_KM_plots[[3]], list_KM_plots[[6]], 
                                list_KM_plots[[1]], list_KM_plots[[4]], list_KM_plots[[7]], 
                                list_KM_plots[[2]], list_KM_plots[[5]], list_KM_plots[[8]]), nrow = 3, ncol = 3),
       height = 30,
       width = 40,
       dpi = 300)

#### Figure S-Excluding Remdesivir ######
SP.AB2.C19 = SP.AB %>% 
  filter(atc != 'J05AB16',
         AB.DETAILED != 'REMDESIVIR',
         zc_admin_route_name  == 'Intravenøs anvendelse') %>%
  AE_infection()

ALL.DISEASES_INF.C19 = ALL.DISEASES %>% 
  filter(patientid %in% cohort_3) %>% 
  left_join(SP.AB2.C19 %>% select(patientid, date_infection, DIFF, n_inf, n_days_IVAB, date_received), 'patientid')

ALL.DISEASES_INF.C19 %>% nrow_npatients()
SP.INF$date_received %>% table() # == date_last_fu**

ALL.DISEASES_INF.C19.1 = ALL.DISEASES_INF.C19 %>% 
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
                   date_truncation = date_max_atc) %>%  ##**
  mutate(n_inf_per_year = n_inf/time_dx_death,
         n_inf_per_year = ifelse(is.na(n_inf_per_year), 0, n_inf_per_year)) %>% 
  filter(time_death_fu_truncated >= 0)

diseases
list_KM_plots.C19 = list()
for (i in 1:length(diseases)) {print(diseases[i])
  list_KM_plots.C19[[i]] = KM_plot(survfit(Surv(time_death_fu_truncated, event) ~ ATC.all.cut, 
                                       ALL.DISEASES_INF.C19.1 %>% filter(Disease == diseases[i])),
                               fun = 'event',
                               title = diseases[i],
                               labs = c('0-3', '4-7', '8-11', '>11'),
                               ylab = '% infected',
                               palette = BLOOD)
}


PLOT.3.C.INFECTIONS.C19 = KM_plot(survfit(Surv(time_death_fu_truncated, event) ~ ATC.all.cut, 
                                      ALL.DISEASES_INF.C19.1),
                              fun = 'event',
                              title = 'Eastern Denmark',
                              labs = c('0-3', '4-7', '8-11', '>11'),
                              ylab = '% infected',
                              palette = BLOOD)


COX_INF.C19 = ALL.DISEASES %>% 
  right_join(ALL.DISEASES_INF.1.ATC %>% select(patientid, time_death_fu_truncated, event)) %>% 
  filter(!is.na(IPI))

utable(~Q(Age)+Sex+IPI+ATC.all.cut+CCI, COX_INF.C19)
PLOT.INF.C19 = ggforest(coxph(Surv(time_death_fu_truncated, event) ~ Age + Sex + IPI + ATC + CCI, 
                          data = COX_INF.C19 %>% 
                            dplyr::rename(ATC = ATC.all.cut) %>% 
                            as.data.frame()), cpositions = c(0.02, 0.15, 0.3))

ggsave(paste0(getwd(), '/Figure_S-remdesivir-sensitivityA.png'),
        arrange_ggsurvplots(list(PLOT.3.C.INFECTIONS.C19)),
        height = 10,
        width = 7,
        dpi = 300)
ggsave(paste0(getwd(), '/Figure_S-remdesivir-sensitivityB.png'),
        PLOT.INF.C19,
        height = 10,
        width = 16,
        dpi = 300)
# gather in image editor as:
ggarrange(PLOT.INF.C19, PLOT.INF.C19)

#### Figure S-Cox Admission ####

PLOT.ADT = ggforest(coxph(Surv(time_golive_ADT_death, ADT) ~ Age + Sex + IPI + ATC + CCI, 
                          data = COX_ADT %>% 
                            dplyr::rename(ATC = ATC.all.cut) %>% 
                            as.data.frame()), cpositions = c(0.02, 0.15, 0.3)) 
# CLL
PLOT.ADT.CLL = ggforest(coxph(Surv(time_golive_ADT_death, ADT) ~ Sex + `CLL-IPI` + ATC + CCI, 
                          data = COX_ADT %>% 
                            filter(Disease == 'CLL') %>% 
                            dplyr::rename(`CLL-IPI` = IPI2, ATC = ATC.all.cut) %>% 
                            as.data.frame()), cpositions = c(0.02, 0.15, 0.3))
#DLBCL
PLOT.ADT.DLBCL = ggforest(coxph(Surv(time_golive_ADT_death, ADT) ~ Sex +`R-IPI` + ATC + CCI, 
                                data = COX_ADT %>% 
                                 filter(Disease == 'DLBCL') %>% 
                                 dplyr::rename(`R-IPI` = IPI, ATC = ATC.all.cut) %>% 
                                 as.data.frame()), cpositions = c(0.02, 0.15, 0.3))

#FL
PLOT.ADT.FL = ggforest(coxph(Surv(time_golive_ADT_death, ADT) ~ Sex + FLIPI2 + ATC + CCI, 
                             data = COX_ADT %>% 
                               filter(Disease == 'FL') %>% 
                           dplyr::rename(FLIPI2 = IPI, ATC = ATC.all.cut)%>% 
                           as.data.frame()), cpositions = c(0.02, 0.15, 0.3))

#HL
PLOT.ADT.HL = ggforest(coxph(Surv(time_golive_ADT_death, ADT) ~ IPS + ATC + CCI, 
                             data = COX_ADT %>% 
                               filter(Disease == 'cHL') %>% 
                           dplyr::rename(ATC = ATC.all.cut)%>% 
                           as.data.frame()), cpositions = c(0.02, 0.15, 0.3)) 

#MCL
PLOT.ADT.MCL = ggforest(coxph(Surv(time_golive_ADT_death, ADT) ~ Sex + MIPI + ATC + CCI, 
                              data = COX_ADT %>% 
                                filter(Disease == 'MCL') %>% 
                            dplyr::rename(MIPI= IPI, ATC = ATC.all.cut)%>% 
                            as.data.frame()), cpositions = c(0.02, 0.15, 0.3)) 

#MZL
PLOT.ADT.MZL = ggforest(coxph(Surv(time_golive_ADT_death, ADT) ~ Sex +  `MALT-IPI` + ATC + CCI, 
                              data = COX_ADT %>% 
                                filter(Disease == 'MZL') %>% 
                            dplyr::rename(`MALT-IPI` = IPI, ATC = ATC.all.cut)%>% 
                            as.data.frame()), cpositions = c(0.02, 0.15, 0.3)) # Try removing AB from ATC?

#LPL
PLOT.ADT.LPL = ggforest(coxph(Surv(time_golive_ADT_death, ADT) ~ Sex + IPSSWM + ATC + CCI,
                              data = COX_ADT %>% 
                                filter(Disease == 'LPL') %>% 
                            dplyr::rename(IPSSWM = IPI, ATC = ATC.all.cut)%>% 
                            as.data.frame()), cpositions = c(0.02, 0.15, 0.3)) 

#MM
PLOT.ADT.MM = ggforest(coxph(Surv(time_golive_ADT_death, ADT) ~ Age + Sex + `R-ISS` + ATC + CCI, 
                             data = COX_ADT %>% 
                               filter(Disease == 'MM') %>% 
                           dplyr::rename(ATC = ATC.all.cut, `R-ISS` = IPI)%>% 
                           as.data.frame()), cpositions = c(0.02, 0.15, 0.3)) 

LABELS = c('All', 'cHL', 'DLBCL', 'FL', 'MZL', 'MCL', 'CLL', 'LPL', 'MM')
LABELS %>% length()

ggsave(paste0(getwd(), '/Figure_S-Cox-ADT.png'),
         ggarrange(PLOT.ADT, PLOT.ADT.HL, PLOT.ADT.DLBCL, PLOT.ADT.FL, PLOT.ADT.MZL, PLOT.ADT.MCL, PLOT.ADT.CLL,  PLOT.ADT.LPL, PLOT.ADT.MM,
                   nrow = 3, ncol = 3, labels = LABELS),
         height = 15,
         width = 22,
         dpi = 300)

#### Figure S-Cox Infection ####
PLOT.INF = ggforest(coxph(Surv(time_death_fu_truncated, event) ~ Age + Sex + IPI + ATC + CCI, 
                          data = COX_INF %>% 
                            dplyr::rename(ATC = ATC.all.cut) %>% 
                            as.data.frame()), cpositions = c(0.02, 0.15, 0.3))

# CLL
PLOT.INF.CLL = ggforest(coxph(Surv(time_death_fu_truncated, event) ~ Sex + `CLL-IPI` + ATC + CCI, 
                              data = COX_INF %>% 
                                filter(Disease == 'CLL') %>% 
                                dplyr::rename(`CLL-IPI` = IPI2, ATC = ATC.all.cut) %>% 
                                as.data.frame()), cpositions = c(0.02, 0.15, 0.3))
#DLBCL
PLOT.INF.DLBCL = ggforest(coxph(Surv(time_death_fu_truncated, event) ~ Sex +`R-IPI` + ATC + CCI, 
                                data = COX_INF %>% 
                                  filter(Disease == 'DLBCL') %>% 
                                  dplyr::rename(`R-IPI` = IPI, ATC = ATC.all.cut) %>% 
                                  as.data.frame()), cpositions = c(0.02, 0.15, 0.3))

#FL
PLOT.INF.FL = ggforest(coxph(Surv(time_death_fu_truncated, event) ~ Sex + FLIPI2 + ATC + CCI, 
                             data = COX_INF %>% 
                               filter(Disease == 'FL') %>% 
                               dplyr::rename(FLIPI2 = IPI, ATC = ATC.all.cut)%>% 
                               as.data.frame()), cpositions = c(0.02, 0.15, 0.3))

#HL
PLOT.INF.HL = ggforest(coxph(Surv(time_death_fu_truncated, event) ~ IPS + ATC + CCI, 
                             data = COX_INF %>% 
                               filter(Disease == 'cHL') %>% 
                               dplyr::rename(ATC = ATC.all.cut)%>% 
                               as.data.frame()), cpositions = c(0.02, 0.15, 0.3)) 

#MCL
PLOT.INF.MCL = ggforest(coxph(Surv(time_death_fu_truncated, event) ~ Sex + MIPI + ATC + CCI, 
                              data = COX_INF %>% 
                                filter(Disease == 'MCL') %>% 
                                dplyr::rename(MIPI= IPI, ATC = ATC.all.cut)%>% 
                                as.data.frame()), cpositions = c(0.02, 0.15, 0.3)) 

#MZL
PLOT.INF.MZL = ggforest(coxph(Surv(time_death_fu_truncated, event) ~ Sex +  `MALT-IPI` + ATC + CCI, 
                              data = COX_INF %>% 
                                filter(Disease == 'MZL') %>% 
                                dplyr::rename(`MALT-IPI` = IPI, ATC = ATC.all.cut)%>% 
                                as.data.frame()), cpositions = c(0.02, 0.15, 0.3))

#LPL
PLOT.INF.LPL = ggforest(coxph(Surv(time_death_fu_truncated, event) ~ Sex + IPSSWM + ATC + CCI,
                              data = COX_INF %>% 
                                filter(Disease == 'LPL') %>% 
                                dplyr::rename(IPSSWM = IPI, ATC = ATC.all.cut)%>% 
                                as.data.frame()), cpositions = c(0.02, 0.15, 0.3)) 

#MM
PLOT.INF.MM = ggforest(coxph(Surv(time_death_fu_truncated, event) ~ Age + Sex + `R-ISS` + ATC + CCI, 
                             data = COX_INF %>% 
                               filter(Disease == 'MM') %>% 
                               dplyr::rename(ATC = ATC.all.cut, `R-ISS` = IPI)%>% 
                               as.data.frame()), cpositions = c(0.02, 0.15, 0.3)) 

LABELS = c('All', 'cHL', 'DLBCL', 'FL', 'MZL', 'MCL', 'CLL', 'LPL', 'MM')
LABELS %>% length()

ggsave(paste0(getwd(), '/Figure_S-Cox-INF.png'),
       ggarrange(PLOT.INF, PLOT.INF.HL, PLOT.INF.DLBCL, PLOT.INF.FL, PLOT.INF.MZL, PLOT.INF.MCL, PLOT.INF.CLL,  PLOT.INF.LPL, PLOT.INF.MM,
                 nrow = 3, ncol = 3, labels = LABELS),
       height = 15,
       width = 22,
       dpi = 300)

##### Figure S-PolyRx Cox #####
PLOT.OS.PolyRX = ggforest(coxph(Surv(time_dx_death, status) ~ Age + Sex + IPI + Polypharmacy + CCI, ALL.DISEASES.COX %>%
                           dplyr::rename(ATC = ATC.all.cut) %>% 
                           as.data.frame()), cpositions = c(0.02, 0.15, 0.3))

PLOT.ADT.PolyRX = ggforest(coxph(Surv(time_golive_ADT_death, ADT) ~ Age + Sex + IPI + Polypharmacy + CCI, data = COX_ADT %>% 
                            dplyr::rename(ATC = ATC.all.cut) %>% 
                            as.data.frame()), cpositions = c(0.02, 0.15, 0.3))


PLOT.INF.PolyRX = ggforest(coxph(Surv(time_death_fu_truncated, event) ~ Age + Sex + IPI + Polypharmacy + CCI, data = COX_INF %>% 
                            dplyr::rename(ATC = ATC.all.cut) %>% 
                            as.data.frame()), cpositions = c(0.02, 0.15, 0.3))

ggsave(paste0(getwd(), '/Figure_S_PolyRX_COX_ABC.png'),
        ggarrange(PLOT.OS.PolyRX, PLOT.ADT.PolyRX, PLOT.INF.PolyRX, nrow = 1, ncol = 3, labels = LETTERS[1:3]),
        height = 10,
        width = 16,
        dpi = 300)

##### Figure S-PolyRX uni #####

PLOT_PolyRX_OS = KM_plot(survfit(Surv(time_dx_death, status) ~ Polypharmacy, 
                                 ALL.DISEASES),
                         title = 'Polypharmacy',
                         labs = c('No', 'Yes'),
                         ylab = '% infected',
                         pval = paste0('HR ', publish(coxph(Surv(time_dx_death, status) ~ Polypharmacy, 
                                                            ALL.DISEASES))$regressionTable$HazardRatio[2], ' ',
                                       publish(coxph(Surv(time_dx_death, status) ~ Polypharmacy, 
                                                     ALL.DISEASES))$regressionTable$CI.95[2]),
                         palette = BLOOD)

PLOT_PolyRX_ADT = KM_plot(survfit(Surv(time_golive_ADT_death, ADT) ~ Polypharmacy, FINAL_COHORT_ADT.1.DURATION),
                          fun = 'event',
                          title = 'Polypharmacy',
                          labs = c('No', 'Yes'),
                          ylab = '% admitted',
                          pval = paste0('HR ', publish(coxph(Surv(time_golive_ADT_death, ADT) ~ Polypharmacy, FINAL_COHORT_ADT.1.DURATION))$regressionTable$HazardRatio[2], ' ',
                                        publish(coxph(Surv(time_golive_ADT_death, ADT) ~ Polypharmacy, FINAL_COHORT_ADT.1.DURATION))$regressionTable$CI.95[2]),
                          palette = BLOOD)

PLOT_PolyRX_INF = KM_plot(survfit(Surv(time_death_fu_truncated, event) ~ Polypharmacy, 
                                  ALL.DISEASES_INF.1),
                          fun = 'event',
                          title = 'Polypharmacy',
                          labs = c('No', 'Yes'),
                          ylab = '% infected',
                          pval = paste0('HR ', publish(coxph(Surv(time_death_fu_truncated, event) ~ Polypharmacy, 
                                                             ALL.DISEASES_INF.1))$regressionTable$HazardRatio[2], ' ',
                                        publish(coxph(Surv(time_death_fu_truncated, event) ~ Polypharmacy, 
                                                      ALL.DISEASES_INF.1))$regressionTable$CI.95[2]),
                          palette = BLOOD)

ggsave(paste0(getwd(), '/Figure_S_polyRx.png'),
       arrange_ggsurvplots(list(PLOT_PolyRX_OS, PLOT_PolyRX_ADT, PLOT_PolyRX_INF), nrow = 1,  ncol = 3),
       height = 18,
       width = 42,
       units = 'cm',
       dpi = 300)

#### Figure S-30 day sensitivity ####
PATIENT_OS2_30 = patient %>% 
  left_join(DX.FIRST %>% select(patientid, date_first_dx = date_diagnosis)) %>% 
  left_join(DX.FIRST.POLY_30 %>% # loaded in prep_data.R
              transmute(patientid, n.ATC, Polypharmacy = factor(Polypharmacy, levels = c('No', 'Yes'))),
            by = 'patientid') %>% 
  left_join(DX.FIRST.CCI %>% transmute(patientid, CCI.score, CCI.2011.update), 'patientid') %>%
  left_join(t_dalycare_last_diagnosis %>% select(patientid, date_last_diagnosis)) %>% 
  mutate(date_death_fu2 = if_else(date_death_fu < date_first_dx, date_last_diagnosis, date_death_fu)) %>% 
  mutate(Time = diff_years(date_first_dx, date_death_fu2))  %>%
  filter(Time >= 0)


#### Figure S-NMI ####
NMI_data = ALL.DISEASES %>% 
  NMI()
NMI_data %>% nrow_npatients()
NMI_data$NMI_score %>% table
ALL.DISEASES.COX_NMI = ALL.DISEASES.COX %>% 
  left_join(NMI_data %>% select(patientid, NMI_score))
COX_ADT_NMI = COX_ADT %>% 
  left_join(NMI_data %>% select(patientid, NMI_score))
COX_INF_NMI = COX_INF %>% 
  left_join(NMI_data %>% select(patientid, NMI_score))

PLOT.OS.NMI = ggforest(coxph(Surv(time_dx_death, status) ~ Age + Sex + IPI + ATC + `NMI score`, ALL.DISEASES.COX_NMI %>% 
                               dplyr::rename(ATC = ATC.all.cut,
                                             `NMI score` = NMI_score) %>% 
                               as.data.frame()), cpositions = c(0.02, 0.15, 0.3))

PLOT.ADT.NMI = ggforest(coxph(Surv(time_golive_ADT_death, ADT) ~ Age + Sex + IPI + ATC + `NMI score`, data = COX_ADT_NMI %>% 
                                dplyr::rename(ATC = ATC.all.cut,
                                              `NMI score` = NMI_score) %>% 
                                as.data.frame()), cpositions = c(0.02, 0.15, 0.3))


PLOT.INF.NMI = ggforest(coxph(Surv(time_death_fu_truncated, event) ~ Age + Sex + IPI + ATC + `NMI score`, data = COX_INF_NMI %>% 
                                dplyr::rename(ATC = ATC.all.cut,
                                              `NMI score` = NMI_score) %>% 
                                as.data.frame()), cpositions = c(0.02, 0.15, 0.3))

ggsave( paste0(getwd(), '/Figure_S_COX_NMI.png'),
        ggarrange(PLOT.OS.NMI, PLOT.ADT.NMI, PLOT.INF.NMI, nrow = 1, ncol = 3, labels = LETTERS[1:3]),
        height = 10,
        width = 16,
        dpi = 300)


#### Figure S-CLL_CI score ####
CLL_CI_data = ALL.DISEASES %>% 
  CLL_CI()

ALL.DISEASES.COX_CI_score = ALL.DISEASES.COX %>% 
  left_join(CLL_CI_data %>% select(patientid, CLL_CI_score))
COX_ADT_CI_score = COX_ADT %>% 
  left_join(CLL_CI_data %>% select(patientid, CLL_CI_score))
COX_INF_CI_score = COX_INF %>% 
  left_join(CLL_CI_data %>% select(patientid, CLL_CI_score))

PLOT.OS.CI = ggforest(coxph(Surv(time_dx_death, status) ~ Age + Sex + IPI + ATC + `CLL CI score`, ALL.DISEASES.COX_CI_score %>% 
                              dplyr::rename(ATC = ATC.all.cut,
                                            `CLL CI score` = CLL_CI_score) %>% 
                              as.data.frame()), cpositions = c(0.02, 0.15, 0.3))

PLOT.ADT.CI = ggforest(coxph(Surv(time_golive_ADT_death, ADT) ~ Age + Sex + IPI + ATC + `CLL CI score`, data = COX_ADT_CI_score %>% 
                               dplyr::rename(ATC = ATC.all.cut,
                                             `CLL CI score` = CLL_CI_score) %>% 
                               as.data.frame()), cpositions = c(0.02, 0.15, 0.3))


PLOT.INF.CI = ggforest(coxph(Surv(time_death_fu_truncated, event) ~ Age + Sex + IPI + ATC + `CLL CI score`, data = COX_INF_CI_score %>% 
                               dplyr::rename(ATC = ATC.all.cut,
                                             `CLL CI score` = CLL_CI_score) %>% 
                               as.data.frame()), cpositions = c(0.02, 0.15, 0.3))

ggsave( paste0(getwd(), '/Figure_S_COX_CI.png'),
        ggarrange(PLOT.OS.CI, PLOT.ADT.CI, PLOT.INF.CI, nrow = 1, ncol = 3, labels = LETTERS[1:3]),
        height = 10,
        width = 16,
        dpi = 300)

#### Table 1 ####
ALL.DISEASES_NMI = ALL.DISEASES %>% 
  left_join(NMI_data %>% select(patientid, NMI_score)) %>% 
  left_join(CLL_CI_data %>% select(patientid, CLL_CI_score))

utable(Disease ~ Q(Age) + Sex + IPI + ATC.all.cut +  CCI.f2 + Q(CCI) + Q(NMI_score) + Q(CLL_CI_score), ALL.DISEASES_NMI,
       "IPI"="Disease-specific IPI",
       "ATC.all.cut" = "No. of medications",
       "CCI.f2" = "CCI intervals",
       "CCI"  = "CCI score",
       "NMI_score" = "NMI score",
       "CLL_CI_score"  = "CLL-CI score") %>% write_utable(table_n = 1)

##### Table 2 #####
write_csv2(utable(ATC.all.cut ~ Q(Age) + Sex + IPI + CCI.f2 + Q(CCI) + Q(NMI_score) + Q(CLL_CI_score), ALL.DISEASES_NMI,
                  "IPI"="Disease-specific IPI",
                  "CCI.f2" = "CCI intervals",
                  "CCI"  = "CCI score",
                  "NMI_score" = "NMI score",
                  "CLL_CI_score"  = "CLL-CI score") %>%  publish %>% as.data.frame(),
           paste0(getwd(), '/Table_2-strata_by_medications.csv') 
)


#### Table 3: adjusting for Age + sex ####
###### Table3A: OS ######
ATC.CLASSES_OS = UNI_OS4$ATC
DX.FIRST.ATC3.OS2 = DX.FIRST.ATC3.OS %>% 
  left_join(ALL.DISEASES %>%  select(patientid, Age, sex)) %>% 
  select(patientid, contains(ATC.CLASSES_OS), time_dx_death, status,  Age, sex)
LIST = purrr::map(ATC.CLASSES_OS, ~coxph(as.formula(paste("Surv(time_dx_death, status) ~  Age + sex + ", .x)), DX.FIRST.ATC3.OS2))
UNI_OS_adjusted = data.table()

for(i in 1:length(LIST)){
  pval = as.numeric(coef(summary(LIST[[i]]))[,5][3])
  UNI_OS_adjusted = publish(LIST[[i]])$regressionTable %>% 
    mutate(n = as.numeric(table(DX.FIRST.ATC3.OS2[,i+1])[2]),
           pval = pval) %>% 
    bind_rows(UNI_OS_adjusted)
}

UNI_OS2_adjusted = UNI_OS_adjusted %>% 
  mutate(Variable = ifelse(Variable =='', NA, Variable)) %>% 
  fill(Variable) %>% 
  filter(CI.95 != '') %>% 
  mutate(HR = as.numeric(HazardRatio),
         lower = as.numeric(gsub('\\[', '', str_split_fixed(CI.95, ';', 2)[,1])),
         upper = as.numeric(gsub('\\]', '', str_split_fixed(CI.95, ';', 2)[,2])),
         n = ifelse(n<5, '<5', n))

UNI_OS3_adjusted = UNI_OS2_adjusted %>% 
  arrange(desc(HR)) %>%
  dplyr::rename(ATC= Variable) %>% 
  left_join(Codes_ATC %>% select(ATC = atc, name))  %>% 
  filter(!is.na(name))

UNI_OS3_adjusted$FDR = p.adjust(UNI_OS3_adjusted$pval, 'fdr')
UNI_OS3_adjusted$bonferroni = p.adjust(UNI_OS3_adjusted$pval, 'bonferroni')

UNI_OS4_adjusted = UNI_OS3_adjusted %>% 
  filter( HR != 0.00,
          !is.na(HR),
          bonferroni < 0.05) %>% 
  mutate(name = ifelse(ATC == 'A07', 'antidiarrheals, intestinal antiinflammatory', name),
         name = ifelse(ATC == 'G03', 'sex hormones and genital modulators', name),
         FDR_round = sprintf('%.3f', round(FDR, 3)),
         FDR_round = ifelse(FDR_round =='0.000', '<0.001', paste(FDR_round)),
         bonferroni_round = sprintf('%.3f', round(FDR, 3)),
         bonferroni_round = ifelse(bonferroni_round =='0.000', '<0.001', paste(bonferroni_round)))

palette(BLOOD)
UNI_OS4_adjusted[,c(1,11,6,14)] %>% head2 #atc, name, n, FDR_round
UNI_OS4_adjusted[,c(1,11,6,15)] %>% head2 #atc, name, n, FDR_round

plotConfidence(x=UNI_OS4_adjusted$HR,
               lower = UNI_OS4_adjusted$lower,
               upper = UNI_OS4_adjusted$upper,
               labels = UNI_OS4_adjusted[,c(1,11,6,15)],
               title.labels= c('ATC', 'Description', 'n', 'Bonferroni'), #15
               title.values=expression(bold(HR (CI[95]))),
               cex=1.5,
               col="black",
               xlab.cex=1.3,
               xlim = c(0.5, 40),
               points.col=c(1:9),
               xlab="Hazard ratio for overall survival",
               plot.log="x")

###### Table 3B: ADT ######
#hospitalization
ATC.CLASSES_ADT = UNI_ADT4$ATC
FINAL_COHORT_ADT.1.ADT2 = FINAL_COHORT_ADT.1.ADT %>% 
  left_join(ALL.DISEASES %>%  select(patientid, Age, sex)) %>% 
  select(patientid, contains(ATC.CLASSES_ADT), time_golive_ADT_death, ADT,  Age, sex)
LIST = purrr::map(ATC.CLASSES_ADT, ~coxph(as.formula(paste("Surv(time_golive_ADT_death, ADT) ~  Age + sex + ", .x)), FINAL_COHORT_ADT.1.ADT2)) 
UNI_ADT_adjusted = data.table()
for(i in 1:length(LIST)){
  pval = as.numeric(coef(summary(LIST[[i]]))[,5][3]) 
  UNI_ADT_adjusted = publish(LIST[[i]])$regressionTable %>% 
    mutate(n = as.numeric(table(FINAL_COHORT_ADT.1.ADT2[,i+1])[2]),
           pval = pval) %>% 
    bind_rows(UNI_ADT_adjusted)
}

UNI_ADT_adjusted2 = UNI_ADT_adjusted %>% 
  mutate(Variable = ifelse(Variable =='', NA, Variable)) %>% 
  fill(Variable) %>% 
  filter(CI.95 != '') %>% 
  mutate(HR = as.numeric(HazardRatio),
         lower = as.numeric(gsub('\\[', '', str_split_fixed(CI.95, ';', 2)[,1])),
         upper = as.numeric(gsub('\\]', '', str_split_fixed(CI.95, ';', 2)[,2])),
         n = ifelse(n<5, '<5', n)) 

UNI_ADT_adjusted3 = UNI_ADT_adjusted2 %>% 
  arrange(desc(HR)) %>% 
  dplyr::rename(ATC= Variable) %>% 
  left_join(Codes_ATC %>% select(ATC = atc, name)) %>% 
  filter(!is.na(name))

UNI_ADT_adjusted3$FDR = p.adjust(UNI_ADT_adjusted3$pval, 'fdr')
UNI_ADT_adjusted3$bonferroni = p.adjust(UNI_ADT_adjusted3$pval, 'bonferroni')

UNI_ADT_adjusted4 = UNI_ADT_adjusted3 %>% 
  filter( HR != 0.00,
          !is.na(HR),
          bonferroni < 0.05) %>% 
  mutate(name = ifelse(ATC == 'A07', 'antidiarrheals, intestinal antiinflammatory', name),
         name = ifelse(ATC == 'G03', 'sex hormones and genital modulators', name),
         FDR_round = sprintf('%.3f', round(FDR, 3)),
         FDR_round = ifelse(FDR_round =='0.000', '<0.001', paste(FDR_round)),
         bonferroni_round = sprintf('%.3f', round(FDR, 3)),
         bonferroni_round = ifelse(bonferroni_round =='0.000', '<0.001', paste(bonferroni_round)))

UNI_ADT_adjusted4 %>% head2
palette(BLOOD)

plotConfidence(x=UNI_ADT_adjusted4$HR,
               lower = UNI_ADT_adjusted4$lower,
               upper = UNI_ADT_adjusted4$upper,
               labels = UNI_ADT_adjusted4[,c(1,11,6,15)],
               title.labels= c('ATC', 'Description', 'n', 'Bonferroni'), #15
               title.values=expression(bold(HR (CI[95]))),
               cex=1.5,
               col="black",
               xlab.cex=1.3,
               xlim = c(1, 3),
               points.col=c(1:9),
               xlab="Hazard ratio for hospitalization",
               plot.log="x")

###### Table 3C: Infection ######
# infections
ATC.CLASSES_INF = UNI_INF4$ATC
ALL.DISEASES_INF.1.ATC2 = ALL.DISEASES_INF.1.ATC %>% 
  select(patientid, contains(ATC.CLASSES_INF), time_death_fu_truncated, event,  Age, sex)
LIST = purrr::map(ATC.CLASSES_INF, ~coxph(as.formula(paste("Surv(time_death_fu_truncated, event) ~  Age + sex + ", .x)), ALL.DISEASES_INF.1.ATC2)) 
UNI_INF_adjusted = data.table()

for(i in 1:length(LIST)){
  pval = as.numeric(coef(summary(LIST[[i]]))[,5][3])
  UNI_INF_adjusted = publish(LIST[[i]])$regressionTable %>% 
    mutate(n = as.numeric(table(ALL.DISEASES_INF.1.ATC2[,i+1])[2]),
           pval = pval) %>% 
    bind_rows(UNI_INF_adjusted)
}

UNI_INF_adjusted2 = UNI_INF_adjusted %>% 
  mutate(Variable = ifelse(Variable =='', NA, Variable)) %>% 
  fill(Variable) %>% 
  filter(CI.95 != '') %>% 
  mutate(HR = as.numeric(HazardRatio),
         lower = as.numeric(gsub('\\[', '', str_split_fixed(CI.95, ';', 2)[,1])),
         upper = as.numeric(gsub('\\]', '', str_split_fixed(CI.95, ';', 2)[,2])),
         n = ifelse(n<5, '<5', n)) 

UNI_INF_adjusted3 = UNI_INF_adjusted2 %>% 
  arrange(desc(HR)) %>% 
  dplyr::rename(ATC= Variable) %>% 
  left_join(Codes_ATC %>% select(ATC = atc, name)) %>% 
  filter(!is.na(name))

UNI_INF_adjusted3$FDR = p.adjust(UNI_INF_adjusted3$pval, 'fdr')
UNI_INF_adjusted3$bonferroni = p.adjust(UNI_INF_adjusted3$pval, 'bonferroni')

UNI_INF_adjusted4 = UNI_INF_adjusted3 %>% 
  filter( HR != 0.00,
          !is.na(HR),
          bonferroni < 0.05) %>% 
  mutate(name = ifelse(ATC == 'A07', 'antidiarrheals, intestinal antiinflammatory', name),
         name = ifelse(ATC == 'G03', 'sex hormones and genital modulators', name),
         FDR_round = sprintf('%.3f', round(FDR, 3)),
         FDR_round = ifelse(FDR_round =='0.000', '<0.001', paste(FDR_round)),
         bonferroni_round = sprintf('%.3f', round(FDR, 3)),
         bonferroni_round = ifelse(bonferroni_round =='0.000', '<0.001', paste(bonferroni_round)))

palette(BLOOD)
plotConfidence(x=UNI_INF_adjusted4$HR,
               lower = UNI_INF_adjusted4$lower,
               upper = UNI_INF_adjusted4$upper,
               labels = UNI_INF_adjusted4[,c(1,11,6,15)],
               title.labels= c('ATC', 'Description', 'n', 'Bonferroni'), #15
               title.values=expression(bold(HR (CI[95]))),
               cex=1.5,
               col="black",
               xlab.cex=1.3,
               xlim = c(0.5, 11),
               points.col=c(1:9),
               xlab="Hazard ratio for severe infection",
               plot.log="x")

###### Gather Table 3 ######
UNI_OS4_adjusted %>% nrow 
UNI_ADT_adjusted4 %>% nrow
UNI_INF_adjusted4 %>%  nrow

UNI_OS4_adjusted$p_adjust %>% summary
UNI_OS4_adjusted$bonferroni_round %>% summary

TABLE_3 = UNI_OS4_adjusted %>% 
  transmute(ATC, n, HR_os = paste(HazardRatio, CI.95), p_val_os = bonferroni_round) %>% 
  full_join(UNI_ADT_adjusted4  %>% 
              transmute(ATC, n_ADT = n, HR_adt = paste(HazardRatio, CI.95), p_val_adt = bonferroni_round), 'ATC') %>% 
  full_join(UNI_INF_adjusted4  %>% 
              transmute(ATC, n_INF = n, HR_inf = paste(HazardRatio, CI.95), p_val_inf = bonferroni_round), 'ATC') %>% 
  left_join(Codes_ATC %>% select(ATC = atc, name)) %>% 
  mutate(n_ADT = ifelse(is.na(n_ADT), n_INF, n_ADT)) %>% 
  select(ATC, name, everything(), -n_INF) 

write_csv2(TABLE_3, paste0(getwd(), '/Table_3.csv'))

#### Table S1 ####
DX.FIRST.CCI2 = DX.FIRST.CCI %>% 
  mutate(across(contains('score'), ~ factor(ifelse(.==0, 'No', 'Yes'), levels = c('No', 'Yes')))) %>% 
  left_join(ALL.DISEASES %>% select(patientid, Disease))

DX.FIRST.CCI2 %>% names %>% paste(collapse = ' + ') #inserted these in utable() below

table_S_CCI = utable(Disease~CCI.CHF.score + CCI.Dementia.score + CCI.CPD.score + CCI.Rheumatic.score + CCI.Liver.Mild.score + CCI.Hemiplegia.score + CCI.Renal.score + CCI.DM.w.compl.score + CCI.Cancer.and.Hem.score + CCI.Liver.Severe.score + CCI.Cancer.Metastatic.score + CCI.AIDS.score, 
                     DX.FIRST.CCI2 %>%  mutate(across(contains('CCI'), ~ factor(., c('Yes', 'No'))))) %>% 
  publish %>% as.data.frame() %>% 
  mutate(V12 = ifelse(str_detect(V12, '0', negate = T), NA, V12)) %>% 
  fill(V12, .direction = 'up') %>% 
  mutate(V12 = ifelse(str_detect(V1, 'Hem.score'), NA, V12)) %>% 
  filter_str_detect(V2, 'Yes') %>% 
  select(-V2)
write_csv2(table_S_CCI, paste0(getwd(), '/table_S1-CCI.csv'))

#### Table S2 ####
# Not adjusted for age + sex

###### Table S2A ######

## UNI ATC analyses
ATC.CLASSES = names(DX.FIRST.ATC3.OS)[2:80] 
ATC.CLASSES
LIST = purrr::map(ATC.CLASSES, ~coxph(as.formula(paste("Surv(time_dx_death, status) ~  ", .x)), DX.FIRST.ATC3.OS))
UNI_OS = data.table()
for(i in 1:length(LIST)){
  pval = coef(summary(LIST[[i]]))[,5]
  UNI_OS = publish(LIST[[i]])$regressionTable %>% 
    mutate(n = as.numeric(table(DX.FIRST.ATC3.OS[,i+1])[2]),
           pval = pval) %>% 
    bind_rows(UNI_OS)
}

UNI_OS2 = UNI_OS %>% 
  mutate(Variable = ifelse(Variable =='', NA, Variable)) %>% 
  fill(Variable) %>% 
  filter(CI.95 != '') %>% 
  mutate(HR = as.numeric(HazardRatio),
         lower = as.numeric(gsub('\\[', '', str_split_fixed(CI.95, ';', 2)[,1])),
         upper = as.numeric(gsub('\\]', '', str_split_fixed(CI.95, ';', 2)[,2])),
         n = ifelse(n<5, '<5', n)) 

UNI_OS3 = UNI_OS2 %>% 
  arrange(desc(HR)) %>% 
  dplyr::rename(ATC= Variable) %>% 
  left_join(Codes_ATC %>% select(ATC = atc, name)) 

UNI_OS3$FDR = p.adjust(UNI_OS3$pval, 'fdr')
UNI_OS3$bonferroni = p.adjust(UNI_OS3$pval, 'bonferroni')

UNI_OS4 = UNI_OS3 %>% 
  filter(bonferroni <0.05) %>% 
  mutate(name = ifelse(ATC == 'A07', 'antidiarrheals, intestinal antiinflammatory', name),
         name = ifelse(ATC == 'G03', 'sex hormones and genital modulators', name),
         FDR_round = sprintf('%.3f', round(FDR, 3)),
         FDR_round = ifelse(FDR_round =='0.000', '<0.001', paste(FDR_round)),
         bonferroni_round = sprintf('%.3f', round(FDR, 3)),
         bonferroni_round = ifelse(bonferroni_round =='0.000', '<0.001', paste(bonferroni_round)))

plotConfidence(x=UNI_OS4$HR,
               lower = UNI_OS4$lower,
               upper = UNI_OS4$upper,
               labels = UNI_OS4[,c(1,13,6,12)],
               title.labels= c('ATC', 'Description', 'n', 'fdr'),
               title.values=expression(bold(HR (CI[95]))),
               cex=1.5,
               col="black",
               xlab.cex=1.3,
               xlim = c(0.1, 20),
               points.col=c(1:9),
               xlab="Hazard ratio for overall survival",
               plot.log="x")

###### Table S2B ######
#hospitalization
FINAL_COHORT_ADT.1.ADT = FINAL_COHORT_ADT.1.DURATION %>% 
  select(patientid, Disease, ADT, time_golive_ADT_death) %>% 
  left_join(DX.FIRST.ATC3.OS) %>% 
  select(patientid, contains(ATC.CLASSES), ADT, time_golive_ADT_death, Disease)
ATC.CLASSES %>% n_distinct()
LIST = purrr::map(ATC.CLASSES, ~coxph(as.formula(paste("Surv(time_golive_ADT_death, ADT) ~  ", .x)), FINAL_COHORT_ADT.1.ADT)) 
UNI_ADT = data.table()
for(i in 1:length(LIST)){
  pval = coef(summary(LIST[[i]]))[,5]
  UNI_ADT = publish(LIST[[i]])$regressionTable %>% 
    mutate(n = as.numeric(table(FINAL_COHORT_ADT.1.ADT[,i+1])[2]),
           pval = pval) %>% 
    bind_rows(UNI_ADT)
}

UNI_ADT2 = UNI_ADT %>% 
  mutate(Variable = ifelse(Variable =='', NA, Variable)) %>% 
  fill(Variable) %>% 
  filter(CI.95 != '') %>% 
  mutate(HR = as.numeric(HazardRatio),
         lower = as.numeric(gsub('\\[', '', str_split_fixed(CI.95, ';', 2)[,1])),
         upper = as.numeric(gsub('\\]', '', str_split_fixed(CI.95, ';', 2)[,2])),
         n = ifelse(n<5, '<5', n)) 

UNI_ADT3 = UNI_ADT2 %>% 
  arrange(desc(HR)) %>% 
  dplyr::rename(ATC= Variable) %>% 
  left_join(Codes_ATC %>% select(ATC = atc, name)) 


UNI_ADT3$FDR = p.adjust(UNI_ADT3$pval, 'fdr')
UNI_ADT3$bonferroni = p.adjust(UNI_ADT3$pval, 'bonferroni')

UNI_ADT4 = UNI_ADT3 %>% 
  filter(bonferroni <0.05) %>% 
  mutate(name = ifelse(ATC == 'A07', 'antidiarrheals, intestinal antiinflammatory', name),
         name = ifelse(ATC == 'G03', 'sex hormones and genital modulators', name),
         FDR_round = sprintf('%.3f', round(FDR, 3)),
         FDR_round = ifelse(FDR_round =='0.000', '<0.001', paste(FDR_round)),
         bonferroni_round = sprintf('%.3f', round(FDR, 3)),
         bonferroni_round = ifelse(bonferroni_round =='0.000', '<0.001', paste(bonferroni_round)))

plotConfidence(x=UNI_ADT4$HR,
               lower = UNI_ADT4$lower,
               upper = UNI_ADT4$upper,
               labels = UNI_ADT4[,c(1,13,6,12)],
               title.labels= c('ATC', 'Description', 'n', 'fdr'),
               title.values=expression(bold(HR (CI[95]))),
               cex=1.5,
               col="black",
               xlab.cex=1.3,
               xlim = c(1, 3),
               points.col=c(1:9),
               xlab="Hazard ratio for hospitalization",
               plot.log="x")

###### Table S2C ######
# infections
ALL.DISEASES_INF.1.ATC = ALL.DISEASES_INF.1 %>% 
  left_join(DX.FIRST.ATC3.OS %>% select(-time_dx_death, -status)) %>% 
  select(patientid, contains(ATC.CLASSES), time_death_fu_truncated, event, Age, sex = sex.x)
LIST = purrr::map(ATC.CLASSES, ~coxph(as.formula(paste("Surv(time_death_fu_truncated, event) ~  ", .x)), ALL.DISEASES_INF.1.ATC)) 
UNI_INF = data.table()
for(i in 1:length(LIST)){
  pval = coef(summary(LIST[[i]]))[,5]
  UNI_INF = publish(LIST[[i]])$regressionTable %>% 
    mutate(n = as.numeric(table(ALL.DISEASES_INF.1.ATC[,i+1])[2]),
           pval = pval) %>% 
    bind_rows(UNI_INF)
}

UNI_INF2 = UNI_INF %>% 
  mutate(Variable = ifelse(Variable =='', NA, Variable)) %>% 
  fill(Variable) %>% 
  filter(CI.95 != '') %>% 
  mutate(HR = as.numeric(HazardRatio),
         lower = as.numeric(gsub('\\[', '', str_split_fixed(CI.95, ';', 2)[,1])),
         upper = as.numeric(gsub('\\]', '', str_split_fixed(CI.95, ';', 2)[,2])),
         n = ifelse(n<5, '<5', n)) 

UNI_INF3 = UNI_INF2 %>% 
  arrange(desc(HR)) %>% 
  dplyr::rename(ATC= Variable) %>% 
  left_join(Codes_ATC %>% select(ATC = atc, name)) 

UNI_INF3$FDR = p.adjust(UNI_INF3$pval, 'fdr')
UNI_INF3$bonferroni = p.adjust(UNI_INF3$pval, 'bonferroni')

UNI_INF4 = UNI_INF3 %>% 
  filter(bonferroni <0.05) %>% 
  mutate(name = ifelse(ATC == 'A07', 'antidiarrheals, intestinal antiinflammatory', name),
         name = ifelse(ATC == 'G03', 'sex hormones and genital modulators', name),
         FDR_round = sprintf('%.3f', round(FDR, 3)),
         FDR_round = ifelse(FDR_round =='0.000', '<0.001', paste(FDR_round)),
         bonferroni_round = sprintf('%.3f', round(FDR, 3)),
         bonferroni_round = ifelse(bonferroni_round =='0.000', '<0.001', paste(bonferroni_round)))

plotConfidence(x=UNI_INF4$HR,
               lower = UNI_INF4$lower,
               upper = UNI_INF4$upper,
               labels = UNI_INF4[,c(1,13,6,12)],
               title.labels= c('ATC', 'Description', 'n', 'fdr'),
               title.values=expression(bold(HR (CI[95]))),
               cex=1.5,
               col="black",
               xlab.cex=1.3,
               xlim = c(0.5, 11),
               points.col=c(1:9),
               xlab="Hazard ratio for severe infection",
               plot.log="x")


###### Gather Table S2 ######
UNI_OS4 %>% nrow()
UNI_ADT4 %>% nrow()
UNI_INF4 %>% nrow()

TABLE_S2 = UNI_OS4 %>% 
  transmute(ATC, n, HR_os = paste(HazardRatio, CI.95), p_val_os = bonferroni_round) %>% 
  full_join(UNI_ADT4  %>% 
              transmute(ATC, n_ADT = n, HR_adt = paste(HazardRatio, CI.95), p_val_adt = bonferroni_round), 'ATC') %>% 
  full_join(UNI_INF4  %>% 
              transmute(ATC, n_INF = n, HR_inf = paste(HazardRatio, CI.95), p_val_inf = bonferroni_round), 'ATC') %>% 
  left_join(Codes_ATC %>% select(ATC = atc, name)) %>% 
  mutate(n_ADT = ifelse(is.na(n_ADT), n_INF, n_ADT)) %>% 
  select(ATC, name, everything(), -n_INF) 

write_csv2(TABLE_S2, paste0(getwd(), '/Table_extra_S2.csv'))

if(SAVE){
save.image(paste0(getwd(), '/After_figures.RData'))
}else{
  cat('\nRData file not saved.')
}
# load('/ngc/projects2/dalyca_r/chribr_r/PolyRx/data/After_figures.RData')

cat(paste('\nAll figures saved to:\n', getwd()))
cat('\nEnd of script.')

# End of script