##################################################
## Project:4th Cycle
## Script purpose: Hazard Ratios
## Datee: 06/28/2020
## Author: Anh-Minh Nguyen
##################################################
rm(list = ls()) #Maybe use with Kaplan Meier 1
library(tidyverse)
library(finalfit)
library(survival)
library(survminer)

medops <- read.csv('MedicationsOperations.csv', stringsAsFactors = F)
blf <- read.csv('BonyLesionandFractureDexamethasoneOperations.csv', stringsAsFactors = F) %>% 
  dplyr::select(-X)
fs <- read.csv('final_stage.csv', stringsAsFactors = F)
labs <- read.csv('Data/Labs2.csv', stringsAsFactors = F)
mye <- read.csv('Data/MyelomaTherapy.csv', stringsAsFactors = F)
pfs <- mye %>% 
  group_by(ID) %>% 
  dplyr::filter(TreatmentPhase == 'Relapse') %>% 
  dplyr::filter(Line == min(Line)) %>% 
  distinct(ID, Line, .keep_all = T) %>% 
  dplyr::select(ID, DaysFromDxStart) %>% 
  ungroup()
colnames(pfs) <- c('ID', 'PFS')
superlist <- unique(c(pfs$ID, labs$ID, medops$ID, blf$ID, fs$ID))
# L Abnormal 'VITD3','EDOL','TEST'
# H Abnormal 'CRE','HSCRP','ESR'
# Both Abnormal 'TSH'
labsmax <- labs %>% 
  dplyr::filter(ObservationId %in% c('CRE','HSCRP','ESR')) %>% 
  dplyr::filter(DaysFromDx >= -90 & DaysFromDx <= 365) %>% 
  group_by(ID, ObservationId) %>% 
  summarise(Count = sum(AbnormalFlags == 'H')) %>% 
  ungroup()
labsmax[labsmax$Count > 1, 'Count'] <- 1
labsmax <- pivot_wider(labsmax, names_from=ObservationId, values_from=Count)

labsmin <- labs %>% 
  dplyr::filter(ObservationId %in% c('VITD3','EDOL','TEST')) %>% 
  dplyr::filter(DaysFromDx >= -90 & DaysFromDx <= 365) %>% 
  group_by(ID, ObservationId) %>% 
  summarise(Count = sum(AbnormalFlags == 'L')) %>% 
  ungroup()
labsmin[labsmin$Count > 1, 'Count'] <- 1
labsmin <- pivot_wider(labsmin, names_from=ObservationId, values_from=Count)

labsany <- labs %>% 
  dplyr::filter(ObservationId %in% c('TSH')) %>% 
  dplyr::filter(DaysFromDx >= -90 & DaysFromDx <= 365) %>% 
  group_by(ID, ObservationId) %>% 
  summarise(Count = sum(AbnormalFlags == 'L' | AbnormalFlags == 'H')) %>% 
  ungroup()
labsany[labsany$Count > 1, 'Count'] <- 1
labsany <- pivot_wider(labsany, names_from=ObservationId, values_from=Count)

### Try to put together the hazard ratio table
letsgo <- merge(blf, pfs, by = 'ID')
letsgo <- left_join(letsgo, fs, by = 'ID') %>% 
  dplyr::select(-StagingSystem)
letsgo$AgeAtDx <- cut(letsgo$AgeAtDx, breaks=c(0,30,70,10000),
                labels=c("<30","30-70",">70"))
grp <- c('MM383', 'MM583', 'MM711')
sum(grp %in% labs$ID)
letsgo <- letsgo[!is.na(letsgo$Stage),] # for now

letsgo$LessPFS <- as.numeric(letsgo$PFS < mean(letsgo$PFS))
l1 <- letsgo[letsgo$BonyLesions == 1 & letsgo$Denosumab == 1,]
l2 <- letsgo[letsgo$BonyLesions == 1 & letsgo$Dexamethasone == 1,]
l21 <- letsgo[letsgo$BonyLesions == 1 & letsgo$Dexamethasone == 0,]
l22 <- letsgo[letsgo$BonyLesions == 0 & letsgo$Dexamethasone == 1,]



# from finalfit
explanatory = c("AgeAtDx", "PatientSex", "Stage", "RacialGroup")
dependent = "Surv(PFS, LessPFS)"
l2 %>%
  hr_plot(dependent, explanatory, dependent_label = "Survival  [Bony Lesions = 1, Dexamethasone = 1]",
          table_text_size=4, title_text_size=14,
          plot_opts=list(xlab("HR, 95% CI"), theme(axis.title = element_text(size=12))))
l22 %>%
  hr_plot(dependent, explanatory, dependent_label = "Survival [Bony Lesions = 0, Dexamethasone = 1]",
          table_text_size=4, title_text_size=14,
          plot_opts=list(xlab("HR, 95% CI"), theme(axis.title = element_text(size=12)))) +
  theme_minimal()
# from survival
bigmodel <- coxph(Surv(PFS, LessPFS) ~ PatientSex + AgeAtDx + Stage + RacialGroup,
        data = letsgo)
ggforest(bigmodel)
