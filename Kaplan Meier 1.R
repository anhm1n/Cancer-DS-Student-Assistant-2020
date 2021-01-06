##################################################
## Project:4th Cycle
## Script purpose: Kaplan Meier
## Date: 06/24/2020
## Author: Anh-Minh Nguyen
##################################################
rm(list = ls())
library(tidyverse)

# Load in datasets/create subsets
blt <- read.csv('Bony & Pain.csv', stringsAsFactors = F)
mye <- read.csv('Data/MyelomaTherapy.csv', stringsAsFactors = F)
pfs <- mye %>% 
  group_by(ID) %>% 
  dplyr::filter(TreatmentPhase == 'Relapse') %>% 
  dplyr::filter(Line == min(Line)) %>% 
  distinct(ID, Line, .keep_all = T) %>% 
  dplyr::select(ID, DaysFromDxStart) %>% 
  ungroup()

# Create new function for filling in missing survivals
fillinSurvival <- function(id, dataset) {
  xx <- mye %>% 
    dplyr::filter(ID == id) %>% 
    dplyr::filter(DaysFromDxStart == max(DaysFromDxStart)) %>% 
    distinct(ID, DaysFromDxStart, .keep_all = T)
  if (is.na(xx[xx$ID == id, 'DaysFromDxStop'])) {
    yy = xx[xx$ID == id, 'DaysFromDxStart']
  } else {
    yy = xx[xx$ID == id, 'DaysFromDxStop']
  }
  dataset[dataset$ID == id, 'OS'] <- yy
  return(dataset)
}

# Replace missing survival
survivaldays <- read.csv('Data/SurvivalDays.csv', stringsAsFactors = F)
survivaldays$OS <- survivaldays$SurvivalDays
survivaldays <- survivaldays[, c('ID', 'OS')]
survivaldays <- fillinSurvival('MM96', survivaldays)
survivaldays <- fillinSurvival('MM501', survivaldays)
survivaldays <- fillinSurvival('MM361', survivaldays)
colnames(pfs) <- c('ID', 'PFS')
blt2 <- merge(blt, pfs, by = 'ID')
blt2 <- merge(blt2, survivaldays, by = 'ID')
# update MM361, MM501, MM96 OS (GO BACK AND RUN)

blt3 <- blt2[blt2$Pains == 1 & blt2$BonyLesions == 1,]

# Add demographics
demographics <- read.csv('Data/Demographics.csv', stringsAsFactors = F)
demographics$AgeAtDx <- cut(demographics$AgeAtDx, breaks=c(0,30,70,10000),
                            labels=c("<30","30-70",">70"))
stages <- read.csv('Data/Stage.csv', stringsAsFactors = F) %>% 
  dplyr::select(-StagingSystem)
stages$Stage <- as.factor(stages$Stage)
blt3 <- left_join(blt3, demographics, by = 'ID')
blt3 <- left_join(blt3, stages, by = 'ID')

blt4 <- pivot_longer(blt3, cols = c(PFS, OS), names_to = 'Measure', values_to = 'DaysFromDx')
ggplot(blt4) +
  geom_histogram(aes(x = DaysFromDx, group = Measure, fill = Measure), 
                 alpha = 0.5, position = 'identity', binwidth = 100) +
  geom_vline(xintercept = c(mean(blt2$PFS), mean(blt2$OS)), lty = c(2,2),
             color = c('red', 'blue')) +
  scale_x_continuous(name ="DaysFromDx") +
  ggtitle('PFS & OS of Patients with Bony Lesions and Reported Pain')

# Add the Labs Tests
labs <- read.csv('Data/Labs2.csv', stringsAsFactors = F)
labsmax <- labs %>% 
  dplyr::filter(ObservationId %in% c('CRE','HSCRP','ESR', 'TSH')) %>% 
  dplyr::filter(DaysFromDx >= -90 & DaysFromDx <= 365) %>% 
  group_by(ID, ObservationId) %>% 
  summarise(Count = sum(AbnormalFlags == 'H')) %>% 
  ungroup()
labsmax[labsmax$Count > 1, 'Count'] <- 1
labsmax <- pivot_wider(labsmax, names_from=ObservationId, values_from=Count)
bonemax <- labsmax[labsmax$ID %in% blt3$ID,]

labsmin <- labs %>% 
  dplyr::filter(ObservationId %in% c('VITD3','EDOL','TEST')) %>% 
  dplyr::filter(DaysFromDx >= -90 & DaysFromDx <= 365) %>% 
  group_by(ID, ObservationId) %>% 
  summarise(Count = sum(AbnormalFlags == 'L')) %>% 
  ungroup()
labsmin[labsmin$Count > 1, 'Count'] <- 1
labsmin <- pivot_wider(labsmin, names_from=ObservationId, values_from=Count)
bonemin <- labsmin[labsmin$ID %in% blt3$ID,]

labsany <- labs %>% 
  dplyr::filter(ObservationId %in% c('TSH')) %>% 
  dplyr::filter(DaysFromDx >= -90 & DaysFromDx <= 365) %>% 
  group_by(ID, ObservationId) %>% 
  summarise(Count = sum(AbnormalFlags == 'L' | AbnormalFlags == 'H')) %>% 
  ungroup()
labsany[labsany$Count > 1, 'Count'] <- 1
labsany <- pivot_wider(labsany, names_from=ObservationId, values_from=Count)
boneany <- labsany[labsany$ID %in% blt3$ID,]


############ Kaplan Meier ##################################################
library(survminer)
library(survival)
blt5 <- blt3
blt5 <- blt5[blt5$OS < 4000,]
blt5 <- blt5[blt5$PFS < 3000,]
blt5$LessPFS <- as.numeric(blt5$PFS < mean(blt5$PFS))
blt5$LessOS <- as.numeric(blt5$OS < mean(blt5$OS))
colnames(blt5)[colnames(blt5) == 'RacialGroup'] <- 'Race'
colnames(blt5)[colnames(blt5) == 'PatientSex'] <- 'Sex'
blt5[blt5$Sex == 'Male', 'Sex'] <- 'M'
blt5[blt5$Sex == 'Female', 'Sex'] <- 'F'
blt6 <- left_join(left_join(blt5,bonemax,by='ID'),bonemin,by='ID')

# Convert to character for quality of life
#blt6[,13:ncol(blt6)][is.na(blt6[,13:ncol(blt6)]) | blt6[,13:ncol(blt6)] == 0] <- '0/NA'
# Created Kaplan Meier Curves
fit1 <- survfit(Surv(PFS, LessPFS) ~ Sex,
                data = blt5)
ggsurvplot(fit1, data = blt5, title = "Less than Average PFS by Sex",
           legend.labs = c(paste('Sex=F:', nrow(blt5[blt5$Sex=='F',])), 
                           paste('Sex=M:', nrow(blt5[blt5$Sex=='M',]))))

fit2 <- survfit(Surv(PFS, LessPFS) ~ AgeAtDx,
                data = blt5)
ggsurvplot(fit2, data = blt5, title = "Less than Average PFS by Age Group",
           legend.labs = c(paste('Age=<30:', nrow(blt5[blt5$AgeAtDx=='<30',])),
                           paste('Age=30-70:', nrow(blt5[blt5$AgeAtDx=='30-70',])),
                           paste('Age=>70:', nrow(blt5[blt5$AgeAtDx=='>70',]))))

fit3 <- survfit(Surv(PFS, LessPFS) ~ Sex + AgeAtDx,
                data = blt5)
ggsurvplot(fit3, data = blt5, title = "Less than Average PFS by Sex & Age",
           legend.labs = c(paste('Sex=F, Age=<30:', nrow(blt5[blt5$AgeAtDx=='<30'&blt5$Sex=='F',])),
                           paste('Sex=F, Age=30-70:', nrow(blt5[blt5$AgeAtDx=='30-70'&blt5$Sex=='F',])),
                           paste('Sex=F, Age=>70:', nrow(blt5[blt5$AgeAtDx=='>70'&blt5$Sex=='F',])),
                           paste('Sex=M, Age=<30:', nrow(blt5[blt5$AgeAtDx=='<30'&blt5$Sex=='M',])),
                           paste('Sex=M, Age=30-70:', nrow(blt5[blt5$AgeAtDx=='30-70'&blt5$Sex=='M',])),
                           paste('Sex=M, Age=>70:', nrow(blt5[blt5$AgeAtDx=='>70'&blt5$Sex=='M',]))))

fit4 <- survfit(Surv(PFS, LessPFS) ~ Race,
                data = blt5)
ggsurvplot(fit4, data = blt5, title = "Less than Average PFS by Race",
           legend.labs = c(paste('Race=Asian:', nrow(blt5[blt5$Race=='Asian',])),
                           paste('Race=Black:', nrow(blt5[blt5$Race=='Black',])),
                           paste('Race=Not reported:', nrow(blt5[blt5$Race=='Not reported',])),
                           paste('Race=Other:', nrow(blt5[blt5$Race=='Other',])),
                           paste('Race=White:', nrow(blt5[blt5$Race=='White',]))))

fit5 <- survfit(Surv(PFS, LessPFS) ~ Race + Sex,
                data = blt5)
ggsurvplot(fit5, data = blt5, title = "Less than Average PFS by Race & Sex",
           legend.labs = c(paste('Race=Asian,Sex=F:', nrow(blt5[blt5$Race=='Asian'&blt5$Sex=='F',])),
                           paste('Race=Asian,Sex=M:', nrow(blt5[blt5$Race=='Asian'&blt5$Sex=='M',])),
                           paste('Race=Black,Sex=F:', nrow(blt5[blt5$Race=='Black'&blt5$Sex=='F',])),
                           paste('Race=Black,Sex=M:', nrow(blt5[blt5$Race=='Black'&blt5$Sex=='M',])),
                           paste('Race=Not reported,Sex=F:', nrow(blt5[blt5$Race=='Not reported'&blt5$Sex=='F',])),
                           paste('Race=Not reported,Sex=M:', nrow(blt5[blt5$Race=='Not reported'&blt5$Sex=='M',])),
                           paste('Race=Other,Sex=F:', nrow(blt5[blt5$Race=='Other'&blt5$Sex=='F',])),
                           paste('Race=Other,Sex=M:', nrow(blt5[blt5$Race=='Other'&blt5$Sex=='M',])),
                           paste('Race=White,Sex=F:', nrow(blt5[blt5$Race=='White'&blt5$Sex=='F',])),
                           paste('Race=White,Sex=M:', nrow(blt5[blt5$Race=='White'&blt5$Sex=='M',]))))

blt7 <- blt6[!(is.na(blt6$VITD3)) & !(is.na(blt6$CRE)),]
fit6 <- survfit(Surv(PFS, LessPFS) ~ VITD3 + CRE,
                data = blt7)
ggsurvplot(fit6, data = blt7, title = "Less than Average PFS by VITD3 & CRE",
          legend.lab = c(paste('VITD3=0,CRE=0:', nrow(blt7[blt7$VITD3==0&blt7$CRE==0,])),
          paste('VITD3=0,CRE=1:', nrow(blt7[blt7$VITD3==0&blt7$CRE==1,])),
          paste('VITD3=1,CRE=0:', nrow(blt7[blt7$VITD3==1&blt7$CRE==0,])),
          paste('VITD3=1,CRE=1:', nrow(blt7[blt7$VITD3==1&blt7$CRE==1,]))))

blt8 <- blt6[!(is.na(blt6$HSCRP)),]
fit7 <- survfit(Surv(PFS, LessPFS) ~ HSCRP,
                data = blt8)
ggsurvplot(fit7, data = blt8, title = "Less than Average PFS by HSCRP",
           legend.labs = c(paste('HSCRP=0:', nrow(blt8[blt8$HSCRP==0,])),
                           paste('HSCRP=1:', nrow(blt8[blt8$HSCRP==1,]))))
