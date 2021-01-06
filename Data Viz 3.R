##################################################
## Project:4th Cycle
## Script purpose: Stem 1st to 1st Relapse Data Visualizations
## Date: 05/25/2020
## Author: Anh-Minh Nguyen
##################################################
rm(list = ls())
library(tidyverse)
library(magrittr)
library(qwraps2)

options(qwraps2_markup = "markdown")
rel1 <- read.csv('PFS Merge.csv', stringsAsFactors = F)
nosct <- read.csv('PFS NO SCT Merge.csv', stringsAsFactors = F)
rel1dem <- read.csv('314_Patients_DemographicInfo.csv', stringsAsFactors = F) %>% 
  dplyr::select(ID, AgeAtDx, PatientSex)
yessct <- rel1[!(rel1$ID %in% nosct$ID),]

yessct <- merge(rel1dem, yessct, by = 'ID')
nosct <- merge(rel1dem, nosct, by = 'ID')

yessct <- yessct[!is.na(yessct$VITD3),]
yessct[yessct$VITD3 == 1, 'VITD3'] <- 'Deficiency'
yessct[yessct$VITD3 == 0, 'VITD3'] <- 'Sufficiency'
nosct <- nosct[!is.na(nosct$VITD3),]
nosct[nosct$VITD3 == 1, 'VITD3'] <- 'Deficiency'
nosct[nosct$VITD3 == 0, 'VITD3'] <- 'Sufficiency'

ggplot(data = yessct, aes(x=VITD3, fill = PatientSex)) +
  geom_bar(position = 'fill') +
  stat_count(geom = "text", 
             aes(label = paste(round((..count..)/sum(..count..)*100), "%")),
             position=position_fill(vjust=0.5), colour="white") +
  ggtitle('Patients with Relapse phase and SCT History') +
  xlab('Vitamin D3 Levels') +
  ylab('Percentages')

ggplot(data = nosct, aes(x=VITD3, fill = PatientSex)) +
  geom_bar(position = 'fill') +
  stat_count(geom = "text", 
             aes(label = paste(round((..count..)/sum(..count..)*100), "%")),
             position=position_fill(vjust=0.5), colour="white") +
  ggtitle('Patients with Relapse phase and NO SCT History') +
  xlab('Vitamin D3 Levels') +
  ylab('Percentages')

yessct2 <- yessct
yessct2$EDOL <- NULL
yessct2$PatientSex <- as.factor(yessct2$PatientSex)
yessct2$VITD3 <- as.factor(yessct2$VITD3)
our_summary1 <- list('Average PFS',
                     list('Mean' = ~ mean(.data$PFS),
                          'Min' = ~ min(.data$PFS))
)


whole <- summary_table(yessct2, our_summary1)
whole

