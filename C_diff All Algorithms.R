##################################################
## Project: C.Diff Cycle
## Script purpose: Infected Analysis Aggregation
## Date: 2020/11/14
## Author: Anh-Minh Nguyen
##################################################
rm(list = ls()) 

############### INSTALL NECESSARY PACKAGES #############
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('caret')) install.packages('caret'); library('caret')
if (!require('glmnet')) install.packages('glmnet'); library('glmnet')
if (!require('naivebayes')) install.packages('naivebayes'); library('naivebayes')
if (!require('e1071')) install.packages('e1071'); library('e1071')
if (!require('formattable')) install.packages('formattable'); library('formattable')

# Read in the data frame
dfc <- read.csv('full_dataset-second.csv') %>% 
  dplyr::select(-X)
y <- read.csv('infections.csv') %>% 
  dplyr::select(-X)


## Aesthetic variable naming, not necessary
names(dfc) <- gsub('..', '_',
                      colnames(dfc), fixed = T)
names(dfc) <- gsub('.', '_',
                      colnames(dfc), fixed = T)


# Turn the ID variable into the row names
# Merge the infected file with the feature dataset
dfc2 <- merge(dfc, y, by = 'ID')
rownames(dfc2) <- dfc2$ID
dfc2$ID <- NULL

# Combine Myeloma and Multiple Myeloma into one feature
dfc2$Multiple_myeloma <- dfc$Multiple_myeloma + dfc$Myeloma
dfc2$Multiple_myeloma[dfc$Multiple_myeloma > 1] <- 1
dfc2$Myeloma <- NULL

# Separate into outcome and features
x <- dfc2[,names(dfc2) != 'Infected']

# for algorithms that require the outcome to be binary
fy <- as.factor(ifelse(dfc2$Infected == 1, 'Yes', 'No'))

# Set cross-validation folds
set.seed(144)
train_control <- trainControl(method="cv", number=5, classProbs = TRUE)
glm1 <- train(x = x, y = fy, method = 'glm', trControl = train_control,
             family = 'binomial')
glm.pred <- predict(glm1, newdata = x)

# Ridge, Lasso
parameters <- c(seq(0.1, 2, by =0.1) ,  seq(2, 5, 0.5) , 
                seq(5, 25, 1))
ridge1 <- train(x = x, y = fy, method = 'glmnet', 
                tuneGrid = expand.grid(alpha = 0, lambda = parameters),
                trControl = train_control,
                family = 'binomial')
ridge.pred <- predict(ridge1, newdata = x)

lasso1 <- train(x = x, y = fy, method = 'glmnet', 
                tuneGrid = expand.grid(alpha = 1, lambda = parameters),
                trControl = train_control,
                family = 'binomial')
lasso.pred <- predict(lasso1, newdata = x)


# Naive Bayes
nb1 <- train(x = x, y = fy, method = 'naive_bayes', 
                trControl = train_control)
nb.pred <- predict(nb1, newdata = x)

# Support Vector Machines with Class Weights
svm1 <- train(x, y = fy, method = 'svmLinearWeights',
              trControl = train_control)
svm.pred <- predict(svm1, newdata = x)
confusionMatrix(svm.pred, fy, positive = 'Yes')

###########################################
# Model Measurements
# Create Table to 
glm.cm <- confusionMatrix(glm.pred, fy)
lasso.cm <- confusionMatrix(lasso.pred, fy)
ridge.cm <- confusionMatrix(ridge.pred, fy)
nb.cm <- confusionMatrix(nb.pred, fy)
svm.cm <- confusionMatrix(svm.pred, fy)

results.df <-
  as.data.frame(round(rbind(c(glm.cm$overall[1], glm.cm$byClass),
  c(lasso.cm$overall[1], lasso.cm$byClass),
  c(ridge.cm$overall[1], ridge.cm$byClass),
  c(nb.cm$overall[1], nb.cm$byClass),
  c(svm.cm$overall[1], svm.cm$byClass)), 4))

rownames(results.df) <- c('GLM', 'Lasso', 'Ridge', 'Naive Bayes',
                          'SVM_Linear')
rf <- as.data.frame(t(results.df))
formattable(rf)
