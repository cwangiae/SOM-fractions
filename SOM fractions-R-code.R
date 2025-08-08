
####----set directory----#####
setwd("~/Documents/SOM fractions")
dir()
getwd()

####----load R package----####
#Data processing
library(tidyverse)
library(dplyr)
library(tidyr)
library(readxl)

#Plotting
library(ggplot2)
library(ggsci)
library(ggpmisc)

#Heatmap
library(ggcorrplot)

#Random forest model
library(randomForest)
library(rfUtilities)
library(e1071)

#Multiple regression models
library(MuMIn)
library(performance)
library(see)

#Structural equation modeling
library(lavaan)

#mean
fun.mean<-function(x){c(m=mean(x),n=length(x),s=sd(x)/sqrt(length(x)))} #mean,se
library(doBy)

####----load data----####

som <- read_excel(path = "POM-MAOM.xlsx", sheet = "data", col_names = T) 
dim(som)

colnames(som)[colnames(som) == "logclay+silt"] <- "logclay.silt"
som$logclay.silt


####----Supplementary Fig. 5--Heatmap----####
som1 <- som %>%dplyr::select(ID, Crop.types,
                             logMGT, logMGP, logNDVI,
                             logpH, logSWC, logNP, logCP,logclay.silt,
                             logMBC,logMBN,logBNC,logFNC, logBNN, logFNN,
                             logBS,
                             logFS,
                             Network_PC1,
                             logPOC,logPON,logMAOC,logMAON) %>%
  rename_with(~ gsub("^log", "", .x), .cols = everything()) %>%
  rename("Soil pH" = pH, 
         "Soil N/P" = NP, "Soil C/P" = CP, Net = Network_PC1,
         "Clay + silt" = clay.silt) 


dim(som1)
str(som1)

datacor<-som1[,c(-1,-2)] 
View(datacor)

which(is.na(datacor))


cormat <- round(cor(datacor),2)
p.mat <- cor_pmat(datacor)

col2 = colorRampPalette(c('#053061', '#2166AC', '#4393C3', '#92C5DE',
                          '#D1E5F0', '#FFFFFF', '#FDDBC7', '#F4A582',
                          '#D6604D', '#B2182B', '#67001F'))

M <- cor(datacor, method = "spearman")


pdf(file = "results/CN-heatmap.pdf", width = 6, height = 6) 
par(mar = c(0, 0, 0, 0)) 

corheatmap <- corrplot(M, col = col2(100),type = "lower", 
                      diag = FALSE,  
                      method = 'square',p.mat = p.mat, insig = "label_sig",
                      sig.level = c(.001, .01, .05),pch.cex = 0.9, tl.cex = 0.8,cl.cex = 0.8, 
                      pch.col = "black",tl.col = "black")
dev.off()

####----Fig.3 Random forest model--POC----####
str(som)

data1 <- som %>%
  dplyr::select(Crop.types,
                logMGT, logMGP, logNDVI,
                logpH, logSWC, logNP, logclay.silt,
                logMBC, 
                logBS,logFS,
                Network_PC1,
                logPOC) 

which(is.na(data1))
str(data1)
data1$Crop.types<-as.factor(data1$Crop.types)


#Parameter estimation
library(e1071)
tune_res <- tune.randomForest(x = data1[,-13], y = data1$logPOC,
                              mtry = c(1:12),
                              ntree = seq(100, 1000, by = 100))
tune_res
summary(tune_res) 

#Re-modeling
set.seed(123)
RF_model_poc <- randomForest(logPOC~logMGT+ logMGP+logNDVI+Crop.types+
                               logpH+ logSWC+ logNP +logclay.silt+
                               logMBC+ 
                               logBS+
                               logFS+
                               Network_PC1,
                             data = data1,
                             ntree=500,
                             mtry =5,
                             importance=TRUE,
                             proximity=TRUE)
summary(RF_model_poc)
varImpPlot(RF_model_poc)

#Model significance testing
set.seed(123)
poc.perm <- rf.significance(RF_model_poc, na.omit(data1), nperm=99, ntree=500) #ntree is consistent with the model
poc.perm

#Importance_scores
importance_scores <- importance(RF_model_poc)
importance_scores
write.csv(importance_scores,"results/RF_POC.csv")


####----Fig.3 Random forest model--MAOC----####

data1 <- som %>%
  dplyr::select(Crop.types,
                logMGT, logMGP, logNDVI,
                logpH, logSWC, logNP, logclay.silt,
                logMBC,logBNC,logFNC, 
                logBS,logFS,
                Network_PC1,
                logMAOC) 

which(is.na(data1))
str(data1)
data1$Crop.types<-as.factor(data1$Crop.types)

#Parameter estimation
library(e1071)
tune_res <- tune.randomForest(x = data1[,-15], y = data1$logMAOC,
                              mtry = c(1:14),
                              ntree = seq(100, 1000, by = 100)
)
tune_res
summary(tune_res)

#Re-modeling
set.seed(123)
RF_model_maoc <- randomForest(logMAOC~logMGT+ logMGP+logNDVI+Crop.types+
                                logpH+ logSWC+ logNP +logclay.silt+
                                logMBC+logBNC+ logFNC+ 
                                logBS+
                                logFS+
                                Network_PC1,
                              data = data1,
                              ntree=1000,
                              mtry =4,
                              importance=TRUE,
                              proximity=TRUE)
RF_model_maoc
summary(RF_model_maoc)
varImpPlot(RF_model_maoc)

#Model significance testing
set.seed(123)
maoc.perm <- rf.significance(RF_model_maoc, na.omit(data1), nperm=99, ntree=1000)
maoc.perm


#Importance_scores
importance_scores <- importance(RF_model_maoc)
importance_scores
write.csv(importance_scores,"results/RF_MAOC.csv")


####----Fig.3 Random forest model--PON----####

data1 <- som %>%
  dplyr::select(Crop.types,
                logMGT, logMGP, logNDVI,
                logpH, logSWC, logCP, logclay.silt,
                logMBN,
                logBS,logFS,
                Network_PC1,
                logPON) 

which(is.na(data1))
str(data1)
data1$Crop.types<-as.factor(data1$Crop.types)


#Parameter estimation
library(e1071)
tune_res <- tune.randomForest(x = data1[,-13], y = data1$logPON,
                              mtry = c(1:12),
                              ntree = seq(100, 1000, by = 100))
tune_res
summary(tune_res)

#Re-modeling
set.seed(123)
RF_model_pon <- randomForest(logPON~logMGT+ logMGP+logNDVI+Crop.types+
                               logpH+ logSWC+ logCP+logclay.silt+
                               logMBN+
                               logBS+
                               logFS+
                               Network_PC1,
                             data = data1,
                             ntree=300,
                             mtry =7,
                             importance=TRUE,
                             proximity=TRUE)
RF_model_pon
summary(RF_model_pon)
varImpPlot(RF_model_pon)

#Model significance testing
set.seed(123)
pon.perm <- rf.significance(RF_model_pon, na.omit(data1), nperm=99, ntree=300)
pon.perm


#Importance_scores
importance_scores <- importance(RF_model_pon)
importance_scores
write.csv(importance_scores,"results/RF_PON.csv")


####----Fig.3 Random forest model--MAON----####

data1 <- som %>%
  dplyr::select(Crop.types,
                logMGT, logMGP, logNDVI,
                logpH, logSWC, logCP, logclay.silt,
                logMBN,logBNN,logFNN, 
                logBS,logFS,
                Network_PC1,
                logMAON) 

which(is.na(data1))
str(data1)
data1$Crop.types<-as.factor(data1$Crop.types)


#Parameter estimation   
library(e1071)
tune_res <- tune.randomForest(x = data1[,-15], y = data1$logMAON,
                              mtry = c(1:14),
                              ntree = seq(100, 1000, by = 100))
tune_res
summary(tune_res)

#Re-modeling  
set.seed(123)
RF_model_maon <- randomForest(logMAON~logMGT+ logMGP+logNDVI+Crop.types+
                                logpH+ logSWC+ logCP+ logclay.silt+
                                logMBN+logBNN+ logFNN+ 
                                logBS+
                                logFS+
                                Network_PC1,
                              data = data1,
                              ntree=300,
                              mtry =3,
                              importance=TRUE,
                              proximity=TRUE)
RF_model_maon
summary(RF_model_maon)
varImpPlot(RF_model_maon)

#Model significance testing  
set.seed(123)
maon.perm <- rf.significance(RF_model_maon, na.omit(data1), nperm=99, ntree=300)
maon.perm

#Importance_scores
importance_scores <- importance(RF_model_maon)
importance_scores
write.csv(importance_scores,"results/RF_MAON.csv")

####----Fig.3 multiple linear regression models--data----####
#data
model<-som %>%
  dplyr::select(Crop.types,
                logMGT, logMGP, logNDVI,
                logpH, logSWC, logNP, logCP,logclay.silt,
                logMBC,logMBN,logBNC,logFNC, logBNN, logFNN,
                logBS,
                logFS,
                Network_PC1,
                logPOC,logPON,logMAOC,logMAON) %>%
  mutate(Crop.types1 = ifelse(Crop.types == "Maize", 1, 0)) %>%  # Create a new column Crop.types1
  mutate(across(-c(Crop.types, Crop.types1), ~ scale(., center = TRUE, scale = TRUE)))  # Standardize all columns except Crop.types and Crop.types1 to mean = 0 and sd = 1

str(model)
names(model)
which(is.na(model))


####----Fig.3 multiple linear regression models--POC----####

model_POC <- lm(logPOC~logMGT+ logMGP+logNDVI+ Crop.types1 +
                  logpH+logSWC+logNP+logclay.silt+
                  logMBC+
                  logBS+
                  logFS+
                  Network_PC1,
                data=model)
summary(model_POC)


##Check collinearity
check_collinearity(model_POC)


##Select best model
dre_model_POC<- dredge(model_POC, 
                       options(na.action = "na.fail"),
                       extra = c("R^2", F = function(x)
                         summary(x)$fstatistic[[1]]))
sub_model_POC <- subset(dre_model_POC, delta < 2) 
summary(sub_model_POC)

##Extract model parameters
mean_R2=mean(sub_model_POC$`R^2`)
mean_AICc=mean(sub_model_POC$`AICc`)
mean_delta_AICc=mean(sub_model_POC$`delta`)
mean_R2
mean_AICc
mean_delta_AICc

##Mean models
aveg_model_POC<-model.avg(dre_model_POC, subset = delta < 2)
summary(aveg_model_POC)
#importance(aveg_model_POC)
sw(aveg_model_POC)

#Extract model coefficients 
a1=summary(aveg_model_POC)
a2 <- as.data.frame(a1$coefmat.full)
a3 <- a2[-1,c(1,2)]
a3$variables <- row.names(a3)
names(a3)[1] <- "coefficients"
names(a3)[2] <- "sd"
a3$proporation <- abs(a3$coefficients)/sum(abs(a3$coefficients))*100
a3
a3$sd <- a3$sd*1.96
a3
str(a3)

##subgroup
subgroup <- read.csv("subgroup.csv", header = T, sep = ",")
comb_date <- a3 %>% left_join (subgroup,
                               by = c('variables' = 'variables1'))
str(comb_date)
write.csv(comb_date,"results/model_POC.csv")

####----Fig.3 multiple linear regression models--MAOC----####

model_MAOC <- lm(logMAOC~logMGT+ logMGP+logNDVI+ Crop.types1 +
                   logpH+logSWC+logNP+logclay.silt+
                   logMBC+logBNC+logFNC+
                   logBS+
                   logFS+
                   Network_PC1,
                 data=model)
summary(model_MAOC)


##Check collinearity 
check_collinearity(model_POC)

##Optimal model selection
dre_model_MAOC<- dredge(model_MAOC, 
                        options(na.action = "na.fail"),
                        extra = c("R^2", F = function(x)
                          summary(x)$fstatistic[[1]]))
sub_model_MAOC <- subset(dre_model_MAOC, delta < 2) 
summary(sub_model_MAOC)

##Extract model parameters  
mean_R2=mean(sub_model_MAOC$`R^2`)
mean_AICc=mean(sub_model_MAOC$`AICc`)
mean_delta_AICc=mean(sub_model_MAOC$`delta`)
mean_R2
mean_AICc
mean_delta_AICc

##Mean models
aveg_model_MAOC<-model.avg(dre_model_MAOC, subset = delta < 2)
summary(aveg_model_MAOC)
sw(aveg_model_MAOC)

#Extract model coefficients 
a1=summary(aveg_model_MAOC)
a2 <- as.data.frame(a1$coefmat.full)
a3 <- a2[-1,c(1,2)]
a3$variables <- row.names(a3)
names(a3)[1] <- "coefficients"
names(a3)[2] <- "sd"
a3$proporation <- abs(a3$coefficients)/sum(abs(a3$coefficients))*100
a3
a3$sd <- a3$sd*1.96
a3
str(a3)

##subgroup

subgroup <- read.csv("subgroup.csv", header = T, sep = ",")
comb_date <- a3 %>% left_join (subgroup,
                               by = c('variables' = 'variables1'))
str(comb_date)

write.csv(comb_date,"results/model_MAOC.csv")


####----Fig.3 multiple linear regression models--PON----####

model_PON <- lm(logPON~logMGT+ logMGP+logNDVI+ Crop.types1 +
                  logpH+logSWC+logCP+logclay.silt+
                  logMBN+
                  logBS+
                  logFS+
                  Network_PC1,
                data=model)
summary(model_PON)


##Check collinearity 
check_collinearity(model_PON)

##Optimal model selection 
dre_model_PON<- dredge(model_PON, 
                       options(na.action = "na.fail"),
                       extra = c("R^2", F = function(x)
                         summary(x)$fstatistic[[1]]))
sub_model_PON <- subset(dre_model_PON, delta < 2) 
summary(sub_model_PON)

##Extract model parameters  
mean_R2=mean(sub_model_PON$`R^2`)
mean_AICc=mean(sub_model_PON$`AICc`)
mean_delta_AICc=mean(sub_model_PON$`delta`)
mean_R2
mean_AICc
mean_delta_AICc

##Mean models
aveg_model_PON<-model.avg(dre_model_PON, subset = delta < 2)
summary(aveg_model_PON)
sw(aveg_model_PON)

#Extract model coefficients 
a1=summary(aveg_model_PON)
a2 <- as.data.frame(a1$coefmat.full)
a3 <- a2[-1,c(1,2)]
a3$variables <- row.names(a3)
names(a3)[1] <- "coefficients"
names(a3)[2] <- "sd"
a3$proporation <- abs(a3$coefficients)/sum(abs(a3$coefficients))*100
a3
a3$sd <- a3$sd*1.96
a3
str(a3)

##subgroup
subgroup <- read.csv("subgroup.csv", header = T, sep = ",")
comb_date <- a3 %>% left_join (subgroup,
                               by = c('variables' = 'variables1'))
str(comb_date)
write.csv(comb_date,"results/model_PON.csv")


####----Fig.3 multiple linear regression models--MAON----####

model_MAON <- lm(logMAON~logMGT+ logMGP+logNDVI+ Crop.types1 +
                   logpH+logSWC+logCP+logclay.silt+
                   logMBN+logBNN+logFNN+
                   logBS+
                   logFS+
                   Network_PC1,
                 data=model)
summary(model_MAON)


##Check collinearity 
check_collinearity(model_POC)

##Optimal model selection 
dre_model_MAON<- dredge(model_MAON, 
                        options(na.action = "na.fail"),
                        extra = c("R^2", F = function(x)
                          summary(x)$fstatistic[[1]]))
sub_model_MAON <- subset(dre_model_MAON, delta < 2) 
summary(sub_model_MAON)

##Extract model parameters  
mean_R2=mean(sub_model_MAON$`R^2`)
mean_AICc=mean(sub_model_MAON$`AICc`)
mean_delta_AICc=mean(sub_model_MAON$`delta`)
mean_R2
mean_AICc
mean_delta_AICc

##Mean models
aveg_model_MAON<-model.avg(dre_model_MAON, subset = delta < 2)
summary(aveg_model_MAON)
sw(aveg_model_MAON)

#Extract model coefficients 
a1=summary(aveg_model_MAON)
a2 <- as.data.frame(a1$coefmat.full)
a3 <- a2[-1,c(1,2)]
a3$variables <- row.names(a3)
names(a3)[1] <- "coefficients"
names(a3)[2] <- "sd"
a3$proporation <- abs(a3$coefficients)/sum(abs(a3$coefficients))*100
a3
a3$sd <- a3$sd*1.96
a3
str(a3)

##subgroup
subgroup <- read.csv("subgroup.csv", header = T, sep = ",")
comb_date <- a3 %>% left_join (subgroup,
                               by = c('variables' = 'variables1'))
str(comb_date)
write.csv(comb_date,"results/model_MAON.csv")

####----Fig.4 SEM--data----####

som2 <- som %>%dplyr::select(ID, Crop.types,
                             logMGT, logMGP, logNDVI,
                             logpH, logSWC, logNP, logCP,logclay.silt,
                             logMBC,logMBN,logBNC,logFNC, logBNN, logFNN,
                             logBS,logFS,Network_PC1,
                             logPOC,logPON,logMAOC,logMAON) %>%
  rename_with(~ gsub("^log", "", .x), .cols = everything()) %>%
  rename(Claysilt = clay.silt, Net = Network_PC1)

dim(som2)

#Create dummy variables using model.matrix
Crop <- model.matrix(~ Crop.types - 1, data = som2)

data <-cbind(som2[,c(1,2)],data.frame(scale(som2[,c(-1,-2)]))) %>% #Standardize all columns to mean = 0 and sd = 1
  cbind(.,Crop) %>%  
  rename(crop=Crop.typesMaize)

which(is.na(data))
str(data)
colnames(data)

####----Fig.4 SEM--C----####

#Define the SEM model
model_C_all <- '
  crop ~ MGT + MGP
  NDVI ~ MGT + MGP + crop
  pH ~ MGT + MGP + NDVI + crop
  SWC ~ MGT + MGP + NDVI + crop
  NP ~ MGT + MGP + NDVI + crop + SWC
  Claysilt ~ MGT + MGP + NDVI + crop + SWC
  Net ~ MGT + MGP + NDVI + crop + pH + SWC + NP + Claysilt
  BS ~ MGT + MGP + NDVI + crop + pH + SWC + NP + Claysilt + Net
  FS ~ MGT + MGP + NDVI + crop + pH + SWC + NP + Claysilt + Net
  MBC ~ MGT + MGP + NDVI + crop + pH + SWC + NP + Claysilt + BS + FS + Net
  BNC ~ MGT + MGP + SWC + Claysilt + MBC + BS + FS + Net
  FNC ~ MGT + MGP + SWC + Claysilt + MBC + BS + FS + Net
  POC ~ MGT + MGP + NDVI + crop + SWC + Claysilt + MBC + BS + FS + Net
  MAOC ~ POC + MGT + MGP + NDVI + crop + SWC + Claysilt + MBC + BNC + FNC + BS + FS + Net
  BNC ~~ FNC
'

#Fit the SEM model
fit_all_C <- sem(model_C_all, data = data, check.gradient = FALSE,estimator="GLS") #check.gradient = FALSE,

fitmeasures(fit_all_C, c("chisq", "df", "pvalue", "gfi", "cfi", "rmr", "srmr", "rmsea"),estimator="GLS")

options(max.print=2000) 
summary(fit_all_C,standardized=TRUE) 
summary(fit_all_C, modindices= T) 
nice_lavaanPlot(fit_all_C)


#Modify the SEM model
model_C <- '
  crop ~ MGT + MGP
  NDVI ~ MGT + MGP + crop + NP
  pH ~ MGT + MGP
  SWC ~ MGT + MGP + crop
  NP ~ MGP  + crop + NDVI + pH + Claysilt
  Claysilt ~ MGT + MGP  + SWC
  Net ~ crop + pH + SWC + NP + Claysilt
  BS ~ MGT + MGP + NDVI + crop + SWC + Net
  FS ~ MGT + crop + NDVI + pH + Claysilt + Net
  MBC ~ MGP + crop + pH + SWC + NP + BS + FS + Net
  BNC ~ MGP + NDVI + crop + pH + NP + SWC + Claysilt + MBC + BS + FS + Net
  FNC ~ MGT + MGP + crop + pH + NP + SWC + Claysilt + MBC + BS + FS + Net
  POC ~ MGP + NDVI + crop + SWC + Claysilt + MBC  + BS + FS + Net
  MAOC ~ POC + MGT + pH + SWC + MBC + BNC + FNC  + BS + FS + Net
  BNC ~~ FNC
  BS ~~ FS
  MGT ~~ MGP
'

#Fit the SEM model
fit_C <- sem(model_C, data = data, check.gradient = FALSE,estimator="GLS") #check.gradient = FALSE,
fitmeasures(fit_C, c("chisq", "df", "pvalue", "gfi", "cfi", "rmr", "srmr", "rmsea"),estimator="GLS")


#pathways
regressions_C <- standardizedSolution(fit_C)[standardizedSolution(fit_C)$op == "~", ] %>%
  mutate(label = if_else(pvalue < 0.001,"***",
                         if_else(pvalue <0.01,"**",
                                 if_else(pvalue <0.05,"*","ns"))))
str(regressions_C)
write.csv(regressions_C,"results/SEM_regressions_C.csv")

#effect
effect_C <- data.frame(coef(fit_C))
write.csv(effect_C,"results/SEM_effect_C.csv")


r2_C <- data.frame(inspect(fit_C, "r2"))
write.csv(r2_C,"results/SEM_r2_C.csv")


####----Fig.4 SEM--N----####

#Define the SEM model
model_N_all <- '
  crop ~ MGT + MGP
  NDVI ~ MGT + MGP + crop
  pH ~ MGT + MGP + NDVI + crop
  SWC ~ MGT + MGP + NDVI + crop
  CP ~ MGT + MGP + NDVI + crop + SWC
  Claysilt ~ MGT + MGP + NDVI + crop + SWC
  Net ~ MGT + MGP + NDVI + crop + pH + SWC + CP + Claysilt
  BS ~ MGT + MGP + NDVI + crop + pH + SWC + CP + Claysilt + Net
  FS ~ MGT + MGP + NDVI + crop + pH + SWC + CP + Claysilt + Net
  MBN ~ MGT + MGP + NDVI + crop + pH + SWC + CP + Claysilt + BS + FS + Net
  BNN ~ MGT + MGP + SWC + Claysilt + MBN + BS + FS + Net
  FNN ~ MGT + MGP + SWC + Claysilt + MBN + BS + FS + Net
  PON ~ MGT + MGP + NDVI + crop + SWC + Claysilt + MBN + BS + FS + Net
  MAON ~ PON + MGT + MGP + NDVI + crop + SWC + Claysilt + MBN + BNN + FNN + BS + FS + Net
  BNN ~~ FNN
'

#Fit the SEM model
fit_all_N <- sem(model_N_all, data = data, check.gradient = FALSE,estimator="GLS") #check.gradient = FALSE,
fitmeasures(fit_all_N, c("chisq", "df", "pvalue", "gfi", "cfi", "rmr", "srmr", "rmsea"),estimator="GLS")

options(max.print=2000)  
summary(fit_all_N,standardized=TRUE) 
summary(fit_all_N, modindices= T)
nice_lavaanPlot(fit_all_N)



#Modify the SEM model
model_N <- '
  crop ~ MGT + MGP
  NDVI ~ MGT + MGP + crop
  pH ~ MGT + MGP
  SWC ~ MGP + crop
  CP ~ MGP + crop + pH
  Claysilt ~ MGT + MGP + SWC
  Net ~ crop + pH + CP + Claysilt + SWC
  BS ~ MGT + MGP + NDVI + crop + SWC + Net
  FS ~ MGT + crop + pH + Claysilt + Net
  MBN ~ MGT + NDVI + crop + pH + SWC  + BS + FS + Net
  BNN ~ MGP + crop + SWC + Claysilt + pH + CP + MBN + BS + FS + Net
  FNN ~ MGP + crop + SWC + Claysilt + pH + CP + MBN + BS + FS + Net
  PON ~ MGP + NDVI + crop + SWC + Claysilt + MBN  + BS + FS + Net
  MAON ~ PON + SWC + Claysilt + CP + MBN + BNN + FNN + BS + FS + Net
  BNN ~~ FNN
  BS ~~ FS
  MGT ~~ MGP
'

#fit SEM model
fit_N <- sem(model_N, data = data, check.gradient = FALSE,estimator="GLS")
fitmeasures(fit_N, c("chisq", "df", "pvalue", "gfi", "cfi", "rmr", "srmr", "rmsea"),estimator="GLS")


#pathways
regressions_N <- standardizedSolution(fit_N)[standardizedSolution(fit_N)$op == "~", ] %>%
  mutate(label = if_else(pvalue < 0.001,"***",
                         if_else(pvalue <0.01,"**",
                                 if_else(pvalue <0.05,"*","ns"))))
str(regressions_N)
write.csv(regressions_N,"results/SEM_regressions_N.csv")

#effect
effect_N <- data.frame(coef(fit_N))
write.csv(effect_N,"results/SEM_effect_N.csv")

r2_N <- data.frame(inspect(fit_N, "r2"))
write.csv(r2_N,"results/SEM_r2_N.csv")



####----Supplementary Fig. 4--Phylum-level relative abundance----#####
ID<-som[,1]

##Bacteria
#bacterial OTU
df.b<-read.csv("OTU_Bac.csv",header = TRUE,row.names = 1, stringsAsFactors = FALSE) %>%
  rownames_to_column(var = "OTUID") %>%  
  gather(key="ID",value = "numb",-OTUID) %>%
  spread(key = "OTUID",value = "numb") 

#Subset samples
df.b1<-left_join(ID, df.b, by=c("ID"="ID")) %>%
  gather(key="OTUID",value = "numb",-ID) %>%
  spread(key = "ID",value = "numb") 

##Filter out OTUs with relative abundance < 0.01 %
rowSums(df.b1[,-1])
colSums(df.b1[,-1])

df.b2 <- df.b1 %>%
  tibble::column_to_rownames(var = names(df.b1)[1]) %>%
  dplyr::mutate(SUM = rowSums(.)) %>%
  dplyr::arrange(desc(SUM)) %>%
  dplyr::mutate(Percent = SUM/sum(SUM))  %>%
  dplyr::filter(Percent > 0.01 * 0.01) %>%  
  dplyr::select(-SUM, -Percent) %>%
  rownames_to_column(var = "OTUID")


taxonomy.b <- read_excel(path = "taxonomy.xlsx",sheet = "bacteria",col_names = T)

tax.b <-as.data.frame( merge(df.b2,taxonomy.b, by = "OTUID")) %>% 
  filter(Kingdom != 'Archaea')  # except Archaea
str(tax.b)

length(unique(tax.b$Phylum))  
colSums(df.b[c(2:445)])  

#phylum
phylum.b<-aggregate(tax.b[c(2:445)],by=list(Phylum=tax.b$Phylum),sum) %>% 
  mutate(across(-1, ~./sum(.)))  


#phyla present at > 400 sites
head(phylum.b)
phylum.b1<-phylum.b[phylum.b %>%
                      tibble::column_to_rownames(var = "Phylum") %>%
                      apply(., 1, function(x){
                        sum(x > 0) > 400
                      }),]

phylum.b2 <- phylum.b1[which(phylum.b1$Phylum != "Unassigned"),] %>% #except unassigned
  gather(key = "ID", value = "per", -Phylum) %>%
  spread(key = "Phylum", value = "per")

phylum.b3 <- phylum.b1 %>% 
  gather(key = "ID", value = "per", -Phylum) %>%
  spread(key = "Phylum", value = "per") %>%
  full_join(.,som[,c(1,2)], by = c("ID"="ID")) %>%
  gather(key = "Phylum", value = "per", -ID, -Crop.types) 
str(phylum.b3)

##Fungi
#fungal OTU
df.f<-read.csv("OTU_Fungi.csv",header = TRUE,row.names = 1, stringsAsFactors = FALSE) %>%
  rownames_to_column(var = "OTUID") %>%  
  gather(key="ID",value = "numb",-OTUID) %>%
  spread(key = "OTUID",value = "numb") 

#Subset samples
df.f1<-left_join(ID, df.f, by=c("ID"="ID")) %>%
  gather(key="OTUID",value = "numb",-ID) %>%
  spread(key = "ID",value = "numb") 


##Filter out OTUs with relative abundance < 0.01 %
rowSums(df.f1[,-1])
colSums(df.f1[,-1])

df.f2 <- df.f1 %>%
  tibble::column_to_rownames(var = names(df.f1)[1]) %>%
  dplyr::mutate(SUM = rowSums(.)) %>%
  dplyr::arrange(desc(SUM)) %>%
  dplyr::mutate(Percent = SUM/sum(SUM))  %>%
  dplyr::filter(Percent > 0.01 * 0.01) %>%  
  dplyr::select(-SUM, -Percent) %>%
  rownames_to_column(var = "OTUID")

taxonomy.f <- read_excel(path = "taxonomy.xlsx",sheet = "fungi",col_names = T)

tax.f <-as.data.frame( merge(df.f2,taxonomy.f, by = "OTUID")) 
str(tax.f)

length(unique(tax.f$Phylum))  
colSums(df.f2[c(2:445)])  

#phylum
phylum.f<-aggregate(tax.f[c(2:445)],by=list(Phylum=tax.f$Phylum),sum) %>% 
  mutate(across(-1, ~./sum(.)))  


#phyla present at > 400 sites
head(phylum.f)

phylum.f1<-phylum.f[phylum.f %>%
                      tibble::column_to_rownames(var = "Phylum") %>%
                      apply(., 1, function(x){
                        sum(x > 0) > 400
                      }),]

phylum.f2 <- phylum.f1[which(phylum.f1$Phylum != "Unassigned"),] %>%  #except unassigned
  gather(key = "ID", value = "per", -Phylum) %>%
  spread(key = "Phylum", value = "per")


phylum.f3 <- phylum.f1 %>% 
  gather(key = "ID", value = "per", -Phylum) %>%
  spread(key = "Phylum", value = "per") %>%
  full_join(.,som[,c(1,2)], by = c("ID"="ID")) %>%
  gather(key = "Phylum", value = "per", -ID, -Crop.types) 
str(phylum.f3)

###Phylum-level relative abundance
#Bateria
str(phylum.b3)
avg.b <- summaryBy(per~Crop.types + Phylum,data = phylum.b3,FUN = fun.mean) %>%
  mutate(mean = per.m*100) %>% 
  arrange(Crop.types, mean) %>% 
  mutate(micro="Bacteria") %>%
  mutate(Phylum = replace(Phylum, Phylum == "Unassigned", "Unassigned1"))

str(avg.b)

str(phylum.f3)
avg.f <- summaryBy(per~Crop.types + Phylum,data = phylum.f3,FUN = fun.mean) %>%
  mutate(mean = per.m*100) %>% 
  arrange(Crop.types, mean) %>% 
  mutate(micro="Fungi")

str(avg.f)

#mixdata
mean.p <- rbind(avg.b, avg.f) %>%
  mutate(micro1=2019)

write.csv(mean.p, "results/phylum.mean.csv")

#plotting
str(mean.p)
mean.p$Phylum <- factor(mean.p$Phylum, levels = unique(mean.p$Phylum)) 

ggplot(data = mean.p,aes(micro1, mean, fill=Phylum))+
  geom_bar(stat="identity", linewidth=5)+
  scale_y_continuous(name = "proporation (%)")+
  facet_grid(Crop.types~micro) + 
  #labs(title = expression("all-bacteria"))+
  labs(x=NULL)+
  scale_color_npg()+
  theme_classic()+
  theme(axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.line.x =element_blank(),
        legend.position="right")

####----Fig.6--heatmap----####
som3<-som2[,c(1,2,20:23,11:16)]
str(som3)

datacor <- full_join(som3,phylum.b2, by=c("ID"="ID")) %>%
  full_join(phylum.f2, by=c("ID"="ID")) 

which(is.na(datacor))
str(datacor)

datacor1 <- datacor[,c(-1,-2)]
which(is.na(datacor1))
str(datacor1)


#Spearman rank correlations 
col2 = colorRampPalette(c('#053061', '#2166AC', '#4393C3', '#92C5DE',
                          '#D1E5F0', '#FFFFFF', '#FDDBC7', '#F4A582',
                          '#D6604D', '#B2182B', '#67001F'))

M <- cor(datacor1, method = "spearman")
p.mat <- cor_pmat(datacor1)


corrplot(M, col = col2(100),type = "lower", 
         diag = FALSE, insig = "blank", 
         method = 'square',p.mat = p.mat)
corrplot(M, col = col2(100),type = "lower", 
         diag = FALSE,  
         method = 'square',p.mat = p.mat, insig = "label_sig",
         sig.level = c(.001, .01, .05),pch.cex = 0.6, tl.cex = 1,cl.cex = 1,
         pch.col = "black")

## r and p
M1<-as.data.frame(M) %>%
  rownames_to_column(var = "phylum") %>%
  slice(11:26) %>% 
  .[, 1:11] %>% 
  gather(key = "Type", value = "value",-phylum)
str(M1)

p.mat1<-as.data.frame(p.mat) %>%
  rownames_to_column(var="phylum") %>%
  slice(11:26) %>%
  .[, 1:11] %>%
  gather(key = "Type", value = "value",-phylum) %>%
  mutate(significance = case_when(
    value >= 0 & value < 0.001 ~ "***",
    value >= 0.001 & value < 0.01 ~ "**",
    value >= 0.01 & value < 0.05 ~ "*",
    value >= 0.05 & value <= 1 ~ " ",
    TRUE ~ NA_character_  
  ))
str(p.mat1)

####----Fig.6--random forest model----####

str(datacor1)

##POC
datacor2 <- datacor1[,c(1,11:26)]
str(datacor2)
dim(datacor2)

#Parameter estimation
tune_res <- tune.randomForest(x = datacor2[,-1], y = datacor2$POC,
                              mtry = c(2:17),
                              ntree = seq(100, 1000, by = 100))
tune_res
#summary(tune_res)

#Build the model
form_reg <- as.formula(paste0("POC ~",paste(colnames(datacor2)[2:17],collapse = "+")))
form_reg
set.seed(123)
POC <- randomForest(form_reg, data=datacor2, 
                    mtry = 2,
                    ntree = 500,
                    importance = TRUE,
                    proximity = TRUE)
POC
varImpPlot(POC)
POC$importance

#Model significance testing
set.seed(123)
poc.perm <- rf.significance(POC, na.omit(datacor2), nperm=99, ntree=500) 
poc.perm


##PON
datacor2 <- datacor1[,c(2,11:26)]
str(datacor2)

#Parameter estimation
tune_res <- tune.randomForest(x = datacor2[,-1], y = datacor2$PON,
                              mtry = c(2:17),
                              ntree = seq(100, 1000, by = 100))
tune_res
#summary(tune_res)

#Build the model
form_reg <- as.formula(paste0("PON ~",paste(colnames(datacor2)[2:17],collapse = "+")))
form_reg
set.seed(123)
PON <- randomForest(form_reg, data=datacor2, 
                    mtry = 5,
                    ntree = 400,
                    importance = TRUE,
                    proximity = TRUE)
PON
varImpPlot(PON)
PON$importance

#Model significance testing
set.seed(123)
pon.perm <- rf.significance(PON, na.omit(datacor2), nperm=99, ntree=400) 
pon.perm

##MAOC
datacor2 <- datacor1[,c(3,11:26)]
str(datacor2)

#Parameter estimation 
tune_res <- tune.randomForest(x = datacor2[,-1], y = datacor2$MAOC,
                              mtry = c(2:17),
                              ntree = seq(100, 1000, by = 100))
tune_res
#summary(tune_res)

# Build the model 
form_reg <- as.formula(paste0("MAOC ~",paste(colnames(datacor2)[2:17],collapse = "+")))
form_reg
set.seed(123)
MAOC <- randomForest(form_reg, data=datacor2, 
                     mtry = 7,
                     ntree = 200,
                     importance = TRUE,
                     proximity = TRUE)
MAOC
varImpPlot(MAOC)
MAOC$importance

# Model significance testing
set.seed(123)
maoc.perm <- rf.significance(MAOC, na.omit(datacor2), nperm=99, ntree=200) 
maoc.perm


##MAON
datacor2 <- datacor1[,c(4,11:26)]
str(datacor2)

#Parameter estimation 
tune_res <- tune.randomForest(x = datacor2[,-1], y = datacor2$MAON,
                              mtry = c(2:17),
                              ntree = seq(100, 1000, by = 100))
tune_res
#summary(tune_res)

#Build the model  
form_reg <- as.formula(paste0("MAON ~",paste(colnames(datacor2)[2:17],collapse = "+")))
form_reg
set.seed(123)
MAON <- randomForest(form_reg, data=datacor2, 
                     mtry = 2,
                     ntree = 400,
                     importance = TRUE,
                     proximity = TRUE)
MAON
varImpPlot(MAON)
MAON$importance

#Model significance testing
set.seed(123)
maon.perm <- rf.significance(MAON, na.omit(datacor2), nperm=99, ntree=400) 
maon.perm


##MBC
datacor2 <- datacor1[,c(5,11:26)]
str(datacor2)

#Parameter estimation 
tune_res <- tune.randomForest(x = datacor2[,-1], y = datacor2$MBC,
                              mtry = c(2:17),
                              ntree = seq(100, 1000, by = 100))
tune_res
#summary(tune_res)

#Build the model  
form_reg <- as.formula(paste0("MBC ~",paste(colnames(datacor2)[2:17],collapse = "+")))
form_reg
set.seed(123)
MBC <- randomForest(form_reg, data=datacor2, 
                    mtry = 2,
                    ntree = 200,
                    importance = TRUE,
                    proximity = TRUE)
MBC
varImpPlot(MBC)

MBC$importance

#Model significance testing
set.seed(123)
mbc.perm <- rf.significance(MBC, na.omit(datacor2), nperm=99, ntree=200) 
mbc.perm


##MBN
datacor2 <- datacor1[,c(6,11:26)]
str(datacor2)

#Parameter estimation 
tune_res <- tune.randomForest(x = datacor2[,-1], y = datacor2$MBN,
                              mtry = c(2:17),
                              ntree = seq(100, 1000, by = 100))
tune_res
#summary(tune_res)

#Build the model  
form_reg <- as.formula(paste0("MBN ~",paste(colnames(datacor2)[2:17],collapse = "+")))
form_reg
set.seed(123)
MBN <- randomForest(form_reg, data=datacor2, 
                    mtry = 2,
                    ntree = 400,
                    importance = TRUE,
                    proximity = TRUE)
MBN
varImpPlot(MBN)

MBN$importance

#Model significance testing
set.seed(123)
mbn.perm <- rf.significance(MBN, na.omit(datacor2), nperm=99, ntree=400) 
mbn.perm


##BNC
datacor2 <- datacor1[,c(7,11:26)]
str(datacor2)

#Parameter estimation 
tune_res <- tune.randomForest(x = datacor2[,-1], y = datacor2$BNC,
                              mtry = c(2:17),
                              ntree = seq(100, 1000, by = 100))
tune_res
#summary(tune_res)

#Build the model  
form_reg <- as.formula(paste0("BNC ~",paste(colnames(datacor2)[2:17],collapse = "+")))
form_reg
set.seed(123)
BNC <- randomForest(form_reg, data=datacor2, 
                    mtry = 6,
                    ntree = 100,
                    importance = TRUE,
                    proximity = TRUE)
BNC
varImpPlot(BNC)
BNC$importance

#Model significance testing
set.seed(123)
bnc.perm <- rf.significance(BNC, na.omit(datacor2), nperm=99, ntree=100) 
bnc.perm


##FNC
datacor2 <- datacor1[,c(8,11:26)]
str(datacor2)

#Parameter estimation 
tune_res <- tune.randomForest(x = datacor2[,-1], y = datacor2$FNC,
                              mtry = c(2:17),
                              ntree = seq(100, 1000, by = 100))
tune_res
#summary(tune_res)

# Build the model  
form_reg <- as.formula(paste0("FNC ~",paste(colnames(datacor2)[2:17],collapse = "+")))
form_reg
set.seed(123)
FNC <- randomForest(form_reg, data=datacor2, 
                    mtry = 6,
                    ntree = 800,
                    importance = TRUE,
                    proximity = TRUE)
FNC
varImpPlot(FNC)
FNC$importance

#Model significance testing
set.seed(123)
fnc.perm <- rf.significance(FNC, na.omit(datacor2), nperm=99, ntree=800) 
fnc.perm


##BNN
datacor2 <- datacor1[,c(9,11:26)]
str(datacor2)

#Parameter estimation 
tune_res <- tune.randomForest(x = datacor2[,-1], y = datacor2$BNN,
                              mtry = c(2:17),
                              ntree = seq(100, 1000, by = 100))
tune_res
#summary(tune_res)

# Build the model  
form_reg <- as.formula(paste0("BNN ~",paste(colnames(datacor2)[2:17],collapse = "+")))
form_reg
set.seed(123)
BNN <- randomForest(form_reg, data=datacor2, 
                    mtry = 6,
                    ntree = 100,
                    importance = TRUE,
                    proximity = TRUE)
BNN
varImpPlot(BNN)
BNN$importance

#Model significance testing
set.seed(123)
bnn.perm <- rf.significance(BNN, na.omit(datacor2), nperm=99, ntree=100) 
bnn.perm


##FNN
datacor2 <- datacor1[,c(10,11:26)]
str(datacor2)

#Parameter estimation 
tune_res <- tune.randomForest(x = datacor2[,-1], y = datacor2$FNN,
                              mtry = c(2:17),
                              ntree = seq(100, 1000, by = 100))
tune_res
#summary(tune_res)

#Build the model  
form_reg <- as.formula(paste0("FNN ~",paste(colnames(datacor2)[2:17],collapse = "+")))
form_reg
set.seed(123)
FNN <- randomForest(form_reg, data=datacor2, 
                    mtry = 4,
                    ntree = 700,
                    importance = TRUE,
                    proximity = TRUE)
FNN
varImpPlot(FNN)
FNN$importance

#Model significance testing
set.seed(123)
fnn.perm <- rf.significance(FNN, na.omit(datacor2), nperm=99, ntree=700) 
fnn.perm


##Importance %IncMSE
df_Importance <-as.data.frame(importance(POC))
df_Importance <- df_Importance[-2]
colnames(df_Importance) <- "POC"
df_Importance$PON <- importance(PON)[1:16]
df_Importance$MAOC <- importance(MAOC)[1:16]
df_Importance$MAON <- importance(MAON)[1:16]
df_Importance$MBC <- importance(MBC)[1:16]
df_Importance$MBN <- importance(MBN)[1:16]
df_Importance$BNC <- importance(BNC)[1:16]
df_Importance$FNC <- importance(FNC)[1:16]
df_Importance$BNN <- importance(BNN)[1:16]
df_Importance$FNN <- importance(FNN)[1:16]

#Data filtering—remove all values ≤ 0 and missing entries
df_Importance[df_Importance<0] <- 0
df_Importance[is.na(df_Importance<0)] <- 0
write.csv(df_Importance,file = "results/RF-importance.csv")

df_Importance1 <- df_Importance %>%
  rownames_to_column(var = "phylum") %>%
  gather(key = "Type", value = "value",-phylum)

max(df_Importance1$value)
str(df_Importance1)

####plotting
M1$phylum <- factor(M1$phylum, levels = unique(M1$phylum)) 
M1$Type <- factor(M1$Type, levels = rev(M1$Type))
df_Importance1$Type <- factor(df_Importance1$Type, order = T,
                              levels = c(FNN,BNN,FNC,BNC,MBN,MBC,MAON,MAOC,PON,POC)) 

picture1 <- ggplot() + 
  geom_tile(data=M1, aes( y = Type, x = phylum, fill = value), 
            color = "black", height = 1, width = 1) +
  geom_text(data = p.mat1, mapping = aes( y = Type, x = phylum, label = significance), 
            size = 5,vjust=-0.2) + 
  scale_fill_gradient2(low = "#2D6DB1", high = "#DC1623", mid = "#FFFFFF", midpoint = 0,
                       breaks=c(-0.4,-0.2,0,0.2,0.4),labels=c(-0.4,-0.2,0,0.2,0.4),limits = c(-0.4,0.4)) +
  
  geom_point(data=df_Importance1[df_Importance1$value>0,], aes( y = Type, x = phylum, size = value), 
             shape = 1, position = position_nudge(y = -0.2))+
  scale_size_continuous(range = c(0.01,9))+ 
  scale_y_discrete(limits = rev)  + 
 
  labs(x = NULL, y = NULL,  size = "Importance(%)", fill = "Correlation \n coefficient")+
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text = element_text(size=10,color = "black",face = "bold"),
        legend.text = element_text(face="bold", size=10),
        legend.title = element_text(face="bold", size=10,color = "black"),
        legend.position = "left",
        panel.grid = element_blank(),
        plot.background = element_rect(linetype = "solid",linewidth =1),
        panel.border = element_rect(color = "black", size = 1, fill = NA) )

picture1
ggsave(file="results/micro-heatmap.pdf", plot = picture1, width = 11,height = 5)