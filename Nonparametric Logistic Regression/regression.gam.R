remove(list=ls())
setwd("C:/Users/E5440/Desktop/esami/Nonparametric statistics/Progetto/Alzheimer")
library(dplyr)
library(mgcv)
library(rgl)
library(splines)
library(pbapply)
library(ISLR2)
library(car)
library(ggplot2)
library(yardstick)

dataset<- read.csv("oasis_longitudinal.csv")
dataset <- filter(dataset, Visit==1)

#EDUC: years of education 
#SES: socio economic status
#MMSE: mini mental state examination (0:30 30=stato mentale perfetto. Anche persone con un iniziale deterioramento cognitivo, ma con un'alta scolarizzazione possono ottenere un punteggio pari a 29 e 30)
#CDR: clinical dementia rating (0:3 3=perdita di memoria grave)
#eTIV: extimated total intracranial volume 
#nWBV: Normalized Whole Brain Volume
#ASF: Atlas Scaling Factor (the volume-scaling factor required to match each individual to the atlas target)

#LABEL = 1 DEMENTED
data <- dataset[,c(6,8,9,11,13,14)]
i.converted <- which(dataset$Group=='Converted')
data <- data[-i.converted,]
data.aux <- dataset[-i.converted,]
data$label <- rep('Nondem', dim(data)[1])
i.dem <- which(data.aux$Group=='Demented')
data$label[i.dem] <- 'Dem'
data$label <- factor(data$label)
data$M <- ifelse(data$M.F=='M',1,0)

dataset2 <- read.table('C:/Users/E5440/Downloads/test set.csv', header = TRUE, sep=';')
dataset2$label <- ifelse(dataset2$CDR > 0, 'Dem', 'Nondem')
dataset2$label <- factor(dataset2$label)
dataset2  <- dataset2[,c(2,4,5,7,9,10,12)]
dataset2$M <- ifelse(dataset2$M.F=='M',1,0)
colnames(dataset2) <- colnames(data)
data <- rbind(data, dataset2) 
data$label <- ifelse(data$label=='Dem', 1, 0)
data$label <- factor(data$label)

train <- data[1:271,]
test <- data[272:dim(data)[1],]


x11()
with(data = data, scatterplotMatrix(data.frame(Age, EDUC, MMSE, eTIV, nWBV), groups = label))

x11()
with(data = data, scatterplotMatrix(data.frame(Age, EDUC, MMSE, eTIV, nWBV), groups = M.F))

#se le distribuzioni sono simili nelle diverse label, la dummy non è significativa, se sono shiftate ma la forma è uguale la dummy influeisce solo nell'intercetta

model_gam=gam(label ~ M + s(EDUC,bs='cr')+ s(I(EDUC*M), bs='cr') + s(nWBV,bs='cr') + s(I(nWBV*M), bs='cr')  + s(Age, bs='cr') + s(I(Age*M), bs='cr') + s(eTIV, bs='cr') + s(I(eTIV*M), bs='cr') , data = train, select = TRUE, family = binomial)
summary(model_gam)
1 - model_gam$deviance/model_gam$null.deviance
model_gam$aic

model_gam2=gam(label ~ M + s(EDUC,bs='cr')+ s(I(EDUC*M), bs='cr') + s(nWBV,bs='cr') + s(Age, bs='cr') + s(eTIV, bs='cr') + s(I(eTIV*M), bs='cr') , data = train, select = TRUE, family = binomial)
summary(model_gam2)
1 - model_gam2$deviance/model_gam2$null.deviance
model_gam2$aic
#sembra che etiv possa essere modellizzata senza la dummy M.F
model_gam3=gam(label ~  M + s(EDUC,bs='cr')+ s(I(EDUC*M), bs='cr') + s(nWBV,bs='cr') + s(Age, bs='cr') + s(eTIV, bs='cr') , data = train, select = TRUE, family = binomial)
summary(model_gam3)#senza educ, con dummy per etiv

anova.gam(model_gam2, model_gam, test = 'Chisq')
anova.gam(model_gam3, model_gam, test = 'Chisq')

#il migliore sembra model_gam2

roc.curve <- function(predicted, test.set){
  p0 <- seq(0,1,by=0.001)
  spec <- NULL
  sens <- NULL
  accuracy <- NULL
  precision <- NULL
  F1.score <- NULL
  for (i in 1:length(p0)) {
    predicted$class.assigned <- rep(0, dim(predicted)[1])
    colnames(predicted) <- c('prob','class.assigned')
    predicted[which(predicted$prob>=p0[i]),2] <- 1
    true.lab <- test.set$label
    i.equal <- which(true.lab==predicted$class.assigned)
    n.equal <- length(i.equal)
    test.equal <- test.set[i.equal,]
    n00 <- length(which(test.equal$label==0))
    n11 <- n.equal - n00
    
    n01 <- length(which(predicted$class.assigned==0&true.lab==1)) #classified,true lab
    n10 <- length(which(predicted$class.assigned==1&true.lab==0)) #classified,true lab
    
    sensitivity <- n11/(n01+n11)
    specificity <- n00/(n00+n10)
    acc <- (n00 + n11)/(n00 + n11 + n01 + n10)
    prec <- n11 / (n11 + n10)
    F1.s <- 2*prec*sensitivity/(prec + sensitivity)
    spec <- c(spec, specificity)
    sens <- c(sens, sensitivity)
    accuracy <- c(accuracy, acc)
    precision <- c(precision, prec)
    F1.score <- c(F1.score, F1.s)
  }
  x.roc <- rep(1,length(spec))-spec
  y.roc <- sens
  data.frame(cbind(x.roc,y.roc,accuracy, precision, F1.score))
}
plot.roc <- function(data.roc, pred.prob){
  ggplot(data.roc,aes(x.roc,y.roc))+
    theme(panel.background = element_rect(fill = 'gray95', colour = 'white'))+
    xlim(c(0,1)) + ylim(c(0,1))+
    geom_path(col='dodgerblue4', lwd=1)+ #geom_path
    geom_line(data = data.frame(cbind(X1=seq(0,1,0.001), X2=seq(0,1,0.001))), aes(X1,X2), col='orangered', lty = 2, lwd=1.3)+
    xlab('1 - specificity') + ylab('sensitivity') + labs(title = 'ROC Curve', subtitle = paste('area under curve is', round(roc_auc_vec(truth = test$label, estimate = pred.prob, event_level = 'second'), digits = 2))) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) 
  
}

pred <- predict(model_gam2, newdata = test, type = 'response')
pred <- as.data.frame(pred)

#ROC curve
#NB è PIù IMPORTANTE MASSIMIZZARE SENSITIVIT
data.roc <- roc.curve(pred, test)
x11()
plot.roc(data.roc = data.roc, pred.prob = pred$pred) + geom_point(aes(x.roc[i1.opt],y.roc[i1.opt]), col = 'purple', cex = 3)

roc_auc_vec(truth = test$label, estimate = pred$pred, event_level = 'second') #da aggiungere nel plot

which(data.roc$y.roc>=0.75 & data.roc$x.roc<=0.4)
data.roc$y.roc[which(data.roc$y.roc>=0.75 & data.roc$x.roc<=0.4)]
data.roc$x.roc[which(data.roc$y.roc>=0.75 & data.roc$x.roc<=0.4)]

i1.opt <- which(data.roc$y.roc>=0.82 & data.roc$x.roc<=0.4)[1]
p1.opt <- seq(0,1,by=0.001)[i1.opt]

sensitivity1.opt <- data.roc$y.roc[i1.opt]
specificity1.opt <- 1 - data.roc$x.roc[i1.opt]
accuracy1.opt <- data.roc$accuracy[i1.opt]
precision1.opt <- data.roc$precision[i1.opt]
F1.score1.opt <- data.roc$F1.score[i1.opt]

#PROVO CON RUBUST
mdl.rob <- glmrob(label ~ M + bs(Age) + bs(nWBV) + bs(eTIV), data = train, method = 'Mqle', family = binomial)
pred.rob <- predict(mdl.rob, newdata = test, type = 'response')
pred.rob <- as.data.frame(pred.rob)
data.roc.rob <- roc.curve(pred.rob, test)
x11()
plot.roc(data.roc = data.roc.rob, pred.prob = pred.rob$pred.rob)



##############################################################################################################################################################
##############################################################################################################################################################

#GLM MODEL
library(tidyverse)
library(caret)
library(leaps)
library(MASS)
library(RcmdrMisc)

model_glm=glm(label ~ M + Age + nWBV + MMSE + eTIV, family = 'binomial', data = train)
summary(model_glm)

1-model_glm$deviance/model_glm$null.deviance


#MODEL SELECTION 
# mod <- stepwise(model_glm, trace = 1,criterion='BIC' ) #criterion=BIC. quello che è considerato come Aic in realtà è Bic, Aik=205.02
# summary(mod)
# 1-mod$deviance/mod$null.deviance
# mod$deviance-model_glm$deviance
# mod$deviance/model_glm$deviance 

pred <- predict(model_glm, newdata = test, type = 'response')
pred <- as.data.frame(pred)
data.roc <- roc.curve(pred,test)
x11()
plot.roc(data.roc = data.roc)
roc_auc_vec(truth = test$label, estimate = pred$pred, event_level = 'second')#da aggiungere nel plot

model_glm2=glm(label ~ Age + nWBV + MMSE, family = 'binomial' , data = train)
summary(model_glm2)

anova(model_glm, model_glm2, test = 'Chisq')

pred <- predict(model_glm2, newdata = test, type = 'response')
pred <- as.data.frame(pred)
data.roc <- roc.curve(pred,test)
x11()
ggplot(data.roc,aes(x.roc,y.roc))+
  theme(panel.background = element_rect(fill = 'gray90', colour = 'white'))+
  xlim(c(0,1)) + ylim(c(0,1))+
  geom_line(col='dodgerblue4', lwd=1.9)+
  geom_ribbon(aes(ymin=0, ymax=y.roc), fill = 'seagreen2')+
  geom_line(data = data.frame(cbind(X1=seq(0,1,0.001), X2=seq(0,1,0.001))), aes(X1,X2), col='orangered', lty = 2, lwd=1.3)+
  xlab('1 - specificity')+ylab('sensitivity')
roc_auc_vec(truth = test$label, estimate = pred$pred, event_level = 'second')#da aggiungere nel plot


######################################################################################################################
######################################################################################################################

