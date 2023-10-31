
setwd("C:/Users/franc/Desktop/non-param")
library(robustbase)
library(psych)
library(MASS)
library(ellipse)
library(here)
library(DescTools)
library(knitr)
library(RobStatTM)
library(mgcv)
library(splines)
library(dplyr)
library(mgcv)
library(rgl)
library(pbapply)
library(ISLR2)
library(car)
library(ggplot2)
library(yardstick)

dataset<-read.table("oasis_longitudinal.csv",header=T,sep=",")
dataset
casestudy=dataset[which(dataset$Group!="Converted"),]
i1=which(casestudy$Group=="Demented")
New=vector(mode="logical",length=336)
New[i1]=1
New[-i1]=0
n=cbind(casestudy,New)

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

fit_glm <- glm( New  ~ Age + nWBV + ASF, data = casestudy,family=binomial )
plot(n)
abline(fit_glm, col="red", lwd=2) #non la plotta
plot(fit_glm)

#228,229,271 sembrano outliers cioè paziente n  OAS2_0108 entrambe le visite e paziente  n 		
# OAS2_0129 terza visita
View(n)

pairs(n[,8:16])


fit_MCD <- covMcd(x = n, alpha = .8)
fit_MCD

plot(fit_MCD, classic=TRUE, labels.id=FALSE, which="distance")
plot(fit_MCD,labels.id=FALSE,which=c("dd"))

data <- read.csv('data.csv', header = TRUE, sep = ',')
data$label <- ifelse(data$label=='Dem', 1, 0)
data$label <- factor(data$label)
train <- data[1:271,]
test <- data[272:371,]
data <- data[,-1]
train <- train[,-1]
test <- test[,-1]


#######modellO con solo MMSE
mod1<-glmrob(label ~ MMSE, data = train, family = binomial)
summary(mod1)
#sia intercetta che MMSE molto significativi

predizione1 <- predict(mod1, newdata = test, type = 'response')
predizione1 <- as.data.frame(predizione1)
data.roc1 <- roc.curve(predizione1, test)
x11()
plot.roc(data.roc = data.roc1, pred.prob = predizione1$predizione1)
roc_auc_vec(truth = test$label, estimate = predizione1$predizione1, event_level = 'second') #non mi funziona il plot

sensitivity=data.roc1$y.roc[575]
sensitivity
specificity=1-data.roc1$x.roc[575]
specificity
data.roc1$accuracy[575]
data.roc1$precision[575]
data.roc1$F1.score[575]

Precision<-array(data=NA, dim= 15)
Recall<- array(data=NA, dim=15)
F1.score.pred1<-array(data=NA,dim=15)

specificity.funct <- function(predicted, test.set){
  p0 <- seq(0,1,by=0.001)
  spec <- NULL
  sens <- NULL
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
    spec <- c(spec, specificity)
    sens <- c(sens, sensitivity)
  }
  data.frame(cbind(spec,sens))
}

spec1<-specificity.funct(predizione1,test)

k=1
s<-seq(from=0.10, to=0.80, by=0.05)
for (j in s){
  label.pred<-ifelse(predizione1>j, 1, 0)
  label.pred <- factor(label.pred)
  True.Positive=0
  True.Negative=0
  False.Positive=0
  False.Negative=0
  for (i in 1:100){
    if(test$label[i]==1 && label.pred[i]==test$label[i])
      True.Positive=True.Positive+1
    if(test$label[i]==0 && label.pred[i]==test$label[i])
      True.Negative=True.Negative+1
    if(test$label[i]==1 && label.pred[i]==0)
      False.Negative=False.Negative+1
    if(test$label[i]==0 && label.pred[i]==1)
      False.Positive=False.Positive+1
  }
  
  Precision[k] = True.Positive / (True.Positive + False.Positive) 
  Recall[k] = True.Positive / (True.Positive + False.Negative)
  F1.score.pred1[k]=2 * (Precision[k] * Recall[k]) / (Precision[k] + Recall[k])
  k=k+1
}
F1.score.pred1

AIC(predizione1)
#introduco splines
mod2<-glmrob(label ~ bs(MMSE), data = train, family = binomial)
summary(mod2)
#perdono tutti di significatività tranne bs(MMSE)2, che è al limite

predizione2 <- predict(mod2, newdata = test, type = 'response')
predizione2 <- as.data.frame(predizione2)
data.roc2 <- roc.curve(predizione2, test)
x11()
plot.roc(data.roc = data.roc2, pred.prob = predizione2$predizione2)
roc_auc_vec(truth = test$label, estimate = predizione2$predizione2, event_level = 'second') #non mi funziona il plot


spec2<-specificity.funct(predizione2,test)

F1.score.pred2<-array(data=NA,dim=15)

k=1
s<-seq(from=0.10, to=0.80, by=0.05)
for (j in s){
  label.pred<-ifelse(predizione2>j, 1, 0)
  label.pred <- factor(label.pred)
  True.Positive=0
  True.Negative=0
  False.Positive=0
  False.Negative=0
  for (i in 1:100){
    if(test$label[i]==1 && label.pred[i]==test$label[i])
      True.Positive=True.Positive+1
    if(test$label[i]==0 && label.pred[i]==test$label[i])
      True.Negative=True.Negative+1
    if(test$label[i]==1 && label.pred[i]==0)
      False.Negative=False.Negative+1
    if(test$label[i]==0 && label.pred[i]==1)
      False.Positive=False.Positive+1
  }
  
  Precision[k] = True.Positive / (True.Positive + False.Positive) 
  Recall[k] = True.Positive / (True.Positive + False.Negative)
  F1.score.pred2[k]=2 * (Precision[k] * Recall[k]) / (Precision[k] + Recall[k])
  k=k+1
}
F1.score.pred2

#introduco tutte le variabili con le splines

mod3<-glmrob(label ~ M + bs(EDUC) + bs(nWBV) + bs(Age)  + bs(MMSE) + bs(eTIV) , data = train, family = binomial)
summary(mod3)
#significative M, e leggermente EDUC e eTIV


predizione3 <- predict(mod3, newdata = test, type = 'response')

predizione3 <- as.data.frame(predizione3)

spec3<-specificity.funct(predizione3,test)

F1.score.pred3<-array(data=NA,dim=15)

k=1
s<-seq(from=0.10, to=0.80, by=0.05)
for (j in s){
  label.pred<-ifelse(predizione3>j, 1, 0)
  label.pred <- factor(label.pred)
  True.Positive=0
  True.Negative=0
  False.Positive=0
  False.Negative=0
  for (i in 1:100){
    if(test$label[i]==1 && label.pred[i]==test$label[i])
      True.Positive=True.Positive+1
    if(test$label[i]==0 && label.pred[i]==test$label[i])
      True.Negative=True.Negative+1
    if(test$label[i]==1 && label.pred[i]==0)
      False.Negative=False.Negative+1
    if(test$label[i]==0 && label.pred[i]==1)
      False.Positive=False.Positive+1
  }
  
  Precision[k] = True.Positive / (True.Positive + False.Positive) 
  Recall[k] = True.Positive / (True.Positive + False.Negative)
  F1.score.pred3[k]=2 * (Precision[k] * Recall[k]) / (Precision[k] + Recall[k])
  k=k+1
}
F1.score.pred3

#provo a diminuire i gradi delle splines

mod4<-glmrob(label ~ M + bs(EDUC) + bs(nWBV,degree=2) + bs(Age,degree=2)  + bs(MMSE,degree=2) + bs(eTIV) , data = train, family = binomial)
summary(mod4)
#modello migliora significativamente
predizione4 <- predict(mod4, newdata = test, type = 'response')
predizione4 <- as.data.frame(predizione4)
data.roc.4<- roc.curve(predizione4, test)
x11()
plot.roc(data.roc = data.roc.4, pred.prob = predizione4$predizione4)
roc_auc_vec(truth = test$label, estimate = predizione4$predizione4, event_level = 'second') 

data.roc.4$F1.score[312]
data.roc.4$accuracy[312]
data.roc.4$precision[312]
sensitivity=data.roc.4$y.roc[312]
specificity=1-data.roc.4$x.roc[312]
sensitivity
specificity
spec4<-specificity.funct(predizione4,test)

F1.score.pred4<-array(data=NA,dim=15)

k=1
s<-seq(from=0.10, to=0.80, by=0.05)
for (j in s){
  label.pred<-ifelse(predizione4>j, 1, 0)
  label.pred <- factor(label.pred)
  True.Positive=0
  True.Negative=0
  False.Positive=0
  False.Negative=0
  for (i in 1:100){
    if(test$label[i]==1 && label.pred[i]==test$label[i])
      True.Positive=True.Positive+1
    if(test$label[i]==0 && label.pred[i]==test$label[i])
      True.Negative=True.Negative+1
    if(test$label[i]==1 && label.pred[i]==0)
      False.Negative=False.Negative+1
    if(test$label[i]==0 && label.pred[i]==1)
      False.Positive=False.Positive+1
  }
  
  Precision[k] = True.Positive / (True.Positive + False.Positive) 
  Recall[k] = True.Positive / (True.Positive + False.Negative)
  F1.score.pred4[k]=2 * (Precision[k] * Recall[k]) / (Precision[k] + Recall[k])
  k=k+1
}
F1.score.pred4
#####
mod5<-glmrob(label ~ M + bs(EDUC,degree=2) + bs(nWBV,degree=2) + bs(Age,degree=2)  + bs(MMSE,degree=2) + bs(eTIV) , data = train, family = binomial)
summary(mod5)
#EDUc non sembra significativo

mod6<-glmrob(label ~ M + bs(nWBV,degree=2) + bs(Age,degree=2)  + bs(MMSE,degree=2) + bs(eTIV) , data = train, family = binomial)
summary(mod6)
#Age non sembra significativio

mod7<-glmrob(label ~ M + bs(nWBV,degree=2) + bs(MMSE,degree=2) + bs(eTIV) , data = train, family = binomial)
summary(mod7)

##nWBV non molto significativo

mod8<-glmrob(label ~ M + bs(MMSE,degree=2) + bs(eTIV) , data = train, family = binomial)
summary(mod8)
predizione8 <- predict(mod8, newdata = test, type = 'response')
predizione8 <- as.data.frame(predizione8)
data.roc.8<- roc.curve(predizione8, test)
x11()
plot.roc(data.roc = data.roc.8, pred.prob = predizione8$predizione8)
roc_auc_vec(truth = test$label, estimate = predizione8$predizione8, event_level = 'second') 

data.roc.8$F1.score[317]
data.roc.8$accuracy[317]
data.roc.8$precision[317]

#questo modello mi piace


########
modello<-glmrob(label ~ M + EDUC + nWBV + Age + MMSE + eTIV , data = data, family = binomial)
summary(modello)

modello1<-glmrob(label ~ M + EDUC + nWBV + Age + MMSE, data = data, family = binomial)
summary(modello1)

modello2<-glmrob(label ~ M + nWBV + Age + MMSE, data = data, family = binomial)
summary(modello2)

modello3<-glmrob(label ~ M + nWBV + MMSE, data = data, family = binomial)
summary(modello3)

modello4<-glmrob(label ~ nWBV + Age + MMSE, data = data, family = binomial)
summary(modello4)

modello5<-glmrob(label ~ nWBV + MMSE, data = data, family = binomial)
summary(modello5)


#Call:  glmrob(formula = label ~ nWBV + MMSE, family = binomial, data = data) 


#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  32.4119     4.2993   7.539 4.74e-14 ***
#  nWBV        -12.7408     4.1353  -3.081  0.00206 ** 
#  MMSE         -0.8341     0.1111  -7.506 6.12e-14 ***
  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#Robustness weights w.r * w.x: 
 # 329 weights are ~= 1. The remaining 42 ones are summarized as
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2529  0.5553  0.6781  0.6894  0.8716  0.9965 

#Number of observations: 371 
#Fitted by method 'Mqle'  (in 5 iterations)


#w.r sono i pesi di robustezza di ogni osservazione (i.e residual*w.r è uguale ai residui di Preason della psi function)
robustness_weights=modello5$w.r
#ci sono 329 pesi che sono circa 1


#divido i dati in train e test e fitto il modello solo su train, poi faccio prediction con i tes data
modello6<-glmrob(label ~ nWBV + MMSE, data = train, family = binomial)
summary(modello6)
#predizione
pred <- predict(modello6, newdata = test, type = 'response')

pred <- as.data.frame(pred)



specificity.funct(pred,test)
#noi abbiamo bisogno di sensitività alta per riconoscere i pazienti malati
k=1
s<-seq(from=0.10, to=0.80, by=0.05)
for (j in s){
label.pred<-ifelse(pred>j, 1, 0)
label.pred <- factor(label.pred)
True.Positive=0
True.Negative=0
False.Positive=0
False.Negative=0
for (i in 1:100){
  if(test$label[i]==1 && label.pred[i]==test$label[i])
    True.Positive=True.Positive+1
  if(test$label[i]==0 && label.pred[i]==test$label[i])
    True.Negative=True.Negative+1
  if(test$label[i]==1 && label.pred[i]==0)
    False.Negative=False.Negative+1
  if(test$label[i]==0 && label.pred[i]==1)
    False.Positive=False.Positive+1
}

Precision[k] = True.Positive / (True.Positive + False.Positive) 
Recall[k] = True.Positive / (True.Positive + False.Negative)
F1.score[k]=2 * (Precision[k] * Recall[k]) / (Precision[k] + Recall[k])
k=k+1
}
F1.score
#modello sembra buono F1 è sempre tra il 70% e l'84% 

modello7<-glmrob(label ~ nWBV + Age + MMSE, data = train, family = binomial)
summary(modello7)
#age non particolarmente significativo ma accettabile

pred2 <- predict(modello7, newdata = test, type = 'response')

pred2 <- as.data.frame(pred2)

k=1
s<-seq(from=0.10, to=0.80, by=0.05)
F1.score.2=array(data=NA,dim=15)
for (j in s){
  label.pred<-ifelse(pred2>j, 1, 0)
  label.pred <- factor(label.pred)
  True.Positive=0
  True.Negative=0
  False.Positive=0
  False.Negative=0
  for (i in 1:100){
    if(test$label[i]==1 && label.pred[i]==test$label[i])
      True.Positive=True.Positive+1
    if(test$label[i]==0 && label.pred[i]==test$label[i])
      True.Negative=True.Negative+1
    if(test$label[i]==1 && label.pred[i]==0)
      False.Negative=False.Negative+1
    if(test$label[i]==0 && label.pred[i]==1)
      False.Positive=False.Positive+1
  }
  
  Precision[k] = True.Positive / (True.Positive + False.Positive) 
  Recall[k] = True.Positive / (True.Positive + False.Negative)
  F1.score.2[k]=2 * (Precision[k] * Recall[k]) / (Precision[k] + Recall[k])
  k=k+1
}
F1.score.2
#anche in questo modello F1 è sempre tra 70% e 83%

#uso dataset train con solo MMSE e MMSE+ splines ( fattore più significativo)

modello8<-glmrob(label ~ MMSE, data = train, family = binomial)
summary(modello8)
#intercetta e MMSE molto significativi, guardo predizione
pred8 <- predict(modello8, newdata = test, type = 'response')

pred8 <- as.data.frame(pred8)

k=1
s<-seq(from=0.10, to=0.80, by=0.05)
F1.score.2=array(data=NA,dim=15)
for (j in s){
  label.pred<-ifelse(pred8>j, 1, 0)
  label.pred <- factor(label.pred)
  True.Positive=0
  True.Negative=0
  False.Positive=0
  False.Negative=0
  for (i in 1:100){
    if(test$label[i]==1 && label.pred[i]==test$label[i])
      True.Positive=True.Positive+1
    if(test$label[i]==0 && label.pred[i]==test$label[i])
      True.Negative=True.Negative+1
    if(test$label[i]==1 && label.pred[i]==0)
      False.Negative=False.Negative+1
    if(test$label[i]==0 && label.pred[i]==1)
      False.Positive=False.Positive+1
  }
  
  Precision[k] = True.Positive / (True.Positive + False.Positive) 
  Recall[k] = True.Positive / (True.Positive + False.Negative)
  F1.score.2[k]=2 * (Precision[k] * Recall[k]) / (Precision[k] + Recall[k])
  k=k+1
}
F1.score.2
#percentuali sempre uguali provo con curva roc



data.roc <- roc.curve(pred8, test)

x11()
plot.roc(data.roc = data.roc)
roc_auc_vec(truth = test$label, estimate = pred8$pred, event_level = 'second') 

###MMSE con splines

modello9<-glmrob(label ~ bs(MMSE), data = train, family = binomial)
summary(modello9)
#intercetta e MMSE perdono significatività con le splines, guardo predizione
pred9 <- predict(modello9, newdata = test, type = 'response')

pred9 <- as.data.frame(pred9)

k=1
s<-seq(from=0.10, to=0.80, by=0.05)
F1.score.2=array(data=NA,dim=15)
for (j in s){
  label.pred<-ifelse(pred8>j, 1, 0)
  label.pred <- factor(label.pred)
  True.Positive=0
  True.Negative=0
  False.Positive=0
  False.Negative=0
  for (i in 1:100){
    if(test$label[i]==1 && label.pred[i]==test$label[i])
      True.Positive=True.Positive+1
    if(test$label[i]==0 && label.pred[i]==test$label[i])
      True.Negative=True.Negative+1
    if(test$label[i]==1 && label.pred[i]==0)
      False.Negative=False.Negative+1
    if(test$label[i]==0 && label.pred[i]==1)
      False.Positive=False.Positive+1
  }
  
  Precision[k] = True.Positive / (True.Positive + False.Positive) 
  Recall[k] = True.Positive / (True.Positive + False.Negative)
  F1.score.2[k]=2 * (Precision[k] * Recall[k]) / (Precision[k] + Recall[k])
  k=k+1
}
F1.score.2
#F1.score leggermente peggiore del precedente, guardo curva roc
data.roc9 <- roc.curve(pred9, test)

x11()
plot.roc(data.roc = data.roc9)
roc_auc_vec(truth = test$label, estimate = pred9$pred, event_level = 'second')

#riparto da modello completo visto che ho cambiato dataset

modello10<-glmrob(label ~ M + EDUC + nWBV + Age + MMSE + eTIV , data = train, family = binomial)
summary(modello10)


predizione10 <- predict(modello10, newdata = test, type = 'response')
predizione10 <- as.data.frame(predizione10)
data.roc.10<- roc.curve(predizione10, test)
x11()
plot.roc(data.roc = data.roc.10, pred.prob = predizione10$predizione10)
roc_auc_vec(truth = test$label, estimate = predizione10$predizione10, event_level = 'second') 

data.roc.10$F1.score[375]
data.roc.10$accuracy[375]
data.roc.10$precision[375]
#EDUC  è il fattore meno significativo

modello11<-glmrob(label ~ M + nWBV + Age + MMSE + eTIV , data = train, family = binomial)
summary(modello11)

#tolgo anche eTIV

modello12<-glmrob(label ~ M + nWBV + Age + MMSE , data = train, family = binomial)
summary(modello12)


modello13<-glmrob(label ~ M + nWBV + MMSE , data = train, family = binomial)
summary(modello13)

modello14<-glmrob(label ~ M + MMSE , data = train, family = binomial)
summary(modello14)


#modello super semplificato: predizione con solo MMSE
modello15<-glmrob(label ~ MMSE , data = train, family = binomial)
summary(modello15)

pred3 <- predict(modello15, newdata = test, type = 'response')

pred3 <- as.data.frame(pred3)

k=1
s<-seq(from=0.10, to=0.80, by=0.05)
F1.score.3=array(data=NA,dim=15)
for (j in s){
  label.pred<-ifelse(pred2>j, 1, 0)
  label.pred <- factor(label.pred)
  True.Positive=0
  True.Negative=0
  False.Positive=0
  False.Negative=0
  for (i in 1:100){
    if(test$label[i]==1 && label.pred[i]==test$label[i])
      True.Positive=True.Positive+1
    if(test$label[i]==0 && label.pred[i]==test$label[i])
      True.Negative=True.Negative+1
    if(test$label[i]==1 && label.pred[i]==0)
      False.Negative=False.Negative+1
    if(test$label[i]==0 && label.pred[i]==1)
      False.Positive=False.Positive+1
  }
  
  Precision[k] = True.Positive / (True.Positive + False.Positive) 
  Recall[k] = True.Positive / (True.Positive + False.Negative)
  F1.score.3[k]=2 * (Precision[k] * Recall[k]) / (Precision[k] + Recall[k])
  k=k+1
}
F1.score.3


modello17<-glmrob(label ~ M + bs(MMSE)  , data = train, family = binomial)
summary(modello17)

pred17 <- predict(modello17, newdata = test, type = 'response')

pred17 <- as.data.frame(pred17)

dataset<-read.table("oasis_longitudinal.csv",header=T,sep=",")
dataset
casestudy=dataset[which(dataset$Group!="Converted"),]
i1=which(casestudy$Group=="Demented")
New=vector(mode="logical",length=336)
New[i1]=1
New[-i1]=0
n=cbind(casestudy,New)

fit_glm <- glm( New  ~ Age + nWBV + ASF, data = casestudy,family=binomial )
plot(n)
abline(fit_glm, col="red", lwd=2) #non la plotta
plot(fit_glm)

#228,229,271 sembrano outliers cioè paziente n  OAS2_0108 entrambe le visite e paziente  n 		
# OAS2_0129 terza visita
View(n)

pairs(data_philips)

fit_MCD <- covMcd(x = n, alpha = .8)
fit_MCD

plot(fit_MCD, classic=TRUE, labels.id=FALSE, which="distance")
plot(fit_MCD,labels.id=FALSE,which=c("dd"))

fit_lms <- lmsreg(New  ~ Age + nWBV + ASF, data= n)
fit_lts <- ltsReg(New  ~ Age + nWBV + ASF ,alpha=.75,mcd=TRUE,data=n)

x11()
plot(n[,1:3])
abline(fit_glm, col="red", lwd=2)
abline(fit_lms, col="darkblue", lwd=2)
abline(fit_lts, col="darkgreen", lwd=2)
legend("bottomleft", c('OLS', 'LMS', 'LTS'), lwd=rep(2,4), col=c("red", "darkblue", "darkgreen"))
#non plotta le linee

plot(fit_lts)

####aggiungo splines
modello10_sp<-glmrob(label ~ M + bs(EDUC) + bs(nWBV) + bs(Age) + bs(MMSE) + bs(eTIV) , data = train, family = binomial)
summary(modello10_sp)


predizione10_sp <- predict(modello10_sp, newdata = test, type = 'response')
predizione10_sp <- as.data.frame(predizione10_sp)
data.roc.10_sp<- roc.curve(predizione10_sp, test)
x11()
plot.roc(data.roc = data.roc.10_sp, pred.prob = predizione10_sp$predizione10_sp)
roc_auc_vec(truth = test$label, estimate = predizione10_sp$predizione10_sp, event_level = 'second') 

data.roc.10_sp$F1.score[375]
data.roc.10_sp$accuracy[375]
data.roc.10_sp$precision[375]

#######
modello<-glmrob(label ~ M + EDUC + nWBV + Age + MMSE + eTIV , data = data, family = binomial)
summary(modello)

modello1<-glmrob(label ~ M + EDUC + nWBV + Age + MMSE, data = data, family = binomial)
summary(modello1)

modello2<-glmrob(label ~ M + nWBV + Age + MMSE, data = data, family = binomial)
summary(modello2)

modello3<-glmrob(label ~ M + nWBV + MMSE, data = data, family = binomial)
summary(modello3)

modello4<-glmrob(label ~ nWBV + Age + MMSE, data = data, family = binomial)
summary(modello4)

modello5<-glmrob(label ~ nWBV + MMSE, data = data, family = binomial)
summary(modello5)


#Call:  glmrob(formula = label ~ nWBV + MMSE, family = binomial, data = data) 


#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  32.4119     4.2993   7.539 4.74e-14 ***
#  nWBV        -12.7408     4.1353  -3.081  0.00206 ** 
#  MMSE         -0.8341     0.1111  -7.506 6.12e-14 ***
  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#Robustness weights w.r * w.x: 
 # 329 weights are ~= 1. The remaining 42 ones are summarized as
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2529  0.5553  0.6781  0.6894  0.8716  0.9965 

#Number of observations: 371 
#Fitted by method 'Mqle'  (in 5 iterations)


#w.r sono i pesi di robustezza di ogni osservazione (i.e residual*w.r è uguale ai residui di Preason della psi function)
robustness_weights=modello5$w.r
#ci sono 329 pesi che sono circa 1


#divido i dati in train e test e fitto il modello solo su train, poi faccio prediction con i tes data
modello6<-glmrob(label ~ nWBV + MMSE, data = train, family = binomial)
summary(modello6)
#predizione
pred <- predict(modello6, newdata = test, type = 'response')

pred <- as.data.frame(pred)


specificity.funct <- function(predicted, test.set){
  p0 <- seq(0,1,by=0.001)
  spec <- NULL
  sens <- NULL
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
    spec <- c(spec, specificity)
    sens <- c(sens, sensitivity)
  }
  data.frame(cbind(spec,sens))
}
specificity.funct(pred,test)
#noi abbiamo bisogno di sensitività alta per riconoscere i pazienti malati
k=1
s<-seq(from=0.10, to=0.80, by=0.05)
for (j in s){
label.pred<-ifelse(pred>j, 1, 0)
label.pred <- factor(label.pred)
True.Positive=0
True.Negative=0
False.Positive=0
False.Negative=0
for (i in 1:100){
  if(test$label[i]==1 && label.pred[i]==test$label[i])
    True.Positive=True.Positive+1
  if(test$label[i]==0 && label.pred[i]==test$label[i])
    True.Negative=True.Negative+1
  if(test$label[i]==1 && label.pred[i]==0)
    False.Negative=False.Negative+1
  if(test$label[i]==0 && label.pred[i]==1)
    False.Positive=False.Positive+1
}

Precision[k] = True.Positive / (True.Positive + False.Positive) 
Recall[k] = True.Positive / (True.Positive + False.Negative)
F1.score[k]=2 * (Precision[k] * Recall[k]) / (Precision[k] + Recall[k])
k=k+1
}
F1.score
#modello sembra buono F1 è sempre tra il 70% e l'84% 

modello7<-glmrob(label ~ nWBV + Age + MMSE, data = train, family = binomial)
summary(modello7)
#age non particolarmente significativo ma accettabile

pred2 <- predict(modello7, newdata = test, type = 'response')

pred2 <- as.data.frame(pred2)

k=1
s<-seq(from=0.10, to=0.80, by=0.05)
F1.score.2=array(data=NA,dim=15)
for (j in s){
  label.pred<-ifelse(pred2>j, 1, 0)
  label.pred <- factor(label.pred)
  True.Positive=0
  True.Negative=0
  False.Positive=0
  False.Negative=0
  for (i in 1:100){
    if(test$label[i]==1 && label.pred[i]==test$label[i])
      True.Positive=True.Positive+1
    if(test$label[i]==0 && label.pred[i]==test$label[i])
      True.Negative=True.Negative+1
    if(test$label[i]==1 && label.pred[i]==0)
      False.Negative=False.Negative+1
    if(test$label[i]==0 && label.pred[i]==1)
      False.Positive=False.Positive+1
  }
  
  Precision[k] = True.Positive / (True.Positive + False.Positive) 
  Recall[k] = True.Positive / (True.Positive + False.Negative)
  F1.score.2[k]=2 * (Precision[k] * Recall[k]) / (Precision[k] + Recall[k])
  k=k+1
}
F1.score.2
#anche in questo modello F1 è sempre tra 70% e 83%

#riparto da modello completo visto che ho cambiato dataset

modello10<-glmrob(label ~ M + EDUC + nWBV + Age + MMSE + eTIV , data = train, family = binomial)
summary(modello10)

#EDUC  è il fattore meno significativo

modello11<-glmrob(label ~ M + nWBV + Age + MMSE + eTIV , data = train, family = binomial)
summary(modello11)

#tolgo anche eTIV

modello12<-glmrob(label ~ M + nWBV + Age + MMSE , data = train, family = binomial)
summary(modello12)


modello13<-glmrob(label ~ M + nWBV + MMSE , data = train, family = binomial)
summary(modello13)

modello14<-glmrob(label ~ M + MMSE , data = train, family = binomial)
summary(modello14)
pred14 <- predict(modello14, newdata = test, type = 'response')
pred14 <- as.data.frame(pred14)
data.roc.14<- roc.curve(pred14, test)
x11()
plot.roc(data.roc = data.roc.14, pred.prob = pred14$pred14)
roc_auc_vec(truth = test$label, estimate = pred14$pred14, event_level = 'second') 
data.roc.14$F1.score[375]
data.roc.14$accuracy[375]
data.roc.14$precision[375]
sensitivity=data.roc.14$y.roc[375]
specificity=1-data.roc.14$x.roc[375]
sensitivity
specificity


#modello super semplificato: predizione con solo MMSE
modello15<-glmrob(label ~ MMSE , data = train, family = binomial)
summary(modello15)

pred15 <- predict(modello15, newdata = test, type = 'response')
pred15 <- as.data.frame(pred15)
data.roc.15<- roc.curve(pred15, test)
x11()
plot.roc(data.roc = data.roc.15, pred.prob = pred15$pred15)
roc_auc_vec(truth = test$label, estimate = pred15$pred15, event_level = 'second') 

data.roc.15$F1.score[375]
data.roc.15$accuracy[375]
data.roc.15$precision[375]


k=1
s<-seq(from=0.10, to=0.80, by=0.05)
F1.score.3=array(data=NA,dim=15)
for (j in s){
  label.pred<-ifelse(pred2>j, 1, 0)
  label.pred <- factor(label.pred)
  True.Positive=0
  True.Negative=0
  False.Positive=0
  False.Negative=0
  for (i in 1:100){
    if(test$label[i]==1 && label.pred[i]==test$label[i])
      True.Positive=True.Positive+1
    if(test$label[i]==0 && label.pred[i]==test$label[i])
      True.Negative=True.Negative+1
    if(test$label[i]==1 && label.pred[i]==0)
      False.Negative=False.Negative+1
    if(test$label[i]==0 && label.pred[i]==1)
      False.Positive=False.Positive+1
  }
  
  Precision[k] = True.Positive / (True.Positive + False.Positive) 
  Recall[k] = True.Positive / (True.Positive + False.Negative)
  F1.score.3[k]=2 * (Precision[k] * Recall[k]) / (Precision[k] + Recall[k])
  k=k+1
}
F1.score.3

######non considero MSSE 
mod11<-glmrob(label ~ M + bs(EDUC) + bs(nWBV) + bs(Age) + bs(eTIV) , data = train, family = binomial)
summary(mod11)

pred11 <- predict(mod11, newdata = test, type = 'response')

pred11 <- as.data.frame(pred11)


###tolgo age
mod12<-glmrob(label ~ M + bs(EDUC) + bs(nWBV) + bs(eTIV) , data = train, family = binomial)
summary(mod12)

pred12 <- predict(mod12, newdata = test, type = 'response')

pred12 <- as.data.frame(pred12)

###tolgo nWBV
mod13<-glmrob(label ~ M + bs(EDUC) + bs(eTIV) , data = train, family = binomial)
summary(mod13)

pred13 <- predict(mod13, newdata = test, type = 'response')

pred13 <- as.data.frame(pred13)

######tolgo eTIV
mod14<-glmrob(label ~ M + bs(EDUC) , data = train, family = binomial)
summary(mod14)

pred14 <- predict(mod14, newdata = test, type = 'response')
pred14 <- as.data.frame(pred14)
data.roc.14<- roc.curve(pred14, test)
x11()
plot.roc(data.roc = data.roc.14, pred.prob = pred14$pred14)
roc_auc_vec(truth = test$label, estimate = pred14$pred14, event_level = 'second') 

data.roc.14$F1.score[375]
data.roc.14$accuracy[375]
data.roc.14$precision[375]


#copio stessi modelli della regressione logistica NON robusta

model_rob=glmrob(label ~ M + bs(EDUC) +bs(nWBV) + bs(Age)  + bs(MMSE) + bs(eTIV) , data = train, family = binomial)
summary(model_rob) # M significativo, rTIV significatico, EDUC leggermente significativo, il resto no

pred_rob <- predict(model_rob, newdata = test, type = 'response')

pred_rob <- as.data.frame(pred_rob)

model_rob2=glmrob(label ~ M + bs(Age) + bs(nWBV) + bs(MMSE) + bs(eTIV), data = train, family = binomial)
summary(model_rob2)#senza educ

#sembra che etiv possa essere modellizzata con la dummy M.F
model_rob3=glmrob(label ~ M + bs(Age) + bs(nWBV) + bs(MMSE), data = train, family = binomial)
summary(model_rob3)#senza educ,  etiv




#TOLGO VARIABILE MMSE
model_red <- glmrob(label ~ M + bs(Age) + bs(nWBV) + bs(eTIV), data = train, family = binomial)
summary(model_red)

pred.red <- predict(model_red, newdata = test, type = 'response')
pred.red <- as.data.frame(pred.red)
data.roc.red <- roc.curve(pred.red, test)
x11()
plot.roc(data.roc = data.roc.red, pred.prob = pred.red$pred.red)
roc_auc_vec(truth = test$label, estimate = pred.red$pred, event_level = 'second') 

data.roc.red$F1.score[426]
data.roc.red$accuracy[426]
data.roc.red$precision[426]
sensitivity=data.roc.red$y.roc[426]
specificity=1-data.roc.red$x.roc[426]
sensitivity
specificity


model_red2 <- glmrob(label ~ M + bs(Age) + bs(nWBV) + eTIV, data = train, family = binomial)
summary(model_red2)# eTIV così aumenta significatività

pred.red2 <- predict(model_red2, newdata = test, type = 'response')
pred.red2 <- as.data.frame(pred.red2)
data.roc.red <- roc.curve(pred.red2, test)
x11()
plot.roc(data.roc = data.roc.red, pred.prob = pred.red$pred.red)
roc_auc_vec(truth = test$label, estimate = pred.red$pred, event_level = 'second')

model_red3 <- glmrob(label ~ M + Age + bs(nWBV) + eTIV, data = train, family = binomial)
summary(model_red3)# Age aumenta significatività e anche NWBV al terzo ordine

pred.red3 <- predict(model_red3, newdata = test, type = 'response')
pred.red3 <- as.data.frame(pred.red3)
data.roc.red <- roc.curve(pred.red3, test)
x11()
plot.roc(data.roc = data.roc.red, pred.prob = pred.red$pred.red)
roc_auc_vec(truth = test$label, estimate = pred.red$pred, event_level = 'second') #non mi funziona il plot

model_red4 <- glmrob(label ~ M + Age + nWBV + eTIV, data = train, family = binomial)
summary(model_red4) #aumenta nWBV diminuiscono Age and eTIV

pred.red4 <- predict(model_red4, newdata = test, type = 'response')
pred.red4 <- as.data.frame(pred.red4)
data.roc.red <- roc.curve(pred.red4, test)
x11()
plot.roc(data.roc = data.roc.red, pred.prob = pred.red$pred.red)
roc_auc_vec(truth = test$label, estimate = pred.red$pred, event_level = 'second') #non mi funziona il plot

model_red5 <- glmrob(label ~ M + Age +  nWBV , data = train, family = binomial)
summary(model_red5) #M diminuisce e il resto non migliora, prossimo modello tolgo Age al posto di togliere eTIV

pred.red5 <- predict(model_red5, newdata = test, type = 'response')
pred.red5 <- as.data.frame(pred.red5)
data.roc.red <- roc.curve(pred.red5, test)
x11()
plot.roc(data.roc = data.roc.red, pred.prob = pred.red5$pred.red5)
roc_auc_vec(truth = test$label, estimate = pred.red5$pred5, event_level = 'second') 

model_red6 <- glmrob(label ~ M + eTIV+  nWBV , data = train, family = binomial)
summary(model_red6)

pred.red6 <- predict(model_red6, newdata = test, type = 'response')
pred.red6 <- as.data.frame(pred.red6)
data.roc.red6 <- roc.curve(pred.red6, test)
x11()
plot.roc(data.roc = data.roc.red6, pred.prob = pred.red6$pred.red6)
roc_auc_vec(truth = test$label, estimate = pred.red6$pred6, event_level = 'second')


data.roc.red6$F1.score[376]
data.roc.red6$accuracy[376]
data.roc.red6$precision[376]
sensitivity=data.roc.red6$y.roc[376]
specificity=1-data.roc.red6$x.roc[376]
sensitivity
specificity

model_red7 <- glmrob(label ~ M +  nWBV , data = train, family = binomial)
summary(model_red7) 

pred.red7 <- predict(model_red7, newdata = test, type = 'response')
pred.red7 <- as.data.frame(pred.red7)
data.roc.red7 <- roc.curve(pred.red7, test)
x11()
plot.roc(data.roc = data.roc.red7, pred.prob = pred.red7$pred.red7)
roc_auc_vec(truth = test$label, estimate = pred.red7$pred7, event_level = 'second') 

data.roc.red7$F1.score[325]
data.roc.red7$accuracy[325]
data.roc.red7$precision[325]
sensitivity=data.roc.red7$y.roc[325]
specificity=1-data.roc.red7$x.roc[325]
sensitivity
specificity

model_red8 <- glmrob(label ~ nWBV , data = train, family = binomial)
summary(model_red8) 

pred.red8 <- predict(model_red8, newdata = test, type = 'response')
pred.red8 <- as.data.frame(pred.red8)
data.roc.red8 <- roc.curve(pred.red8, test)
x11()
plot.roc(data.roc = data.roc.red8, pred.prob = pred.red8$pred.red8)
roc_auc_vec(truth = test$label, estimate = pred.red8$pred8, event_level = 'second') 


data.roc.red8$F1.score[376]
data.roc.red8$accuracy[376]
data.roc.red8$precision[376]
sensitivity=data.roc.red8$y.roc[376]
specificity=1-data.roc.red8$x.roc[376]
sensitivity
specificity


