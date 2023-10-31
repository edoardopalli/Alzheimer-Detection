setwd("C:/Users/User/Downloads")
library(robustbase)
library(psych)
library(MASS)
library(ellipse)
library(here)
library(DescTools)
library(knitr)
library(RobStatTM)
library(mgcv)
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

data <- read.csv('data.csv', header = TRUE, sep = ',')
data$label <- ifelse(data$label=='Dem', 1, 0)
data$label <- factor(data$label)
train <- data[1:271,]
test <- data[272:371,]
data <- data[,-1]
train <- train[,-1]
test <- test[,-1]


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