setwd("C:/Users/User/Downloads")
library("dplyr")
alz<-read.table("oasis_cross-sectional.csv",header=T,sep=",")
head(alz)
alz2<-read.table("oasis_longitudinal.csv",header=T,sep=",")
head(alz2)
ND2<-alz2[which(alz2$Group=="Nondemented"),]
ND2
D2<-alz2[which(alz2$Group=="Demented"),]
D2
C2<-alz2[which(alz2$Group=="Converted"),]
C2

load("C:/Users/User/OneDrive/Desktop/ERICA/POLI/APPLIED STATISTICS/LAB_5/mcshapiro.test.RData")

NDrid<-ND2[,4:15]
NDrid<-NDrid[,-3]
NDrid<-NDrid[,-3]
NDrid
#controllo gaussianità
mcshapiro.test ( NDrid )
mcshapiro.test(NDrid)$pvalue
#pvalue=0
shapiro.test(NDrid$Age)
shapiro.test(NDrid$EDUC)
shapiro.test(NDrid$SES)
shapiro.test(NDrid$MMSE)
shapiro.test(NDrid$eTIV)
shapiro.test(NDrid$nWBV) #ok normalità per questo, per il resto no
shapiro.test(NDrid$CDR)
shapiro.test(NDrid$ASF)

n = nrow(NDrid) 
p = ncol(NDrid)
sample.mean = sapply ( NDrid, mean )
sample.var = var ( NDrid )
alpha <- 0.01
cfr.fisher = ( (n - 1) * p / (n - p) ) * qf ( 1 - alpha, p, n - p )


casestudy=alz2[which(alz2$Group!="Converted"),]
i1=which(casestudy$Group=="Demented")
New=vector(mode="logical",length=336)
New[i1]=1
New[-i1]=0
n=cbind(casestudy,New)
head(n)
fit = glm( New  ~ Age + EDUC  + eTIV + nWBV + CDR+ ASF, data = n,family=binomial )
summary ( fit )
fit$coefficients

fit1 = glm( New  ~ Age + EDUC  + nWBV + CDR+ ASF, data = casestudy,family=binomial )
summary ( fit1 )
fit1$coefficients

fit2 = glm( New  ~ Age + nWBV + CDR+ ASF, data = casestudy,family=binomial )
summary ( fit2 )
fit2$coefficients

fit3 = glm( New  ~ Age + nWBV , data = casestudy,family=binomial )
summary ( fit3 )
fit3$coefficients

fit4 = glm( New  ~ Age + nWBV + ASF, data = casestudy,family=binomial )
summary ( fit4 )
predict.glm(fit4,newdata=C2,type="response")
fit4$coefficients

###forse il 4 è il modello migliore

fit5 = glm( New  ~ Age + nWBV + ASF + SES, data = casestudy,family=binomial )
summary ( fit5 )
fit5$coefficients

predict.glm(fit5,newdata=C2,type="response")
#modello sembra ok ma non mi convince la prediction

fit6 = glm( New  ~ Age + nWBV  + MMSE, data = casestudy,family=binomial )
summary ( fit6 )
fit6$coefficients
predict.glm(fit6,newdata=C2,type="response")

fit7 = glm( New  ~ Age + nWBV + ASF + MMSE , data = casestudy,family=binomial )
summary ( fit7 )
fit7$coefficients
#prob 0 or 1 occurred
sex<-as.factor(casestudy$M.F)

fit8 = glm( New  ~ M.F + Age + nWBV + eTIV +MMSE, data = casestudy,family=binomial )
summary ( fit8 )
fit8$coefficients

pred<-predict.glm(fit8,newdata=C2,type="response")
#modello sembra buono
conv<-2*C2$CDR
l<-length(pred)
precison<-sum(round(pred)==conv)/l



fit9 = glm( New  ~ Age + nWBV  + SES, data = casestudy,family=binomial )
summary ( fit9 )

predict.glm(fit9,newdata=C2,type="response")

library(rgl)
library(ISLR2)
library(car)
library(np)
library(splines)
library(fda)

m_loc = npreg(Group  ~ Age + EDUC + SES + MMSE  + eTIV + nWBV + CDR+ ASF,
              ckertype = 'uniform',
              bws = 3200, # bandwidth
              data = nd2)
typeof(ND2)
typeof(alz2)

#####svarioni a caso

income_newdata=data.frame(income=with(Prestige, seq(range(income)[1],range(income)[2],by=0.5)))
preds=predict(m_loc,newdata=income_newdata,se=T)
se.bands=cbind(preds$fit +2* preds$se.fit ,preds$fit -2* preds$se.fit)
with(
  Prestige,
  plot(
    income ,
    prestige ,
    xlim = range(income_newdata$income) ,
    cex = .5,
    col = " darkgrey ",
    main = 'Local Averaging - bws3200 - Uniform kernel'
  )
)
lines(income_newdata$income,preds$fit ,lwd =2, col =" blue")
matlines(income_newdata$income,se.bands ,lwd =1, col =" blue",lty =3)





with(ND2, scatterplotMatrix(data.frame(Visit, MR.Delay, Age ,EDUC, SES, MMSE, CDR, eTIV,nWBV, ASF)))
model_lm=lm(~ education + income, data=alz2)
summary(model_lm)

#####plot 
library(ISLR2)
library(car)
library(mgcv)
library(rgl)
library(splines)
library(pbapply)
with(alz2, scatterplotMatrix(data.frame(Age, eTIV, ASF,nWBV,MMSE,SES)))

pazienti=vector(mode="integer",length=373)
pazienti[1]=1
j=0
for(k in 1:372){
  if (alz2[k,1]==alz2[k+1,1]){
     pazienti[k+1]=j
     }
  else
    j=j+1
  pazienti[k]=j+1
}
pazienti


matplot(alz2$Subject.ID,alz2$MMSE,)




########prova con test set 
 
 cs<-n[-1,]
 cs<-cs[-1,]
 cs<-cs[-4,]
 cs<-cs[-4,]
 cs<-cs[-4,]
 cs<-cs[-4,]
 cs<-cs[-4,]
 cs<-cs[-7,]
 cs<-cs[-7,]
 dt<-ND2[1:9,]

 fit_t = glm( New  ~ Age + nWBV  + MMSE, data = cs,family=binomial )
 summary ( fit_t )
 fit_t$coefficients
 predict.glm(fit_t,newdata=dt,type="response")
 
 fit_t1 = glm( New  ~ Age + nWBV  + MMSE + ASF, data = cs,family=binomial )
 summary ( fit_t1 )
 fit_t1$coefficients
 predict.glm(fit_t1,newdata=dt,type="response")
 
 