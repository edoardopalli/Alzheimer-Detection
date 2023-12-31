### PCA ALZEHEIMER

library(dplyr)
setwd('C:/Users/franc/Desktop/NONPA/PROGETTO/alzheimer')
dataset_completo<- read.csv("oasis_longitudinal.csv")

visit1 <- filter(dataset_completo, Visit==1)
library(rgl)
library(mgcv)
library(rgl)
library(splines)
library(pbapply)



DATA.label <- visit1[,3]
DATA <- visit1[,c(8,9,11,13,14,15)] 
n <- dim(DATA)[1] 
p <- dim(DATA)[2] 


x11()
pairs(DATA, pch=19)
boxplot(DATA)
matplot(t(DATA), type='l', xlab='DATA', ylab= '', lty=1, col=rainbow(33), las=1)
M <- sapply(DATA,mean)
M
S <- cov(DATA)
S
R <- cor(DATA)
R
#CHECK VARIABILITIES TO KNOW IF STANDARDIZE OR NOT  
x11()
boxplot(DATA,col=3)
data<-DATA[,-2]
data<-data[,-5]
DATA.sd <- scale(data)
DATA.sd <- data.frame(DATA.sd)

head(DATA.sd)
sapply(DATA,mean)
sapply(DATA.sd,mean)
sapply(DATA,sd) 
sapply(DATA.sd,sd) #now sd are all =1
cov(DATA.sd) 
x11()
pairs(DATA[,-(1:2)])
pairs(DATA.sd)


#PCA: --------------------------------------------------------------------
pc.DATA <- princomp(DATA.sd, scores=T)
pc.DATA
summary(pc.DATA)
load<-pc.DATA$loadings
x11()
par(mar=c(1,4,0,2),mfrow=c(ncol(DATA.sd),1))
for(i in 1:ncol(DATA.sd))barplot(load[,i],ylim=c(-1,1))
pc<-data*load[,1:3]
pc<-matrix(data=0,nrow=nrow(data),ncol=3)
for(j in 1:nrow(data)){
  x=crossprod(as.numeric(data[j,]),load[,1:3])
  pc[j,]=t(x)
  
}

casestudy=visit1[which(visit1$Group!="Converted"),]
i1=which(casestudy$Group=="Demented")

New=vector(mode="logical",length=136)
New[i1]=1
New[-i1]=0


i2=which(alz2$Group=="Converted")
newdataset=pc[-i2,]
newdataset <- data.frame(newdataset)

model_gam=gam(New ~ s(X1,bs='cr') + s(X2,bs='cr')+ s(X3,bs='cr'), data = newdataset,family= binomial)
summary(model_gam)

hist(model_gam$residuals)
qqnorm(model_gam$residuals)

pred_gam=predict(model_gam,newdata=grid)
shapiro.test(model_gam$residuals)


model_gam_ns <-
  lm(New ~ ns(X1, df = 3) + ns(X2, df = 3) + ns(X3,df=3), data = newdataset)

plot(model_gam_ns$residuals,model_gam$residuals)
cor(model_gam_ns$residuals,model_gam$residuals)

model_gam_inter = gam(New~ s(X1, bs = 'cr') + 
                        s(X2, bs ='cr') + s(X3,bs='cr') +
                        s(I(X1 * X2), bs='cr')+s(I(X1 * X3), bs = 'cr')+s(I(X3 * X2), bs = 'cr'),
                      data = newdataset)
summary(model_gam_inter)

model_gam_inter_2 = gam(New~ s(X1, bs = 'cr') + 
                        s(X2, bs ='cr') + s(X3,bs='cr') +
                        s(I(X1 * X2), bs='cr')+s(I(X1 * X3), bs = 'cr'),
                      data = newdataset)
summary(model_gam_inter_2)

hist(model_gam_inter_2$residuals)
qqnorm(model_gam_inter_2$residuals)

# EXPLAINED VARIANCE: -----------------------------------------------------

y_range=c(0,  max(pc.DATA$sdev^2))
y_range2=c(0,  max(sapply(DATA.sd,sd)^2))
x11()
layout(matrix(c(2,3,1,3),2,byrow=T))
barplot(pc.DATA$sdev^2, las=2, main='Principal Components', ylim=y_range, ylab='Variances')
barplot(sapply(DATA.sd,sd)^2, las=2, main='Original Variables', ylim=y_range2, ylab='Variances')
#screeplot:
plot(cumsum(pc.DATA$sdev^2)/sum(pc.DATA$sde^2), type='b', axes=F, xlab='number of components',
     ylab='contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(DATA.sd),labels=1:ncol(DATA.sd),las=2)



# SCORES ------------------------------------------------------------------
#i-th score: y_i=e_i'*(x-mu) where e_i is the i-th PC (vector in R^p) and x is the matrix of observations
scores.DATA <- pc.DATA$scores
scores.DATA


# LOADINGS ----------------------------------------------------------------

TOT_N_COMP <- dim(pc.DATA$loadings)[1]
  NEEDED_K <-3

load.DATA    <- pc.DATA$loadings
load.DATA

x11()#attenta a quanti loadings devi plottare
par(mar = c(1,4,0,2), mfrow = c(3,2)) #for ex tot_n= -> a=4 b=2
for(i in 1: TOT_N_COMP)barplot(load.DATA[,i], ylim = c(-1, 1))

# let's plot only the most significant loadings
x11()
par(mar = c(1,4,0,2), mfrow = c(3,1))
for(i in 1:NEEDED_K) barplot(ifelse(abs(load.DATA[,i]) < 0.3, 0, load.DATA[,i]) , ylim = c(-1, 1));abline(h=0)
for(i in 1:NEEDED_K) {print(paste('loading e_',i,'=(' ))
  print(load.DATA[,i])
  }


# PROJECTION: -------------------------------------------------------------
scores.DATA <- pc.DATA$scores
scores.DATA
# Projection on the space generated by the k-th principal component
x11(width=21, height=7)
par(mfrow=c(2,5))
matplot(t(DATA.sd), type='l', main = 'Data', ylim=range(DATA.sd))

meanF <- colMeans(DATA.sd)
matplot(meanF, type='l', main = '0 PC', lwd=2, ylim=range(DATA.sd))
for(i in 1:TOT_N_COMP)
{
  projection <- matrix(meanF, dim(DATA.sd)[[1]], dim(DATA.sd)[[2]], byrow=T) + scores.DATA[,i] %*% t(load.DATA[,i])
  matplot(t(projection), type='l', main = paste(i, 'PC'), ylim=range(DATA.sd))
  matplot(meanF, type='l', lwd=2, add=T)
}

# Projection on the space generated by the first k principal components (cumulative)
x11(width=21, height=7)
par(mfrow=c(2,5))
matplot(t(DATA.sd), type='l', main = 'Data', ylim=range(DATA.sd))
meanF <- colMeans(DATA.sd)
matplot(meanF, type='l', main = 'First 0 PCs', lwd=2, ylim=range(DATA.sd))
projection <- matrix(meanF, dim(DATA.sd)[[1]], dim(DATA.sd)[[2]], byrow=T)
for(i in 1:TOT_N_COMP)
{
  projection <- projection + scores.DATA[,i] %*% t(load.DATA[,i])
  matplot(t(projection), type='l', main = paste('First', i, 'PCs'), ylim=range(DATA.sd))
  matplot(meanF, type='l', lwd=2, add=T)
}

#after NN PCs we see no significant differences in the plot so NN PCs are sufficient to explain majority of the variability















############################################################################
# TWO WAYS ANOVA
Y <- DATA.sd
#Y <- visit1[,8]
group1 <- factor(visit1$Group)
group2 <- factor(visit1$M.F,c("M","F"))
tt <- paste(group1,group2) #AR FF COLUMNS NAME OF THE COLUMNS WITH THE GROUPS
#AR FF COLUMNS NAME OF THE COLUMNS WITH THE GROUPS
group_12 <- factor( tt)
#ExAd <- paste(Ex, Ad)

g <- length(levels(group1))
b <- length(levels(group2))
g
b

M             <- colMeans(Y)
Mgroup1     <- tapply(Y, group1,mean)
Mgroup2         <- tapply(Y, group2, mean)
Minteraction <- tapply(Y, group_12, mean)


plot(group_12, Y, col=rainbow(5)[2:5], ylim=c(0,24))
# the group2s look different from the plot


# Parametric test:-------
attach(DATA.sd)
Y<-cbind(Age,EDUC,MMSE,nWBV,ASF)
detach(DATA.sd)

summary.manova(manova(Y~ group1+ group2+ group1:group2 ,data=visit1))
#interaction is not significant: pval=0.3828 
# Without interaction
summary.manova(manova(Y ~ group1 + group2))
#BOTH SIGNIFICANT


# PERMUTATIONAL SETTING: -----


#TEST STATISTICS
fit <- manova(Y~ group1+ group2+ group1:group2 ,data=DATA.sd)
T0_group1_group2 <- summary.manova(fit, test="Wilks")$stats[3,3]

# permutational distribution:

manova.H0group1_group2 <- manova(Y ~ group1 + group2)
manova.H0group1_group2
residuals.H0group1_group2 <- manova.H0group1_group2$residuals
n <- 150

B=1000
seed = 26111992
set.seed(seed)
T_group1_group2 <- numeric(B)
for(perm in 1:B){
  permutation <- sample(n)
  residuals.H0group1_group2 <- residuals.H0group1_group2[permutation]
  Y.perm.H0group1_group2 <- manova.H0group1_group2$fitted + residuals.H0group1_group2
  T_group1_group2[perm] <- summary.manova(manova(Y.perm.H0group1_group2 ~ group1 + group2 + group1:group2))$stats[3,3]
}

#pvalue
sum(T_group1_group2 >= T0_group1_group2)/B # 0.32

#if pvalue high -> Not significant, reduce the model and then perform my test on my main effects.

#H_0:alpha=0 vs $H_1:alpha!=0, I am assuming under $H_0$ the following model $Y = \mu + \alpha_i + \epsilon$, while for $H_0:\alpha=0$ vs $H_1:\alpha\neq0$ I am assuming $Y = \mu + \beta_j + \epsilon$

T0_group1 <- 15.054
# residuals under H0:
# Y = mu + beta*group2
manova.H0group1 <- manova(Y ~ group2)
residuals.H0group1 <- manova.H0group1$residuals

# TEST OF FACTOR group2   (H0: beta=0)
T0_group2 <- 12.944      
# residuals under H0:
# Y = mu + alpha*group1
manova.H0group2 <- manova(Y ~ group1)
residuals.H0group2 <- manova.H0group2$residuals

# permutational distribution and P-value

# TEST OF FACTOR group1 ANF TEST OF FACTOR group2: 

# p-values
B <- 1000
T_group2 <- T_group1 <- numeric(B)
for(perm in 1:B){
  permutation <- sample(n)
  
  Y.perm.H0group1 <- manova.H0group1$fitted + residuals.H0group1[permutation]
  T_group1[perm] <- summary.manova(manova(Y.perm.H0group1 ~ group1 + group2))[[1]][1,4]
  
  Y.perm.H0group2 <- manova.H0group2$fitted + residuals.H0group2[permutation]
  T_group2[perm] <- summary.manova(manova(Y.perm.H0group2 ~ group1 + group2))[[1]][2,4]
}

#pvalue
sum(T_group1 >= T0_group1)/B #0 significant
sum(T_group2 >= T0_group2)/B #0 significant


