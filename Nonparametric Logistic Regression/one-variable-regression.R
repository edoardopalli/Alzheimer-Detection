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
library(np)
library(locfit)

#EDUC: years of education 
#SES: socio economic status
#MMSE: mini mental state examination (0:30 30=stato mentale sano. Anche persone con un iniziale deterioramento cognitivo, ma con un'alta scolarizzazione possono ottenere un punteggio pari a 29 e 30)
#CDR: clinical dementia rating (0:3 3=perdita di memoria grave)
#eTIV: extimated total intracranial volume 
#nWBV: Normalize Whole Brain Volume
#ASF: Atlas Scaling Factor (the volume-scaling factor required to match each individual to the atlas target)
#LABEL:Demented=1

data <- read.csv('data_unavar.csv', header = TRUE, sep = ',')
data <- data[,-1]
data$label <- ifelse(data$label=='Dem', 1, 0)
data$label <- factor(data$label)
train <- data[1:271,]
test <- data[272:371,]



freq.mmse <- NULL
for(i in 1:(range(train$MMSE)[2]-range(train$MMSE)[1]+1)){
  freq.mmse <- c(freq.mmse, length(which(train$MMSE==(range(train$MMSE)[1]+i-1)))/dim(train)[1])
}
hist(train$MMSE)

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
    geom_path(col='blue', lwd=1)+ #geom_path
    geom_line(data = data.frame(cbind(X1=seq(0,1,0.001), X2=seq(0,1,0.001))), aes(X1,X2), col='orangered', lty = 2, lwd=1.3)+
    xlab('1 - Specificity') + ylab('Sensitivity') + labs(title = 'ROC Curve', subtitle = paste('Area Under Curve is', round(roc_auc_vec(truth = test$label, estimate = pred.prob, event_level = 'second'), digits = 2))) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) 
  
}


#GLM WITH 1 COVARIATE
mdl <- glm(label ~ MMSE, data = train, family = binomial)
summary(mdl)
1 - mdl$deviance/mdl$null.deviance
mdl$aic

x <- data.frame(MMSE = seq(range(train$MMSE)[1], range(train$MMSE)[2], by=0.25))
p.hat.param <- predict(mdl, newdata = x, type = 'response') #1/(1+exp(-mdl$coefficients[1]-mdl$coefficients[2]*x))
x11() #qua come scrivere le frequenze per ciascun valore di MMSE
ggplot()+
  geom_line(aes(x$MMSE,p.hat.param), col='red')+
  geom_point(aes(train$MMSE, as.numeric(train$label==1)))
pred.glm <- predict(mdl, newdata = test, type = 'response')
pred.glm <- as.data.frame(pred.glm)
data.roc.glm <- roc.curve(pred.glm, test.set = test)
x11()
plot.roc(data.roc.glm, pred.glm$pred.glm)
roc_auc_vec(truth = test$label, estimate = rep(1,dim(test)[1]) - pred.glm$pred.glm)

# Global polynomial regression
mdl_poly=lapply(1:7,function(degree){glm(label ~ poly(MMSE,degree=degree), data = train, family = binomial)})
do.call(what = anova, args = mdl_poly)
anova(mdl_poly[[3]], mdl_poly[[4]], test = 'Chisq')
anova(mdl_poly[[4]], mdl_poly[[5]], test = 'Chisq')
anova(mdl_poly[[4]], mdl_poly[[6]], test = 'Chisq')
anova(mdl_poly[[4]], mdl_poly[[7]], test = 'Chisq')

aic <- NULL
for (i in 1:7) {
  aic <- c(aic, mdl_poly[[i]]$aic)
}

plot(1:7, aic, main = 'AIC', xlab= 'Degree of polynomial', ylab = ' ', col='blue', pch=16)
lines(1:7, aic, col='blue')
mdl_poly[[4]]$aic
1 - mdl_poly[[4]]$deviance/mdl_poly[[4]]$null.deviance

auc_poly <- NULL
for (i in 1:7) {
  p_ <- predict(mdl_poly[[i]], newdata = test, typse = 'response')
  auc_poly <- c(auc_poly,roc_auc_vec(truth = test$label, estimate = rep(1,dim(test)[1])-p_))
}

plot(1:7, auc_poly, main = 'Area under ROC curve', xlab= 'Degree of polynomial', ylab = ' ', col='blue', pch=16)
lines(1:7, auc_poly, col='blue')
#I choose model 4 
p.hat.poly <- predict(mdl_poly[[4]], newdata = x, type = 'response') #1/(1+exp(-mdl$coefficients[1]-mdl$coefficients[2]*x))
x11()
ggplot() +
  geom_line(aes(x$MMSE,p.hat.poly), col='Blue') +
  geom_point(aes(train$MMSE, as.numeric(train$label==1))) + 
  ylab('Probability to be Demented') +xlab('MMSE') + labs(title = 'Logistic Curve') + theme(plot.title = element_text(hjust = 0.5))

pred.poly <- predict(mdl_poly[[4]], newdata = test, type = 'response')
pred.poly <- as.data.frame(pred.poly)
data.roc.poly <- roc.curve(pred.poly, test.set = test)
x11()
plot.roc(data.roc.poly, pred.poly$pred.poly) + geom_point(aes(x.roc[i.opt],y.roc[i.opt]), col = 'purple', cex = 3)

i.opt <- which(data.roc.poly$y.roc>=0.75 & data.roc.poly$x.roc<=0.25)[322]
p0.opt <- seq(0,1, length.out = 1000)[322]
sens.poly <- data.roc.poly$y.roc[i.opt]
spec.poly <- 1 - data.roc.poly$x.roc[i.opt]
F1.poly <- data.roc.poly$F1.score[i.opt]
accuracy.poly <- data.roc.poly$accuracy[i.opt]
precision.poly <- data.roc.poly$precision[i.opt]

#Local Likelihood 
fit_locfit <- locfit(as.numeric(label==1) ~ locfit::lp(MMSE, deg = 1), data = train, family = "binomial", kern = "gauss")
summary(fit_locfit)

x11()
plot(fit_locfit, main = 'Logistic Curve', xlab = 'MMSE', ylab = 'Probability to be Demented', col='blue', ylim=c(0,1))
points(train$MMSE, as.numeric(train$label==1), pch = 16)
grid()


gcv.stat <- NULL
nn <- c(1/3,1/2,0.7,0.8,1,2,3,4,5,6)
for(i in 1:10){
  fit <- locfit(as.numeric(label==1) ~ locfit::lp(MMSE, deg = 1, nn=nn[i]), data = train,
                family = "binomial", kern = "gauss")
  gcv.stat <- c(gcv.stat, gcv(fit)[4])
}
nn.min <- which(gcv.stat==min(gcv.stat))

gcv.locfit <- locfit(as.numeric(label==1) ~ locfit::lp(MMSE, deg = 1, nn=nn[nn.min]), data = train,
                     family = "binomial", kern = "gauss")

pred.locfit <- predict(fit_locfit, newdata = test, type = 'response')
pred.locfit <- as.data.frame(pred.locfit)
data.roc.locfit <- roc.curve(pred.locfit, test)
x11()
plot.roc(data.roc.locfit, pred.locfit$pred.locfit) + geom_point(aes(x.roc[i.opt],y.roc[i.opt]), col = 'purple', cex = 3)
roc_auc_vec(truth = test$label, estimate = pred.locfit$pred.locfit, event_level = 'second') #da aggiungere nel plot

i.opt <- which(data.roc.locfit$y.roc>=0.75 & data.roc.locfit$x.roc<=0.25)[317]
p0.opt <- seq(0,1, length.out = 1000)[317]
sens.locfit <- data.roc.locfit$y.roc[i.opt]
spec.locfit <- 1 - data.roc.locfit$x.roc[i.opt]
F1.locfit <- data.roc.locfit$F1.score[i.opt]
accuracy.locfit <- data.roc.locfit$accuracy[i.opt]
precision.locfit <- data.roc.locfit$precision[i.opt]


p.hat.locfit <- predict(fit_locfit, newdata = x, type = 'response') #1/(1+exp(-mdl$coefficients[1]-mdl$coefficients[2]*x))
x11()
ggplot() +
  geom_line(aes(x$MMSE, p.hat.locfit), col='blue') +
  geom_point(aes(train$MMSE, as.numeric(train$label==1))) + ylab('Probability to be Demented') +xlab('MMSE') + labs(title = 'Logistic Curve') + theme(plot.title = element_text(hjust = 0.5))

#Spline
mdl_spline <- glm(label ~ bs(MMSE, degree = 1), data = train, family = binomial )
summary(mdl_spline)
1- mdl_spline$deviance/mdl_spline$null.deviance
mdl_spline$aic
pred.spline <- predict(mdl_spline, newdata = test, type = 'response')
pred.spline <- as.data.frame(pred.spline)
data.roc <- roc.curve(pred.spline, test)
roc_auc_vec(truth = test$label, estimate = pred.spline$pred.spline, event_level = 'second') #da aggiungere nel plot

knots <- quantile(train$MMSE, probs = c(0.5))
mdl_spline2 <- glm(label ~ bs(MMSE, degree = 1, knots = knots), data = train, family = binomial )
summary(mdl_spline2)
1- mdl_spline2$deviance/mdl_spline2$null.deviance
mdl_spline2$aic
pred.spline2 <- predict(mdl_spline2, newdata = test, type = 'response')
pred.spline2 <- as.data.frame(pred.spline2)
data.roc.spline <- roc.curve(pred.spline2, test)
x11()
plot.roc(data.roc.spline, pred.spline2$pred.spline2) + geom_point(aes(x.roc[i.opt],y.roc[i.opt]), col = 'purple', cex = 3)
roc_auc_vec(truth = test$label, estimate = pred.spline2$pred.spline2, event_level = 'second') #da aggiungere nel plot

i.opt <- which(data.roc.spline$y.roc>=0.75 & data.roc.spline$x.roc<=0.25)[295]
p0.opt <- seq(0,1, length.out = 1000)[295]
sens.spline <- data.roc.spline$y.roc[i.opt]
spec.spline <- 1 - data.roc.spline$x.roc[i.opt]
F1.spline <- data.roc.spline$F1.score[i.opt]
accuracy.spline <- data.roc.spline$accuracy[i.opt]
precision.spline <- data.roc.spline$precision[i.opt]
# knots <- quantile(train$MMSE, probs = c(0.18, 0.26, 0.4, 0.5, 0.7))
# boundary.knots <- quantile(train$MMSE, probs = c(0.05, 0.8))
# mdl_n_spline <- glm(label ~ ns(MMSE, knots = knots, Boundary.knots = boundary.knots), data = train, family = binomial )
# summary(mdl_n_spline)
# pred.n.spline <- predict(mdl_n_spline, newdata = test, type = 'response')
# pred.n.spline <- as.data.frame(pred.n.spline)
# data.roc <- roc.curve(pred.n.spline, test)
# roc_auc_vec(truth = test$label, estimate = pred.n.spline$pred.n.spline, event_level = 'second') #da aggiungere nel plot

#AIC model spline = 217, AIC model Natural Spline=220
1-mdl_spline$deviance/mdl_spline$null.deviance
1-mdl_spline2$deviance/mdl_spline2$null.deviance
#scelgo secondo modello
p.hat.spline <- predict(mdl_spline2, newdata = x, type = 'response')
x11()
ggplot() +
  geom_line(aes(x$MMSE, p.hat.spline), col='blue', lwd = 1.3) +
  geom_point(aes(train$MMSE, as.numeric(train$label==1))) + ylab('Probability to be Demented') +xlab('MMSE') + labs(title = 'Logistic Curve') + theme(plot.title = element_text(hjust = 0.5))


#General Plot
colors <- c("Parametric" = "red", "Polynomial" = "blue", "Loc likelihood" = "green", "Spline" = "purple")
x11()
ggplot()+
  geom_line(aes(x$MMSE,p.hat.param, col='Parametric'), size = 0.8)+
  geom_line(aes(x$MMSE,p.hat.poly, col='Polynomial'), size = 0.8) +
  geom_line(aes(x$MMSE, p.hat.locfit, col='Loc Likelihood'), size = 0.8) +
  #geom_line(aes(x[,1],p.hat.test), col='orange') +
  geom_line(aes(x$MMSE, p.hat.spline, col='Spline'), size = 0.8) +
  geom_point(aes(train$MMSE, as.numeric(train$label==1)), size = 0.8) +
  xlab('MMSE') + ylab('Probability to be Demented') + labs(color = "Legend") +
  scale_color_manual(values = colors)






















