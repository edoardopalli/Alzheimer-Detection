setwd("C:/Users/franc/Desktop/NONPA/PROGETTO")

library(survival)
library(survminer)
library(dplyr) 
library(ggplot2)
library(knitr)
library(broom)

library(stringr)
dataset_old <- read.csv("oasis_longitudinal.csv")
head(dataset_old)
colnames(dataset_old)

## DEMENTI  
dataset_old$ID = str_sub(dataset_old$Subject.ID,-3) 
dementi_prima <- dataset_old[dataset_old$Group=='Demented'& dataset_old$Visit==1,]
dementi_prima$status <- 1
dementi_prima$group_new <- 'Demented'
##NON DEMENTI
nondementi_last <- dataset_old[dataset_old$Group=='Nondemented',]
nondementi_last <-nondementi_last%>%group_by(ID)%>% mutate(max_visit=max(Visit))%>%ungroup()%>% filter(Visit==max_visit)
length(unique(nondementi_last$ID))
nondementi_last$status<-0
nondementi_last$group_new <- 'Nondemented'
## CONVERTED

converted_first <- dataset_old[dataset_old$Group=='Converted'&dataset_old$CDR >0,]
converted_first <-converted_first%>%group_by(ID)%>% mutate(min_visit=min(Visit))%>%ungroup()%>% filter(Visit==min_visit)
converted_first$group_new <- 'Demented'
converted_first$status <- 1
#ALL TOGETHER 
dataset <- rbind(nondementi_last[,-17],dementi_prima,converted_first[,-17])

dataset$time_y <- dataset$MR.Delay / 365
dataset$status_fact <- factor(dataset$Group)
old <- dataset[dataset$Age>90,]
sesantenni<- dataset[dataset$Age>60&dataset$Age<70,]
numero_174 <- dataset_old[dataset_old$Subject.ID=='OAS2_0174',]

x11()
ggplot(data=dataset,aes(x=ID,y=time_y)) + 
  geom_bar(stat='identity',width=0.2) +
  geom_point(aes(color=status_fact,shape=status_fact),size=6) +
  coord_flip()

#without demented at first

dataset_nodem <- rbind(nondementi_last[,-17],converted_first[,-17])

dataset_nodem$time_y <- dataset_nodem$MR.Delay / 365
dataset_nodem$status_fact <- factor(dataset_nodem$Group)

x11()
ggplot(data=dataset_nodem,aes(x=ID,y=time_y)) + 
  geom_bar(stat='identity',width=0.2) +
  geom_point(aes(color=status_fact,shape=status_fact),size=6) +
  coord_flip()

## Survival Object

head(Surv(dataset$MR.Delay, dataset$group_new=='Demented')) #Surv(time, event) -> create a surv obj

# Kaplan-Meier estimator for survival curve --------------

fit <- survfit(Surv(dataset$MR.Delay, status==1) ~ 1, data = dataset)
# status==1 -> when the event occurred

summary(fit)

kable(head(tidy(fit),20))
surv_median(fit)# median survival times: time at which the survival probability, S(t), is 0.5.


# plots
x11()
plot(fit, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', col='red',
     main="Kaplan-Meier Curve for Lung Cancer Survival")

ggsurvplot(fit,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=90,
           title="Kaplan-Meier Curve ")

dev.off()

cumulative_incidence <- 1 - fit$surv # cumulative failure probability CFP(t) = P(T<t)

# adding fun= 'event'
x11()
ggsurvplot(fit,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=90,
           fun='event',  #with event you plot 1-y
           title="Cumulative Incidence Curve for Alzheimer Survival")


dev.off()

 

H <- fit$cumhaz # Nelson-Aalen cumulative hazard rate estimator

x11()
ggsurvplot(fit,
           risk.table = TRUE, # Add risk table
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=90,
           fun='cumhaz', #plots the cumulative hazard function
           title="Cumulative Hazard Curve for Alzheimer Survival")

#######################

# Kaplan-Meier Curves using AGE

#investigate if there is a difference in terms of survival among the two groups.
# works best with cathegorical variable
min(dementi_prima$Age) #61
max(dementi_prima$Age) #96
min(converted_first$Age) #65
max(converted_first$Age) #92
min(nondementi_last$Age) #62
max(nondementi_last$Age) #97
min(dataset$Age) #61
max(dataset$Age) #97
fit.age <- survfit(Surv(MR.Delay, status) ~ Age, data=dataset)
summary(fit.age)
x11()
ggsurvplot(fit.age, conf.int = F, risk.table.col = "strata", legend='none')

# WRONG: this will get a separate curve for every unique value of age

# distribution of age
x11()
hist(dataset$Age, xlab='Age [years]', main='Histogram of age in Lung Cancer Data',freq=FALSE)
summary(dataset$Age)


#cut point at 70 ->'young' subjects with age less or equal than 70 
# ->'old' subject aged more than 70.

dataset$cutpoint <- cut(dataset$Age, breaks=c(60,70,80,90 , Inf), labels=c("60", "70","80","90"))


#KM plot

fit.age <- survfit(Surv(MR.Delay, status) ~ cutpoint, data=dataset)
x11()
ggsurvplot(fit.age, conf.int = T,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=90,
           legend.labs=c("60", "70","80","90"), legend.title="Age class",  
           palette=c("darkblue","cyan3","red","green"), 
           title="Kaplan-Meier Curves by age class for Lung Cancer Survival")
#COMMENT:
# It looks like there's some differences in the curves between "old" and
# "young" patients, with older patients having slightly worse survival
# odds. Is there statistical evidence for that difference?
dev.off()

log_rank_test <- survdiff(Surv(MR.Delay, status) ~ cutpoint, data=dataset)
log_rank_test  #0.09

#Pvalue small->the difference in survival between those groups is significant

#???????????????????POSSO FARLO CON PIù GRUPPI????????????????????
#gruppi 1 e 2
hazard_ratio <- (log_rank_test$obs[1]/log_rank_test$exp[1])/(log_rank_test$obs[2]/log_rank_test$exp[2])
hazard_ratio #1.30735
hazard_ratio_14 <- (log_rank_test$obs[1]/log_rank_test$exp[1])/(log_rank_test$obs[4]/log_rank_test$exp[4])
hazard_ratio_14 #2.4392
                        # HR< 1 -> group 1 protective factor ( less risk)
                        # HR= 1 -> similar risk
                        # HR> 1 -> group 1 risk factor (more risk)

#group1= younger than 70
#HR = 0.674 < 1$ indicating that the risk of deaths in younger than 70
#years old is 0.674 times the risk in older than 70: being young is a
#protective factor.


# univariate Cox regression model ------


cox.age <- coxph(Surv(MR.Delay, status) ~ Age, data = dataset)

summary(cox.age) #pvalue =0.0116 small ->we can say that age is significant


# baseline survival function S_0(t)
x11()
plot(survfit(cox.age, data=dataset), 
     col="darkorange2", lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Baseline estimated survival probability')
grid()


# BUILD A NEW DATAFRAME ( TEST DATASET FOR PREDICTION)
age_df <- with(dataset,
               data.frame(Age = c(50,65,80) )#consider ages equal to 50, 65 and 80
)

#estimate survival:


fit.age <- survfit(cox.age, newdata = age_df)
fit.age
# HIGHER THE AGE LOWER THE SURVIVAL PROB


plot(fit.age, conf.int=T,
     col=c("dodgerblue2","navy","darkmagenta"), lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Adjusted Survival Probability Plot')
grid()
legend('topright', c("Age = 50", "Age = 65", "Age = 80"),
       lty=c(1,1,1), lwd=c(2,2,2), col=c("dodgerblue2","navy","darkmagenta"))

## Multivariate Cox regression -----


glimpse(dataset)

dataset$sex <- as.factor(dataset$M.F)

# Cox's regression model:
mod.cox <- coxph(Surv(MR.Delay, status) ~ Age + sex + EDUC + MMSE, data =  dataset)
summary(mod.cox)

# COMMENT:-----
# The p-values for all three overall tests (likelihood, Wald, and score)
# are extremely small, indicating that the model is significant. These
# tests evaluate the omnibus null hypothesis that all of the $\beta$s are
# 0. In the above example, the test statistics are in close agreement, and
# the omnibus null hypothesis is soundly rejected.
# 
# In the multivariate Cox analysis, the covariates sex and ph.karno are
# significant (p \< 0.05). However, the covariates age and wt.loss fail to
# be significant.
# 
# The HR for sex is exp(coef) = exp(0.514) = 1.67 with 95% CI = [1.19;
# 2.35]. The hazard ratios of covariates are interpretable as
# multiplicative effects on the hazard. For example, holding the other
# covariates constant, being a male increases the hazard by a factor of
# 1.67, or 67%. We conclude that, being male is associated with bad
# prognostic.
# 
# The HR for ph.karno is exp(coef) = exp(-0.013) = 0.987 with 95% CI =
# [0.975;0.999], indicating a strong relationship between the ph.karno
# value and decreased risk of death. Holding the other covariates
# constant, a higher value of ph.karno is associated with a better
# survival.
# 
# The hazard ratio HR of age is exp(coef) = 1.01, with a 95% CI =
# [0.996;1.035]. Because the confidence interval for HR includes 1, these
# results indicate that age makes a smaller contribution to the difference
# in the HR after adjusting for the other covariates.
# 
# Similarly, the hazard ratio HR of wt.loss is exp(coef) = 0.998, with a
# 95% CI = [0.985;1.010]. Because the confidence interval for HR includes
# 1, these results indicate that wt.loss makes a smaller contribution to
# the difference in the HR after adjusting for the other covariates.
################################################################################
### Visualizing Hazard ratios----

#Hr and its CIs 
ggforest(mod.cox, data=dataset)
#standard errors represented by bar


#baseline survival function $S_0(t)$

plot(survfit(mod.cox, data=dataset), 
     col="darkorange2", lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Baseline estimated survival probability')
grid()


### Cox Model Assumptions and Goodness of fit

#graphical check the goodness of fit: Martingale Residuals------

plot(predict(mod.cox), residuals(mod.cox, type='martingale'),
     xlab='Fitted values', ylab='Martingale residuals', main='Residual Plot', las=1)
# Add a line for residual=0
abline(h=0, col='red')
# Fit a smoother for the points
lines(smooth.spline(predict(mod.cox), residuals(mod.cox, type='martingale')), col='blue')
# you want the residual around 0

ggcoxdiagnostics(mod.cox, type = "martingale")

# hard coded:-----
cox_model_matrix <- model.matrix(mod.cox)
centered_model_matrix <- sweep(cox_model_matrix,MARGIN = 2,STATS = mod.cox$means,FUN = "-")
manual_linear_pred <- c(centered_model_matrix%*%mod.cox$coefficients)
bese_haz_est <- basehaz(mod.cox)

# Martingale residuals
basehaz_for_i <- Vectorize(function(t){
  pos_haz <- which(t<=bese_haz_est$time)[1]
  bese_haz_est$hazard[pos_haz]
})

obs_NA <- which(is.na(dataset$wt.loss))
dataset_no_NA <- dataset[-obs_NA,]
martingale_res_manual <- (dataset_no_NA$status-1)- basehaz_for_i(t=dataset_no_NA$time)*exp(manual_linear_pred)
######

# deviance residual (normalized transform of the martingale residual)----


#These residuals should be roughly symmetrically distributed about zero
#with a standard deviation of 1.
# 
# -   Positive values correspond to individuals that "died too soon"
#     compared to expected survival times.
# -   Negative values correspond to individual that "lived too long".
# -   Very large or small values are outliers, which are poorly predicted
#     by the model.


ggcoxdiagnostics(mod.cox, type = "deviance")

# hard code
deviance_res_manual <- sign(martingale_res_manual)*(-2*(martingale_res_manual+(dataset_no_NA$status-1)*log((dataset_no_NA$status-1)-martingale_res_manual)))^(1/2)


# Schoenfeld residuals

ggcoxdiagnostics(mod.cox, type = "schoenfeld")

#hard code:-----
df_4_schoenfeld <- dataset %>%
  tidyr::drop_na(age, sex, ph.karno, wt.loss) %>%
  arrange(time)

df_4_schoenfeld_sub <- df_4_schoenfeld %>%
  select(age, sex, ph.karno, wt.loss) %>% 
  mutate(sex=if_else(sex=="Male",1,0))

cox_coef <- mod.cox$coefficients
X_denom <- exp(data.matrix(df_4_schoenfeld_sub) %*% cox_coef)

X_num <-
  sweep(
    x = df_4_schoenfeld_sub,
    MARGIN = 1,
    STATS = X_denom,
    FUN = "*"
  )

X_bar <- matrix(nrow = nrow(X_num), ncol = ncol(X_num))
X_den <- numeric(nrow(X_num))

for(i in 1:nrow(X_num)){
  time_instant <- df_4_schoenfeld$time[i]
  risk_set <- which(time_instant<=df_4_schoenfeld$time)
  X_bar[i,] <- colSums(X_num[risk_set,])
  X_den[i] <- sum(X_denom[risk_set])
}

schoenfeld_calc <- df_4_schoenfeld_sub-sweep(X_bar,MARGIN = 1,STATS = X_den,FUN = "/")
schoenfeld_manual <- schoenfeld_calc[df_4_schoenfeld$status==2,]
############
####  plot log(-log(KM(t))) vs. t or log(t) ####
# look for parallelism. This can be done only for categorical covariates.

sex.km <- survfit(Surv(time, status) ~ sex, data = dataset[!is.na(dataset$wt.loss) & !is.na(dataset$ph.karno),])

plot(sex.km, fun='cloglog', 
     col=c("deeppink2","dodgerblue2"), lwd=2, lty=1,
     ylab="log(-log(Survival Probability))")
grid()
legend('topleft', c("Female", "Male"),
       lty=c(1,1), lwd=c(2,2), col=c("deeppink2","dodgerblue2"))

#Curves seem to be parallel -\> PH assumption seems to be satisfied for gender.



#Test for PH using scaled Schoenfeld test for PH
# H0: Hazards are proportional vs  H1: Hazards are NOT proportional


test.ph <- cox.zph(mod.cox)
test.ph

# pvalue small:  the global test is statistically significant.
# we can not assume the proportional hazards
# check if for one of the characteristic is particularly

#scaled schoenfeld residuals:

par(mfrow=c(2,2))
for(i in 1:4){
  plot(test.ph[i])
  abline(h=0, col='red')
}

#or
ggcoxdiagnostics(mod.cox, type = "scaledsch")

#############################################
## Stratified Cox Model


mod.cox.strata <- coxph(Surv(time, status) ~ age + sex + strata(ph.karno) + wt.loss, data =  dataset)
summary(mod.cox.strata)

# PH ASSUMPTIONS
test.ph.strata <- cox.zph(mod.cox.strata)
test.ph.strata

# PH assumptions are satisfied for all variables and for the global model.
# if pvalue large