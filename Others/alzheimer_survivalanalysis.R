#library(sf)
library(dplyr)
#library(gstat)
#library(maps)
#library(maptools)
#library(rgdal)
#library(sp)
#library(shp2graph)
#library(lubridate)

library(survival)
library(survminer)
library(dplyr) 
library(ggplot2)
library(knitr)
library(broom)


# controllo ingressi in area c
setwd('C:/Users/Elena/Desktop/Elena/Polimi/MAGISTRALE/Nonparametric statistics/Progetto/github repository/ALZHEIMER_prognonpa/Elena')
dataset_xsectional <- read.csv("oasis_cross-sectional.csv", header = T)
dataset_longitudinal <- read.csv("oasis_longitudinal.csv", header = T)


# voglio tenere solo i nondemented (censored) e i converted (event)

dataset_surv <- dataset_longitudinal

dataset_surv <- dataset_surv %>%  filter(Group != 'Demented')

# v <- vector()
# 
# for (i in 1:length(unique(dataset_surv$Subject.ID))){
#   print(i)
#   for(j in unique(dataset_surv$Subject.ID)){
#     print(j)
#     for (k in 1:227){
#       print(k)
#       if (dataset_surv$Subject.ID[k] == j){
#         v[i] <- k
#         print(dataset_surv$Visit[k])
#       }
#     }
#   }
# }
# 
# 
# 
# v <- vector()
# id <- 'OAS2_0001'
# 
# for (i in 1:length(unique(dataset_surv$Subject.ID))){
#   print(i)
#   for(k in 1:227){
#     print(k)
#     for(j in unique(dataset_surv$Subject.ID)){
#       if (dataset_surv$Subject.ID[k] == id){
#         v[i] <- dataset_surv$Visit[k]
#       }
#     }
#   }
# }





#v <- vector()
v <- rep(0, 227)
#w <- vector()
#v[1] <- 1

for (i in 2:226){
  #for (j in 2:86){
  if(dataset_surv$Group[i] != 'Converted'){
    if(dataset_surv$Subject.ID[i] == dataset_surv$Subject.ID[i-1] & dataset_surv$Subject.ID[i] != dataset_surv$Subject.ID[i+1]){
      #v[j] <- i
      v[i] <- 1
      #append(w, i)
      #break
    }
    #break
  }
  #}
}

#v[86] <- 227
v[227] <- 1
#append(w, 227)

v[22] <- 1
v[24] <- 1
v[38] <- 1
v[51] <- 1
v[69] <- 1
v[119] <- 1
v[136] <- 1
v[147] <- 1
v[159] <- 1
v[169] <- 0  # 168-169 è outlier? paziente OAS2_0131
v[171] <- 1
v[184] <- 1
v[186] <- 1
v[215] <- 1





#dataset_surv %>% slice(dataset_surv, v)
dataset_surv <- filter(dataset_surv, v!=0)






dataset_surv$ID <- factor(dataset_surv$Subject.ID , labels = seq(1:85))  # 85 invece di 86

#dataset_longitudinal$ID <- factor(seq(1:length(unique(dataset_longitudinal$Subject.ID))))
dataset_surv$time_y <- dataset_surv$MR.Delay
#dataset_longitudinal$status_fact <- ifelse(dataset_longitudinal$Group == 'Nondemented',1, 2)
dataset_surv$status_fact <- factor(dataset_surv$Group , labels = (c('Event', 'Censor')))
subs <- dataset_surv[5:10,]

x11()
ggplot(data=subs,aes(x=ID,y=time_y)) + 
  geom_bar(stat='identity',width=0.2) +
  geom_point(aes(color=status_fact,shape=status_fact),size=6) +
  coord_flip()


##################################






head(Surv(dataset_surv$time_y, dataset_surv$status_fact=='Event'))


# kaplan - meier estimator for the survival of our group of patients

fit <- survfit(Surv(time_y, status_fact=='Event') ~ 1, data = dataset_surv)

summary(fit)

kable(head(tidy(fit),20))

median_St <- fit$time[fit$surv<=0.6][1]  # non va sotto 0.5923960, dunque non posso trovare la median perchè dovrei scendere sotto 0.5!

#surv_median(fit) 

x11()
plot(fit, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', col='red',
     main="Kaplan-Meier Curve for Dementia Survival")

x11()  
ggsurvplot(fit,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           #surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=120,
           title="Kaplan-Meier Curve for Dementia Survival")




# 
# # riduco a 13 e 13
# 
# data_surv_ridotto <- slice(dataset_surv, 1:16, 18, 26, 44, 51, 55, 60, 63, 69, 70, 81)
# 
# fit2 <- survfit(Surv(time_y, status_fact=='Event') ~ 1, data = data_surv_ridotto)
# 
# summary(fit2)
# 
# kable(head(tidy(fit2),20))
# 
# median_St2 <- fit2$time[fit2$surv<=0.5][1]  
# 
# #surv_median(fit) 
# 
# x11()
# plot(fit2, conf.int = T, xlab='Time [days]', ylab = 'Survival Probability', col='red',
#      main="Kaplan-Meier Curve for Dementia Survival")
# 
# x11()  
# ggsurvplot(fit2,
#            risk.table = TRUE, # Add risk table
#            risk.table.col = "strata", # Change risk table color by groups
#            surv.median.line = "hv", # Specify median survival
#            ggtheme = theme_bw(), # Change ggplot2 theme
#            break.time.by=120,
#            title="Kaplan-Meier Curve for Dementia Survival")
# 











cumulative_incidence <- 1 - fit$surv

x11()
ggsurvplot(fit,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=120,
           fun='event',
           title="Cumulative Incidence Curve for Dementia Survival")



# 
# cumulative_incidence2 <- 1 - fit2$surv
# 
# x11()
# ggsurvplot(fit2,
#            risk.table = TRUE, # Add risk table
#            risk.table.col = "strata", # Change risk table color by groups
#            surv.median.line = "hv", # Specify median survival
#            ggtheme = theme_bw(), # Change ggplot2 theme
#            break.time.by=120,
#            fun='event',
#            title="Cumulative Incidence Curve for Dementia Survival")
# 





H <- fit$cumhaz

x11()
ggsurvplot(fit,
           risk.table = TRUE, # Add risk table
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=120,
           fun='cumhaz',
           title="Cumulative Hazard Curve for Dementia Survival")










fit.age <- survfit(Surv(time_y, status_fact=='Event') ~ Age, data=dataset_surv)
x11()
ggsurvplot(fit.age, conf.int = F, risk.table.col = "strata", legend='none')


x11()
hist(dataset_surv$Age, xlab='Age [years]', main='Histogram of age in Dementia Data')


summary(dataset_surv$Age)







fit.educ <- survfit(Surv(time_y, status_fact=='Event') ~ EDUC, data=dataset_surv)
x11()
ggsurvplot(fit.educ, conf.int = F, risk.table.col = "strata", legend='none')


x11()
hist(dataset_surv$EDUC, xlab='EDUC [levels]', main='Histogram of education in Dementia Data')


summary(dataset_surv$EDUC)




fit.MMSE <- survfit(Surv(time_y, status_fact=='Event') ~ MMSE, data=dataset_surv)
x11()
ggsurvplot(fit.MMSE, conf.int = F, risk.table.col = "strata", legend='none')


x11()
hist(dataset_surv$MMSE, xlab='MMSE', main='Histogram of MMSE in Dementia Data')


summary(dataset_surv$MMSE)




fit.etiv <- survfit(Surv(time_y, status_fact=='Event') ~ eTIV, data=dataset_surv)
x11()
ggsurvplot(fit.etiv, conf.int = F, risk.table.col = "strata", legend='none')


x11()
hist(dataset_surv$eTIV, xlab='eTIV', main='Histogram of eTIV in Dementia Data')


summary(dataset_surv$eTIV)
















lung$agecat70 <- cut(lung$age, breaks=c(0, 70, Inf), labels=c("young", "old"))


fit.age <- survfit(Surv(time, status) ~ agecat70, data=lung)
ggsurvplot(fit.age, conf.int = T,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by=90,
           legend.labs=c("Young (<= 70)","Old (> 70)"), legend.title="Age class",  
           palette=c("darkblue","cyan3"), 
           title="Kaplan-Meier Curves by age class for Lung Cancer Survival")


log_rank_test <- survdiff(Surv(time, status) ~ agecat70, data=lung)
log_rank_test








glimpse(dataset_surv)

mod.cox <- coxph(Surv(time_y, status_fact) ~ Age + M.F + EDUC + SES, MMSE, CDR, eTIV, nWBV, ASF, data =  dataset_surv)
summary(mod.cox)

















