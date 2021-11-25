library(readxl)
library(survival)
library(survminer)
library(dplyr)
library(Hmisc)
library(splines)
library(naniar)
library(zoo)
library(mice)
setwd('C:/Users/derck/surfdrive/AMC/H3. Quantib mortaliteit')

df <- read_excel("complete_data - Correct2.xlsx")

df<-subset(df, quality==1)
df<-subset(df, observatietijd>-1)

df <- replace_with_na(df,replace=list(length = -9))
df <- replace_with_na(df,replace=list(gewicht = -1))
df <- replace_with_na(df,replace=list(NYHA = -1))
df$NYHA <-ifelse(df$NYHA>2,1,2)

#LVF
df <- replace_with_na(df,replace=list(LVF = -1))
df$LVF<-ifelse(df$LVF>45,0,df$LVF)
df$LVF<-ifelse(df$LVF==1,0,df$LVF)
df$LVF<-ifelse(df$LVF<=45&df$LVF>5 ,1,df$LVF)
df$LVF<-ifelse(df$LVF>1 & df$LVF<5,1,df$LVF)

#eGFR
df$eGFR <-ifelse(df$geslacht==1, 175*(df$kreatine/88.4)^-1.154 * (df$Age)^-0.203, (175*(df$kreatine/88.4)^-1.154 * (df$Age)^-0.203)* 0.742)
df$CKD<- ifelse(df$eGFR<60,1,0)

#Imputeren 
notimputing<-df[c(1:4,13:33)]
imputing<-select(df, c(1,5:12,34,38))

x<-mice(data=imputing, m=1, maxit=10, seed=2)
complete<- complete(x, action=1)

df<-inner_join(notimputing,complete, by=c("Pseudo"))

rm(x, complete, imputing, notimputing)

#Prepare dataset
df$musclearea<-(df$Psoas_muscle+df$Abdominal_muscle+df$Longspine_muscle)
df$SMI<-df$musclearea/((df$length/100)*(df$length/100))

df$musclearea_thres<-(df$Psoas_muscle_thresholded+df$Abdominal_muscle_thresholded+df$Longspine_muscle_thresholded)
df$SMI_thres<-df$musclearea_thres/((df$length/100)*(df$length/100))

df$HUmuscle<-(df$HU_Psoas_muscle*df$Psoas_muscle+df$HU_Abdominal_muscle*df$Abdominal_muscle+df$HU_Longspine_muscle*df$Longspine_muscle)/df$musclearea

df$BMI<-df$gewicht/((df$length/100)*(df$length/100))
df$BSA<-sqrt((df$length*df$gewicht)/3600)

df$IMAT<-df$musclearea-df$musclearea_thres
df$IMATper<-df$IMAT/df$musclearea*100

#tijden instellen
df$overleden_30d<-ifelse(df$Overleden==1 & df$observatietijd<=30, 1, 0)
df$overleden_30d<-as.factor(df$overleden_30d)

df$observatietijd30<-df$observatietijd
df$observatietijd30[df$observatietijd>30]<- 31

df$overleden_1jaar<-ifelse(df$Overleden==1 & df$observatietijd<=365, 1, 0)
df$overleden_1jaar<-as.factor(df$overleden_1jaar)

df$observatietijd365<-df$observatietijd
df$observatietijd365[df$observatietijd>365]<- 366

#opdelen in 3 groepen
male<-subset(df, df$geslacht==1)
female<-subset(df, df$geslacht==2)

cut2(male$SMI_thres, g=3,onlycuts=TRUE)
cut2(female$SMI_thres, g=3,onlycuts=TRUE)

male$muscleareaclass <- as.factor(cut2(male$SMI_thres,g=3))
levels(male$muscleareaclass) <- c('low','mid', 'high')

female$muscleareaclass <- as.factor(cut2(female$SMI_thres,g=3))
levels(female$muscleareaclass) <- c('low','mid', 'high')

cut2(male$HUmuscle, g=3,onlycuts=TRUE)
cut2(female$HUmuscle, g=3,onlycuts=TRUE)

male$HUclass <- as.factor(cut2(male$HUmuscle,g=3))
levels(male$HUclass) <- c('low','mid', 'high')

female$HUclass <- as.factor(cut2(female$HUmuscle,g=3))
levels(female$HUclass) <- c('low','mid', 'high')

cut2(male$IMATper, g=3,onlycuts=TRUE)
cut2(female$IMATper, g=3,onlycuts=TRUE)

male$IMATclass <- as.factor(cut2(-male$IMATper,g=3))
levels(male$IMATclass) <- c('low','mid', 'high')

female$IMATclass <- as.factor(cut2(-female$IMATper,g=3))
levels(female$IMATclass) <- c('low','mid', 'high')

df<- as.data.frame(rbind(female, male))

#theme
theme <- theme(axis.line = element_line(colour = "black", size=0.4),
               panel.grid.major = element_line(colour = "grey90", size=0.4),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_blank(),
               axis.text = element_text(color = "black",size=12, face = 'plain'),
               text=element_text(color = "black",size=12, face = 'plain')
) 

######
#SMI
#Kaplan Meier
fit <- survfit(Surv(observatietijd, Overleden) ~ muscleareaclass, data = df)
fit
summary(fit)$table

p<-ggsurvplot(fit,xlab="Time (Years)",xscale="d_y", size=0.7,axes.offset=FALSE,
           break.time.by=(365.25),pval=F, test.for.trend=FALSE,
           conf.int=FALSE,label.curves=TRUE, ggtheme=theme_minimal()+theme,
           tables.theme=theme_cleantable(), xlim=c(0,1700),expand_limits=c(0,0),
           risk.table=T, risk.table.fontsize=4, risk.table.y.text = F, censor=F, palette = c("red",'black','blue'), ylim = c(0.4, 1.0))
print(p)
#export => Save as PDF => 5x6 inch

#log rank
surv_diff <- survdiff(Surv(observatietijd, Overleden) ~ muscleareaclass, data = df)
surv_diff
res <- pairwise_survdiff(Surv(observatietijd, Overleden) ~ muscleareaclass, data = df)
res
#Cox regression 
df$observatietijd<-ifelse(df$observatietijd==0, 1, df$observatietijd)

#catagorical
df$muscleareaclass2[df$muscleareaclass== 'low'] <- 1
df$muscleareaclass2[df$muscleareaclass== 'mid'|df$muscleareaclass== 'high'] <- 2
df$muscleareaclass2<-factor(df$muscleareaclass2, labels=c('low', 'high'))
df$muscleareaclass2<-relevel(df$muscleareaclass2, ref = 'high')

#crude model
#1jaar
res.cox <- coxph(Surv(observatietijd365, overleden_1jaar==1) ~ muscleareaclass2 , data = df)
summary(res.cox)
#5jaar
res.cox <- coxph(Surv(observatietijd, Overleden==1) ~ muscleareaclass2 , data = df)
summary(res.cox)

#adjusted model
#1jaar
res.cox <- coxph(Surv(observatietijd365, overleden_1jaar==1) ~ muscleareaclass2+Age+ geslacht+EURO2+NYHA+COPD+LVF+Accessroute+CKD , data = df)
summary(res.cox)
#5jaar
res.cox <- coxph(Surv(observatietijd, Overleden==1) ~ muscleareaclass2+Age+ geslacht+EURO2+NYHA+COPD+LVF+Accessroute+CKD , data = df)
summary(res.cox)

###### 
#SMD
#Kaplan Meier
fit <- survfit(Surv(observatietijd, Overleden) ~ HUclass, data = df)
fit
summary(fit)$table

p<-ggsurvplot(fit,xlab="Time (Years)",xscale="d_y", size=0.7,axes.offset=FALSE,
              break.time.by=(365.25),pval=F, test.for.trend=FALSE,
              conf.int=FALSE,label.curves=TRUE, ggtheme=theme_minimal()+theme,
              tables.theme=theme_cleantable(), xlim=c(0,1700),expand_limits=c(0,0),
              risk.table=T, risk.table.fontsize=4,  risk.table.y.text = F, censor=F, palette = c("red",'black','blue'), ylim = c(0.4, 1.0))
print(p)
#export => Save as PDF => 5x6 inch

#log rank
surv_diff <- survdiff(Surv(observatietijd, Overleden) ~ HUclass, data = df)
surv_diff
res <- pairwise_survdiff(Surv(observatietijd, Overleden) ~ HUclass, data = df)
res
#Cox regression 
df$observatietijd<-ifelse(df$observatietijd==0, 1, df$observatietijd)

#categorical
df$HUclass2[df$HUclass== 'low'] <- 1
df$HUclass2[df$HUclass== 'mid'|df$HUclass== 'high'] <- 2
df$HUclass2<-factor(df$HUclass2, labels=c('low', 'high'))
df$HUclass2<-relevel(df$HUclass2, ref = 'high')

#crude model
#1jaar
res.cox <- coxph(Surv(observatietijd365, overleden_1jaar==1) ~ HUclass2 , data = df)
summary(res.cox)
#5jaar
res.cox <- coxph(Surv(observatietijd, Overleden==1) ~ HUclass2 , data = df)
summary(res.cox)

#adjusted model
#1jaar
res.cox <- coxph(Surv(observatietijd365, overleden_1jaar==1) ~ HUclass2+Age+ geslacht+EURO2+NYHA+COPD+LVF+Accessroute+CKD , data = df)
summary(res.cox)
#5jaar
res.cox <- coxph(Surv(observatietijd, Overleden==1) ~ HUclass2+Age+ geslacht+EURO2+NYHA+COPD+LVF+Accessroute+CKD , data = df)
summary(res.cox)

###### 
#IMAT
#Kaplan Meier
fit <- survfit(Surv(observatietijd, Overleden) ~ IMATclass, data = df)
fit
summary(fit)$table

p<-ggsurvplot(fit,xlab="Time (Years)",xscale="d_y", size=0.7,axes.offset=FALSE,
           break.time.by=(365.25),pval=F, test.for.trend=FALSE,
           conf.int=FALSE,label.curves=TRUE, ggtheme=theme_minimal()+theme,
           tables.theme=theme_cleantable(), xlim=c(0,1700),expand_limits=c(0,0),
           risk.table=T, risk.table.fontsize=4, censor=F, risk.table.y.text = F, palette = c("red",'black', "blue"), ylim = c(0.4, 1.0))
print(p)
#export => Save as PDF => 5x6 inch

#log rank
surv_diff <- survdiff(Surv(observatietijd, Overleden) ~ IMATclass, data = df)
surv_diff
res <- pairwise_survdiff(Surv(observatietijd, Overleden) ~ muscleareaclass, data = df)
res
#Cox regression 
df$observatietijd<-ifelse(df$observatietijd==0, 1, df$observatietijd)

#catagorical
df$IMATclass2[df$IMATclass== 'low'] <- 1
df$IMATclass2[df$IMATclass== 'mid'|df$IMATclass== 'high'] <- 2
df$IMATclass2<-factor(df$IMATclass2, labels=c('low', 'high'))
df$IMATclass2<-relevel(df$IMATclass2, ref = 'high')

#crude model
#1jaar
res.cox <- coxph(Surv(observatietijd365, overleden_1jaar==1) ~ IMATclass2 , data = df)
summary(res.cox)
#5jaar
res.cox <- coxph(Surv(observatietijd, Overleden==1) ~ IMATclass2 , data = df)
summary(res.cox)

#adjusted model  
#1jaar
res.cox <- coxph(Surv(observatietijd365, overleden_1jaar==1) ~ IMATclass2+Age+ geslacht+EURO2+NYHA+COPD+LVF+Accessroute+CKD , data = df)
summary(res.cox)
#5jaar
res.cox <- coxph(Surv(observatietijd, Overleden==1) ~ IMATclass2+Age+ geslacht+EURO2+NYHA+COPD+LVF+Accessroute+CKD , data = df)
summary(res.cox)

######
#combi
df$combi[df$muscleareaclass== 'low' & df$HUclass=='low'] <- 1
df$combi[] <- 2
df$combi[df$muscleareaclass== 'low' & df$HUclass=='low'] <- 1

df$combi2[df$muscleareaclass== 'low' & df$IMATclass=='low'] <- 1
df$combi2[] <- 2
df$combi2[df$muscleareaclass== 'low' & df$IMATclass=='low'] <- 1

#KM
fit <- survfit(Surv(observatietijd, Overleden) ~ df$combi2, data = df)
fit
summary(fit)$table

p<-ggsurvplot(fit,xlab="Time (Years)",xscale="d_y", size=0.7,axes.offset=FALSE,
           break.time.by=(365.25),pval=F, test.for.trend=FALSE,
           conf.int=FALSE,label.curves=TRUE, ggtheme=theme_minimal()+theme,
           tables.theme=theme_cleantable(), xlim=c(0,1700),expand_limits=c(0,0),
           risk.table=T, risk.table.fontsize=5, censor=F, risk.table.y.text = F, palette = c("red",'blue'), ylim = c(0.4, 1.0))
print(p)
#export => Save as PDF => 6.5x6.5 inch

#log rank
surv_diff <- survdiff(Surv(observatietijd, Overleden) ~ combi, data = df)
surv_diff

df$combi<-factor(df$combi, labels=c('lowlow','rest'))
df$combi<-relevel(df$combi, ref = 'rest')

df$combi2<-factor(df$combi2, labels=c('lowlow','rest'))
df$combi2<-relevel(df$combi2, ref = 'rest')

#combi SMI-SMD
#Crude
#1jaar
res.cox <- coxph(Surv(observatietijd365, overleden_1jaar==1) ~ df$combi, data = df)
summary(res.cox)
#4jaar
res.cox <- coxph(Surv(observatietijd, Overleden==1) ~ combi , data = df)
summary(res.cox)

#Adjusted model        
#1jaar
res.cox <- coxph(Surv(observatietijd365, overleden_1jaar==1) ~ combi+Age+ geslacht+EURO2+NYHA+COPD+LVF+Accessroute+CKD , data = df)
summary(res.cox)
#5jaar
res.cox <- coxph(Surv(observatietijd, Overleden==1) ~ combi+Age+ geslacht+EURO2+NYHA+COPD+LVF+Accessroute+CKD , data = df)
summary(res.cox)

#combi2 SMI-IMAT
#Crude
#1jaar
res.cox <- coxph(Surv(observatietijd365, overleden_1jaar==1) ~ df$combi2, data = df)
summary(res.cox)
#4 jaar
res.cox <- coxph(Surv(observatietijd, Overleden==1) ~ combi2 , data = df)
summary(res.cox)

#Adjusted model        
#1jaar
res.cox <- coxph(Surv(observatietijd365, overleden_1jaar==1) ~ combi2+Age+ geslacht+EURO2+NYHA+COPD+LVF+Accessroute+CKD , data = df)
summary(res.cox)
#4-year
res.cox <- coxph(Surv(observatietijd, Overleden==1) ~ combi2+Age+ geslacht+EURO2+NYHA+COPD+LVF+Accessroute+CKD , data = df)
summary(res.cox)
