library(Rmisc)
library(magrittr)
library(tidyverse)
library(nlme)
library(emmeans)
library(multcomp)
library(car)
library(lme4)
library(predictmeans)
library(purrr)
library(performance)

dat<-read.csv("data/canola_data_2026.csv", na.strings="na")
str(dat)
dat$trial<-as.factor(dat$trial)
dat$site<-as.factor(dat$site)
dat$plot<-as.factor(dat$plot)
dat$block<-as.factor(dat$block)
dat$treatment<-as.factor(dat$treatment)
dat$start_bloomDate<-as.Date(dat$start_bloomDate, format="%m/%d/%y")
dat<-dat %>% 
  mutate(yld=seed_yield_lbA*((100-seed_moisture)/(100-8.5)))
pps = position_dodge(width = .75)

#Outlier test ----
#Run these one-by-one instead of an automated script to visually assess the result
vdat<-dat %>% 
  filter(trial=="v_can_26")
#yield
outlier_x<-check_outliers(lmer(yld~treatment+(1|block), vdat), method="zscore")
cleandat<-vdat[-which(outlier_x),]
#moisture
check_outliers(lmer(Moisture~treatment+(1|block), vdat), method="zscore")
#protein
check_outliers(lmer(Protein.Dry.basis~treatment+(1|block), vdat), method="zscore")
#oil
check_outliers(lmer(Oil.Dry.basis~treatment+(1|block), vdat), method="zscore")
#test weight
check_outliers(lmer(seed_test_weight_lb.bu~treatment+(1|block), vdat), method="zscore") #same row of data as yield outlier

#Variety trial analysis ----
v_table<- cleandat %>% 
  group_by(treatment, treatment_2) %>% 
  summarise("Seed yield" = mean(yld, na.rm=T),
            "Protein" = mean(Protein.Dry.basis, na.rm=T),
            "Oil" = mean(Oil.Dry.basis, na.rm=T),
            "Test weight" = mean(seed_test_weight_lb.bu, na.rm=T),
            "Moisture at harvest" = mean(Moisture, na.rm=T))
v_table #table of means

#Calculate LSDs for each trait, one LSD for TF and one for LL
LSD_table<-data.frame("traits"=c("TF", "LL"),
                      "yield LSD"=NA,
                      "yield CV"=NA,
                      "protein LSD"=NA,
                      "protein CV"=NA,
                      "oil LSD"=NA,
                      "oil CV"=NA,
                      "tw LSD"=NA,
                      "tw CV"=NA)

LSD_table$yield.LSD[1]<-cleandat %>% 
  filter(traits=="TF") %>% 
  mutate(yld_lsd_tf=as.numeric(predictmeans(
    model=lmer(yld~treatment + (1|block), data=.),
    modelterm="treatment",
    pairwise=TRUE,
    adj="none")$LSD[3]))


zzz<-as.numeric(predictmeans(
  model=lmer(yld~treatment + (1|block), data=cleandat),
  modelterm="treatment",
  pairwise=TRUE,
  adj="none")$LSD[3])

mod1<-lmer(yld~treatment + (1|block), data=clean_v_dat)
anova(mod1)
emm<-emmeans(mod1, ~treatment)
summary(emm)
pm<-predictmeans(
  model=mod1,
  modelterm="treatment",
  pairwise=TRUE,
  adj="none")
lsd_val<-pm$LSD


clean_v_dat %>% 
  group_by(treatment, treatment_2) %>% 
  summarise(sd_grain = sd(yld, na.rm=T), 
            n_grain = n(),
            mean_grain = mean(yld, na.rm=T)) %>% #order matters if you want to keep the g_ry column name, must do mean after sd
  mutate(se_grain=sd_grain/sqrt(n_grain)) %>%
  ggplot(aes(x=fct_reorder(treatment, mean_grain), y=mean_grain, color=treatment_2)) + 
  geom_point(size=1.5, position = pps)+ 
  #facet_grid(~fyear) +
  #geom_vline(xintercept = 80, color="black")+
  geom_errorbar(aes(ymin=mean_grain-se_grain, ymax=mean_grain+se_grain), 
                width=0.2, position = pps) +
  xlab("Variety") +
  ylab("Seed yield (lbs/acre)")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.border=element_blank(),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.y=element_text(size=10, color="black"),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x=element_text(size=10,color='black'))




