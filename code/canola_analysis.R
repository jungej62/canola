library(readxl)
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

#dat<-read.csv("data/canola_data_2026.csv", na.strings="na")
dat<-read_xlsx("data/canola_data_2026.xlsx", na="na") #I started with a csv that was extracted from the field book, but then manually edited it and saved it as an excel so I could make a metadata tab. I should have edited this in the google sheet or using code, but alas... next time.
str(dat)
dat$trial<-as.factor(dat$trial)
dat$site<-as.factor(dat$site)
dat$plot<-as.factor(dat$plot)
dat$block<-as.factor(dat$block)
dat$treatment<-as.factor(dat$treatment)
#dat$start_bloomDate<-as.Date(dat$start_bloomDate, format="%m/%d/%y")
dat<-dat %>% 
  mutate(yld=seed_yield_lb_a*((100-seed_moisture)/(100-8.5)),
         plant_counts=(plant_count_1+plant_count_2)/2)
pps = position_dodge(width = .75)

#Outlier test ----
#Here are the response variables we'll analyze
response_vars<-set_names(c("ESV", "groundcover", "bloom_start", "height", "lodging", "shattering", "test_weight_lb_a", "protein_dry", "oil_dry", "yld", "plant_counts"))
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

#Help with data anlysis workflow
response_vars<-set_names(c("ESV", "groundcover", "height", "lodging", "test_weight_lb_a", "protein_dry", "oil_dry", "yld"))
treatment_col<-"treatment"
block_col<-"block"

results <- response_vars %>% 
  map(function(resp) {
    
    # 1. Estimate Population Mean and Population CV (%)
    pop_mean <- mean(cleandat[[resp]], na.rm = TRUE)
    pop_sd   <- sd(cleandat[[resp]], na.rm = TRUE)
    pop_cv   <- (pop_sd / pop_mean) * 100
    
    # 2. Fit Mixed-Effects Model: resp ~ treatment + (1 | block)
    formula_obj <- as.formula(paste(resp, "~", treatment_col, "+ (1 |", block_col, ")"))
    model_fit   <- lmer(formula_obj, data = cleandat)
    
    # 3. Compute Estimated Marginal Means per treatment level
    emm <- emmeans(model_fit, specs = treatment_col)
    
    trt_means <- emm %>% 
      as_tibble() %>% 
      dplyr::select(all_of(treatment_col), emmean, SE)
    
    # 4. Extract Pairwise Contrasts to Calculate Exact 5% LSD Value
    # Fishers LSD uses adjust = "none"
    pairwise_diffs <- contrast(emm, method = "pairwise", adjust = "none") %>% 
      as_tibble()
    
    # Calculate LSD = t_crit * SE_diff
    # (Since design is balanced, SE_diff and df are identical across pairs)
    se_diff <- pairwise_diffs$SE[1]
    df_res  <- pairwise_diffs$df[1]
    
    lsd_5pct <- qt(1 - 0.05 / 2, df = df_res) * se_diff
    
    # Return structured results
    list(
      trt_means = trt_means,
      pop_mean  = pop_mean,
      pop_cv    = pop_cv,
      lsd_5pct  = lsd_5pct,
      df        = df_res
    )
  })

#Table creation
wide_trt_means <- map_dfr(results, ~ .x$trt_means, .id = "Response") %>% 
  dplyr::select(Response, treatment, emmean) %>% 
  pivot_wider(names_from = Response, values_from = emmean)

# Extract summary metrics (Pop Mean, Pop CV, LSD) into matching format
summary_rows <- tibble(
  treatment = c("Mean", "CV (%)", "LSD (0.05)"),
  ESV = c(results$ESV$pop_mean, results$ESV$pop_cv, results$ESV$lsd_5pct),
  groundcover = c(results$groundcover$pop_mean, results$groundcover$pop_cv, results$groundcover$lsd_5pct),
  height = c(results$height$pop_mean, results$height$pop_cv, results$height$lsd_5pct),
  lodging = c(results$lodging$pop_mean, results$lodging$pop_cv, results$lodging$lsd_5pct),
  test_weight_lb_a = c(results$test_weight_lb_a$pop_mean, results$test_weight_lb_a$pop_cv, results$test_weight_lb_a$lsd_5pct),
  protein_dry = c(results$protein_dry$pop_mean, results$protein_dry$pop_cv, results$protein_dry$lsd_5pct),
  oil_dry = c(results$oil_dry$pop_mean, results$oil_dry$pop_cv, results$oil_dry$lsd_5pct),
  yld = c(results$yld$pop_mean, results$yld$pop_cv, results$yld$lsd_5pct)
)


response_vars<-set_names(c("ESV", "groundcover", "height", "lodging", "test_weight_lb_a", "protein_dry", "oil_dry", "yld"))
# Combine into single final presentation table

final_summary_table <- bind_rows(wide_trt_means, summary_rows)
final_summary_table$YieldPercent<-final_summary_table$yld/results$yld$pop_mean

write.csv(final_summary_table, file="data/variety_table2.csv")

#Jake's cumbersome method of variety trial analysis ----
v_table<- cleandat %>% 
  group_by(treatment, treatment_2, traits) %>% 
  summarise("Seed yield" = mean(yld, na.rm=T),
            "Protein" = mean(Protein.Dry.basis, na.rm=T),
            "Oil" = mean(Oil.Dry.basis, na.rm=T),
            "Test weight" = mean(seed_test_weight_lb.bu, na.rm=T),
            "Moisture at harvest" = mean(Moisture, na.rm=T))
v_table #table of means
write.csv(v_table, file="data/variety_table.csv")

#Calculate LSDs for each trait, one LSD for TF and one for LL
# LSD_table<-data.frame("traits"=c("TF", "LL"),
#                       "yield LSD"=NA,
#                       "yield CV"=NA,
#                       "protein LSD"=NA,
#                       "protein CV"=NA,
#                       "oil LSD"=NA,
#                       "oil CV"=NA,
#                       "tw LSD"=NA,
#                       "tw CV"=NA)

as.numeric(predictmeans(
    model=lmer(seed_test_weight_lb.bu~treatment + (1|block),
               data=cleandat),
#            data=subset(cleandat, traits=="TF")),
    modelterm="treatment",
    pairwise=TRUE,
    adj="none")$LSD[3])

stats_out<-cleandat %>% 
  mutate(mean_yld = mean(yld, na.rm=T),
         CV_yld = sd(yld, na.rm=T)/mean(yld, na.rm=T),
         mean_protein = mean(Protein.Dry.basis, na.rm=T),
         CV_protein = sd(Protein.Dry.basis, na.rm=T)/mean(Protein.Dry.basis, na.rm=T),
         mean_oil = mean(Oil.Dry.basis, na.rm=T),
         CV_oil = sd(Oil.Dry.basis, na.rm=T)/mean(Oil.Dry.basis, na.rm=T),
         mean_tw = mean(seed_test_weight_lb.bu, na.rm=T),
         CV_tw = sd(seed_test_weight_lb.bu, na.rm=T)/mean(seed_test_weight_lb.bu, na.rm=T)) %>% 
  dplyr::select(mean_yld, CV_yld, mean_protein, CV_protein, mean_oil, CV_oil, mean_tw, CV_tw) %>%
  slice_head(n = 1)
write.csv(stats_out, file="data/stats_table.csv")


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




