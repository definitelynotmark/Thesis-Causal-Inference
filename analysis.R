## @knitr analysis_setup
library(tidyverse)
library(magrittr)
library(survival)
library(mice)
library(ahw)
library(data.table)
library(ggfortify)
theme_set(theme_bw())
options(ggplot2.discrete.colour = function(...) scale_colour_brewer(..., palette = "Set2"))
options(ggplot2.discrete.fill = function(...) scale_fill_brewer(..., palette = "Set2"))

## @knitr analysis_data
eso_stent <- read_csv("../Data/data_MH.csv")

# Change some variables
eso_stent <- eso_stent %>%
  mutate(dysphagia_grade = case_when(dysphagia == 0 ~ 0, 
                                     dysphagia_grade < 6 ~ dysphagia_grade - 1),
         dysphagia_treatment_type = ifelse(is.na(dysphagia_treatment_type), 
                                           "None", dysphagia_treatment_type),
         sex = ifelse(sex == 0, "M", "F"),
         scope_passable_tumor = case_when(scope_passable_tumor < 2 ~ scope_passable_tumor),
         n_stage = case_when(n_stage != 4 ~ n_stage),
         m_stage = case_when(m_stage != 2 ~ m_stage))

# Recode discrete variables as factors
eso_stent <- eso_stent %>%
  mutate(across(c(sex, dysphagia, dysphagia_grade, scope_passable_tumor, 
                  histological_type, sys_chemo, radiotherapy, performance_score, 
                  ps_tertiary, asa, t_stage, t_stage_binary, n_stage, m_stage), factor))

## @knitr analysis_missing_data
eso_stent_baseline <- eso_stent %>%
  select(sex, age_diagnosis, asa, histological_type, dysphagia_grade, 
         scope_passable_tumor, cci, ps_tertiary, t_stage, n_stage, m_stage)

eso_stent_timevar <- eso_stent %>%
  select(record_id, dead, follow_up_time, titdd, dysphagia_treatment_type_binary, 
         time_to_chemo, time_to_radiotherapy)

# Imputing missing values with mice package
set.seed(1)
eso_stent_baseline_imp <- mice(eso_stent_baseline, m = 1, seed = 1, 
                               maxit = 5, printFlag = FALSE) %>% complete()
eso_stent_analysis <- cbind(eso_stent_baseline_imp, eso_stent_timevar)

## @knitr analysis_missing_plot
md <- md.pattern(eso_stent_baseline, rotate.names = TRUE)

## @knitr analysis_long
# Create a long_format data.table
eso_stent_analysis_long <- eso_stent_analysis %>%
  tmerge(data2 = eso_stent_analysis, 
         id = record_id, 
         tstop = follow_up_time,
         from.state = tdc(titdd, ifelse(!is.na(dysphagia_treatment_type_binary), 
                                        "Treated", "Not treated"), 
                          "Not treated"),
         to.state = event(titdd, ifelse(!is.na(dysphagia_treatment_type_binary), 
                                        "Treated", "Not treated")),
         to.state = event(follow_up_time, ifelse(dead, "Death", "Censored")),
         intervention = tdc(titdd),
         chemo = tdc(time_to_chemo),
         radio = tdc(time_to_radiotherapy),
         options = list(tstartname = "from", tstopname = "to", idname = "id")) %>%
  data.table() %>%
  relocate(id, from, to, from.state, to.state, intervention, chemo, radio)

# Add noise to event times
set.seed(1) # Errors will ensue when fitting K-M if using wrong seed ......
eso_stent_analysis_long <- eso_stent_analysis_long %>%
  addNoiseAtEventTimes()

## @knitr analysis_treatment_hazards
# Data before treatment is received
eso_stent_analysis_long_treatment <- eso_stent_analysis_long %>% filter(intervention == 0)

# The observed hazard of receiving treatment
fFit <- aalen(Surv(from, to, to.state == "Treated") ~ 1 + sex + age_diagnosis + 
                asa + histological_type + dysphagia_grade + scope_passable_tumor + 
                cci + ps_tertiary + t_stage + n_stage + m_stage + chemo + radio, 
              data = eso_stent_analysis_long_treatment, 
              id = eso_stent_analysis_long_treatment$id)

# The counterfactual/hypothetical hazard of receiving treatment, does not depend on L
cfFit <- aalen(Surv(from, to, to.state == "Treated") ~ 1, 
               data = eso_stent_analysis_long_treatment, id = eso_stent_analysis_long_treatment$id)

## @knitr analysis_weights
# Estimating likelihood ratio process
eso_stent_analysis_long_weights <- makeContWeights(
  fFit, # factual hazard
  cfFit, # hypothetical hazard
  eso_stent_analysis_long, # data
  "Not treated", # at-risk of receiving stent
  "Treated", # receives stent at end of interval
  "from", # name of variable containing start-times
  "to", # name of variable containing stop-times
  "from.state", # name of variable with starting state 
  "to.state", # name of variable with end state
  "id", # name of variable with id
  b = 1, # smoothing bandwidth
  weightRange = c(0.01, 5), # truncates weights outside interval
  willPlotWeights = FALSE
) %>%
  filter(!is.na(weights)) %>%
  relocate(id, from, to, from.state, to.state, intervention, chemo, radio, weights)

## @knitr analysis_KM
KM_MSM_cont <- survfit(Surv(from, to, to.state == "Death") ~ intervention, 
                       data = eso_stent_analysis_long_weights, 
                       weights = eso_stent_analysis_long_weights$weights, 
                       id = id)
KM_naive <- survfit(Surv(from, to, to.state == "Death") ~ intervention, 
                    data = eso_stent_analysis_long, 
                    id = id)

## @knitr analysis_survival_plot
adjusted_data <- KM_MSM_cont %>% 
  fortify() %>% 
  as.data.frame() %>% 
  mutate(Method = "Weighted Kaplan-Meier")
unadjusted_data <- KM_naive %>% 
  fortify() %>% 
  as.data.frame() %>% 
  mutate(Method = "Unweighted Kaplan-Meier")

ggplot(mapping = aes(x = time, y = surv)) +
  geom_step(aes(color = strata, linetype = Method), data = adjusted_data) +
  ggfortify:::geom_confint(aes(ymin = lower, ymax = upper, fill = strata), 
                           data = adjusted_data, alpha = 0.2) +
  geom_step(aes(color = strata, linetype = Method), data = unadjusted_data) +
  labs(fill = "Treatment", linetype = NULL) + xlab("Time since diagnosis (days)") + 
  scale_y_continuous(name = "Probability of survival", labels = scales::percent) +
  scale_color_discrete(name = NULL, labels = c("Untreated", "Treated")) +
  scale_fill_discrete(name = NULL, labels = c("Untreated", "Treated"))

## @knitr analysis_weights_plot
set.seed(1)
eso_stent_analysis_long_weights %>%
  pivot_longer(cols = c(from, to)) %>%
  filter(id %in% sample(601, 20)) %>%
  ggplot(aes(x = value, y = weights, group = id, color = factor(intervention))) +
  geom_step() +
  geom_point(aes(shape = to.state), 
             data = . %>% filter(to.state == "Death" | to.state == "Censored", name == "to")) +
  ylab("Weight") + xlab("Time since diagnosis (days)") + labs(shape = NULL) +
  scale_color_discrete(name = NULL, labels = c("Untreated", "Treated"))
