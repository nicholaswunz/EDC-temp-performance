# Author: Nicholas Wu (nicholas.wu.nz@gmail.com)
# Date: xx/07/2020
# R version: 3.5.1 -- "Feather Spray"
# Paper ID: EDC/temp performance paper
# Description: Data analysis & figure production

## R PACKAGES ##----------------------------------------------------------------------------------
# install and load packages
library(ggplot2)
library(dplyr) 
library(lmerTest)
library(emmeans)
library(tidyr)

# functions
mytheme <- function() {
  theme_bw() + 
    theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.8), # Border around plotting area.
                         panel.grid.major = element_blank(), # Major grid lines blank
                         panel.grid.minor = element_blank(), # Minor grid lines blank
                         axis.line = element_blank(), # axis line size
                         axis.ticks = element_line(colour = "black", size = 0.8),
                         axis.text = element_text(size = 10, colour = "black"), # axis text size
                         axis.title = element_text(size = 10), #axis title size
                         panel.background = element_rect(fill = "transparent"), # bg of the panel
                         plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                         legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                         legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel bg)
} # set up plot theme

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m - (sd(x)/sqrt(length(x)))
  ymax <- m + (sd(x)/sqrt(length(x)))
  n <- length(x)
  return(c(y = m, ymin = ymin, ymax = ymax, num = n))
} # calculate mean +/- SE

#create clean data function

# set directory
setwd('xxxx')

## 1. UCRIT ##-----------------------------------------------------------------------------------
# load dataset
ucrit_data <- read.csv("ucrit_all.csv")

# clean dataset for analysis
clean_ucrit_data <- ucrit_data %>% 
  dplyr::mutate(temp = as.factor(temp),
         temp.acute = as.factor(temp.acute),
         treatment = factor(treatment, levels = c("18-C", "18-expose", "28-C", "28-expose")),
         exposure = factor(exposure, levels = c("control", "expose"))
         )
str(clean_ucrit_data)

## 1a. DATA SUMMARY ##-------------------------------------------------------
ucrit_sum <- clean_ucrit_data %>% 
  dplyr::group_by(bisphenol, temp, exposure, temp.acute) %>% 
  dplyr::summarise(ucrit_mean = mean(ucrit2_post[!is.na(ucrit2_post)]),
                   ucrit_n = length(ucrit2_post[!is.na(ucrit2_post)]),
                   body_length = mean(length2_post[!is.na(length2_post)]))

# mean length and mass
clean_ucrit_data %>% 
  dplyr::group_by(bisphenol, temp, exposure) %>% 
  dplyr::summarise(length_mean = mean(length_post[!is.na(length_post)]),
                   length_sd = sd(length_post[!is.na(length_post)]),
                   mass_mean = mean(body.mass_post[!is.na(body.mass_post)]),
                   mass_sd = sd(body.mass_post[!is.na(body.mass_post)]))

## 1b. ANALYSIS ##---------------------------------------------
# three-way interaction with length as covariate and ID as random effect
ucrit_BPA_data <- subset(clean_ucrit_data, bisphenol == "BPA")
ucrit_BPF_data <- subset(clean_ucrit_data, bisphenol == "BPF")
ucrit_BPS_data <- subset(clean_ucrit_data, bisphenol == "BPS")

ucrit_BPA_model <- lmerTest::lmer(
  ucrit_post ~ ucrit_pre + temp*exposure*temp.acute + length_post + (1|ID), 
  data = ucrit_BPA_data
  )

ucrit_BPF_model <- lmerTest::lmer(
  ucrit_post ~ ucrit_pre + temp*exposure*temp.acute + length_post + (1|ID), 
  data = ucrit_BPF_data
  )

ucrit_BPS_model <- lmerTest::lmer(
  ucrit_post ~ ucrit_pre + temp*exposure*temp.acute + length_post + (1|ID), 
  data = ucrit_BPS_data
  )

anova(ucrit_BPA_model)
anova(ucrit_BPS_model)
anova(ucrit_BPF_model)

## 1c. FIGURES ##----------------------------------------------------
# Extract  estimated means
ucrit_BPA_es <- ucrit_BPA_model %>% 
  emmeans::emmeans(pairwise ~ temp*exposure*temp.acute, adjust = "tukey")
ucrit_BPF_es <- ucrit_BPS_model %>% 
  emmeans::emmeans(pairwise ~ temp*exposure*temp.acute, adjust = "tukey")
ucrit_BPS_es <- ucrit_BPF_model %>% 
  emmeans::emmeans(pairwise ~ temp*exposure*temp.acute, adjust = "tukey")

ucrit_BPA_emmeans <- data.frame(ucrit_BPA_es$emmeans, bisphenol = "BPA")
ucrit_BPF_emmeans <- data.frame(ucrit_BPF_es$emmeans, bisphenol = "BPF")
ucrit_BPS_emmeans <- data.frame(ucrit_BPS_es$emmeans, bisphenol = "BPS")

# unite all bisphenols
emmeans_all <- ucrit_BPA_emmeans %>% 
  dplyr::union(ucrit_BPF_emmeans) %>% 
  dplyr::union(ucrit_BPS_emmeans) %>% 
  tidyr::unite("temp.treat", bisphenol,temp.acute, remove = FALSE) # for distinct colours by bisphenols

# convert m s to BL s from ucrit sum 
emmeans_BL <- emmeans_all %>% 
  dplyr::mutate(emmean2 = emmean / ucrit_sum$body_length,
                lower.CL2 = lower.CL / ucrit_sum$body_length,
                upper.CL2 = upper.CL / ucrit_sum$body_length)

# Plot Fig 1 - Plot estimated marginal means while accounting for pre-Ucrit, BL, and random effect ID
emmeans_BL %>% 
  ggplot(aes(x = temp.acute, y = emmean2, colour = temp.treat, shape = exposure)) +
  geom_point(data = clean_ucrit_data, aes(x = temp.acute, y = ucrit2_post, group = exposure), 
             size = 3, colour = "grey", position = position_dodge(0.5), alpha = 0.3) +
  geom_errorbar(aes(ymin = lower.CL2, ymax = upper.CL2), size = 0.8, width = 0.1, 
                position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  xlab("Acute test temperature (°C)") +
  ylab(expression(italic("U")["crit"]~"(BL s"^"-1"*")")) + 
  scale_color_manual(values = c("#879FDB", "#002F70",       # blue
                                "#CB6CA2", "#6D1C68",       # purple
                                "#81CA9F", "#14505C")) +    # green
  mytheme() + facet_grid(bisphenol ~ temp)


## 2. ENZYME ASSAY ##---------------------------------------------------------------------------------
# load dataset
enzyme_data <- read.csv("enzyme_all.csv")

clean_enzyme_data <- enzyme_data %>% 
  dplyr::mutate(temp = as.factor(temp),
                temp.acute = as.factor(temp.acute),
                treatment = factor(treatment, levels = c("18-C", "18-expose", "28-C", "28-expose")),
                exposure = factor(exposure, levels = c("control", "expose"))
  )

## 2a DATA SUMMARY ##----------------------------------------------------------------------------------------
enzyme_sum <- clean_enzyme_data %>% 
  dplyr::group_by(bisphenol, temp, exposure, temp.acute) %>% 
  dplyr::summarise(CS_mean = mean(CS.activity[!is.na(CS.activity)]),
                   CS_SD = sd(CS.activity[!is.na(CS.activity)]),
                   CS_n = length(CS.activity[!is.na(CS.activity)]),
                   LDH_mean = mean(LDH.activity[!is.na(LDH.activity)]),
                   LDH_SD = sd(LDH.activity[!is.na(LDH.activity)]),
                   LDH_n = length(LDH.activity[!is.na(LDH.activity)]))

## 2b ANALYSIS ##---------------------------------------------
# three-way interaction with ID as random effect
enzyme_BPA_data <- subset(clean_enzyme_data, bisphenol == "BPA")
enzyme_BPF_data <- subset(clean_enzyme_data, bisphenol == "BPF")
enzyme_BPS_data <- subset(clean_enzyme_data, bisphenol == "BPS")

# citrate synthase
CS_BPA_model <- lmerTest::lmer(
  CS.activity ~ temp*exposure*temp.acute + (1|ID), 
  data = enzyme_BPA_data
  )

CS_BPF_model <- lmerTest::lmer(
  CS.activity ~ temp*exposure*temp.acute + (1|ID), 
  data = enzyme_BPF_data
  )

CS_BPS_model <- lmerTest::lmer(
  CS.activity ~ temp*exposure*temp.acute + (1|ID), 
  data = enzyme_BPS_data
  )

anova(CS_BPA_model)
anova(CS_BPF_model)
anova(CS_BPS_model)

# Lactate dehydrogenase
LDH_BPA_model <- lmerTest::lmer(
  LDH.activity ~ temp*exposure*temp.acute + (1|ID), 
  data = enzyme_BPA_data
  )

LDH_BPF_model <- lmerTest::lmer(
  LDH.activity ~ temp*exposure*temp.acute + (1|ID), 
  data = enzyme_BPF_data
  )

LDH_BPS_model <- lmerTest::lmer(
  LDH.activity ~ temp*exposure*temp.acute + (1|ID), 
  data = enzyme_BPS_data
  )

anova(LDH_BPA_model)
anova(LDH_BPF_model)
anova(LDH_BPS_model)

## 2c FIGURES ##----------------------------------------------------------------------------------------
# Citrate synthase
# Extract estimated marginal means + 95% CI
CS_BPA_es <- CS_BPA_model %>% 
  emmeans::emmeans(pairwise ~ temp*exposure*temp.acute, adjust = "tukey")
CS_BPF_es <- CS_BPF_model %>% 
  emmeans::emmeans(pairwise ~ temp*exposure*temp.acute, adjust = "tukey")
CS_BPS_es <- CS_BPS_model %>% 
  emmeans::emmeans(pairwise ~ temp*exposure*temp.acute, adjust = "tukey")

CS_BPA_emmeans <- data.frame(CS_BPA_es$emmeans, bisphenol = "BPA")
CS_BPF_emmeans <- data.frame(CS_BPF_es$emmeans, bisphenol = "BPF")
CS_BPS_emmeans <- data.frame(CS_BPS_es$emmeans, bisphenol = "BPS")

# unite all bisphenols
CS_emmeans_all <- CS_BPA_emmeans %>% 
  dplyr::union(CS_BPF_emmeans) %>% 
  dplyr::union(CS_BPS_emmeans) %>% 
  tidyr::unite("temp.treat", bisphenol,temp.acute, remove = FALSE) # for distinct colours by bisphenols

# Plot Fig 2 - Plot estimated marginal means while accounting for random effect ID
CS_emmeans_all %>% 
  ggplot(aes(x = temp.acute, y = emmean, colour = temp.treat, shape = exposure)) +
  geom_point(data = clean_enzyme_data, aes(x = temp.acute, y = CS.activity, group = exposure), 
             size = 3, colour = "grey", position = position_dodge(0.5), alpha = 0.3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 0.8, width = 0.1, 
                position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  xlab("Acute test temperature (°C)") +
  ylab(expression("CS activity (µmol min"^"-1"*" g"^"-1"*")")) +
  scale_color_manual(values = c("#879FDB", "#002F70",       # blue
                                "#CB6CA2", "#6D1C68",       # purple
                                "#81CA9F", "#14505C")) +    # green
  mytheme() + facet_grid(bisphenol ~ temp)

# Lactate dehydrogenase figure
# Extract estimated marginal means + 95% CI
LDH_BPA_es <- LDH_BPA_model %>% 
  emmeans::emmeans(pairwise ~ temp*exposure*temp.acute, adjust = "tukey")
LDH_BPF_es <- LDH_BPF_model %>% 
  emmeans::emmeans(pairwise ~ temp*exposure*temp.acute, adjust = "tukey")
LDH_BPS_es <- LDH_BPS_model %>% 
  emmeans::emmeans(pairwise ~ temp*exposure*temp.acute, adjust = "tukey")

LDH_BPA_emmeans <- data.frame(LDH_BPA_es$emmeans, bisphenol = "BPA")
LDH_BPF_emmeans <- data.frame(LDH_BPF_es$emmeans, bisphenol = "BPF")
LDH_BPS_emmeans <- data.frame(LDH_BPS_es$emmeans, bisphenol = "BPS")

# unite all bisphenols
LDH_emmeans_all <- LDH_BPA_emmeans %>% 
  dplyr::union(LDH_BPF_emmeans) %>% 
  dplyr::union(LDH_BPS_emmeans) %>% 
  tidyr::unite("temp.treat", bisphenol,temp.acute, remove = FALSE) # for distinct colours by bisphenols

# Plot Fig 3 - Plot estimated marginal means while accounting for random effect ID
LDH_emmeans_all %>% 
  ggplot(aes(x = temp.acute, y = emmean, colour = temp.treat, shape = exposure)) +
  geom_point(data = clean_enzyme_data, aes(x = temp.acute, y = LDH.activity, group = exposure), 
             size = 3, colour = "grey", position = position_dodge(0.5), alpha = 0.3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), size = 0.8, width = 0.1, 
                position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  xlab("Acute test temperature (°C)") +
  ylab(expression("LDH activity (µmol min"^"-1"*" g"^"-1"*")")) +
  scale_color_manual(values = c("#879FDB", "#002F70",       # blue
                                "#CB6CA2", "#6D1C68",       # purple
                                "#81CA9F", "#14505C")) +    # green
  mytheme() + facet_grid(bisphenol ~ temp)



## RAW DATA MEAN +/- SE ##--------------------------------------------------------
#ucrit
clean_ucrit_data %>% 
  ggplot(aes(x = temp.acute, y = ucrit2_post, colour = temp, shape = treatment)) +
  geom_point(size = 3, colour = "grey", position = position_dodge(0.5), alpha = 0.3) +
  stat_summary(fun.data = data_summary, geom = "errorbar", 
               position = position_dodge(0.5), size = 1, width = 0.2) +
  stat_summary(fun.y = mean, geom = "point", 
               position = position_dodge(0.5), size = 3) +
  xlab("Test temperature (°C)") +
  ylab(expression(italic("U")["crit"]~"(BL s"^"-1"*")")) + 
  scale_shape_manual(values = c(16, 16, 17, 17)) +
  scale_color_manual(values = c("#023FA5", "#8E063B")) +
  mytheme() + facet_grid(bisphenol ~ temp)


clean_enzyme_data

# CS activity
clean_enzyme_data %>% 
  ggplot(aes(x = temp.acute, y = CS.activity, colour = temp, shape = treatment)) +
  geom_point(size = 3, colour = "grey", position = position_dodge(0.5), alpha = 0.3) +
  stat_summary(fun.data = data_summary, geom = "errorbar", 
               position = position_dodge(0.5), size = 1, width = 0.2) +
  stat_summary(fun.y = mean, geom = "point", 
               position = position_dodge(0.5), size = 3) +
  xlab("Test temperature (°C)") +
  ylab(expression("CS activity (µmol min"^"-1"*" g"^"-1"*")")) +
  scale_shape_manual(values = c(16, 16, 17, 17)) +
  scale_color_manual(values = c("#023FA5", "#8E063B")) +
  mytheme() + facet_grid(bisphenol ~ temp)


# LDH activity
clean_enzyme_data %>% 
  ggplot(aes(x = temp.acute, y = LDH.activity, colour = temp, shape = treatment)) +
  geom_point(size = 3, colour = "grey", position = position_dodge(0.5), alpha = 0.3) +
  stat_summary(fun.data = data_summary, geom = "errorbar", 
               position = position_dodge(0.5), size = 1, width = 0.2) +
  stat_summary(fun.y = mean, geom = "point", 
               position = position_dodge(0.5), size = 3) +
  xlab("Test temperature (°C)") +
  ylab(expression("LDH activity (µmol min"^"-1"*" g"^"-1"*")")) +
  scale_shape_manual(values = c(16, 16, 17, 17)) +
  scale_color_manual(values = c("#023FA5", "#8E063B")) +
  mytheme() + facet_grid(bisphenol ~ temp)

## SUPPLEMENTARY INFO ##-------------------------------------------------------------
# Ucrit % change
clean_ucrit_data <- clean_ucrit_data %>%
  unite("temp.treat", bisphenol,temp.acute, remove = FALSE)

clean_ucrit_data %>%
  ggplot(aes(x = temp.acute, y = ucrit_change, colour = temp.treat, shape = exposure)) +
  geom_jitter(size = 3, position = position_dodge(0.5), colour = "grey", alpha = 0.3) +
  stat_summary(fun.data = data_summary, geom = "errorbar", position = position_dodge(0.5), size = 1, width = 0.1) +
  stat_summary(fun.y = mean, geom = "point", position = position_dodge(0.5), size = 3) +
  xlab("Acute test temperature (°C)") +
  ylab(expression("Percentage change in"~italic("U")["crit"]~"(%)")) + 
  scale_color_manual(values = c("#879FDB", "#002F70",       # blue
                                "#CB6CA2", "#6D1C68",       # purple
                                "#81CA9F", "#14505C")) +    # green
  geom_hline(yintercept = 0, linetype = "dashed") + 
  mytheme() + facet_grid(bisphenol ~ temp)


 