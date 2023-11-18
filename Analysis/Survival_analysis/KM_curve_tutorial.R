# ---------------------- #
# KM_curve_tutorial.R
# ---------------------- #
# Author: Laura Twomey
# Date of creation: 12.10.2023
# Version: 1.0

# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Set seed
set.seed(42)

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(survival)
library(ggsurvfit)
library(survminer)
#library(ranger)
#library(ggfortify)

# Theme ====================================================================
biostatsquid_theme <- theme(plot.title = element_text(size = rel(2)),
                            panel.grid.major.y = element_line(colour = 'gray'),
                            panel.grid.minor.y = element_line(colour = 'gray'),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            plot.background = element_rect(fill = NULL, colour = 'white'),
                            panel.background = element_rect(fill = 'white'),
                            # Axis stuff
                            axis.line = element_line(colour = 'black', linewidth = 1),
                            axis.text = element_text(colour = "black", face = 'bold'),
                            axis.text.x = element_text(size = rel(1)),
                            axis.text.y = element_text(size = rel(1)),
                            axis.title = element_text(size = rel(1.2)),
                            axis.ticks = element_line(colour = 'black', linewidth = 1.2),
                            # Legend stuff
                            legend.position = "bottom",
                            legend.margin = margin(6, 6, 6, 6),
                            legend.title = element_text(face = 'bold'),
                            legend.background = element_blank(),
                            legend.box.background = element_rect(colour = "black"))

# Load the data ============================================================
df <- survival::rotterdam
head(df)
colnames(df)

# Key columns for survival analysis
# 1. Censoring status: 1 = event happened, 0 = censored (or TRUE and FALSE)
table(df$death)
df <- df %>% mutate(status = death)

# 2. Time-to-event (we can use either rtime or dtime)
df$rtime %>% head()
# In years
df <- df %>% mutate(rtime_yrs = rtime/365.25)

# # If you need to calculate time-to-event from dates
# df2 <- data.frame(surgery_date = sample(seq(as.Date('1999/01/01'), as.Date('2000/01/01'), by = "day"), 12),
#                   day_of_last_followup = sample(seq(as.Date('2000/01/01'), as.Date('2020/01/01'), by = "day"), 12))
# head(df2)
# class(df2$surgery_date)
# df2 <- df2 %>% mutate(os_yrs = as.numeric(as.duration(day_of_last_followup - surgery_date), 'years'),
#                       os_months = as.numeric(as.duration(day_of_last_followup - surgery_date), 'months'),
#                       os_days = as.numeric(as.duration(day_of_last_followup - surgery_date), 'days'))
# df2
# 

# Create a survival object
surv_obj <- Surv(df$rtime_yrs, df$status)
head(surv_obj)
# Create survival curve
s1 <- survfit(surv_obj ~ 1, data = df)

# Kaplan-Meier plots ======================================================
## Plot -------------------------
km_curve <- ggsurvfit(s1, linewidth = 1) +
  labs(x = "Years",
       y = "Overall survival probability") +
  add_confidence_interval() +
  add_risktable() +
  #add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  biostatsquid_theme

# Plot customisation - it's ggplot!
# limit plot to show 8 years and less
km_curve + 
  coord_cartesian(xlim = c(0, 8))

## Estimating x-year survival ----------------------
# What is the probability of surviving beyond a certain number of years, x?
summary(s1, times = 1)$surv
summary(s1, times = c(0, 0.5, 2, 5, 7)) # we can specify whatever times we like
km_curve +
  geom_vline(xintercept = 1, linetype = 'dashed', colour = 'red', size = 1) +
  geom_hline(yintercept = summary(s1, times = 5)$surv, linetype = 'dashed', colour = 'red', size = 1)

# Ignoring censoring leads to an overestimate of the overall survival probability!
df$nocensoring <- 1
surv_obj_nocensoring <- Surv(df$rtime_yrs, df$nocensoring)
s1_nocensoring <- survfit(surv_obj_nocensoring ~ 1, data = df)
km_curve_nc <- ggsurvfit(s1_nocensoring, linewidth = 1) +
  labs(x = "Years",
       y = "Overall survival probability") +
  add_confidence_interval() +
  add_risktable() +
  #add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  biostatsquid_theme

km_curve_nc + km_curve
summary(s1_nocensoring, times = 1)$surv

# Estimating median survival rtime --------------------------------------------
survfit(Surv(rtime, status) ~ 1, data = df)

# What happens if you use a “naive” estimate? Here “naive” means that you exclude the censored patients from the calculation entirely to estimate median survival rtime among the patients who have had the event.
# 
# Summarize the median survival rtime among the 165 patients who died:
#   
#   df %>% 
#   filter(status == 1) %>% 
#   summarize(median_surv = median(rtime))
# 
# Ignoring censoring will lead to an underestimate of median survival rtime because the follow-up rtime that censored patients contribute is excluded (blue line). The true survival curve accounting for censoring in the df data is shown in yellow for comparison.


survfit(Surv(rtime, status) ~ 1, data = df) %>% 
  tbl_survfit(
    probs = 0.5,
    label_header = "**Median survival (95% CI)**"
  )


# Or do it with survminer
# https://rpkgs.datanovia.com/survminer/#ggsurvplot-drawing-survival-curves
km_curve <- ggsurvplot(s1, data = df,
                       size = 1,                 # change line size
                       palette = c("#E7B800", "#2E9FDF"),# custom color palettes
                       censor.shape="|", censor.size = 4,
                       conf.int = TRUE,          # Add confidence interval
                       pval = TRUE,              # Add p-value
                       risk.table = TRUE,        # Add risk table
                       risk.table.col = "strata",# Risk table color by groups
                       legend.labs =
                         c("Male", "Female"),    # Change legend labels
                       risk.table.height = 0.25, # Useful to change when you have multiple groups
                       ggtheme = theme_bw()      # Change ggplot2 theme
)

# Comparing survival rtimes between groups
# aka Log rank test
survdiff(Surv(rtime, status) ~ sex, data = df)

# Cox regression
coxph(Surv(rtime, status) ~ sex, data = df)

coxph(Surv(rtime, status) ~ sex, data = df) %>% 
  tbl_regression(exp = TRUE) 



# https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html#The_df_dataset
# https://yuzar-blog.netlify.app/posts/2021-01-03-survival-analysis-1-a-gentle-introduction-into-kaplan-meier-curves/