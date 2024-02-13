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
table(df$status)
# 2. Time-to-event (we can use either rtime or dtime)
df$rtime %>% head()
# In years
df <- df %>% mutate(rtime_yrs = rtime/365.25)

# # If you need to calculate time-to-event from dates
# df2 <- data.frame(surgery_date = sample(seq(as.Date('1999/01/01'), as.Date('2000/01/01'), by = "day"), 12),
#                   day_of_last_followup = sample(seq(as.Date('2000/01/01'), as.Date('2020/01/01'), by = "day"), 12))
# head(df2)
# class(df2$surgery_date)
df2 <- df2 %>% mutate(os_yrs = as.numeric(as.duration(day_of_last_followup - surgery_date), 'years'),
                      os_months = as.numeric(as.duration(day_of_last_followup - surgery_date), 'months'),
                      os_days = as.numeric(as.duration(day_of_last_followup - surgery_date), 'days'))
df2


# Create a survival object
surv_obj <- Surv(time = df$rtime_yrs, event = df$status)
head(surv_obj) # a vector of follow-up time, with “+” to represent if an observation was right-censored. 
# Create survival curve
s1 <- survfit(surv_obj ~ 1, data = df)
summary(s1)

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
km_curve <- km_curve + 
  coord_cartesian(xlim = c(0, 8))
km_curve

## Estimating x-year survival ----------------------
# What is the probability of surviving beyond a certain number of years, x?
summary(s1, times = 1)$surv
summary(s1, times = c(0, 0.5, 2, 5, 7)) # we can specify whatever times we like
km_curve +
  geom_vline(xintercept = 1, linetype = 'dashed', colour = 'red', size = 1) +
  geom_hline(yintercept = summary(s1, times = 5)$surv, linetype = 'dashed', colour = 'red', size = 1) + 
  coord_cartesian(xlim = c(0, 8))

# Ignoring censoring leads to an overestimate of the overall survival probability!
df$rtime_yrs_nc <- ifelse(df$status == 0, max(df$rtime_yrs), df$rtime_yrs)
surv_obj_nocensoring <- Surv(df$rtime_yrs_nc, df$status)
s1_nocensoring <- survfit(surv_obj_nocensoring ~ 1, data = df)
km_curve_nc <- ggsurvfit(s1_nocensoring, linewidth = 1) +
  labs(x = "Years",
       y = "Overall survival probability") +
  add_confidence_interval() +
  add_risktable() +
  #add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() + 
  coord_cartesian(xlim = c(0, 8)) +
  biostatsquid_theme 

# Create a plot with both survival curves
plot(s1_nocensoring, col = "blue", lty = 1, main = "Survival Curves with and without Censoring")
lines(s1, col = "red", lty = 2, add = TRUE)
# add a legend to the plot
legend("topright", legend = c("without censoring", "with censoring"), 
       lty = c(1, 1), col = c('blue', 'red'), cex = .85, bty = "n")

# Estimating median survival rtime --------------------------------------------
s1 
# With no censoring we get an overestimate of the median survival time
s1_nocensoring # in this case we cannot even calculate it because the survival probability of 0.5 was not reached.

# Or do it with survminer
# https://rpkgs.datanovia.com/survminer/#ggsurvplot-drawing-survival-curves

# Log rank test --------------------------------------------

unique(df$chemo) # if they had chemotherapy or not
# Create the new survfit object based on chemo
s2 <- survfit(surv_obj  ~ chemo, data = df)
# This time let's use the ggsurvplot function
km_curve <- ggsurvplot(s2, data = df, # survminer functions require that you specify the survival object and again specify the data used to fit the survival object. Remember to do this to avoid non-specific error messages. 
                       size = 1,                 # change line size
                       palette = c("#E7B800", "#2E9FDF"),# custom color palettes
                       censor.shape = "|", censor.size = 4,
                       conf.int = TRUE,          # Add confidence interval
                       pval = TRUE,              # Add p-value
                       risk.table = TRUE,        # Add risk table
                       risk.table.col = "strata",# Risk table color by groups
                       legend.labs = list('0' = "No chemo", '1' = "Chemo"),    # Change legend labels
                       risk.table.height = 0.25, # Useful to change when you have multiple groups
                       ggtheme = theme_bw()      # Change ggplot2 theme
)
km_curve

# Now we can compute the test of the difference between the survival curves using survdiff
# Comparing survival rtimes between groups aka Log rank test
logrankres_chemo <- survdiff(Surv(rtime_yrs, status) ~ chemo, data = df)
logrankres_chemo
# p > 0.05 - non significant difference between getting chemo or not.

# What about hormonal treatment?
logrankres_hormon <- survdiff(Surv(rtime_yrs, status) ~ hormon, data = df)
logrankres_hormon
# p < 0.05, seems like the difference in survival probabilities between patients who received hormonal data or not is significant!

# Cox regression --------------------------------------------
# Fit the model
cox_res <- coxph(Surv(rtime_yrs, status) ~ hormon + chemo + size + er + pgr + nodes + meno + grade + age, data = df)
cox_res
table(df$hormon)
unique(df$size) # < 20 group is the reference group for size 

# If we have a binary categorical column
df$treatment <- rep(c('drugA', 'placebo'), nrow(df)/2)
head(df)
unique(df$treatment)
cox_res <- coxph(Surv(rtime_yrs, status) ~ hormon + treatment, data = df)
cox_res
# If we wanted placebo to be the reference group:
df$treatment <- factor(df$treatment, levels = c('drugA', 'placebo'))
df$treatment <- relevel(df$treatment, ref = 'placebo')
unique(df$treatment)
cox_res <- coxph(Surv(rtime_yrs, status) ~ hormon + treatment, data = df)
cox_res
df$treatment <- NULL

# Get a summary of the model - we get the CI of the estimated HR and the different test scores
cox_res <- coxph(Surv(rtime_yrs, status) ~ hormon + chemo + size + er + pgr + nodes + meno + grade + age, data = df)
cox_res
# Before interpreting the results, verify whether the proportional hazards assumptions are respected
# Cox model assumes that the HR between any two groups remains constant over time
test <- survival::cox.zph(cox_res) # this tests for proportional hazards based on Schoenfeld residuals
test 
# Interpretation: 
# - p-val < 0.05 --> evidence against proportional hazards assumption, HRs are not constant over time
# - chisq ++ --> stronger violation of assumptions

# Plot the Schoenfeld residuals over time for each covariate
survminer::ggcoxzph(test, point.size = 0.1) #  if the residuals show a clear pattern over time, it may indicate a violation of the proportional hazards assumption.
# Interpretation: 
# - No Pattern (Constant Residuals): If the residuals appear randomly scattered around zero, 
# with no clear trend or pattern, it suggests that the proportional hazards assumption is reasonable:)
# - Linear Trend: A linear trend (increasing or decreasing) in the residuals over time 
# might suggest a violation of the proportional hazards assumption. 
# For example, if the residuals are consistently positive or negative across time, 
# it indicates a time-dependent effect.
# - Nonlinear Pattern: If the residuals exhibit a non-linear pattern or specific shapes 
# (e.g., U-shape, V-shape), it may indicate deviations from proportional hazards.
# Parallelism: In addition to looking for patterns, you should check for parallelism. 
# Parallelism means that the spread and distribution of residuals are relatively 
# constant across time. If the residuals widen or narrow over time, it may suggest a violation of the assumption.

# Get a summary of our model
summary(cox_res)
# Hormonal treatment: HR = 0.79 -> if we hold all the other covariates constant, 
# patients who get hormonal treatment have a decrease in the risk of death of 21% 
# Size: positive association between the size of the tumour and mortality: 2.33x if it's > 50 mm

# Forest plots ================================================================
# Visualise your Cox model results
forest_plot <- ggforest(cox_res, data = df) 
forest_plot


# --------------------------------------------------------------------------









# Extract and format the global p-value
summary(cox_res)[["logtest"]][["pvalue"]]
formatted_p_val <- format(summary(cox_res)[["logtest"]][["pvalue"]], scientific = TRUE)
forest_plot + 
  scale_y_continuous(labels = scales::scientific_format(scale = 1e-3, decimal.mark = ".", big.mark = ",")) +
  theme_minimal()


#https://epirhandbook.com/en/survival-analysis.html
# https://rpkgs.datanovia.com/survminer/#ggsurvplot-drawing-survival-curves 
# this to edit

# https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html#The_df_dataset
# https://yuzar-blog.netlify.app/posts/2021-01-03-survival-analysis-1-a-gentle-introduction-into-kaplan-meier-curves/
# https://epirhandbook.com/en/survival-analysis.html