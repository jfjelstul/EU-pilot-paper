###########################################################################
# Joshua C. Fjelstul
# Replication code for:
# Can International Institutions Help States to Comply with International Law?
# Encouraging Evidence from the European Union's Pilot Program
# by Sivaram Cheruvu and Joshua C. Fjelstul
###########################################################################

# libraries
library(lubridate)
library(stringr)
library(sandwich)
library(lmtest)
library(rworldmap)
library(dplyr)
library(tidyr)
library(stargazer)
library(ggplot2)
library(gridExtra)

# set working directory
setwd("~/Documents/EU-pilot-paper")

###########################################################################
###########################################################################
# data preparation
###########################################################################
###########################################################################

# read in data
cases <- read.csv("data/infringement-cases-duration.csv", stringsAsFactors = FALSE)

# filter
cases <- filter(cases, case_year <= 2019)

##################################################
# resolution
##################################################

# first stage
cases$resolved_LFN258_stage <- as.numeric(cases$stage_LFN258 == 1 & cases$stage_RO258 == 0)

# second stage
cases$resolved_RO258_stage <- as.numeric(cases$stage_RO258 == 1 & cases$stage_RF258 == 0)

# pre-litigation phase (first or second stage)
cases$resolved_prelitigation <- as.numeric(cases$stage_RF258 == 0)

# litigation phase
cases$resolved_postlitigation <- as.numeric(cases$stage_RF258 == 1)

#################################################
# workload
#################################################

# workload
cases$workload <- NA
for(i in 1:nrow(cases)) {
  cases$workload[i] <- sum(cases$date_LFN258 <= cases$date_LFN258[i] & cases$date_closing >= cases$date_LFN258[i] & cases$member_state_code == cases$member_state_code[i] & cases$directorate_general_code == cases$directorate_general_code[i])
}

#################################################
# pilot data
#################################################

# participants
participants <- c("Austria", "Czech Republic", "Denmark", "Germany", "Finland", "Hungary", "Ireland", "Italy", "Lithuania", "Netherlands", "Portugal", "Slovenia", "Sweden", "Spain", "United Kingdom")

# date range
pilot <- filter(cases, date_LFN258 > ymd("2004-05-01") & date_closing < ymd("2010-03-01"))

# treatment variable
pilot$treated <- as.numeric(pilot$member_state %in% participants)

# pre/post variable
pilot$post_treatment <- as.numeric(pilot$date_LFN258 > ymd("2008-04-15"))

# group variable (for fixed effects)
pilot$group <- as.factor(str_c(pilot$member_state_code, "-", pilot$directorate_general_code))

###########################################################################
###########################################################################
# analysis
###########################################################################
###########################################################################

##################################################
# function to run models
##################################################

run <- function(f, dat) {
  mod <- lm(f, dat)
  vcov <- vcovHC(mod, type = "HC0")
  table <- coeftest(mod, vcov = vcov)
  se <- table[,2]
  p <- table[,4]
  out <- list(mod = mod, table = table, se = se, p = p)
  return(out)
} 

##################################################
# DiD analysis
##################################################

# model 1
f <- duration_LFN258 ~ treated * post_treatment
mod1 <- run(f, pilot)
# mod1$table

# model 2
f <- duration_LFN258 ~ treated * post_treatment + workload
mod2 <- run(f, pilot)
# mod2$table

# model 3
f <- log(duration_LFN258) ~ treated * post_treatment + workload + factor(member_state_code) + factor(directorate_general_code)
mod3 <- run(f, pilot)
# mod3$table

# model 4
f <- duration_LFN258 ~ treated * post_treatment + factor(group)
mod4 <- run(f, pilot)
# mod4$table

# model 5
f <- duration_LFN258 ~ treated * post_treatment + workload + factor(member_state_code) + factor(directorate_general_code) + factor(case_year)
mod5 <- run(f, pilot)
# mod5$table

# save
stargazer(mod1$mod, mod2$mod, mod3$mod, mod4$mod, mod5$mod,
          out = "tables/pilot-DiD.tex", type = "latex", style = "default", digits = 3, no.space = TRUE, 
          align = TRUE, omit = "factor",
          order = c("^treated$", "^post_treatment$", "^treated:post_treatment$", "^workload_subject$"),
          se = list(mod1$se, mod2$se, mod3$se, mod4$se, mod5$se),
          p = list(mod1$p, mod2$p, mod3$p, mod4$p, mod5$p),
          covariate.labels = c("\\textsc{pilot}", "\\textsc{period}", "\\textsc{pilot $\\times$ period}", "\\textsc{workload}", "\\textit{Constant}"),
          keep.stat = c("n", "rsq"))

##################################################
# DiDiD analysis
##################################################

# model 1
f <- duration_LFN258 ~ treated * post_treatment * type_nonconformity
mod1 <- run(f, pilot)
# mod1$table

# model 2
f <- duration_LFN258 ~ treated * post_treatment * type_nonconformity + workload
mod2 <- run(f, pilot)
# mod2$table

# model 3
f <- duration_LFN258 ~ treated * post_treatment * type_nonconformity + workload + factor(member_state_code) + factor(directorate_general_code)
mod3 <- run(f, pilot)
# mod3$table

# model 4
f <- duration_LFN258 ~ treated * post_treatment * type_nonconformity + factor(group)
mod4 <- run(f, pilot)
# mod4$table

# model 5
f <- duration_LFN258 ~ treated * post_treatment * type_nonconformity + workload + factor(member_state_code) + factor(directorate_general_code) + factor(case_year)
mod5 <- run(f, pilot)
# mod5$table

# save
stargazer(mod1$mod, mod2$mod, mod3$mod, mod4$mod, mod5$mod,
          out = "tables/pilot-DiDiD.tex", type = "latex", style = "default", digits = 3, no.space = TRUE, 
          align = TRUE, omit = "factor",
          order = c("^treated$", "^post_treatment$", "^incorrect$", "^treated:post_treatment$", "^treated:incorrect$", "^post_treatment:incorrect$", "^treated:post_treatment:incorrect$", "^workload_subject$"),
          se = list(mod1$se, mod2$se, mod3$se, mod4$se, mod5$se),
          p = list(mod1$p, mod2$p, mod3$p, mod4$p, mod5$p),
          covariate.labels = c("\\textsc{pilot}", "\\textsc{period}","\\textsc{non-conformity}", "\\textsc{pilot $\\times$ period}","\\textsc{pilot $\\times$ non-conformity}", "\\textsc{period $\\times$ non-conformity}", "\\textsc{pilot $\\times$ period $\\times$ non-conformity}", "\\textsc{workload}", "\\textit{Constant}"),
          keep.stat = c("n", "rsq"))

##################################################
# causal mechanism
##################################################

# model 1
f <- resolved_LFN258_stage ~ treated * post_treatment * type_nonconformity + factor(member_state_code) + factor(directorate_general_code)
mod1 <- run(f, pilot)
# mod1$table

# model 2
f <- resolved_RO258_stage ~ treated * post_treatment * type_nonconformity + factor(member_state_code) + factor(directorate_general_code)
mod2 <- run(f, pilot)
# mod2$table

# model 3
f <- resolved_prelitigation ~ treated * post_treatment * type_nonconformity + factor(member_state_code) + factor(directorate_general_code)
mod3 <- run(f, pilot)
# mod3$table

# save
stargazer(mod1$mod, mod2$mod, mod3$mod,
          out = "tables/causal-mechanism.tex", type = "latex", style = "default", digits = 3, no.space = TRUE, 
          align = TRUE, omit = "factor",
          order = c("^treated$", "^post_treatment$", "^incorrect$", "^treated:post_treatment$", "^treated:incorrect$", "^post_treatment:incorrect$", "^treated:post_treatment:incorrect$"),
          se = list(mod1$se, mod2$se, mod3$se),
          p = list(mod1$p, mod2$p, mod3$p),
          covariate.labels = c("\\textsc{pilot}", "\\textsc{period}","\\textsc{non-conformity}", "\\textsc{pilot $\\times$ period}","\\textsc{pilot $\\times$ non-conformity}", "\\textsc{period $\\times$ non-conformity}", "\\textsc{pilot $\\times$ period $\\times$ non-conformity}", "\\textit{Constant}"),
          keep.stat = c("n", "rsq"))

##################################################
# by hand
##################################################

# DiD estimate
# first subscript is treated vs not treated
# second subscript is pre-treatment vs post-treatment
y00 <- mean(pilot$duration_LFN258[pilot$treated == 0 & pilot$post_treatment == 0])
y10 <- mean(pilot$duration_LFN258[pilot$treated == 1 & pilot$post_treatment == 0])
y01 <- mean(pilot$duration_LFN258[pilot$treated == 0 & pilot$post_treatment == 1])
y11 <- mean(pilot$duration_LFN258[pilot$treated == 1 & pilot$post_treatment == 1])
# first difference is pre-treatment vs post-treatment
# second difference is treated vs not treated
DiD <- (y11 - y10) - (y01 - y00)
# -37.77206

# DiD estimate by regression
f <- duration_LFN258 ~ treated * post_treatment
mod <- run(f, pilot)
DiD <- mod$table[4,1]
# -37.77206

# DiDiD estimate
# first subscript is treated vs not treated
# second subscript is pre-treatment vs post-treatment
# third subscript is nonconformity vs noncommunication
y000 <- mean(pilot$duration_LFN258[pilot$treated == 0 & pilot$post_treatment == 0 & pilot$type_nonconformity == 0])
y001 <- mean(pilot$duration_LFN258[pilot$treated == 0 & pilot$post_treatment == 0 & pilot$type_nonconformity == 1])
y100 <- mean(pilot$duration_LFN258[pilot$treated == 1 & pilot$post_treatment == 0 & pilot$type_nonconformity == 0])
y101 <- mean(pilot$duration_LFN258[pilot$treated == 1 & pilot$post_treatment == 0 & pilot$type_nonconformity == 1])
y010 <- mean(pilot$duration_LFN258[pilot$treated == 0 & pilot$post_treatment == 1 & pilot$type_nonconformity == 0])
y011 <- mean(pilot$duration_LFN258[pilot$treated == 0 & pilot$post_treatment == 1 & pilot$type_nonconformity == 1])
y110 <- mean(pilot$duration_LFN258[pilot$treated == 1 & pilot$post_treatment == 1 & pilot$type_nonconformity == 0])
y111 <- mean(pilot$duration_LFN258[pilot$treated == 1 & pilot$post_treatment == 1 & pilot$type_nonconformity == 1])
# first difference is pre-treatment vs post-treatment
# second difference is treated vs not treated
# third difference is nonconformity vs noncommunication
DiDiD <- (y111 - y101) - (y011 - y001) - ((y110 - y100) - (y010 - y000))
# -72.3254

# DiD estimate by regression
f <- duration_LFN258 ~ treated * post_treatment * type_nonconformity
mod <- run(f, pilot)
DiDiD <- mod$table[8,1]
# -72.3254

###########################################################################
###########################################################################
# figures
###########################################################################
###########################################################################

#################################################
# plot theme functions
#################################################

custom_theme <- function (title = 10, text = 8, margin.b = 15, angle = 0) {
  out <- theme_minimal()
  out <- out + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(color = "black", hjust = 0.5, margin = margin(t = 0, r = 0, b = 7, l = 0), size = title),
    axis.title.y = element_text(color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0), size = title),
    axis.title.x = element_text(color = "black", margin = margin(t = 10, r = 0, b = 0, l = 0), size = title),
    axis.text.x = element_text(color = "black", margin = margin(t = 7, r = 0, b = 0, l = 0), size = text),
    axis.text.y = element_text(color = "black", margin = margin(t = 0, r = 7, b = 0, l = 0), size = text),
    strip.background = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5, color = "black"),
    axis.ticks.length = unit(5, "pt"),
    strip.text = element_text(color = "black", size = title),
    plot.margin = margin(t = 5, r = 5, b = margin.b, l = 5)
  )
  return(out)
}

custom_titles <- function(x = NULL, y = NULL, main = NULL) {
  list(xlab(x), ylab(y), ggtitle(main))
}

#################################################
# backlog
#################################################

# make data frame
backlog <- data.frame(date =  seq.Date(from = ymd("2005-01-01"), to = ymd("2019-12-30"), by = "day"))

# read in data
master <- read.csv("data/infringement-cases.csv", stringsAsFactors = FALSE)

# earliest and latest dates
master$earliest_date <- str_extract(master$case_history, "[0-9]{4}-[0-9]{2}-[0-9]{2}")
master$latest_date <- ifelse(master$stage_closing == 1 | master$stage_withdrawal == 1, str_extract(master$case_history, "[0-9]{4}-[0-9]{2}-[0-9]{2}$"), "2019-12-30")
master$earliest_date <- ymd(master$earliest_date)
master$latest_date <- ymd(master$latest_date)

# count cases
for(i in 1:nrow(backlog)) {
  backlog$count[i] <- sum(backlog$date[i] > master$earliest_date & backlog$date[i] < master$latest_date)
}

# breaks
breaks <-  c(ymd("2005-01-01"), ymd("2006-01-01"), ymd("2007-01-01"),
             ymd("2008-01-01"), ymd("2009-01-01"), ymd("2010-01-01"),
             ymd("2011-01-01"), ymd("2012-01-01"), ymd("2013-01-01"),
             ymd("2014-01-01"), ymd("2015-01-01"), ymd("2016-01-01"),
             ymd("2017-01-01"), ymd("2018-01-01"), ymd("2019-01-01"),
             ymd("2020-01-01"))

# make plot
plot <- ggplot() +
  geom_line(data = backlog, aes(x = date, y = count), color = "black", size = 0.5) +
  scale_x_date(breaks = breaks, labels = 2005:2020) +
  scale_y_continuous(limits = c(1200, 3200), breaks = seq(1200, 3200, 200)) +
  custom_theme() +
  custom_titles(x = "Date", y = "Backlog", main = "Backlog of Infringement Cases")

# save plot
pdf(file = "plots/backlog.pdf", width = 8, height = 6)
plot
dev.off()

#################################################
# maps
#################################################

# map
world_map <- fortify(spTransform(getMap(resolution = "low"), CRS("+proj=wintri")))

# pilot
world_map$pilot <- as.numeric(world_map$id %in% participants)
EU <- unique(cases$member_state[cases$date_LFN258 < ymd("2008-04-15")])
world_map$EU <- as.numeric(world_map$id %in% EU)
world_map$pilot[world_map$EU == 0] <- 2
world_map$pilot <- factor(world_map$pilot, levels = c(0, 1, 2), labels = c("Not an Original Participant", "Original Participant", "Not an EU Member State"))

# make map
plot <- ggplot() +
  geom_map(data = world_map, map = world_map, mapping = aes(map_id = id, fill = pilot), color = "black", size = 0.2) +
  scale_fill_manual(name = NULL, values = c("gray50", "gray30", "gray80")) +
  coord_equal() +
  ylim(3900000, 7700000) +
  xlim(-950000, 2500000) +
  custom_theme() +
  theme(panel.border = element_rect(fill = NA, size = 1),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.key = element_rect(fill = NA, size = 0.5))

# save plot
pdf(file = "plots/pilot-map.pdf", width = 8, height = 6)
plot
dev.off()

#################################################
# DiD plot
#################################################

# calculate means
# first subscript is treated vs not treated
# second subscript is pre-treatment vs post-treatment
y00 <- mean(pilot$duration_LFN258[pilot$treated == 0 & pilot$post_treatment == 0])
y10 <- mean(pilot$duration_LFN258[pilot$treated == 1 & pilot$post_treatment == 0])
y01 <- mean(pilot$duration_LFN258[pilot$treated == 0 & pilot$post_treatment == 1])
y11 <- mean(pilot$duration_LFN258[pilot$treated == 1 & pilot$post_treatment == 1])

# treatment effect
DiD <- (y11 - y10) - (y01 - y00)

# make plot
plot <- ggplot() +
  
  # participant line
  geom_segment(aes(x = 0, y = y10, xend = 1, yend = y11), linetype = "solid") +
  geom_point(aes(0, y10)) +
  geom_point(aes(1, y11)) +
  
  # non-participant line
  geom_segment(aes(x = 0, y = y00, xend = 1, yend = y01), linetype = "solid") +
  geom_point(aes(0, y00)) +
  geom_point(aes(1, y01)) +
  
  # counterfactual
  geom_segment(aes(x = 0, y = y10, xend = 1, yend = y11 - DiD), linetype = "dashed") +
  geom_point(aes(0, y10)) +
  geom_point(aes(1, y11 - DiD)) +
  
  # treatment size
  geom_segment(aes(x = 1 + 0.05, y = y11, xend = 1 + 0.05, yend = y11 - DiD), linetype = "solid") +

  # treatment size segment caps
  geom_segment(aes(x = 1 + 0.05 - 0.015, y = y11, xend = 1 + 0.05 + 0.015, yend = y11), linetype = "solid") +
  geom_segment(aes(x = 1 + 0.05 - 0.015, y = y11 - DiD, xend = 1 + 0.05 + 0.015, yend = y11 - DiD), linetype = "solid") +
  
  # labels
  annotate(geom = "text", label = "Counterfactual", x = 0.5, y = 250, hjust = 0, color = "black", size = 3.5) +
  annotate(geom = "text", label = "Participant", x = 0.15, y = 252, hjust = 1, color = "black", size = 3.5) +
  annotate(geom = "text", label = "Non-Participant", x = 0.15, y = 222, hjust = 1, color = "black", size = 3.5) +
  annotate(geom = "text", label = "Effect", x = 1.015, y = 206, hjust = 1, color = "black", size = 3.5) +
  
  # scale
  scale_x_continuous(expand = c(0.2, 0.1), breaks = c(0, 1), labels = c("Pre-Treatment", "Post-Treatment")) +
  scale_y_continuous(limits = c(175, 280)) +
  
  # labels
  xlab(NULL) +
  ylab("Duration of LFN Stage (Days)") +
  
  # theme
  custom_theme() +
  custom_titles(y = "Duration of LFN Stage")

# save plot
pdf(file = "plots/pilot-DiD.pdf", width = 8, height = 6)
plot
dev.off()

#################################################
# end R script
#################################################
