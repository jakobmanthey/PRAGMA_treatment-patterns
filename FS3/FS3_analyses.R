# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# PROJECT TITLE:  PRAGMA
# CODE AUTHOR:    JM
# DATE STARTED:   251106

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 0) ESSENTIALS
# ______________________________________________________________________________________________________________________

# clean workspace
rm(list=ls())

# input path
inpath <- paste0("input/")

# output path
outpath <- paste0("FS3/out")

# load libraries
library(tidyverse) 
library( nnet )
library( kableExtra )
library( data.table )
library( ggplot2 )
library( ggthemes )
library( tidyr )
library( stringr )
library( lubridate )
library( here )
library( sjPlot )

here::i_am("2_sample and sequence description.R")
# themes and options
theme_set( theme_gdocs() )
options(scipen = 999)

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

##  this code builds on the LCA to identify treatment pathways to answer FS3:
#   "Unterscheiden sich Morbidität und Mortalität zwischen Patient*innen mit verschiedenen Versorgungswegen?"

#   compare treatment types with regard to inpatient/outpatient and mortality?

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 1) LOAD DATA
# ______________________________________________________________________________________________________________________

##  1) PERSON LEVEL DATA WITH CLASS INFORMATION
# -------------------------------------------------------

filename <- paste0(inpath,"input data_person level.RDS")
person.dat <- readRDS(filename)[order(pragmaid)]
rm(filename)

##  2) ELIXHAUSER DETAILED DATA
# -------------------------------------------------------

filename <- paste0(inpath,"input data_elixhauser detailed data.RDS")
elix.dat <- readRDS(filename)
rm(filename)

##  3) INTERVENTION LEVEL DATA
# -------------------------------------------------------

filename <- paste0(inpath,"input data_intervention level.RDS")
interv.dat <- readRDS(filename)
rm(filename)

##  4) OUTCOME: INFORMATION ON INPATIENT AND OUTPATIENT INFO IN 12M SUBSEQUENT TO 24M TREATMENT PERIOD
# -------------------------------------------------------

# all diagnoses
filename <- paste0(inpath,"1_data_all diagnoses_","2024-05-28",".rds")
diag.dat <- readRDS(filename)
rm(filename)

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 2) PREPARE
# ______________________________________________________________________________________________________________________

##  define inpatient and outpatient outcomes
# -------------------------------------------------------

# merge (only classic inpatient and outpatient)
outc.dat <- merge(person.dat[,.(gkv, pragmaid, start = date.aud, end = date.period.end-1)], 
                  diag.dat[setting %like% "inpatient|outpatient$",.(gkv,pragmaid,setting,case.id,icd,icd_type,icd.alc,date.diag.start,date.diag.end)],
                  by = c("gkv","pragmaid"), all.x = T)

# which interventions?
outc.dat[, table(setting,icd_type)]
outc.dat <- outc.dat[icd_type %in% c("admission","confirmed","primary","secondary")]
nrow(outc.dat) # 2849448

# keep only interventions with start or end date within 24M period
outc.dat <- outc.dat[(date.diag.start %between% list(start,end) |
                       date.diag.end %between% list(start,end))]
outc.dat[, date.diag.start := fifelse(date.diag.start<start,start,date.diag.start)]
outc.dat[, date.diag.end := fifelse(date.diag.end>end,end,date.diag.end)]

# INPATIENT
outc.inpat <- outc.dat[setting == "inpatient"]
nrow(outc.inpat) # 132913

##  for each case, determine whether alc-related:
outc.inpat <- unique(outc.inpat[, .(date.diag.start,
                                    date.diag.end,
                                    alc.related = any(icd.alc), 
                                    icds = paste0(unique(icd),collapse = ";")), by = .(gkv,pragmaid,case.id)])

##  aggregate
outc.inpat.agg <- outc.inpat[, .(
  inpat_n = uniqueN(case.id),
  inpat_los = as.numeric(sum(date.diag.end - date.diag.start + 1)),
  inpat_n.alc = uniqueN(case.id[alc.related == TRUE]),
  inpat_los.alc = as.numeric(sum((date.diag.end - date.diag.start + 1)[alc.related == TRUE])),
  inpat_n.noalc = uniqueN(case.id[alc.related == F]),
  inpat_los.noalc = as.numeric(sum((date.diag.end - date.diag.start + 1)[alc.related == F]))
), by = .(gkv, pragmaid)]

summary(outc.inpat.agg$inpat_los)
nrow(outc.inpat.agg) # 5208 persons with at least one hospitalisation during the entire time

##  correct: number of alcohol hospitalisations and days in intervention (independent variable) should not be part of outcome (dependent variable)
corrdat <- interv.dat[interv.type %like% "qwt|inpat" & !is.na(interv_id),.(pragmaid,start = date.aud, end = date.period.end-1,interv.type,interv_id,date.interv.start,date.interv.end)]
corrdat[, date.interv.start := fifelse(date.interv.start<start,start,date.interv.start)]
corrdat[, date.interv.end := fifelse(date.interv.end>end,end,date.interv.end)]
corrdat <- corrdat[, .(corr_n = uniqueN(interv_id),
                       corr_los = as.numeric(sum(date.interv.end - date.interv.start + 1))), by = pragmaid]

outc.inpat.agg <- merge(outc.inpat.agg,
                        corrdat,
                        by = c("pragmaid"),
                        all.x = T)

outc.inpat.agg[!is.na(corr_n), ':=' (
  inpat_n.alc = ifelse(inpat_n.alc-corr_n<0,0,inpat_n.alc-corr_n),
  inpat_los.alc = ifelse(inpat_los.alc-corr_los<0,0,inpat_los.alc-corr_los)
  )]

outc.inpat.agg[, inpat_n := rowSums(.SD), .SDcols = c("inpat_n.alc", "inpat_n.noalc")]
outc.inpat.agg[, inpat_los := rowSums(.SD), .SDcols = c("inpat_los.alc", "inpat_los.noalc")]
outc.inpat.agg[inpat_n.alc + inpat_n.noalc != inpat_n] # none
outc.inpat.agg[inpat_los.alc + inpat_los.noalc != inpat_los] # none
outc.inpat.agg$corr_n <- outc.inpat.agg$corr_los <- NULL

# OUTPATIENT
outc.outpat <- outc.dat[setting == "outpatient"]
nrow(outc.outpat) # 874616

##  for each case, determine whether alc-related:
outc.outpat <- unique(outc.outpat[, .(date.diag.start,
                                      date.diag.end,
                                      alc.related = any(icd.alc), 
                                      icds = paste0(unique(icd),collapse = ";")), by = .(gkv,pragmaid,case.id)])

##  aggregate
outc.outpat.agg <- outc.outpat[, .(
  outpat_n = uniqueN(case.id),
  outpat_n.alc = uniqueN(case.id[alc.related == TRUE]),
  outpat_n.noalc = uniqueN(case.id[alc.related == F])
  ), by = .(gkv, pragmaid)]
nrow(outc.outpat.agg) # 9471 persons with at least one outpatient contact during the entire time



##  Get analytic data set
# -------------------------------------------------------

person.dat[date.death < date.period.start] # none
person.dat[date.death < date.period.end & date.death > date.period.start]

# combine all relevant data
data <- copy(person.dat[,.(gkv,pragmaid,sex,age,agegroup,emp.type,
                           elix_sum_nomental,
                           date.period.start,date.aud,date.period.end,
                           predclass,class_lab2,
                           died,date.death)])

data[, died.during := ifelse(died ==F, F, date.death < date.period.end & date.death > date.aud)]
data[, table(died.during)]

data <- merge(data,
              outc.inpat.agg,
              by = c("gkv","pragmaid"), all.x = T)
data <- merge(data,
              outc.outpat.agg,
              by = c("gkv","pragmaid"), all.x = T)

data[, ':=' (inpat_n = ifelse(is.na(inpat_n),0,inpat_n),
             inpat_los = ifelse(is.na(inpat_los),0,inpat_los),
             inpat_n.alc = ifelse(is.na(inpat_n.alc),0,inpat_n.alc),
             inpat_los.alc = ifelse(is.na(inpat_los.alc),0,inpat_los.alc),
             inpat_n.noalc = ifelse(is.na(inpat_n.noalc),0,inpat_n.noalc),
             inpat_los.noalc = ifelse(is.na(inpat_los.noalc),0,inpat_los.noalc),
             outpat_n = ifelse(is.na(outpat_n),0,outpat_n),
             outpat_n.alc = ifelse(is.na(outpat_n.alc),0,outpat_n.alc),
             outpat_n.noalc = ifelse(is.na(outpat_n.noalc),0,outpat_n.noalc)
             )]
data[, ':=' (inpat_any = inpat_n > 0,
             inpat_any.noalc = inpat_n.noalc > 0,
             outpat_any = outpat_n > 0,
             outpat_any.noalc = outpat_n.noalc > 0)]

# class labels
data$predclass <- factor(data$predclass)
data$class_lab2 <- factor(data$class_lab2)

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 2) ANALYSES
# ______________________________________________________________________________________________________________________

##  1) Sample description
# -------------------------------------------------------

# sex and age
data[, table(sex)]
data[, table(agegroup)]

# Elixhauser
data[, mean(elix_sum_nomental == 0)]
data[, mean(elix_sum_nomental %in% c(1,2))]
data[, mean(elix_sum_nomental >=3)]


##  2) Inpatient utilization
# -------------------------------------------------------

# main: 
data[, mean(inpat_any), by = predclass][order(predclass)] # lowest for class 0

# secondary 1:
data[, mean(inpat_any.noalc), by = predclass][order(predclass)] # lowest for class 3, 4 + 6

# secondary 2:
data[, mean(inpat_los), by = predclass][order(predclass)] # lowest for class 0 + 5
hist(data$inpat_los, breaks = 100)

# secondary 3:
data[, mean(inpat_los.noalc), by = predclass][order(predclass)] # lowest for class 3, 5 + 6
hist(data$inpat_los.noalc, breaks = 100)

# all together:
tab1 <- data[, {
  main_test <- prop.test(sum(inpat_any), .N, conf.level = 0.95, correct = FALSE)
  sec1_test <- prop.test(sum(inpat_any.noalc), .N, conf.level = 0.95, correct = FALSE)
  .(
    main_mean = mean(inpat_any),
    main_ci_low = main_test$conf.int[1],
    main_ci_high = main_test$conf.int[2],
    secondary1_mean = mean(inpat_any.noalc),
    secondary1_ci_low = sec1_test$conf.int[1],
    secondary1_ci_high = sec1_test$conf.int[2],
    secondary2_mean = mean(inpat_los),
    secondary2_iqr_low = quantile(inpat_los, 0.25),
    secondary2_iqr_high = quantile(inpat_los, 0.75),
    secondary3_mean = mean(inpat_los.noalc),
    secondary3_iqr_low = quantile(inpat_los.noalc, 0.25),
    secondary3_iqr_high = quantile(inpat_los.noalc, 0.75)
  )
}, by = predclass][order(predclass)]

openxlsx::write.xlsx(tab1, paste0("FS3/out/table1_inpatient_descriptives_",Sys.Date(),".xlsx"))

##  3) Outpatient utilization
# -------------------------------------------------------

# main: 
data[, mean(outpat_n), by = predclass][order(predclass)] # highest for brief psych care
hist(data$outpat_n, breaks = 100)

# secondary 1:
data[, mean(outpat_n.noalc), by = predclass][order(predclass)] # highest for brief psych care
hist(data$outpat_n.noalc, breaks = 100)

# all together:
tab4 <- data[, .(
  main_mean = mean(outpat_n), 
  main_iqr_low = quantile(outpat_n, 0.25), 
  main_iqr_high = quantile(outpat_n, 0.75),
  secondary1_mean = mean(outpat_n.noalc), 
  secondary1_iqr_low = quantile(outpat_n.noalc, 0.25), 
  secondary1_iqr_high = quantile(outpat_n.noalc, 0.75)
  ), by = predclass][order(predclass)] # 

openxlsx::write.xlsx(tab4, paste0("FS3/out/table4_outpatient_descriptives_",Sys.Date(),".xlsx"))


##  4) Mortality
# -------------------------------------------------------

data[, mean(died.during), by = predclass][order(predclass)]


##  3) MODELS for INPAT
# -------------------------------------------------------

# main: any inpat
summary(glm(inpat_any ~ sex + agegroup + elix_sum_nomental + emp.type, data, family = "binomial"))
summary(glm(inpat_any ~ sex + agegroup + elix_sum_nomental + emp.type + predclass, data, family = "binomial"))

# secondary 1: noalc inpat
summary(glm(inpat_any.noalc ~ sex + agegroup + elix_sum_nomental + emp.type, data, family = "binomial"))
summary(glm(inpat_any.noalc ~ sex + agegroup + elix_sum_nomental + emp.type + predclass, data, family = "binomial"))

# MAIN: number inpat
poistest <- glm(inpat_n ~ sex + agegroup + elix_sum_nomental + emp.type + predclass, data, family = poisson)
# Calculate dispersion
dispersion <- sum(residuals(poistest, type="pearson")^2) / df.residual(poistest)
# If dispersion > 1.5-2, fit Negative Binomial:
dispersion > 1.5 # T - 3.5
mod.inpat.main <- MASS::glm.nb(inpat_n ~ sex + agegroup + elix_sum_nomental + emp.type + predclass, data)
summary(mod.inpat.main)
rm(poistest, dispersion)


# SEC1: number inpat noalc
poistest <- glm(inpat_n.noalc ~ sex + agegroup + elix_sum_nomental + emp.type + predclass, data, family = poisson)
# Calculate dispersion
dispersion <- sum(residuals(poistest, type="pearson")^2) / df.residual(poistest)
# If dispersion > 1.5-2, fit Negative Binomial:
dispersion > 1.5 # T - 2.8
mod.inpat.sec1 <- MASS::glm.nb(inpat_n.noalc ~ sex + agegroup + elix_sum_nomental + emp.type + predclass, data)
summary(mod.inpat.sec1)
rm(poistest, dispersion)


# SEC2: los inpat - ZINB model
poistest <- glm(inpat_los ~ sex + agegroup + elix_sum_nomental + emp.type + predclass, data, family = poisson)
# Calculate dispersion
dispersion <- sum(residuals(poistest, type="pearson")^2) / df.residual(poistest)
# If dispersion > 1.5-2, fit Negative Binomial:
dispersion > 1.5 # T - 120
mod.inpat.sec2 <- MASS::glm.nb(inpat_los ~ sex + agegroup + elix_sum_nomental + emp.type + predclass, data)
summary(mod.inpat.sec2)

mod.inpat.sec2.zinb <- pscl::zeroinfl(
  inpat_los ~ sex + agegroup + elix_sum_nomental + emp.type + predclass | sex + agegroup + emp.type + predclass,
  data = data,
  dist = "negbin"
)
summary(mod.inpat.sec2.zinb)
rm(poistest, dispersion)


# SEC3: los inpat noalc
poistest <- glm(inpat_los.noalc ~ sex + agegroup + elix_sum_nomental + emp.type + predclass, data, family = poisson)
# Calculate dispersion
dispersion <- sum(residuals(poistest, type="pearson")^2) / df.residual(poistest)
# If dispersion > 1.5-2, fit Negative Binomial:
dispersion > 1.5 # T - 78
mod.inpat.sec3 <- MASS::glm.nb(inpat_los.noalc ~ sex + agegroup + elix_sum_nomental + emp.type + predclass, data)
summary(mod.inpat.sec3)

mod.inpat.sec3.zinb <- pscl::zeroinfl(
  inpat_los.noalc ~ sex + agegroup + elix_sum_nomental + emp.type + predclass | sex + agegroup + emp.type + predclass,
  data = data,
  dist = "negbin"
)
summary(mod.inpat.sec3.zinb)
rm(poistest, dispersion)


#tab_model(mod.inpat.main, mod.inpat.sec1, mod.inpat.sec2, mod.inpat.sec3, 
#          dv.labels = c("Main (N all-cause)",
#                        "Secondary1 (N ohne alc)",
#                        "Secondary2 (LOS)",
#                        "Secondary3 (LOS ohne alc)"),
#         file = paste0("FS3/out/table2_models inpat",Sys.Date(),".html"))

tab_model(mod.inpat.main, mod.inpat.sec1, 
          dv.labels = c("Main (N all-cause)",
                        "Secondary1 (N ohne alc)"),
          file = paste0("FS3/out/table2_models inpat MAIN und SEC1_",Sys.Date(),".html"))
          
tab_model(mod.inpat.sec2.zinb, mod.inpat.sec3.zinb, 
          dv.labels = c("Secondary2 (LOS) - ZINB",
                        "Secondary3 (LOS ohne alc) - ZINB"),
          file = paste0("FS3/out/table3_models inpat SEC2 und SEC3_",Sys.Date(),".html"))






##  4) MODELS for OUTPAT
# -------------------------------------------------------

# main: number outpat all-cause
poistest <- glm(outpat_n ~ sex + agegroup + elix_sum_nomental + emp.type + predclass, data, family = poisson)
# Calculate dispersion
dispersion <- sum(residuals(poistest, type="pearson")^2) / df.residual(poistest)
# If dispersion > 1.5-2, fit Negative Binomial:
dispersion > 1.5 # T - 5.6
mod.outpat.main <- MASS::glm.nb(outpat_n ~ sex + agegroup + elix_sum_nomental + emp.type + predclass, data)
summary(mod.outpat.main)
rm(poistest, dispersion)

# secondary 1: number outpat noalc
poistest <- glm(outpat_n.noalc ~ sex + agegroup + elix_sum_nomental + emp.type + predclass, data, family = poisson)
# Calculate dispersion
dispersion <- sum(residuals(poistest, type="pearson")^2) / df.residual(poistest)
# If dispersion > 1.5-2, fit Negative Binomial:
dispersion > 1.5 # T - 8.0
mod.outpat.sec1 <- MASS::glm.nb(outpat_n.noalc ~ sex + agegroup + elix_sum_nomental + emp.type + predclass, data)
summary(mod.outpat.sec1)
rm(poistest, dispersion)

tab_model(mod.outpat.main, mod.outpat.sec1,
          dv.labels = c("Main (N all-cause)",
                        "Secondary1 (N ohne alc)"),
          file = paste0("FS3/out/table5_models inpat_",Sys.Date(),".html"))



##  5) MODELS for MORTALITÄT
# -------------------------------------------------------

sum(data$died.during)
data[, .(
  n = sum(died.during),
  mean = mean(died.during)
  ), by = predclass][order(predclass)] # 


# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 3) TABLES
# ______________________________________________________________________________________________________________________

##  1) TABLE 1
#   .............................................


# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 4) FIGURES
# ______________________________________________________________________________________________________________________

##  Fig 1) INPAT descriptives N
#   .............................................

pdat <- copy(data[,.SD, .SDcols = names(data)[names(data) %like% "pragmaid|class_lab2|inpat_n$|inpat_n.noalc|inpat_los$|inpat_los.noalc"]])
pdat <- melt(pdat, id.vars = c("pragmaid","class_lab2"))[order(variable)]
pdat[, outcome := dplyr::recode(variable,
  "inpat_n" = "Main (N all-cause)",
  "inpat_n.noalc" = "Secondary1 (N ohne alc)", 
  "inpat_los" = "Secondary2 (LOS)",
  "inpat_los.noalc" = "Secondary3 (LOS ohne alc)")]
pdat$outcome <- factor(pdat$outcome, levels = c("Main (N all-cause)","Secondary1 (N ohne alc)", "Secondary2 (LOS)", "Secondary3 (LOS ohne alc)"))
pdat[, table(variable,outcome)]
pdat[, mean := mean(value), by = .(outcome,class_lab2)]

pdat[, group := ifelse(outcome %like% "LOS", "LOS", "N")]

##  main und sec1
ggplot(pdat[group == "N"], aes(x = class_lab2, y = value, fill = class_lab2)) +
  ggtitle("Anzahl Hospitalisierungen",
          "all-cause (main) vs. ohne Alkoholbezug (secondary1); Dreieck = Mittelwert\nzur besseren Darstellungen werden nur Werte bis 5 dargestellt") +
  facet_wrap(outcome ~ .) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(aes(y = mean), size = 3, shape = 25, fill = "black") +
  scale_fill_viridis_d("") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 3, nrow = 3)) + 
  scale_x_discrete("") + 
  scale_y_continuous("Anzahl Hospitalisierungen") +
  coord_flip(ylim = c(0,5))

ggsave(paste0(outpath,"/Fig 1A_inpat n dist_",Sys.Date(),".png"), width = 12, height = 6)
ggsave(paste0(outpath,"/Fig 1A_inpat n dist_",Sys.Date(),".svg"), width = 12, height = 6)

##  sec2 und sec3
ggplot(pdat[group == "LOS"], aes(x = class_lab2, y = value, fill = class_lab2)) +
  ggtitle("Anzahl Tage im Krankenhaus (length of stay/LOS)",
          "all-cause (secondary2) vs. ohne Alkoholbezug (secondary3); Dreieck = Mittelwert\nzur besseren Darstellungen werden nur Werte bis 50 dargestellt") +
  facet_wrap(outcome ~ .) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(aes(y = mean), size = 3, shape = 25, fill = "black") +
  scale_fill_viridis_d("") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 3, nrow = 3)) + 
  scale_x_discrete("") + 
  scale_y_continuous("Anzahl Hospitalisierungen") +
  coord_flip(ylim = c(0,50))

ggsave(paste0(outpath,"/Fig 1B_inpat los dist_",Sys.Date(),".png"), width = 12, height = 6)
ggsave(paste0(outpath,"/Fig 1B_inpat los dist_",Sys.Date(),".svg"), width = 12, height = 6)

rm(pdat)


##  Fig 2) OUTPAT descriptives N
#   .............................................

pdat <- copy(data[,.SD, .SDcols = names(data)[names(data) %like% "pragmaid|class_lab2|outpat_n$|outpat_n.noalc"]])
pdat <- melt(pdat, id.vars = c("pragmaid","class_lab2"))[order(variable)]
pdat[, outcome := dplyr::recode(variable,
                                "outpat_n" = "Main (N all-cause)",
                                "outpat_n.noalc" = "Secondary1 (N ohne alc)")]
pdat$outcome <- factor(pdat$outcome, levels = c("Main (N all-cause)","Secondary1 (N ohne alc)"))
pdat[, table(variable,outcome)]
pdat[, mean := mean(value), by = .(outcome,class_lab2)]

##  main und sec1
ggplot(pdat, aes(x = class_lab2, y = value, fill = class_lab2)) +
  ggtitle("Anzahl ambulante Kontakte",
          "all-cause (main) vs. ohne Alkoholbezug (secondary1); Dreieck = Mittelwert\nzur besseren Darstellungen werden nur Werte bis 30 dargestellt") +
  facet_wrap(outcome ~ .) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(aes(y = mean), size = 3, shape = 25, fill = "black") +
  scale_fill_viridis_d("") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 3, nrow = 3)) + 
  scale_x_discrete("") + 
  scale_y_continuous("Anzahl ambulante Kontakte") +
  coord_flip(ylim = c(0,30))

ggsave(paste0(outpath,"/Fig 2_outpat n dist_",Sys.Date(),".png"), width = 12, height = 6)
ggsave(paste0(outpath,"/Fig 2_outpat n dist_",Sys.Date(),".svg"), width = 12, height = 6)

rm(pdat)

