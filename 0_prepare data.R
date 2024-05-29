# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# PROJECT TITLE:  PRAGMA
# CODE AUTHOR:    JM
# DATE STARTED:   240429

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
outpath <- paste0("output/")

# load libraries
library( data.table )
library( ggplot2 )
library( ggthemes )
library( tidyr )
library( stringr )
library( lubridate )

# themes and options
theme_set( theme_gdocs() )
options(scipen = 999)

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 1) LOAD DATA
# ______________________________________________________________________________________________________________________

DATE <- "2024-05-28"

##  1) PRAGMA-ID
# -------------------------------------------------------

filename <- paste0(inpath,"0_pragma_id_GKV with Stammdata_",DATE,".rds")
id.dat <- readRDS(filename)[order(pragmaid)]
rm(filename)

##  2) INTERVENTIONS
# -------------------------------------------------------

filename <- paste0(inpath,"4_data_AUD diagnosis and intervention sequences_",DATE,".rds")
interv.dat <- readRDS(filename)
interv.dat <- interv.dat[period == "24m"]
rm(filename)

##  3) DIAGNOSES
# -------------------------------------------------------

filename <- paste0(inpath,"1_data_all diagnoses_",DATE,".rds")
diag.dat <- readRDS(filename)
aud.dat <- diag.dat[icd.alc == T]
rm(filename)

##  4) EMPLOYMENT
# -------------------------------------------------------

filename <- paste0(inpath,"1_data_employment periods_",DATE,".rds")
emp.dat <- readRDS(filename)
rm(filename)

##  5) INCOME
# -------------------------------------------------------

filename <- paste0(inpath,"1_data_income data_",DATE,".rds")
inc.dat <- readRDS(filename)
rm(filename)

##  6) SEQUENCE GROUPINGS
# -------------------------------------------------------

filename <- paste0(outpath,"1_sequence_data_24m_monthly_lca_results_",DATE,".csv")
class.dat <- data.table(read.csv(filename))
rm(filename)

##  7) INTERVENTIONS LONG
# -------------------------------------------------------

filename <- paste0("output/treatment_series_at_resolution_cache/treatment_series_24m_quarter_followup_1_aud_only_intervs_",DATE,".csv")
interv.long.dat <- data.table(read.csv(filename))
rm(filename)

##  X) Temporary data
# -------------------------------------------------------

filename <- paste0(inpath,"1_data_insurance periods_",DATE,".rds")
temp.dat <- readRDS(filename)
rm(filename)

temp.dat[pragmaid == "i704ZVu6sF"]
rm(temp.dat)



# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 2) Prepare data
# ______________________________________________________________________________________________________________________

##  1) SOCIODEMOGRAPHICS AND EMPLOYMENT
#   .............................................

data <- copy(id.dat[,.(pragmaid,gkv,sex,yob,nationality,died,date.death)])

# dimensions
nrow(data) # 24313
length(unique(data$pragmaid)) # 24313

# sex and nationality
table(data$sex)
table(data$nationality)
data$nationality <- factor(data$nationality, levels = c("deutsch","nicht deutsch","unbekannt"))

# only relevant sample
data <- data[pragmaid %in% unique(interv.dat$pragmaid)]

# add employment during 12m before diagnosis
add <- copy(emp.dat)

  ##  keep only periods starting/ending in relevant period
  add <- merge(unique(interv.dat[,.(pragmaid,lookbehind.start = date.period.start,lookbehind.end = date.aud)]),
               add,
               by = "pragmaid", all.x = T)
  
  ##  calculate days of period between start and end to select relevant periods
  add[, days_overlap := pmin(lookbehind.end, date.emp.end) - pmax(lookbehind.start, date.emp.start) +1]
  add <- add[days_overlap >=0]
  
  ##  determine dominant period
  add[, days_all := sum(days_overlap), by = pragmaid]
  add[, table(days_all)] # should be at around 365 days --> 313-367 = ok
  add[, days_prop := as.numeric(days_overlap) / as.numeric(days_all)]
  add[, select := max(days_prop), by = pragmaid]
  add[, select := days_prop == select]
  add <- add[select == T,.(pragmaid,emp.type,days_prop)]
  
  
    ### more than 1 period: employed > unemployed > retired > other
    select <- add[, .N, by = pragmaid][N>1]$pragmaid
    add[pragmaid %in% select, check := ifelse(any(emp.type == "employed"), "employed",
                                              ifelse(any(emp.type == "unemployed"), "unemployed",
                                                     ifelse(any(emp.type == "retired"), "retired", "other"))), 
        by = pragmaid]
    add <- add[is.na(check) | check == emp.type,.(pragmaid,emp.type)]
    
data <- merge(data,add,by = "pragmaid")
rm(select,add)

# emp.type
data$emp.type <- factor(data$emp.type, 
                        levels = c("employed","unemployed","retired","other"))

# dimensions
nrow(data) # 9491
length(unique(data$pragmaid)) # 9491


##  2) ADD COMORBIDITY
#   .............................................

comorb.dat <- diag.dat[pragmaid %in% unique(data$pragmaid)]

##  keep all (or only certain?) diagnoses
comorb.dat[, table(setting, icd_type)]
comorb.dat <- comorb.dat[icd_type %in% c("confirmed","primary","any")]
nrow(comorb.dat) # 2700719

##  3 periods: 1) baseline // 2) 1st year post AUD // 3) 2nd year post AUD 

sel.dat.period1 <- unique(interv.dat[,.(pragmaid,
                                        start = date.period.start,
                                        end = date.aud)])
sel.dat.period2 <- unique(interv.dat[,.(pragmaid,
                                        start = date.aud,
                                        end = date.aud+years(1))])
sel.dat.period3 <- unique(interv.dat[,.(pragmaid,
                                        start = date.aud+years(1),
                                        end = date.period.end)])

### period 1:
comorb.period1 <- merge(comorb.dat,
                        sel.dat.period1, 
                        by = c("pragmaid"), all.x = T)
comorb.period1 <- unique(comorb.period1[date.diag.start %between% list(start,end),
                                        .(pragmaid,icd)])
nrow(comorb.period1) # 199086

### period 2:
comorb.period2 <- merge(comorb.dat,
                        sel.dat.period2, 
                        by = c("pragmaid"), all.x = T)
comorb.period2 <- unique(comorb.period2[date.diag.start %between% list(start,end),
                                 .(pragmaid,icd)])
nrow(comorb.period2) # 209187

### period 3:
comorb.period3 <- merge(comorb.dat,
                        sel.dat.period3, 
                        by = c("pragmaid"), all.x = T)
comorb.period3 <- unique(comorb.period3[date.diag.start %between% list(start,end),
                                 .(pragmaid,icd)])
nrow(comorb.period3) # 200495

rm(sel.dat.period1, sel.dat.period2, sel.dat.period3)

##  any non-matches?
comorb.period1[is.na(pragmaid)]
comorb.period2[is.na(pragmaid)]
comorb.period3[is.na(pragmaid)]

##  get comorbidity indices
### period 1
elix.dat.period1 <- comorbidity::comorbidity(x = comorb.period1,
                                             id = "pragmaid",
                                             code = "icd",
                                             map = "elixhauser_icd10_quan",
                                             assign0 = T,
                                             tidy.codes = T)

names(elix.dat.period1)[!names(elix.dat.period1) %like% "pragmaid"] <- paste0("elix_",names(elix.dat.period1)[!names(elix.dat.period1) %like% "pragmaid"])

### period 2
elix.dat.period2 <- comorbidity::comorbidity(x = comorb.period2,
                                             id = "pragmaid",
                                             code = "icd",
                                             map = "elixhauser_icd10_quan",
                                             assign0 = T,
                                             tidy.codes = T)

names(elix.dat.period2)[!names(elix.dat.period2) %like% "pragmaid"] <- paste0("elix_",names(elix.dat.period2)[!names(elix.dat.period2) %like% "pragmaid"])

### period 3
elix.dat.period3 <- comorbidity::comorbidity(x = comorb.period3,
                                             id = "pragmaid",
                                             code = "icd",
                                             map = "elixhauser_icd10_quan",
                                             assign0 = T,
                                             tidy.codes = T)

names(elix.dat.period3)[!names(elix.dat.period3) %like% "pragmaid"] <- paste0("elix_",names(elix.dat.period3)[!names(elix.dat.period3) %like% "pragmaid"])

##  get sum score
elix.dat.period1$period <- 1
elix.dat.period2$period <- 2
elix.dat.period3$period <- 3
elix.dat <- rbind(elix.dat.period1,
                  elix.dat.period2,
                  elix.dat.period3)

elix.dat_sum <- data.table(elix.dat)

  ### A) all categories
  elix.dat_sum$elix_sum_all <- rowSums(elix.dat_sum[,.SD, .SDcols = names(elix.dat_sum)[names(elix.dat_sum) %like% "^elix"]])
  ### B) all categories but not alcohol, drug, psycho, depre
  elix.dat_sum$elix_sum_nomental <- rowSums(elix.dat_sum[,.SD, .SDcols = names(elix.dat_sum)[names(elix.dat_sum) %like% "^elix" &
                                                                              !names(elix.dat_sum) %like% "alcohol$|drug$|psycho$|depre$" & 
                                                                              !names(elix.dat_sum) %like% "^elix_sum"]])

  elix.dat_sum <- dcast(elix.dat_sum[,.(pragmaid,period,elix_sum_all,elix_sum_nomental)], 
                      pragmaid ~ period, value.var = c("elix_sum_all","elix_sum_nomental"))
  
data <- merge(data,
              elix.dat_sum, 
              by = "pragmaid", all.x = T)

# dimensions
nrow(data) # 9491
length(unique(data$pragmaid)) # 9491

# missing elix data:
data[is.na(elix_sum_nomental_1)] # none
data[is.na(elix_sum_nomental_2)] # none
data[is.na(elix_sum_nomental_3)] # 255 --> presumably 0 diagnoses in that period...

  ##  example:
  comorb.period3[pragmaid == elix.dat_sum[elix_sum_nomental_3 == 0, pragmaid][2]] ## -> 0 assigned only if any diagnoses are present

data[is.na(elix_sum_nomental_3), table(died)] # 21 died
data[is.na(elix_sum_nomental_3), elix_sum_nomental_3 := 0]

##  3) ADD OTHER RELEVANT TIME-INVARYING VARIABLES
#   .............................................

add <- unique(interv.dat[,.(pragmaid,date.aud,date.period.start,date.period.end,interv.any)])

data <- merge(data,add,by = "pragmaid")
rm(add)

# dimensions
nrow(data) # 9491
length(unique(data$pragmaid)) # 9491

# any deaths ? --> remove before LCA!
data[date.death %between% list(date.period.start,date.period.end),.(gkv,date.death,date.period.start,date.period.end)]


##  4) DEFINE AGEGROUP
#   .............................................

data[, age := year(date.aud) - yob]
data[, summary(age)]
data[is.na(age)] ## none

# define agegroups
data[, agegroup := ifelse(age <= 24, "18-24",
                          ifelse(age <= 34, "25-34",
                                 ifelse(age <= 44, "35-44",
                                        ifelse(age <= 54, "45-54",
                                               ifelse(age <= 64, "55-64",
                                                      ifelse(age <= 74, "65-74", "75-96"))))))]
data[, prop.table(table(agegroup))]
data[, table(age, agegroup)]
data$agegroup <- factor(data$agegroup)

# dimensions
nrow(data) # 9491
length(unique(data$pragmaid)) # 9491

##  5) ADD INCOME FOR EMPLOYED
#   .............................................

add <- copy(inc.dat[pragmaid %in% data[emp.type == "employed"]$pragmaid])

# get relevant year:
  ##  keep only periods starting/ending in relevant period
  add <- merge(unique(interv.dat[,.(pragmaid,
                                    lookbehind.start = date.period.start,
                                    lookbehind.start.year = year(date.period.start),
                                    lookbehind.end = date.aud,
                                    lookbehind.end.year = year(date.aud))]),
               add,
               by = "pragmaid", all.y = T)[order(pragmaid,income.year)]
  
  
  add <- add[income.year == lookbehind.start.year | income.year == lookbehind.end.year]
  add[, year_n := 1:.N, by = pragmaid]
  add[, table(year_n)] # 1 or 2 only
  add[]
  
  ##  calculate weights for >1 year
  add[, year1 := as.numeric(as.Date(paste0(year(lookbehind.start),"-12-31"))- lookbehind.start)]
  add[, year2 := as.numeric(lookbehind.end - as.Date(paste0(year(lookbehind.end),"-01-01"))) ]
  add[, total := year1 + year2]
  add[, year1 := year1/total]
  add[, year2 := year2/total]
  add <- add[,.(pragmaid,income,wt = ifelse(year_n == 1,year1,year2))]
  add <- add[,.(income = Hmisc::wtd.mean(income,wt, na.rm = T)), by = pragmaid]
  
  
data <- merge(data,add,by = "pragmaid",all.x = T)
rm(add)

data[, table(emp.type, is.na(income))] # missing for 909 employed people
data[gkv == "aok", table(emp.type, is.na(income))] # missing for 55 employed people
data[gkv == "dak", table(emp.type, is.na(income))] # missing for 854 employed people

# dimensions
nrow(data) # 9491
length(unique(data$pragmaid)) # 9491


##  6) ADD CLASSES FROM LCA
#   .............................................

# A) SEQUENCE CLASSES
data <- merge(data,
              class.dat, 
              by = "pragmaid", all.x = T)

# missing class
data[is.na(predclass), predclass := 0]
data[, table(predclass)]
round(data[, prop.table(table(predclass))],3)

# class labels:
data[, class := dplyr::recode(predclass,
                              "0" = "0: no intervention (n=6631)",
                              "3" = "1: brief psychiatric care (n=1267)",
                              "4" = "2: inpatient standard only (n=255)",
                              "5" = "3: inpatient intense (n=597)",
                              "1" = "4: rehabilitation (n=366)",
                              "6" = "5: counselling (n=267)",
                              "2" = "6: mixed (n=108)")]

data[, class_lab := dplyr::recode(predclass,
                              "0" = "0\nno intervention\n(69%)",
                              "3" = "1\nbrief psychiatric care\n(13%)",
                              "4" = "2\ninpatient standard only\n(3%)",
                              "5" = "3\ninpatient intense\n(6%)",
                              "1" = "4\nrehabilitation\n(4%)",
                              "6" = "5\ncounselling\n(3%)",
                              "2" = "6\nmixed\n(1%)")]

data[, class_lab2 := dplyr::recode(predclass,
                                  "0" = "0:no intervention",
                                  "3" = "1:brief psychiatric care",
                                  "4" = "2:inpatient standard only",
                                  "5" = "3:inpatient intense",
                                  "1" = "4:rehabilitation",
                                  "6" = "5:counselling",
                                  "2" = "6:mixed")]

data$class <- as.factor(data$class)
data$class_lab <- as.factor(data$class_lab)
data$class_lab2 <- as.factor(data$class_lab2)
data[, table(class, useNA = "always")]
data[, table(class_lab, useNA = "always")]
#data$predclass <- NULL

data[, table(class,interv.any)] # should fit 100%

# dimensions
nrow(data) # 9491
length(unique(data$pragmaid)) # 9491


##  7) ADD INTERVENTIONS
#   .............................................

add <- unique(interv.dat[,.(pragmaid,interv.type,interv_id)])
add <- unique(add[, .(check = any(!is.na(interv_id))), by = .(pragmaid,interv.type)])
add <- dcast(add, pragmaid ~ interv.type, value.var = "check")

data <- merge(data,
              add, 
              by = "pragmaid", all.x = T)
rm(add)

# get number of intervention types:
vars <- names(data)[names(data) %like% "bado|inpat|medi|psych_full|psych_short|qwt|reha"]
data[, interv_types_n := rowSums(.SD),by = pragmaid, .SDcols = vars]
rm(vars)

# compare with interv.any:
data[, table(interv_types_n,interv.any)] # consistent
data[interv.any == T] # 2860

# dimensions
nrow(data) # 9491
length(unique(data$pragmaid)) # 9491



# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 3) OUTPUT
# ______________________________________________________________________________________________________________________

saveRDS(data, paste0(inpath,"input data_person level.RDS"))
saveRDS(interv.dat, paste0(inpath,"input data_intervention level.RDS"))
saveRDS(elix.dat, paste0(inpath,"input data_elixhauser detailed data.RDS"))

