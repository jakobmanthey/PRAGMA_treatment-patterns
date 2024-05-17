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

##  1) PRAGMA-ID
# -------------------------------------------------------

filename <- paste0(inpath,"0_pragma_id_GKV with Stammdata_2024-05-16.rds")
id.dat <- readRDS(filename)[order(pragmaid)]
rm(filename)

##  2) INTERVENTIONS
# -------------------------------------------------------

filename <- paste0(inpath,"4_data_AUD diagnosis and intervention sequences_2024-05-16.rds")
interv.dat <- readRDS(filename)
interv.dat <- interv.dat[period == "24m"]
rm(filename)

##  3) DIAGNOSES
# -------------------------------------------------------

filename <- paste0(inpath,"1_data_all diagnoses_2024-05-16.rds")
diag.dat <- readRDS(filename)
aud.dat <- diag.dat[icd.alc == T]
rm(filename)

##  4) EMPLOYMENT
# -------------------------------------------------------

filename <- paste0(inpath,"1_data_employment periods_2024-05-16.rds")
emp.dat <- readRDS(filename)
rm(filename)

##  5) INCOME
# -------------------------------------------------------

filename <- paste0(inpath,"1_data_income data_2024-05-16.rds")
inc.dat <- readRDS(filename)
rm(filename)

##  6) SEQUENCE GROUPINGS
# -------------------------------------------------------

filename <- paste0(outpath,"sequence_data_24m_monthly_lca_results.csv")
class.dat <- data.table(read.csv(filename))
rm(filename)

##  7) Temporary data
# -------------------------------------------------------

filename <- paste0(inpath,"1_data_insurance periods_2024-05-16.rds")
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

data <- copy(id.dat[,.(pragmaid,gkv,sex,yob,nationality)])

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

# dimensions
nrow(data) # 9496
length(unique(data$pragmaid)) # 9496


##  2) ADD COMORBIDITY
#   .............................................

comorb.dat <- diag.dat[pragmaid %in% unique(data$pragmaid)]

##  keep all (or only certain?) diagnoses
comorb.dat[, table(setting, icd_type)]
comorb.dat <- comorb.dat[icd_type %in% c("confirmed","primary","any")]
nrow(comorb.dat) # 2695273

##  add boundaries of look-behind period
sel.dat <- unique(interv.dat[,.(pragmaid,
                   start = date.period.start,
                   end = date.aud)])

comorb.dat <- merge(comorb.dat,
                    sel.dat, 
                    by = c("pragmaid"), all.x = T)

nrow(comorb.dat) # 2695273
rm(sel.dat)

##  restrict to 12m period prior to analyses
comorb.dat <- comorb.dat[date.diag.start %between% list(start,end)]
nrow(comorb.dat) # 526427

##  remove people not matched
comorb.dat <- comorb.dat[!(is.na(start) | is.na(end))]
nrow(comorb.dat) # 526427 -> no removals

##  keep only relevant information:
comorb.dat <- unique(comorb.dat[,.(pragmaid,icd)])

##  get comorbidity indices
elix.dat <- comorbidity::comorbidity(x = comorb.dat,
                                     id = "pragmaid",
                                     code = "icd",
                                     map = "elixhauser_icd10_quan",
                                     assign0 = T,
                                     tidy.codes = T)

names(elix.dat)[!names(elix.dat) %like% "pragmaid"] <- paste0("elix_",names(elix.dat)[!names(elix.dat) %like% "pragmaid"])

##  get sum score
elix.dat_sum <- data.table(elix.dat)

  ### A) all categories
  elix.dat_sum$elix_sum_all <- rowSums(elix.dat_sum[,.SD, .SDcols = names(elix.dat_sum)[names(elix.dat_sum) %like% "^elix"]])
  ### B) all categories but not alcohol, drug, psycho, depre
  elix.dat_sum$elix_sum_nomental <- rowSums(elix.dat_sum[,.SD, .SDcols = names(elix.dat_sum)[names(elix.dat_sum) %like% "^elix" &
                                                                              !names(elix.dat_sum) %like% "alcohol$|drug$|psycho$|depre$" & 
                                                                                !names(elix.dat_sum) %like% "^elix_sum"]])

data <- merge(data,
              elix.dat_sum[,.(pragmaid,elix_sum_all,elix_sum_nomental)], 
              by = "pragmaid", all.x = T)

# dimensions
nrow(data) # 9496
length(unique(data$pragmaid)) # 9496

##  3) ADD OTHER RELEVANT TIME-INVARYING VARIABLES
#   .............................................

add <- unique(interv.dat[,.(pragmaid,date.aud,date.period.start,date.period.end,interv.any)])

data <- merge(data,add,by = "pragmaid")
rm(add)

# dimensions
nrow(data) # 9496
length(unique(data$pragmaid)) # 9496

##  4) DEFINE AGEGROUP
#   .............................................

data[, age := year(date.aud) - yob]
data[, summary(age)]
data[is.na(age)] ## none

# remove people aged >18 at time of AUD diagnoses --> remove before LCA
data <- data[age>=18]

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
nrow(data) # 9488
length(unique(data$pragmaid)) # 9488

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
nrow(data) # 9488
length(unique(data$pragmaid)) # 9488


##  6) ADD GROUPINGS
#   .............................................

# A) SEQUENCE CLASSES
data <- merge(data,
              class.dat, 
              by = "pragmaid", all.x = T)

# missing class
data[is.na(predclass), predclass := 0]
data[, table(predclass)]

# class labels:
data[, class := dplyr::recode(predclass,
                              "0" = "0: no intervention (n=6651)",
                              "4" = "4: psychiatric/psychological care (n=1256)",
                              "2" = "2: counselling (n=158)",
                              "1" = "1: inpatient contacts (n=535)",
                              "3" = "3: qwt and other (n=433)",
                              "5" = "5: rehabilitation (n=360)",
                              "6" = "6: mixed (n=99)")]
data[, class_de := dplyr::recode(predclass,
                              "0" = "0: Keine Intervention",
                              "4" = "1: Vorwiegend psych. Kurzkontakte",
                              "2" = "5: Nur Suchtberatung",
                              "1" = "2: Vorwiegend stationäre Behandlung",
                              "3" = "3: Vorwiegend qualifizierter Entzug",
                              "5" = "4: Rehabilitation+",
                              "6" = "6: Alles")]

data$class <- as.factor(data$class)
data$class_de <- factor(data$class_de, levels = unique(data$class_de)[order(unique(data$class_de))])
data[, table(class, useNA = "always")]
data[, table(class_de, useNA = "always")]
data[, table(class_de,class)]
data$predclass <- NULL

data[, table(class,interv.any)] # should fit 100%

# dimensions
nrow(data) # 9488
length(unique(data$pragmaid)) # 9488


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
data[, interv_n := rowSums(.SD),by = pragmaid, .SDcols = vars]
rm(vars)

# compare with interv.any:
data[, table(interv_n,interv.any)] # consistent
data[interv.any == T] # 2837

# dimensions
nrow(data) # 9481
length(unique(data$pragmaid)) # 9481






# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 3) ANALYSES
# ______________________________________________________________________________________________________________________

# proportion with any intervention:
data[interv.any == T] # 2837
data[, mean(interv.any)] # 29.9%
data[, mean(interv.any), by = sex] # more female
data[, mean(interv.any), by = agegroup] # more younger 
data[, mean(interv.any), by = age <= 64] # more younger
data[, mean(interv.any), by = emp.type] # more unemployed

# number of interventions:
interv.dat[!is.na(date.interv.start), sum(!is.na(date.interv.start)), by = .(pragmaid,interv.type)][,mean(V1), by = interv.type]

# % intervention
data[, mean(bado)] # 3.6%
data[, mean(inpat)] # 7%
data[, mean(medi)] # 1%
data[, mean(psych_full)] # 0,7%
data[, mean(psych_short)] # 17%
data[, mean(qwt)] # 9%
data[, mean(reha)] # 4%

# Rehabilitation and QWT
data[, table(qwt,reha)]
data[, prop.table(table(qwt,reha),1)]


# sex:
summary(glm(sex == "female" ~  age, data, family = "binomial"))
summary(glm(sex == "female" ~  age + class_de, data, family = "binomial"))



# unemp:
summary(glm(emp.type == "unemployed" ~  age, data, family = "binomial"))
summary(glm(emp.type == "unemployed" ~  age + sex, data, family = "binomial"))
summary(glm(emp.type == "unemployed" ~  age + sex + class_de, data, family = "binomial"))


# comorbidity:
data[, .(elix_all = round(mean(elix_sum_all),1),
         elix_nomental = round(mean(elix_sum_nomental),1)), by = class][order(class)]

summary(glm(elix_sum_nomental ~  age, data, family = "gaussian"))
summary(glm(elix_sum_nomental ~  age, data, family = "poisson"))

summary(glm(elix_sum_nomental ~  age + sex, data, family = "poisson"))
summary(glm(elix_sum_nomental ~  age + sex + class_de, data, family = "poisson"))



# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 4) TABLES
# ______________________________________________________________________________________________________________________

##  1) TABLE 1
#   .............................................

vars <- names(data)[names(data) %like% "sex|age$|emp.type|elix_sum_nomental"]

tab1 <- rbind(data[,c(list(group = "Alle"), .SD),.SDcols = vars],
              data[interv.any == T,c(list(group = "Irgendeine Intervention"), .SD),.SDcols = vars],
              data[bado == T,c(list(group = "Suchtberatung"), .SD),.SDcols = vars],
              data[inpat == T,c(list(group = "Stationär"), .SD),.SDcols = vars],
              data[qwt == T,c(list(group = "Qualifizierter Entzug"), .SD),.SDcols = vars],
              data[reha == T,c(list(group = "Reha"), .SD),.SDcols = vars],
              data[medi == T,c(list(group = "Pharma"), .SD),.SDcols = vars],
              data[psych_full == T,c(list(group = "Psychotherapie"), .SD),.SDcols = vars],
              data[psych_short == T,c(list(group = "Psych. Kurzkontakt"), .SD),.SDcols = vars])

tab1$group <- factor(tab1$group, levels = unique(tab1$group))
tab1[, sex := ifelse(sex == "female", "weiblich", "männlich")]
tab1$sex <- factor(tab1$sex)
tab1[, emp.type := ifelse(emp.type == "employed", "in Arbeit",
                          ifelse(emp.type == "unemployed", "arbeitslos",
                                 ifelse(emp.type == "retired", "in Rente", "sonstiges")))]
tab1$emp.type <- factor(tab1$emp.type)

tab1out <- tab1[, .(.N, 
                    female_prop = sum(sex == "weiblich")/.N,
                    age_mean = mean(age),
                    age_low = quantile(age, 0.25),
                    age_high = quantile(age, 0.75),
                    empl_prop = sum(emp.type == "in Arbeit")/.N,
                    elix_mean = mean(elix_sum_nomental),
                    elix_low = quantile(elix_sum_nomental, 0.25),
                    elix_high = quantile(elix_sum_nomental, 0.75)), by = group]



# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 5) FIGURES INTERVENTIONS (GERMAN!)
# ______________________________________________________________________________________________________________________

##  1) Intervention
#   .............................................

pdat <- copy(tab1[!group %like% "Alle|Intervention"])
colors <- ggthemes::ggthemes_data$gdocs$colors$value[3:9]

ggplot(pdat, aes(x = group, fill = group)) +
  #ggtitle("Anzahl Personen mit mind. 1 Intervention nach der ersten Alkohol-Diagnose") + 
  geom_bar() + 
  scale_fill_manual("", values = colors, guide = F) +
  scale_x_discrete("") +
  scale_y_continuous("")

ggsave(paste0(outpath,"figures/DEU_groups_count_",Sys.Date(),".png"), width = 11, height = 6)


##  2) sex and age
#   .............................................

pdat <- copy(tab1)

pdat$group_rev <- factor(pdat$group, levels = rev(levels(pdat$group)))
levels(pdat$group)
levels(pdat$group_rev)

pdat_prop <- unique(pdat[, .(n_female = sum(sex == "weiblich"),n_group = .N), by = .(group_rev)])
pdat_prop[, prop_class_female := round(n_female/n_group,2)]

colors <- rev(ggthemes::ggthemes_data$gdocs$colors$value[1:9])

ggplot(pdat, aes(x = group_rev, y = age)) +
  geom_boxplot(aes(fill = group_rev)) + 
  geom_label(data = pdat_prop, aes(x = group_rev,y = 80, label = scales::percent(prop_class_female), group = group_rev), 
             position = position_dodge(width = 0.75), fill = "#FFDF00") +
  scale_fill_manual(values = colors, guide = F) +
  scale_x_discrete("") +
  scale_y_continuous("Alter in Jahren") +
  theme(legend.position = "bottom") +
  coord_flip()
  
ggsave(paste0(outpath,"figures/DEU_groups_sex and age",Sys.Date(),".png"), width = 10, height = 5)


##  3) employment
#   .............................................

pdat <- copy(tab1)

pdat$group_rev <- factor(pdat$group, levels = rev(levels(pdat$group)))
levels(pdat$group)
levels(pdat$group_rev)

ggplot(pdat, aes(x = group_rev, fill = emp.type)) +
  geom_bar(position = "fill") +
  scale_fill_viridis_d("") +
  coord_flip() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 4, nrow = 1, reverse = TRUE)) + 
  scale_x_discrete("") + 
  scale_y_continuous("", label = scales::percent) 

ggsave(paste0(outpath,"figures/DEU_groups_employment_",Sys.Date(),".png"), width = 10, height = 5)

rm(pdat)

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 5) FIGURES CLASSES (GERMAN!)
# ______________________________________________________________________________________________________________________

##  1) age and sex
#   .............................................

pdat <- copy(data[,.SD, .SDcols = names(data)[names(data) %like% "pragmaid|class_de|sex|age"]])

pdat$class_rev <- factor(pdat$class_de, levels = rev(levels(pdat$class_de)))

pdat_prop <- unique(pdat[, .(n_female = sum(sex == "female"),n_group = .N), by = .(class_rev)])
pdat_prop[, prop_class_female := round(n_female/n_group,2)]

ggplot(pdat, aes(x = class_rev, y = age)) +
  geom_boxplot(aes(fill = class_rev)) + 
  geom_label(data = pdat_prop, aes(x = class_rev,y = 80, label = scales::percent(prop_class_female), group = class_rev), 
             position = position_dodge(width = 0.75), fill = "#FFDF00") +
  scale_fill_viridis_d(guide = F) +
  scale_x_discrete("") +
  scale_y_continuous("Alter in Jahren") +
  theme(legend.position = "bottom") +
  coord_flip()

ggsave(paste0(outpath,"figures/DEU_classes_sex and age_",Sys.Date(),".png"), width = 10, height = 5)

rm(pdat, pdat_prop)


##  2) employment
#   .............................................

pdat <- copy(data[,.SD, .SDcols = names(data)[names(data) %like% "pragmaid|class_de|emp.type|sex"]])
pdat$class_rev <- factor(pdat$class_de, levels = rev(levels(pdat$class_de)))

pdat[, emp.type := ifelse(emp.type == "employed", "in Arbeit",
                          ifelse(emp.type == "unemployed", "arbeitslos",
                                 ifelse(emp.type == "retired", "in Rente", "sonstiges")))]

ggplot(pdat, aes(x = class_rev, fill = emp.type)) +
  geom_bar(position = "fill") +
  scale_fill_viridis_d("") +
  coord_flip() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 4, nrow = 1, reverse = TRUE)) + 
  scale_x_discrete("") + 
  scale_y_continuous("", label = scales::percent)

ggsave(paste0(outpath,"figures/DEU_classes_employment_",Sys.Date(),".png"), width = 10, height = 5)

rm(pdat)


##  3) comorbidity
#   .............................................

pdat <- copy(data[,.SD, .SDcols = names(data)[names(data) %like% "pragmaid|class_de|elix|sex"]])
pdat$class_rev <- factor(pdat$class_de, levels = rev(levels(pdat$class_de)))
pdat[, mean := mean(elix_sum_nomental), by = class_rev]

ggplot(pdat, aes(x = class_rev, y = elix_sum_nomental, fill = class_rev)) +
  #geom_boxplot() +
  geom_violin() +
  geom_point(aes(y = mean), size = 3, shape = 25, fill = "black") +
  scale_fill_viridis_d("", guide = F) +
  coord_flip() +
  theme(legend.position = "bottom") +
  scale_x_discrete("") + 
  scale_y_continuous("Summenwert Elixhauser Komorbiditätsindex (1-27)")

ggsave(paste0(outpath,"figures/DEU_classes_comorbidity_",Sys.Date(),".png"), width = 10, height = 5)

rm(pdat)



##  4) heatmap elixhauser
#   .............................................

pdat <- merge(data[,.(pragmaid,class_de)],
              elix.dat,by = "pragmaid", all.x = T)
pdat <- melt(pdat, id.vars = c("pragmaid", "class_de"))
pdat[, n_class := length(unique(pragmaid)), by = class_de]

pdat2 <- unique(pdat[, .(percentage = sum(value == 1)/n_class,n_class), by = .(class_de,variable) ])

pdat2$label <- pdat2$variable
levels(pdat2$label) <- attr(elix.dat, "variable.labels")[2:32]

ggplot(pdat2, aes(x = class_de, y = label, fill = percentage)) +
  geom_tile(show.legend = F) + 
  geom_text(aes(label = scales::percent(percentage, accuracy = 1)), color = "black") +
  scale_fill_gradient(low = "#FFDAB9", high = "#FF0000") +
  scale_x_discrete("") +
  scale_y_discrete("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(paste0(outpath,"figures/DEU_classes_heatmap comorbidity_",Sys.Date(),".png"), width = 10, height = 10)

rm(pdat, pdat2)


# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 5) FIGURES CLASSES
# ______________________________________________________________________________________________________________________

##  2) age and sex
#   .............................................

pdat <- copy(data[,.SD, .SDcols = names(data)[names(data) %like% "pragmaid|class|sex|age"]])
pdat_prop <- unique(pdat[, .(n_male = sum(sex == "male"),n_class = .N), by = .(class)])
pdat_prop[, prop_class_sex := ifelse(sex=="male", n_male/n_class,  (n_class-n_male)/n_class)]

ggplot(pdat, aes(x = sex, y = age, fill = class_de)) +
  ggtitle("boxplots of age and % of class female/male sex") + 
  geom_boxplot() + 
  geom_label(data = pdat_prop, aes(x = sex, y = 80, label = scales::percent(prop_class_sex)), 
             position = position_dodge(width = 0.75)) +
  scale_fill_viridis_d() +
  scale_x_discrete("") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 4, nrow = 2))

ggsave(paste0(outpath,"figures/classes_sex and age.png"), width = 12, height = 8)

rm(pdat, pdat_prop)


##  3) employment
#   .............................................

pdat <- copy(data[,.SD, .SDcols = names(data)[names(data) %like% "pragmaid|class|emp.type|sex"]])
levels(pdat$class) <- rev(levels(pdat$class))

ggplot(pdat, aes(x = class, fill = emp.type)) +
  geom_bar(position = "fill") +
  facet_grid(. ~ sex) +
  scale_fill_viridis_d("") +
  coord_flip() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 2, nrow = 2)) + 
  scale_x_discrete("") + 
  scale_y_continuous("", label = scales::percent)

ggsave(paste0(outpath,"figures/classes_employment.png"), width = 10, height = 5)

rm(pdat)





##  1) heatmap elixhauser
#   .............................................

pdat <- merge(data[,.(pragmaid,class)],
              elix.dat,by = "pragmaid", all.x = T)
pdat <- melt(pdat, id.vars = c("pragmaid", "class"))
pdat[, n_class := length(unique(pragmaid)), by = class]

pdat2 <- unique(pdat[, .(percentage = sum(value == 1)/n_class,n_class), by = .(class,variable) ])

pdat2$label <- pdat2$variable
levels(pdat2$label) <- attr(elix.dat, "variable.labels")[2:32]

ggplot(pdat2, aes(x = class, y = label, fill = percentage)) +
  geom_tile(show.legend = F) + 
  geom_text(aes(label = scales::percent(percentage, accuracy = 1)), color = "black") +
  scale_fill_gradient(low = "#FFDAB9", high = "#FF0000") +
  scale_x_discrete("") +
  scale_y_discrete("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(paste0(outpath,"figures/fig_classes_heatmap comorbidity.png"), width = 10, height = 10)

rm(pdat, pdat2)