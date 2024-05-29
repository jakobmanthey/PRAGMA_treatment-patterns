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

##  1) PERSON LEVEL DATA (no intervention-time information)
# -------------------------------------------------------

filename <- paste0(inpath,"input data_person level.RDS")
data <- readRDS(filename)[order(pragmaid)]
rm(filename)

##  2) INTERVENTION LEVEL DATA
# -------------------------------------------------------

filename <- paste0(inpath,"input data_intervention level.RDS")
interv.dat <- readRDS(filename)
rm(filename)

DATE <- "2024-05-28"
filename <- paste0("output/treatment_series_at_resolution_cache/treatment_series_24m_quarter_followup_1_aud_only_intervs_",DATE,".csv")
interv.long.dat <- data.table(read.csv(filename))
rm(filename)

##  3) ELIXHAUSER DETAILED DATA
# -------------------------------------------------------

filename <- paste0(inpath,"input data_elixhauser detailed data.RDS")
elix.dat <- readRDS(filename)
rm(filename)


# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 2) ANALYSES
# ______________________________________________________________________________________________________________________

##  1) Sample description
# -------------------------------------------------------


##  2) Treatment utilization
# -------------------------------------------------------

# proportion with any intervention:
data[interv.any == T] # 2860
data[, mean(interv.any)] # 30.1%
data[, mean(interv.any), by = sex] # more female
data[, mean(interv.any), by = agegroup] # more younger 
data[, mean(interv.any), by = age <= 64] # more younger
data[, mean(interv.any), by = emp.type] # more unemployed

# number of intervention types:
data[, summary(interv_types_n)]
data[interv_types_n>0,.(p1 = sum(interv_types_n==1)/.N,
                        p2 = sum(interv_types_n==2)/.N,
                        p3 = sum(interv_types_n==3)/.N,
                        p4p = sum(interv_types_n>=4)/.N)]


##  3) Treatment utilization patterns
# -------------------------------------------------------

# number of interventions:
interv.dat[!is.na(date.interv.start), sum(!is.na(date.interv.start)), by = .(pragmaid,interv.type)][,mean(V1), by = interv.type]

# % intervention
data[, mean(bado)] # 4%
data[, mean(inpat)] # 7%
data[, mean(medi)] # 1%
data[, mean(psych_full)] # 0.7%
data[, mean(psych_short)] # 17%
data[, mean(qwt)] # 9%
data[, mean(reha)] # 4%

# Rehabilitation and QWT
data[, table(qwt,reha)]
data[, prop.table(table(qwt,reha),1)]

##  CLASS PROPORTIONS:
data[, prop.table(table(class))]
data[interv.any == T, prop.table(table(class))]

##  INTERVENTIONS WITHIN EACH CLASS:
data[, mean(bado), by = class] # class5: 100%
data[, mean(inpat), by = class] # class2: 100%
data[, mean(medi), by = class] # class6: 92%
data[, mean(psych_full), by = class] # class6: 10%
data[, mean(psych_short), by = class] # class1: 100% 
data[, mean(qwt), by = class] # class4: 100%
data[, mean(reha), by = class] # class3: 100%

### better presentation:
temp <- data[,.(class,bado,inpat,medi,psych_full,psych_short,qwt,reha)]
temp <- melt(temp, id.vars = "class")
temp[, mean(value), by = .(class,variable)][order(class,variable)]

##  4) Treatment utilization patterns and sociodemographics
# -------------------------------------------------------

# sex:
data[, median(age), by = class]
summary(glm(sex == "female" ~  age + emp.type, data, family = "binomial"))
summary(glm(sex == "female" ~  age + emp.type + class, data, family = "binomial"))
# --> sig diff for class 1 and 6 (p 0.01)

# age:
summary(glm(age ~ sex + emp.type, data, family = "gaussian"))
summary(glm(age ~ sex + emp.type + class, data, family = "gaussian"))
# --> sig diff for class 1-6 (p 0.01)

# unemp:
summary(glm(emp.type == "unemployed" ~  age + sex, data, family = "binomial"))
summary(glm(emp.type == "unemployed" ~  age + sex + class, data, family = "binomial"))
# --> sig diff for class 1 (p 0.01): exp(0.258032) = 1.29

##  5) Treatment utilization patterns and comorbidity
# -------------------------------------------------------

# comorbidity baseline:
data[, .(elix_all = round(mean(elix_sum_all_1),1),
         elix_nomental = round(mean(elix_sum_nomental_1),1)), by = class][order(class)]

summary(glm(elix_sum_nomental_1 ~  age + sex + emp.type, data, family = "poisson"))
summary(glm(elix_sum_nomental_1 ~  age + sex + emp.type + class, data, family = "poisson"))
# --> sig diff for classes 1-4



# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 3) TABLES
# ______________________________________________________________________________________________________________________

##  1) TABLE 1
#   .............................................

# description of interventions

##  2) TABLE 2
#   .............................................

vars <- names(data)[names(data) %like% "sex|age$|emp.type|elix_sum_nomental"]

tab1 <- rbind(data[,c(list(group = "All"), .SD),.SDcols = vars],
              data[interv.any == T,c(list(group = "Any intervention"), .SD),.SDcols = vars],
              data[bado == T,c(list(group = "Counselling"), .SD),.SDcols = vars],
              data[inpat == T,c(list(group = "Inpatient (standard)"), .SD),.SDcols = vars],
              data[qwt == T,c(list(group = "Inpatient (intense)"), .SD),.SDcols = vars],
              data[reha == T,c(list(group = "Rehabilitation"), .SD),.SDcols = vars],
              data[medi == T,c(list(group = "Pharmacological"), .SD),.SDcols = vars],
              data[psych_full == T,c(list(group = "Psychotherapy"), .SD),.SDcols = vars],
              data[psych_short == T,c(list(group = "Psych. brief contact"), .SD),.SDcols = vars])

tab1$group <- factor(tab1$group, levels = unique(tab1$group))
tab1$sex <- factor(tab1$sex)
tab1$emp.type <- factor(tab1$emp.type, 
                        levels = c("employed","unemployed","retired","other"))

tab1out <- tab1[, .(.N, 
                    female_prop = format(sum(sex == "female")/.N * 100, digits = 3, nsmall = 1),
                    age_mean = format(mean(age), digits = 3, nsmall = 1),
                    age_iqr_low = format(quantile(age, 0.25), digits = 3, nsmall = 1),
                    age_iqr_high = format(quantile(age, 0.75), digits = 3, nsmall = 1),
                    empl_prop = format(sum(emp.type == "employed")/.N * 100, digits = 3, nsmall = 1),
                    elix_mean = format(mean(elix_sum_nomental_1), digits = 2, nsmall = 1),
                    elix_iqr_low = format(quantile(elix_sum_nomental_1, 0.25), digits = 3, nsmall = 1),
                    elix_iqr_high = format(quantile(elix_sum_nomental_1, 0.75), digits = 3, nsmall = 1)), 
                by = group]

tab1out <- tab1out[,.(group,N,female_prop,
                      age = paste0(age_mean," (",age_iqr_low," to ",age_iqr_high,")"),
                      empl_prop,
                      elix = paste0(elix_mean," (",elix_iqr_low," to ",elix_iqr_high,")"))]

write.csv(tab1out,
          paste0(outpath,"tables/TAB1_descriptives_",Sys.Date(),".csv"),
          row.names = F)



# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 4) FIGURES
# ______________________________________________________________________________________________________________________

##  Fig 1) class comparison over time
#   .............................................

pdat <- merge(data[,.(pragmaid,class_lab)],
              interv.long.dat, by = "pragmaid")[!class_lab %like% "no intervention"]
pdat[, n_class := length(unique(pragmaid)), by = class_lab]
vars <- names(pdat)[names(pdat) %like% "bado|inpat|medi|psych_full|psych_short|qwt|reha"]
pdat <- melt(pdat, id.vars = c("pragmaid","class_lab","n_class","rel_time"), measure.vars = vars)[order(variable)]
pdat[,dv := sum(value)/n_class, by = .(class_lab,rel_time,variable)] 

pdat[, table(class_lab)]

pdat$intervention <- factor(pdat$variable, levels = c("bado","inpat","qwt","reha","medi","psych_full","psych_short"))
levels(pdat$intervention) <- c("Counselling",
                               "Inpatient (standard)",
                               "Inpatient (intense)",
                               "Rehabilitation",
                               "Pharmacological",
                               "Psychotherapy",
                               "Psych. brief contact")
pdat[, table(intervention,variable)]
colors <- ggthemes::ggthemes_data$gdocs$colors$value[3:9]

ggplot(pdat, aes(x = rel_time, y = dv, col = intervention)) +
  facet_grid(.~class_lab) + 
  geom_line(linewidth = 1.2) +
  scale_color_manual("",values = colors) +
  scale_x_continuous("Number of quarters after AUD diagnosis") +
  scale_y_continuous("[%] Persons in class with treatment type at time X", label = scales::percent) + 
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 4, nrow = 2))

ggsave(paste0(outpath,"figures/Fig 1_classes_overview over time_",Sys.Date(),".png"), width = 12, height = 6)

rm(pdat)

##  Fig 2) comorbidity
#   .............................................

pdat <- copy(data[,.SD, .SDcols = names(data)[names(data) %like% "pragmaid|class_lab2|elix_sum_nomental_1"]])
pdat$class_rev <- factor(pdat$class_lab2, levels = rev(levels(pdat$class_lab2)))

pdat[, mean := mean(elix_sum_nomental_1), by = class_rev]

ggplot(pdat, aes(x = class_rev, y = elix_sum_nomental_1, fill = class_rev)) +
  geom_violin() +
  geom_point(aes(y = mean), size = 3, shape = 25, fill = "black") +
  scale_fill_viridis_d("", guide = F) +
  coord_flip() +
  theme(legend.position = "bottom") +
  scale_x_discrete("") + 
  scale_y_continuous("Sum Score Elixhauser Comorbidity Index (0-27)")

ggsave(paste0(outpath,"figures/Fig 3_classes_comorbidity_",Sys.Date(),".png"), width = 9, height = 5)

rm(pdat)



# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 5) SUPPLEMENTARY FIGURES
# ______________________________________________________________________________________________________________________

##  Supplementary Figure 1) classes by age and sex
#   .............................................

pdat <- copy(data[,.SD, .SDcols = names(data)[names(data) %like% "pragmaid|class_lab2|sex|age"]])

pdat$class_rev <- factor(pdat$class_lab2, levels = rev(levels(pdat$class_lab2)))

pdat_prop <- unique(pdat[, .(n_female = sum(sex == "female"),n_group = .N), by = .(class_rev)])
pdat_prop[, prop_class_female := round(n_female/n_group,2)]

ggplot(pdat, aes(x = class_rev, y = age)) +
  geom_boxplot(aes(fill = class_rev)) + 
  geom_label(data = pdat_prop, aes(x = class_rev,y = 80, label = scales::percent(prop_class_female), group = class_rev), 
             position = position_dodge(width = 0.75), fill = "#FFDF00") +
  scale_fill_viridis_d(guide = F) +
  scale_x_discrete("") +
  scale_y_continuous("Age in years") +
  theme(legend.position = "bottom") +
  coord_flip()

ggsave(paste0(outpath,"figures/Suppl Fig 1_classes_sex and age_",Sys.Date(),".png"), width = 10, height = 5)

rm(pdat, pdat_prop)


##  Supplementary Figure 2) employment
#   .............................................

pdat <- copy(data[,.SD, .SDcols = names(data)[names(data) %like% "pragmaid|class_lab2|emp.type|sex"]])
pdat$class_rev <- factor(pdat$class_lab2, levels = rev(levels(pdat$class_lab2)))

pdat$emp.type   <- factor(pdat$emp.type, levels = rev(levels(pdat$emp.type)))

ggplot(pdat, aes(x = class_rev, fill = emp.type)) +
  geom_bar(position = "fill") +
  scale_fill_viridis_d("") +
  coord_flip() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 4, nrow = 1, reverse = TRUE)) + 
  scale_x_discrete("") + 
  scale_y_continuous("", label = scales::percent)

ggsave(paste0(outpath,"figures/Suppl Fig 2_classes_employment_",Sys.Date(),".png"), width = 10, height = 5)

rm(pdat)


##  Supplementary Figure 3) heatmap elixhauser
#   .............................................

pdat <- merge(data[,.(pragmaid,class_lab2)],
              elix.dat[elix.dat$period ==1,],by = "pragmaid", all.x = T)
pdat <- melt(pdat, id.vars = c("pragmaid","period", "class_lab2"))
pdat[, n_class := length(unique(pragmaid)), by = class_lab2]

pdat2 <- unique(pdat[, .(percentage = sum(value == 1)/n_class,n_class), by = .(class_lab2,variable)])

pdat2$label <- pdat2$variable
levels(pdat2$label) <- attr(elix.dat, "variable.labels")[2:32]

ggplot(pdat2, aes(x = class_lab2, y = label, fill = percentage)) +
  geom_tile(show.legend = F) + 
  geom_text(aes(label = scales::percent(percentage, accuracy = 1)), color = "black") +
  scale_fill_gradient(low = "#FFDAB9", high = "#FF0000") +
  scale_x_discrete("") +
  scale_y_discrete("") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) 

ggsave(paste0(outpath,"figures/Suppl Fig 3_classes_heatmap comorbidity_",Sys.Date(),".png"), width = 10, height = 10)

rm(pdat, pdat2)


##  Supplementary Figure 4) comorbidity over time
#   .............................................

pdat <- copy(data[,.SD, .SDcols = names(data)[names(data) %like% "pragmaid|class_lab2|elix_sum_nomental"]])
pdat <- melt(pdat, id.vars = c("pragmaid","class_lab2"), variable.name = "period")
pdat[, period := substr(period,19,19)]
pdat[, period_num := as.numeric(period)]

pdat[, mean := mean(value), by = .(class_lab2,period)]

#pdat$class_rev <- factor(pdat$class_lab2, levels = rev(levels(pdat$class_lab2)))

ggplot(pdat, aes(x = period_num, y = mean, fill = class_lab2, color = class_lab2)) +
  geom_point(size = 3, shape = 25) +
  geom_line() +
  scale_fill_viridis_d("") +
  scale_color_viridis_d("") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 4, nrow = 2)) + 
  scale_x_continuous("Year", breaks = c(1,2,3)) + 
  scale_y_continuous("Sum Score\nElixhauser Comorbidity Index (0-27)")

ggsave(paste0(outpath,"figures/Suppl Fig 4_classes_comorbidity over time_",Sys.Date(),".png"), width = 9, height = 5)

rm(pdat)

