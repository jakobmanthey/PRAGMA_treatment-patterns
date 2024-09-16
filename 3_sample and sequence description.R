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


here::i_am("3_sample and sequence description.R")
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

# Elixhauser
data[, mean(elix_sum_nomental == 0)]
data[, mean(elix_sum_nomental %in% c(1,2))]
data[, mean(elix_sum_nomental >=3)]

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
data[, table(class)]
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

##  4) Treatment utilization patterns and sociodemographics --> section to be added by Kilian, including export of Suppl Table 2
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
data[, .(elix_all = round(mean(elix_sum_all),1),
         elix_nomental = round(mean(elix_sum_nomental),1)), by = class][order(class)]
data[, nat_deutsch := nationality == "deutsch"]
data[, age_z := (age-mean(age))/sd(age)]

# poisson (not reported)
summary(glm(elix_sum_nomental ~  age_z + sex + emp.type + nat_deutsch, data, family = "poisson"))
summary(glm(elix_sum_nomental ~  age_z + sex + emp.type + nat_deutsch + class, data, family = "poisson"))
# --> sig diff for classes 1-4

mod.out <- glm(elix_sum_nomental ~  age_z + sex + emp.type + nat_deutsch + class, data, family = "poisson")
stargazer::stargazer(mod.out, type = "html", ci = T, apply.coef = exp, p.auto = F, #report = c("v","c","s","p"),
                     out = paste0("output/tables/SUPP_TAB1_",Sys.Date(),".html"))

# zero-inflated

mod.out <- pscl::zeroinfl(elix_sum_nomental ~ age_z + sex + emp.type + nat_deutsch + class | age_z + sex, data)
summary(mod.out)
out <- data.table(var = names(coef(mod.out)),
                  coef = format(exp(coef(mod.out)), digits = 2, nsmall = 1),
                  cil = format(exp(confint(mod.out))[,1], digits = 2, nsmall = 1),
                  ciu = format(exp(confint(mod.out))[,2], digits = 2, nsmall = 1))
out <- out[,.(var, out = paste0(coef, " (",cil, " to ", ciu, ")"))]

write.csv(out,
          paste0(outpath,"tables/SUPP_TAB1_ZINB_",Sys.Date(),".csv"),
          row.names = F)

## 6) Multinomial Regression Analysis of Covariate influence on Latent Class Memberships
# -------------------------------------------------------


local({
covariates = c("pragmaid","class_lab2", "predclass", "emp.type", "age", "sex", "nationality", "income")
dat4 = readr::read_rds(here("input", "sequence_data_240528","input data_person level.RDS"))

dat4 = dat4 %>%
  select(all_of(covariates)) %>%
  mutate(sex = as.factor(sex),
         emp.type = as.factor(emp.type),
         nationality = as.factor(nationality),
         nationality = fct_collapse(nationality, "nicht deutsch" = c("nicht deutsch", "unbekannt")) , 
         age = case_when( age <= 34 ~ "18-34",
                          age > 34 & age <= 54 ~ "35-54",
                          age > 54 & age <= 64 ~ "55-64", 
                          age >= 65 ~ "65+") %>% as.factor() %>% fct_rev()
  )

set.seed(924)

mod_fit = multinom(class_lab2 ~ emp.type + age + sex + nationality,
                   data = dat4)
sm = summary(mod_fit)
sm


z = summary(mod_fit)$coefficients/summary(mod_fit)$standard.errors
p = (1 - pnorm(abs(z), 0, 1)) * 2
p %>% round(digits = 3) %>% print() #%>% knitr::kable() %>% print() #%>% kableExtra::save_kable(here("output", "multinom_regression_p_values.html"))

# Risk Ratios
rs = exp(coef(mod_fit)) %>% as.data.frame() 
rs %>% mutate(LC = rownames(rs), .before = "(Intercept)") %>% 
  pivot_longer(cols =  !c("(Intercept)", "LC"),
               names_to = "Independent Variable", 
               values_to = "Exp(coeffitient") %>%
   print() #%>% kableExtra::save_kable(here("output", "multinom_regression_exp_coefficients.html"))

round(exp(confint(mod_fit)),3) %>% print() #%>% knitr::kable() %>% print()# %>% kableExtra::save_kable(here("output", "multinom_regression_confidence_intervals.html"))
})



# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 3) TABLES
# ______________________________________________________________________________________________________________________

##  1) TABLE 1
#   .............................................

# description of interventions

# number of interventions:
interv.dat[!is.na(date.interv.start), sum(!is.na(date.interv.start)), by = .(pragmaid,interv.type)][,mean(V1), by = interv.type]

##  PSYCH-BRIEF
interv.dat[interv.type == "psych_short" & !is.na(date.interv.start),.(n = length(unique(date.interv.start))), by = pragmaid][,summary(n)]
##  PSYCH-LONG
interv.dat[interv.type == "psych_full" & !is.na(date.interv.start),.(n = length(unique(date.interv.start))), by = pragmaid][,summary(n)]
##  PHARMA
interv.dat[interv.type == "medi" & !is.na(date.interv.start),.(n = length(unique(date.interv.start))), by = pragmaid][,summary(n)]
##  INPAT-STAND
interv.dat[interv.type == "inpat" & !is.na(date.interv.start),.(n = length(unique(date.interv.start))), by = pragmaid][,summary(n)]
interv.dat[interv.type == "inpat" & !is.na(date.interv.start),.(n.days = date.interv.end - date.interv.start), by = pragmaid][,summary(as.numeric(n.days+1))]
##  REHA
interv.dat[interv.type == "reha" & !is.na(date.interv.start),.(n = length(unique(date.interv.start))), by = pragmaid][,summary(n)]
interv.dat[interv.type == "reha" & !is.na(date.interv.start),.(n.days = date.interv.end - date.interv.start), by = pragmaid][,summary(as.numeric(n.days+1))]

##  2) TABLE 2
#   .............................................

vars <- names(data)[names(data) %like% "sex|age$|emp.type|elix_sum_nomental|nationality"]

tab1 <- rbind(data[,c(list(group = "All"), .SD),.SDcols = vars],
              data[interv.any == T,c(list(group = "Any intervention"), .SD),.SDcols = vars],
              data[bado == T,c(list(group = "COUNSEL"), .SD),.SDcols = vars],
              data[inpat == T,c(list(group = "INPAT-STAND"), .SD),.SDcols = vars],
              data[qwt == T,c(list(group = "INPAT-INTENSE"), .SD),.SDcols = vars],
              data[reha == T,c(list(group = "REHA"), .SD),.SDcols = vars],
              data[medi == T,c(list(group = "PHARMA"), .SD),.SDcols = vars],
              data[psych_full == T,c(list(group = "PSYCH-FULL"), .SD),.SDcols = vars],
              data[psych_short == T,c(list(group = "PSYCH-BRIEF"), .SD),.SDcols = vars])

tab1[, table(nationality)] # regroup unbekannt --> 18 cases only
tab1[, nationality := dplyr::recode(nationality,"unbekannt" = "nicht deutsch")]

tab1$group <- factor(tab1$group, levels = unique(tab1$group))
tab1$sex <- factor(tab1$sex)
tab1$emp.type <- factor(tab1$emp.type, 
                        levels = c("employed","unemployed","retired","other"))

tab1out <- tab1[, .(.N, 
                    female_prop = format(sum(sex == "female")/.N * 100, digits = 3, nsmall = 1),
                    age_mean = format(mean(age), digits = 3, nsmall = 1),
                    age_iqr_low = format(quantile(age, 0.25), digits = 3, nsmall = 1),
                    age_iqr_high = format(quantile(age, 0.75), digits = 3, nsmall = 1),
                    nat_prop = format(sum(nationality == "deutsch")/.N * 100, digits = 3, nsmall = 1),
                    empl_prop = format(sum(emp.type == "employed")/.N * 100, digits = 3, nsmall = 1),
                    elix_mean = format(mean(elix_sum_nomental), digits = 2, nsmall = 1),
                    elix_iqr_low = format(quantile(elix_sum_nomental, 0.25), digits = 3, nsmall = 1),
                    elix_iqr_high = format(quantile(elix_sum_nomental, 0.75), digits = 3, nsmall = 1)), 
                by = group][order(N, decreasing = T)]

tab1out <- tab1out[,.(group,N,female_prop,
                      age = paste0(age_mean," (",age_iqr_low," to ",age_iqr_high,")"),
                      empl_prop,nat_prop,
                      elix = paste0(elix_mean," (",elix_iqr_low," to ",elix_iqr_high,")"))]

write.csv(tab1out,
          paste0(outpath,"tables/TAB1_descriptives_",Sys.Date(),".csv"),
          row.names = F)

##  test any vs no intervention:
temp <- data[,c(list(group = "Any intervention"), .SD),.SDcols = c(vars,"interv.any")]

  # sex difference
  temp[, prop.table(table(sex,interv.any),2)]
  temp[, chisq.test(interv.any,sex)] # p <.001

  # age difference
  temp[, t.test(age ~ interv.any)] # p <.001
  
  # employment difference
  temp[, prop.table(table(emp.type == "employed",interv.any),2)]
  temp[, chisq.test(emp.type == "employed",interv.any)] # p <.001

  # nationality difference
  temp[, prop.table(table(nationality == "deutsch",interv.any),2)]
  temp[, chisq.test(interv.any,nationality == "deutsch")] # p <.001
  
  # Elixhauser
  temp[, t.test(elix_sum_nomental ~ interv.any)] # p <.001
  

##
rm(vars, tab1, tab1out, temp)


##  3) Supplementary TABLE 1
#   .............................................

local({
m2 = readRDS(here("input", "LCA_model_2.rds"))
m3 = readRDS(here("input", "LCA_model_3.rds"))
m4 = readRDS(here("input", "LCA_model_4.rds"))
m5 = readRDS(here("input", "LCA_model_5.rds"))
m6 = readRDS(here("input", "LCA_model_6.rds"))
m7 = readRDS(here("input", "LCA_model_7.rds"))


source(here("upsample_daily_utils.R"))

sumtab = data.frame(Models = str_c(c(2,3,4,5,6,7), " Class"),
                    LL = c(m2$llik, m3$llik, m4$llik, m5$llik, m6$llik, m7$llik),
                    BIC = c(m2$bic, m3$bic, m4$bic, m5$bic, m6$bic, m7$bic), 
                    "Smallest Class (n)" = c(smallest_class(m2$predclass, "n"),
                                             smallest_class(m3$predclass, "n"),
                                             smallest_class(m4$predclass, "n"),
                                             smallest_class(m5$predclass, "n"),
                                             smallest_class(m6$predclass, "n"),
                                             smallest_class(m7$predclass, "n")),
                    
                    "Smallest Class (%)" = c(smallest_class(m2$predclass, "rel"),
                                             smallest_class(m3$predclass, "rel"),
                                             smallest_class(m4$predclass, "rel"),
                                             smallest_class(m5$predclass, "rel"),
                                             smallest_class(m6$predclass, "rel"),
                                             smallest_class(m7$predclass, "rel")),
                    Entropy = c(get_entropy(m2), 
                                get_entropy(m3),
                                get_entropy(m4),
                                get_entropy(m5),
                                get_entropy(m6),
                                get_entropy(m7)),
                    
                    "ALCPP min" = c(smallest_lcprob(m2), 
                                    smallest_lcprob(m3),
                                    smallest_lcprob(m4),
                                    smallest_lcprob(m5),
                                    smallest_lcprob(m6), 
                                    smallest_lcprob(m7))
)

sumtab %>% write_csv("output", paste0("Supp Tab 1_Model_Selection", Sys.Date(), ".csv"))

})
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
                               "Inpatient (intensive)",
                               "Rehabilitation",
                               "Pharmacological",
                               "Psychotherapy",
                               "Psych. brief contact")
pdat[, table(intervention,variable)]
colors <- ggthemes::ggthemes_data$gdocs$colors$value[3:9]
colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2")

ggplot(pdat, aes(x = rel_time, y = dv, col = intervention)) +
  facet_grid(.~class_lab) + 
  geom_line(linewidth = 1.2) +
  scale_color_manual("",values = colors) +
  scale_x_continuous("Number of quarters after AUD diagnosis") +
  scale_y_continuous("[%] Persons in class with treatment type at time X", label = scales::percent) + 
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 4, nrow = 2))

ggsave(paste0(outpath,"figures/Fig 1_classes_overview over time_",Sys.Date(),".png"), width = 12, height = 6)
ggsave(paste0(outpath,"figures/Fig 1_classes_overview over time_",Sys.Date(),".tiff"), width = 12, height = 6)

rm(pdat)

##  Fig 2) to be added by Kilian
#   .............................................
dat_tmp = read_csv(here("input", "data_for_fig2.csv"))
get_categorical_interv_num = function(num_intervs){
  if(num_intervs == 1) return("1")
  if(num_intervs > 1 & num_intervs <= 4) return("2 - 4")
  if(num_intervs > 4) return("5+")
}

get_categorical_interv = Vectorize(get_categorical_interv_num)

repeated_intervs = dat_tmp %>% 
  group_by(predclass) %>% 
  mutate(n_class = n_distinct(pragmaid)) %>% 
  ungroup() %>% 
  group_by(pragmaid, predclass, interv.type) %>%
  mutate(num_intervs = n(), more_than_one_interv = get_categorical_interv(num_intervs)) %>% 
  group_by(predclass, interv.type) %>% 
  mutate(n_had_interv = n_distinct(pragmaid)) %>% 
  group_by(predclass, interv.type, more_than_one_interv) %>% 
  summarize(n = n_distinct(pragmaid), frac = n/mean(n_had_interv), n_class = unique(n_class)) %>% 
  group_by(predclass, interv.type) %>% 
  mutate(class_frac = round(sum(n)/n_class, 2)* 100) %>% 
  ungroup()


repeated_intervs %>% 
  ggplot(aes(x = interv.type, y = frac * 100)) +
  geom_col(aes(fill = more_than_one_interv)) +
  geom_text(aes(x = interv.type, y = 90,label = paste0(class_frac, "%")), data = repeated_intervs %>% select(predclass, interv.type, class_frac) %>% distinct() ) +
  facet_wrap(~ predclass) +
  labs(x = "Intervention type", y = "Patients within intervention [%]", fill = "Intervention repetitions") + 
  coord_flip() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(size = 12), 
        axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15))


ggsave(paste0("output/","figures/","Figure 2_fig_interv_repetitions_barplot",Sys.Date(),".png"), width = 12, height = 6)
ggsave(paste0("output/","figures/","Figure 2_fig_interv_repetitions_barplot",Sys.Date(),".tiff"), width = 12, height = 6)

rm(dat_tmp)
rm(repeated_intervs)

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 5) SUPPLEMENTARY FIGURES
# ______________________________________________________________________________________________________________________

##  Supplementary Figure 1) heatmap elixhauser
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

##  Supplementary Figure 2) heatmap cross utilisation --> to be added by Kilian
#   .............................................

dat_tmp = read_csv(here("input", "data_for_supp_fig_2.csv"))

lookup = c("REHAB" = "reha",
           "PSYCH-BRIEF" = "psych_short", 
           "INPAT-STANDARD" = "inpat",
           "COUNSEL" = "bado",
           "PSYCH-LONG" = "psych_full",
           "PHARMA" = "medi",
           "INPAT-INTENSIVE" = "qwt"
)

# helper function  
get_overlap = function(dat,
                       cls,
                       var1,
                       var2,
                       perspective = c("class", "intervention")){
  
  
  perspective = match.arg(perspective)
  if(perspective == "class"){
    overlaps = dat %>%
      ungroup() %>% 
      filter(predclass == cls) %>% 
      select(pragmaid, all_of(c(var1,var2))) %>% 
      group_by(pragmaid) %>% 
      summarize(overlap = all(pick(var1),  pick(var2)))
  }
  else{
    dat = dat[c(dat[var1] == TRUE), ]
    overlaps = dat %>%
      ungroup() %>% 
      filter(predclass == cls) %>% 
      select(pragmaid, all_of(c(var1,var2))) %>% 
      group_by(pragmaid) %>% 
      summarize(overlap = all(pick(var1),  pick(var2)))
  }
  
  overlap = sum(overlaps$overlap)/nrow(overlaps)
  
  # res = data.frame(var = var1, val = overlap) 
  # names(res)[2] = var2
  return(overlap)
}
# helper function end

# helper function
get_inter_variable_overlaps = function(dat,
                                       cls,
                                       analysis_vars,
                                       perspective =  c("class", "intervention")){
  
  perspective = match.arg(perspective)
  
  
  resresres = data.frame()
  for(v1 in analysis_vars){
    resres = data.frame(var = v1)
    
    for(v2 in analysis_vars){
      res = data.frame(val = get_overlap(dat,
                                         cls,
                                         v1,
                                         v2,
                                         perspective = perspective))
      names(res) = v2
      resres = cbind(resres, res)
    }
    
    resresres = rbind(resresres, resres)
  }
  resresres$predclass = cls
  return(resresres)
}




overlaps_class_persp = map_dfr(unique(dat_tmp$predclass),
                               ~ get_inter_variable_overlaps(
                                 dat_tmp,
                                 cls =.x, 
                                 analysis_vars = names(lookup),
                                 perspective = "class"))


overlaps_interv_persp = map_dfr(unique(dat_tmp$predclass), 
                                ~ get_inter_variable_overlaps(
                                  dat_tmp,
                                  cls =.x, 
                                  analysis_vars = names(lookup),
                                  perspective = "intervention"))


remove_zero_perc = function(x){
  if_else(x == "0%", "", x)
}

heatmap_pal = viridis::viridis_pal(begin = 0.2, direction = -1)
cols = c("white", heatmap_pal(99))


overlaps_class_persp %>%
  pivot_longer(cols = all_of(names(lookup)),
               names_to = "var2", 
               values_to = "val") %>% 
  mutate(var = as.factor(var), var2 = as.factor(var2)) %>% 
  ggplot(aes(x = reorder(var, desc(var)), y = var2, fill = val*100)) +
  geom_tile() +
  facet_wrap(~ predclass, scales = "free") +
  geom_text(aes(label = paste0(round(val * 100, 0), "%") %>% remove_zero_perc()), size = 2) +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .4, size = 5),
        axis.text.y = element_text(size = 5),
        legend.position = "bottom") +
  scale_x_discrete(limits = rev(levels(var))) +
  labs(x = "Further Interventions",
       y = "Intervention", 
       fill = "Fraction of Class[%]") +
  scale_fill_gradientn(guide = guide_colorbar(
    theme = theme(
      legend.key.width = unit(10, "cm"),
      legend.key.height = unit(.2, "cm"), 
      legend.title = element_text(size = 10, vjust = 1.1), 
    )),
    colours  = cols)

ggsave(paste0("output/","figures/","Suppl Fig 2_fig_heatmap_class_perspektive",Sys.Date(),".png"), width = 10, height = 6)

overlaps_interv_persp %>%
  pivot_longer(cols = all_of(names(lookup)),
               names_to = "var2", 
               values_to = "val") %>% 
  mutate(var = as.factor(var), var2 = as.factor(var2), 
         val = ifelse(is.nan(val), 0, val)) %>% 
  ggplot(aes(x = reorder(var, desc(var)), y = var2, fill = val*100)) +
  geom_tile() +
  facet_wrap(~predclass) +
  geom_text(aes(label = paste0(round(val * 100, 0), "%") %>% remove_zero_perc()), size = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .4, size = 7),
        axis.text.y = element_text(size = 7), 
        legend.position = "bottom", 
        strip.text = element_text(size = 10)) +
  scale_x_discrete(limits = rev(levels(var))) +
  labs(x = "Primary Intervention",
       y = "Further Interventions",
       fill = "Persons with intervention X who also used Y [%]") +
  scale_fill_gradientn(guide = guide_colorbar(
    theme = theme(
      legend.key.width = unit(6, "cm"),
      legend.key.height = unit(.2, "cm"), 
      legend.title = element_text(size = 10, vjust = 1.1), 
    )),
    colours  = cols)


ggsave(paste0("output/","figures/","Suppl Fig 2_fig_heatmap_intervention_perspektive",Sys.Date(),"png"), width = 10, height = 6)


rm(dat_tmp)
rm(overlaps_class_persp)
rm(overlaps_interv_persp)
rm(lookup)

##  Supplementary Figure 3) distribution age and sex by class
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

ggsave(paste0(outpath,"figures/Suppl Fig 3_classes_sex and age_",Sys.Date(),".png"), width = 10, height = 5)

rm(pdat, pdat_prop)


##  Supplementary Figure 4) distribution employment by class
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

ggsave(paste0(outpath,"figures/Suppl Fig 4_classes_employment_",Sys.Date(),".png"), width = 10, height = 5)

rm(pdat)


##  Supplementary Figure 5) distribution comorbidity by class
#   .............................................

pdat <- copy(data[,.SD, .SDcols = names(data)[names(data) %like% "pragmaid|class_lab2|elix_sum_nomental"]])
pdat$class_rev <- factor(pdat$class_lab2, levels = rev(levels(pdat$class_lab2)))

pdat[, mean := mean(elix_sum_nomental), by = class_rev]

ggplot(pdat, aes(x = class_rev, y = elix_sum_nomental, fill = class_rev)) +
  geom_violin() +
  geom_point(aes(y = mean), size = 3, shape = 25, fill = "black") +
  scale_fill_viridis_d("", guide = "none") +
  coord_flip() +
  theme(legend.position = "bottom") +
  scale_x_discrete("") + 
  scale_y_continuous("Sum Score Elixhauser Comorbidity Index (0-27)")

ggsave(paste0(outpath,"figures/Suppl Fig 5_classes_comorbidity_",Sys.Date(),".png"), width = 9, height = 5)

rm(pdat)













