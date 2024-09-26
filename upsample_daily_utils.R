
####################TIDY AND PREPARE DATA FUNCTIONS##########################

################################ helper functions
rename_point_treatment = function(varname){
  treatname = str_extract(varname, "(?<=_).*" )
  return(paste0("date.", treatname))
}
rename_point_treatments = function(vars){
  return(unique(vars %>% map_chr(~rename_point_treatment(.x))))
}

rename_span_treatment = function(varname){
  treatname = str_extract(varname, "(?<=_).*" )
  if(str_detect(varname, "start"))  {return(paste0("date.", treatname, ".start"))}
  return(paste0("date.", treatname, ".end"))
}
rename_span_treatments = function(vars){
  return(unique(vars %>% map_chr(~rename_span_treatment(.x))))
}

# returns a vector of interventions which last only for one day within the data (such as medi pickup)
get_point_vars = function(dat){
  pointvars = dat %>%
      group_by(interv.type) %>% 
      dplyr::select(date.interv.start, date.interv.end) %>%
      drop_na() %>%
      summarize(is_pointvar = all(date.interv.start == date.interv.end)) %>%
      filter(is_pointvar == T) %>% 
      dplyr::select(interv.type) %>% 
      pull()
  return(pointvars)
}

# returns a vector of interventions which last for a span of time (such as reha or qwt)
get_span_vars = function(dat){
  spanvars = dat %>%
    group_by(interv.type) %>% 
  dplyr::select(date.interv.start, date.interv.end) %>%
    drop_na() %>%
    summarize(is_pointvar = all(date.interv.start == date.interv.end)) %>%
    filter(is_pointvar == F) %>% 
    dplyr::select(interv.type) %>% 
    pull()
  return(spanvars)
}
##################################### main function

# Pivots and tidys up sequence data provided by Jakob into a wide which can be up sampled by 
# upsample_daily_by_id(). Further standardises names and deletes unnecessary columns
widen_sequence_data = function(dat, # data as in sequence_data_for_kilian_xyz (with aud)
                               period_duration = "24m"){ # what observation period? must be one of the values contaied in dat$period
  
  if(!(period_duration %in% dat$period)) stop("provide a period duration which also exists in the column dat$period")
  variables_point = get_point_vars(dat) 
  variables_span = get_span_vars(dat)
  varnames_span = c(map_chr(variables_span, ~paste0("date.interv.start_", .x)),
                    map_chr(variables_span, ~paste0("date.interv.end_", .x))
  ) 
  varnames_point = map_chr(variables_point, ~paste0("date.interv.start_", .x))
  varnames_exclude = map_chr(variables_point, ~paste0("date.interv.end_", .x))             
  
  dat = dat %>%
    filter(period == period_duration) %>%  # Set the desired period length here!
    pivot_wider(id_cols = c(pragmaid, # TODO:this still could be automated
                            date.aud,
                            period,
                            date.period.start,
                            date.period.end,
                            interv.any,
                            interv_id),
                names_from = interv.type,
                values_from = c(date.interv.start, date.interv.end)) %>% 
    dplyr::select(-{{varnames_exclude}})  %>% 
    group_by(pragmaid) %>% 
    mutate(interv.any.overall = any(interv.any)) %>% 
    rename_with(rename_span_treatments, all_of(varnames_span)) %>% 
    rename_with(rename_point_treatments, all_of(varnames_point))
  return(dat)  
}

###############UPSAMPLING and DOWNSAMPLING FUNCTIONS##################################################

#################### helper functions #############
# join all dataframes within  contained in a list of "joinables"
join_it = function(joinables){
  suppressMessages({
    if(length(joinables) == 1) return(joinables[[1]])
    if(is_empty(joinables)) return(data.frame(pragmaid = character()))
    res = full_join(joinables[[1]], joinables[[2]])#, multiple = "any")
    
    if(length(joinables)>2){
      
      for(i in 3:length(joinables)){
        res = full_join(res, joinables[[i]])#, multiple = "any")
        
      }
    }  
  })
  return(res)
}



# dat is the data, variable_name is a string with the name of the respective variable
upsample_variable_timespan <- function(dat, variable_name) {
  
  # Column name adjustments
  start_name <- paste0("date.", variable_name, ".start")
  end_name <- paste0("date.", variable_name, ".end")
  lookup = c()
  
  # Pre-allocate result data frame
  res <- data.frame(
    pragmaid = character(),
    date = Date(),
    variable_name = logical(),  #  insert variable name later
    variable_name_id = numeric()
  )
  
  pragmaid <- dat$pragmaid %>% unique()
  tmp_y <- dat
  
  for (obs in 1:nrow(tmp_y)) {
    tmp_y_obs <- tmp_y[obs, ]
    
    if (!is.na(tmp_y_obs[[start_name]]) & !is.na(tmp_y_obs[[end_name]])) {
      res_tmp <- data.frame(
        pragmaid = pragmaid,
        date = seq(tmp_y_obs[[start_name]],
                   tmp_y_obs[[end_name]],
                   by = 1),
        variable_name = TRUE, 
        variable_name_id = tmp_y_obs$interv_id
      )
      res <- rbind(res, res_tmp)
    }
  }
  names(res) = c("pragmaid", "date", variable_name, paste0(variable_name, "_id"))
  return(res)
}


# dat is the data, variable_name is a string with the name of the respective variable
upsample_variable_timepoint = function(dat, variable_name){
  
  date_name <- paste0("date.", variable_name)
  
  
  res = data.frame(pragmaid = character(),
                   date = Date(),
                   variable_name = logical(),
                   variable_id = numeric()
  )
  
  pragmaid = dat$pragmaid %>% unique()
  tmp_y = dat
  
  for(obs in 1:nrow(tmp_y)){
    
    tmp_y_obs = tmp_y[obs, ]
    
    if(!is.na(tmp_y_obs[date_name])){
      
      res_tmp = data.frame(pragmaid = pragmaid,
                           date = tmp_y_obs[date_name],
                           variable_name = T,
                           variable_id = tmp_y_obs$interv_id
      )
      
      res = rbind(res, res_tmp)
    }
  }
  names(res) = c("pragmaid", "date", variable_name, paste0(variable_name, "_id"))
  return(res)
}


upsample_vars = function(spl,variable_name, point = F){
  if(point){return(upsample_variable_timepoint(spl, variable_name))}
  return(upsample_variable_timespan(spl, variable_name))
}



# Now wrap up above functions into one function which upsamples a data split for one pragmaid

# takes a split of data split by pragmaid as created by split(dat, .f = dat$pragmaid)[[x]]
upsample_daily_by_id = function(spl, variables_span, variables_point){
  res = join_it(list(
                      join_it(future_map(variables_span, ~ upsample_vars(spl = spl, variable_name = .x, point = F))),
                      join_it(future_map(variables_point, ~ upsample_vars(spl = spl, variable_name = .x, point = T)))
    )
  )
  return(res)
}


############ tests for upsample_daily_by_id
# variables_span = c("reha", "qwt", "psych_full", "inpat", "psych_short")
# variables_point = c("aud","medi")
# 
# # 90 rows, but first one is empty (no interv), thus expect upsample_daily_by_id
# # to output 89 rows since 89 psych_short treatments were applied for one day each
# tst_id = "OKWBEX8ZJ2" # had lots of psych treatments
# 
# # 25 rows but long treatment spans
# tst_id = "kXcIe1YPsx" # had all treatments
# tst = dat %>% filter(pragmaid == tst_id) 
# 
# tst_yearly = tibble(yearly %>% filter(pragmaid == tst_id))
# #test function
# tst_res = upsample_daily_by_id(tst, variables_span, variables_point)
# tst_resres = left_join(tst_yearly, tst_res)
# # expect 89 and 497 respectively
# tst_resres %>% dplyr::select(all_of(variables_span), medi) %>% sum(na.rm=T)


# this function creates a dataframe with the according timesteps and relative time scale grouping variables for the subsequent
# by time summarisation of treatments. such that its result can be joined with the result of upsample_data_to_daily. 
get_person_timesteps = function(person,
                                  period_duration = "24m",
                                  resol = c("day","week", "month", "bimonth", "quarter", "halfyear", "year")){
  
  resol = match.arg(resol)
  #TODO: Implement event handling so not only aud but also any intervention works as event
  aud_date = person %>% filter(period == period_duration) %>% 
    dplyr::select(date.aud, date.period.end) %>%
    distinct()
  id = person %>% dplyr::select(pragmaid) %>% distinct() %>% pull()
  
  
  
  resols = list("day" = days(1),
                "week" = weeks(1),
                "month" = months(1),
                "bimonth" = months(2),
                "quarter" = months(3),
                "halfyear" = months(6),
                "year" = years(1))
  
  timestep = resols[[resol]]
  date_cnt = aud_date$date.aud
  period_num = as.numeric(str_extract(period_duration, "\\d+"))
  
  timesteps = map_dfr(0:(months(period_num)/timestep), ~ data.frame(date_timestep = aud_date$date.aud %m+% (.x*timestep)))
  num_steps = length(0:(months(period_num)/timestep)) -1
  ts.x = timesteps[1:nrow(timesteps) - 1, ]
  ts.y = timesteps[2:nrow(timesteps), ]
  
  map_dfr(1:num_steps,
          ~data.frame(pragmaid = id,
                      date_timestep = ts.x[.x],
                      date = seq.Date(ts.x[.x],
                                      ts.y[.x],
                                      by = 1)[1:length(seq.Date(ts.x[.x],ts.y[.x], by = 1)) -1], 
                      rel_time = .x
          )
  ) %>% 
    return() 
  
}


  
  
  
upsample_data_to_daily = function(dat,
                                  variables_span,
                                  variables_point,
                                  period_duration = "24m",
                                  resol = c("day","week", "month", "bimonth", "quarter", "halfyear", "year")
                                  ){
  
 # Since I want to downsample to days and then upsample to month or quarter year,
#  I create a Dataframe containing all possible days of the year once for each ID
indivs = dat$pragmaid %>% unique()
num_indivs = length(indivs)
daterange = seq(min(dat$date.period.start), max(dat$date.period.end), by = 1)
num_days = length(daterange)

# yearly = tibble(pragmaid = rep(indivs, each = num_days),
#                 date = rep(daterange, times = num_indivs)) %>%
#   left_join(dat %>% dplyr::select(pragmaid),
#             by = c("pragmaid"="pragmaid"),
#             multiple = "any")
# Split data for core distribution and run helper function over it
# This step is computationally pretty expensive. thus hsouldnt be run too often. 
# Set eval to true on top of the cell in order to run this chunk
suppressMessages({
plan(multisession, workers = parallel::detectCores()-2)
res = future_map_dfr(split(dat, dat$pragmaid),
                     ~left_join(get_person_timesteps(.x, period_duration = period_duration, resol = resol),
                                upsample_daily_by_id(.x, variables_span, variables_point),
                                by = c("pragmaid", "date")
                               ),
                     .progress =T)
})

return(res)
}

daily_to_custom_resolution = function(dat, 
                                      analysis_vars
                                      ){

  res = dat %>% 
    group_by(pragmaid, date_timestep, rel_time) %>%
    summarize(across(all_of(analysis_vars), ~any(!is.na(.x))))
  return(res)
}

################Main Function for changing data resolution###########
change_resolution = function(dat,
                             period_duration = "24m",
                             resol = c("day","week", "month", "bimonth", "quarter", "halfyear", "year"), 
                             return_raw_daily = F, 
                             analysis_vars = NULL){
  #TODO: Implement a way to specify which variables should be summarized
  #TODO: Implement a way to automaticly rejoin data with covariates if available
  resol = match.arg(resol)
  variables_span = get_span_vars(dat)
  variables_point = c(get_point_vars(dat), "aud") # need to add aud here since it is not coded as an interv in sequence data
  timestep_lookup = c(day = 365, week = 364/7, month = 12, bimonth = 6, quarter = 4, halfyear = 2, year = 1)
  allvars = c(get_span_vars(dat), get_point_vars(dat)) # all variables which are interventions
  if(is_null(analysis_vars)) analysis_vars = allvars
  
  if(return_raw_daily) {return(
    dat %>%
      widen_sequence_data(period_duration = period_duration) %>%
      upsample_data_to_daily(variables_span = variables_span,
                             variables_point = variables_point,
                             period_duration = period_duration, 
                             resol = resol))}
 
  res = dat %>%
    widen_sequence_data(period_duration = period_duration) %>% 
    upsample_data_to_daily(variables_span = variables_span,
                           variables_point = variables_point,
                           period_duration = period_duration, 
                           resol = resol) %>% 
    daily_to_custom_resolution(analysis_vars = analysis_vars) %>% 
    mutate(
        year = year(date_timestep),
        .after = "date_timestep") %>%
    group_by(pragmaid, rel_time) %>%
    mutate(interv.any = any(pick(all_of(allvars)))) %>%
    group_by(pragmaid) %>%
    mutate(interv.any.allquarters = any(interv.any)) %>% 
    rename_with(~str_replace(.x, "timestep", resol) ,.cols = all_of(c("date_timestep")))
  
   return(res)
}


###############################Treatment Series Extraction##########################################


################### helper functions ################
# a problem with clustering curve shapes is, that line distances between curves are computed using time coordinate. Since the temporal differences aren't of interest here, 
# patient timelines need to be lagged, such that they start with first intervention or diagnosis.
# helper function to extract only quarters after first indicator event -> or other treatment event, depending on what variable comes after quarter$var in function
# modified function to only give data followup = x quarters after indicator event
extract_treatment = function(person,
                             followup = 8, # this is dependant on resol. If resol is quarter, 8 = 2 years. if day, 8 = 8 days
                             resol = c("day","week", "month", "bimonth", "quarter", "halfyear", "year"), 
                             event = c("aud", "interv.any")){
  resol = match.arg(resol)
  varsym = rlang::sym(resol)
  event = match.arg(event)
  num_timesteps = nrow(person)
  timesteps = person[resol] %>% pull()
  for(i in 1:num_timesteps){
    if(i > num_timesteps - followup) return(data.frame())
    timestep = person %>% filter(!!varsym == timesteps[i])
    if(timestep[event] %>% pull()) return(person[i:(i+followup), ]) # i:num_quarters returns whole treatment from first qwt, i:i+followup returns followup timepoints after first indicator event
  }
}

# function needs to know the resolution of the data which was chosen in change_resolution()
extract_treatment_series = function(dat,
                                    followup = 8, # this is dependant on resol. If resol is quarter, 8 = 2 years. if day, 8 = 8 days
                                    resol = c("day","week", "month", "bimonth", "quarter", "halfyear", "year"), 
                                    event = c("aud", "interv.any")){
  resol = match.arg(resol)
  event = match.arg(event)
  ids = dat %>%
    dplyr::select(pragmaid) %>%
    unique() %>%
    pull() 
  
  plan(multisession, workers = parallel::detectCores()-2)
  res = future_map_dfr(dat %>% split(dat$pragmaid),
                       ~extract_treatment(.x,
                                           followup = followup,
                                           resol = resol,
                                           event = event), 
                       .progress =T)
  
  # add relative time
  res = res %>% 
    group_by(pragmaid) %>% 
    mutate(rel_time = 1:n())
  
  return(tibble(res))
}



################################################MASTER FUNCTION ##########################################
get_treatment_series_at_resolution = function(dat, # data as in sequence_data_for_kilian_xyz (must contain aud)
                                              period_duration = "24m", # what observation period should be taken (in months, must exist in dat$period)? 
                                              resol = c("day","week", "month", "bimonth", "quarter", "halfyear", "year"), # what timestep?
                                              event = c("aud", "interv.any"), # what event should be used to start time series?
                                              followup = 1,# how many timesteps should be included in time series after event? must be lower than period. has the unit of resol, thus must be converted into months to match with perio # 8 quarters = 24months
                                              save = c(T,F), # saves result to current directory
                                              load_if_possible = c(T,F), # looks if function can find a saved version with similar configurations to load instead of computing it
                                              render = c(T,F), # when rendering user input cant be taken. If True and finds saved data version , it loads it
                                              only_intervs = c(T,F), # only includes individuals who had an intervention within observational period.
                                              analysis_vars = NULL) { 
  resol = match.arg(resol)                                   
  event = match.arg(event)
  
  
    if(save | load_if_possible){
    save_load_dir = here("output", "treatment_series_at_resolution_cache")
    fs::dir_create(save_load_dir)
    save_name = paste0("treatment_series_",
                       period_duration,"_",
                       resol,"_",
                       "followup_",followup,"_",
                       event,"_",
                       ifelse(only_intervs, yes = "only_intervs_", no = ""),
                       today(),
                       ".csv")
    
    look_name = paste0("treatment_series_",
                       period_duration,"_",
                       resol,"_",
                       "followup_",
                       followup,"_",
                       event,"_",
                       ifelse(only_intervs, yes = "only_intervs_", no = "")
                       )
    
    match_files = fs::dir_ls(here(save_load_dir))
    if(!is_empty(match_files)) match_files = max(match_files[str_detect(match_files, look_name)])
    match_files =  na.omit(match_files)
    if(load_if_possible & !is_empty( match_files)){
      
      if(!render) inp = readline(prompt = paste0("Found stored precomputed file with same config:\n ", match_files[1], "\nshould the file be loaded instead? (Y/N)"))
      else inp = "y"
      
      if(!(tolower(inp) %in%c("y","n"))) stop("please provide a valid answer (Y -Yes | N -No)")
      if(tolower(inp) == "y") return(vroom::vroom(match_files))
    }
  }

  period_num = as.numeric(str_extract(period_duration, "\\d+"))
  period_lookup = c(day = 30.5,
                    week = 4,
                    month = 1,
                    bimonth = 1/2,
                    quarter = 1/3,
                    halfyear = 1/6,
                    year = 1/12)
  if(followup > period_lookup[resol] * period_num) stop("Followup must encode for less or equal time as contained in period")
  
  
  if(only_intervs) dat = dat %>%  group_by(pragmaid) %>% mutate(interv.any.overall = any(interv.any)) %>% filter(interv.any.overall ==T)
  res = change_resolution(dat, period_duration = period_duration, resol = resol, analysis_vars = analysis_vars) #%>% 
    #extract_treatment_series(resol = resol, event = event, followup = followup)
  if(save){vroom::vroom_write(res, here(save_load_dir, save_name), delim = ",")}
  return(res)
  
}


###################################Venn Diagram Plotting Function###########################

# Takes a dataframe "dat" with the structure of dat_qwt but with an added column "predclass" containing class membership predictions for every individual. 
# varname is the variable (one of analysis_vars) which should be plotted.
# class is the predclass which should be plotted
# Watch it, this function RECURSES SLOWLY! When using a lot of rel_time steps (high time resolution) this function might take forever since it computes all intersections of all variable combinations (up to third degree of crossing). 
# Should only be used with less than 4 rel_time steps. 
# Alternatively define "limit" to be the maximum rel_time which should be included in the calculation and plot.
# showplot = T makes the function plot. Anyways a plottable (venneuler(obj)) object will be returned
venngardium_leviosa = function(dat, varname, class, limit = NULL, showplot = F){
  
  get_vecs = function(dat, varname, class){
    res = list()
    var = rlang::sym(varname)
    for(step in unique(dat$rel_time)){
      res[[paste0(varname, "_",step)]] = dat %>%
        filter(rel_time == step, !!var == T, predclass == class) %>% 
        dplyr::select(pragmaid) %>%
        pull()
    }
    return(res)
  }
  
  get_vec_intersect = function(vecs, varname){
    e = environment()
    pe = parent.env(e)
    lookup = get("lookup", envir = pe)
    res = data.frame(var = varname)
    current = vecs[[1]]
    current_name = names(vecs)[1]
    vecs = vecs[-1]
    vecnames = names(vecs)
    res[current_name] = length(current)
    lookup[[current_name]] = current
    if(!is_empty(vecs)){
      for(i in 1:length(vecs)){
        #     print(current)
        #      print(vecs[[i]])
        #      print(intersect(current, vecs[[i]]))
        resname = paste0(current_name, "&", vecnames[i])
        res[resname] = length(intersect(current, vecs[[i]]))
        
        lookup[[resname]] = intersect(current, vecs[[i]])
        
      }
    }
    
    
    if(!is_empty(vecs)) {
      assign("lookup", c(get("lookup", envir = pe),lookup), envir = pe)
      res = suppressMessages(left_join(res, Recall(vecs, varname)))
    }
    
    else {
      assign("lookup", c(get("lookup", envir = pe),lookup), envir = pe)
      return(res)
    }
  }
  lookup = list()
  if(!is_null(limit)) dat = dat %>% filter(rel_time <= limit)
  vecs = get_vecs(dat, varname, class)
  ress = get_vec_intersect(vecs, varname)
  lookup = lookup[!duplicated(lookup)]
  nams = names(ress)[-1]
  nams_orig = nams[nams %in% names(vecs)]
  for(i in 1:length(nams_orig)){
    current_name = nams_orig[i]   
    # print(paste("Current", current_name))
    
    nams_missing_comb = nams[str_detect(nams, "&") &
                               !str_detect(nams, current_name)]
    
    #print(paste("Missing", nams_missing_comb))
    if(!is_empty(nams_missing_comb)){
      for(j in 1:length(nams_missing_comb)){
        missing_val = intersect(lookup[[nams_missing_comb[j] ]],
                                lookup[[current_name]])
        # print(paste("VAL", missing_val))
        lookup[[paste0(nams_missing_comb[j],"&",current_name)]] = missing_val
        ress[paste0(nams_missing_comb[j],"&",current_name)] = length(missing_val)
        
      }
    }
    
  }
  
  res = as.numeric(ress[1,-1])
  names(res) = names(ress)[-1]
  if(showplot) plot(venneuler(res), main = paste0(varname, " Class ", class ))
  return(res)
}



# MODEL REPORTING CONVENIENCE FUNCTIONS 
get_entropy = function(m, digits = 2){
  error_prior <- entropy(m$P) # Class proportions
  error_post <- mean(apply(m$posterior, 1, entropy))
  R2_entropy <- (error_prior - error_post) / error_prior
  return(round(R2_entropy, 2))
}

smallest_class = function(predclass, type = c("n", "rel")){
  type = match.arg(type)
  res = data.frame(predclass = predclass) %>%
    group_by(predclass) %>% 
    summarize(n = n()) %>% 
    ungroup() %>% 
    mutate(rel = n/sum(n)) %>% 
    dplyr::select(all_of(type)) %>% min()
  if(type == "rel") res = round(res*100, 1)
  return(res)
}

smallest_lcprob = function(m){
  
  res = data.frame(predclass = m$predclass)
  lc = n_distinct(res$predclass)
  for(i in 1:lc){ 
    resres = data.frame(x = m$posterior[ , i])
    names(resres) = paste0("LC", i)
    res = cbind(res, resres)
  }
  
  res %>%
    group_by(predclass) %>% 
    summarize(across(starts_with("LC"), mean)) %>% 
    pivot_longer(cols = starts_with("LC"),
                 names_to = "Class",
                 values_to = "prob") %>% 
    group_by(predclass) %>% 
    summarize(maxprob = max(prob)) %>% 
    ungroup() %>% 
    mutate(maxprob = round(maxprob*100, 1)) %>% 
    filter(maxprob == min(maxprob)) %>% 
    dplyr::select(maxprob) %>% pull()
  
}

mean_lcprob = function(m){
  res = data.frame(predclass = m$predclass)
  lc = n_distinct(res$predclass)
  for(i in 1:lc){ 
    resres = data.frame(x = m$posterior[ , i])
    names(resres) = paste0("LC", i)
    res = cbind(res, resres)
  }
  
  res %>%
    group_by(predclass) %>% 
    summarize(across(starts_with("LC"), mean)) %>% 
    pivot_longer(cols = starts_with("LC"),
                 names_to = "Class",
                 values_to = "prob") %>% 
    group_by(predclass) %>% 
    summarize(maxprob = max(prob)) %>% 
    ungroup() %>% 
    mutate(maxprob = round(maxprob*100, 1)) %>% 
    summarize(meanprob = mean(maxprob)) %>% 
    dplyr::select(meanprob) %>% 
    pull()
  
}


get_entropy = function(m){
  error_prior <- entropy(m$P) # Class proportions
  error_post <- mean(apply(m$posterior, 1, entropy))
  R2_entropy <- (error_prior - error_post) / error_prior
  return(round(R2_entropy, 3))
}
entropy <- function(p) sum(-p * log(p), na.rm =T) 

