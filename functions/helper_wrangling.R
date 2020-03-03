
###'######################################################################
###'
###' Helper functions for data analysis
###' 
###' 
###' 20180901 JoonHo Lee
###' 
###' 

### Package dependency
library(tidyverse)
library(scales)



###'######################################################################
###'
###' tabdf(): Tabulate frequencies
###' 
###' Similar to the Stata command "tab"
###'
###'

tabdf <- function(df, 
                  variable){
  
  ### Enquote variables
  x <- enquo(variable)
  
  ### Generate table
  tibble_tbl <- df %>%
    group_by(!!x) %>%
    summarise(Freq = n()) %>%
    ungroup() %>%
    mutate(total_n = sum(Freq, na.rm = TRUE),
           Percent = round((Freq/total_n)*100, 1), 
           CumFreq = cumsum(Freq), 
           CumPercent = round((CumFreq/total_n)*100, 1)) %>%
    dplyr::select(-total_n) 
  
  ### Display table as data.frame format
  data.frame(tibble_tbl)
}



###'######################################################################
###'
###' classmode(): Check classes and modes of selected variables
###'
###'

classmode <- function(df, ...){
  
  ### Enquote variables
  vars <- quos(...)  # any rules for dplyr::select() works. ex) everything(), ends_with(), etc.
  
  ### Select variables
  df_select <- df %>% 
    dplyr::select(!!!vars)
  
  ### Return classes and modes
  mat <- cbind(sapply(df_select, class), 
               sapply(df_select, mode))
  
  ### Convert to data.frame format
  df_mat <- data.frame(rownames(mat), mat)
  rownames(df_mat) <- NULL
  names(df_mat) <- c("variable", "class", "mode")
  
  return(df_mat)
}



###'######################################################################
###'
###' listvars(): list selected variables
###'
###'

listvars <- function(df, ..., nrow = 100){
  
  ### Enquote variables
  vars <- quos(...)  # any rules for dplyr::select() works. ex) everything(), ends_with(), etc.
  
  ### Select variables
  df_select <- df %>% 
    dplyr::select(!!!vars)
  
  ### Return rows for the selected variables
  df_select[1:nrow, ]

}



###'######################################################################
###'
###' empty_as_na(): Convert empty strings to NA
###' 
###' : it's important to include "trimws" to remove blank spaces 
###'
###'

empty_as_na <- function(x){
  
  if("factor" %in% class(x)) x <- as.character(x) 
  
  ## since ifelse wont work with factors
  
  ifelse(trimws(as.character(x)) !="", x, NA)
}



###'######################################################################
###'
###' operation14(): Remove districts with insufficient years of data
###' 
###' => Analyze only traditional schools in elementary, high, and unified 
###'    school districts that have been in continuous operation (14 years) 
###'    in California from 2003 through 2017

setwd(work_dir)
load("~/SACS/processed_data/years_of_operation.rda")   # Data dependency: years_of_operation.csv

operation14 <- function(df){
  df %>%
    left_join(years_of_operation[, !names(years_of_operation) %in% c("Dname", "Dtype")], 
              by = c("Ccode", "Dcode")) %>%
    filter(opr_years == 14)
}



###'######################################################################
###'
###' get_weighted_mean(): Get weighted district averages
###' 
###' 

get_weighted_mean <- function(df, 
                              x = Fiscalyear, 
                              y = sum_value_PP_16, 
                              weight = K12ADA_C, 
                              ...){
  
  ### Enquote variables
  x <- enquo(x)
  y <- enquo(y)
  weight <- enquo(weight)
  group_var <- quos(...)
  
  df %>% 
    group_by(!!x, !!!group_var) %>%
    summarise(mean_value = round(weighted.mean(!!y, !!weight, na.rm = TRUE), 2))
}



###'######################################################################
###'
###' Calculate percentages based on groups
###'
###'

group_percent <- function(df, 
                          value = mean_value, 
                          ...){
  
  # Enquote variables
  value <- enquo(value)
  group <- quos(...)
  
  
  #' (1) Calculate the percentages
  #' (2) Format the labels and calculate their positions
  df %>%
    group_by(!!!group) %>%
    mutate(group_sum = sum(!!value, na.rm = TRUE), 
           percent = !!value/group_sum * 100, 
           # don't need to calculate the label positions from ggplot 2.1.0 
           # position = cumsum(amount) - 0.5 * amount,  
           label_text = paste0(sprintf("%.1f", percent), "%")) -> df
  return(df)
}



###'######################################################################
###'
###' is.numeric.elementwise(): 
###' Check whether each element is numeric
###'
###'

is.numeric.elementwise <- function(vector){
  
  lvector <- c()
  
  for (i in seq_along(vector)){
    
    element <- vector[i]
    
    elem_TF <- is.na(as.numeric(element))
    
    lvector <- c(lvector, elem_TF)
    
  }
}



###'######################################################################
###'
###' Get a dataframe of regression estimates 
###'
###'

get_lm_est_df <- function(lm_fit){
  
  summary <- summary(lm_fit)
  
  df <- data.frame(summary$coefficients)
  
  names(df) <- c("estimate", "std_error", "t_value", "p_value")
  
  df <- round(df, 3)
  
  df$variable <- rownames(df)
  rownames(df) <- NULL
  
  df <- df %>% 
    dplyr::select(variable, everything())
  
  df
}


###'######################################################################
###'
###' school_composition_count()
###' 
###' => Generate table for summarizing school-level compositions
###' => Based on "COUNTING" variable
###'
###'

school_composition_count <- function(df, 
                                     var_to_count, 
                                     factor, 
                                     levels_to_replace = NULL, 
                                     table_name = "Staff Type", 
                                     year = NULL){
  
  ### Enquote variables
  var_to_count <- enquo(var_to_count)
  factor <- enquo(factor)
  
  ### Calculate Subtotal
  df_subtotal <- df %>%
    group_by(CountyCode, DistrictCode, SchoolCode) %>%
    summarise(subtotal = n_distinct(!!var_to_count)) 
  
  ### Calculate Subgroup Counts
  df_by_subgroup <- df %>% 
    group_by(CountyCode, DistrictCode, SchoolCode, 
             !!factor) %>% 
    summarise(N = n_distinct(!!var_to_count)) 
  
  ### Calculate Subgroup Percentages
  df_by_subgroup <- df_by_subgroup %>%
    left_join(df_subtotal, 
              by = c("CountyCode", "DistrictCode", "SchoolCode")) %>%
    mutate(PCT = 100*(N/subtotal))
  
  ###' Reshape from long to wide data format 
  ###' (1) Long data format to spread multiple variables
  df_temp <- df_by_subgroup %>%
    gather(stat, value, N, PCT)  
  
  ###' (2) Recode factor levels to brief names
  df_temp <- df_temp %>%
    rename(factor = !!factor)

  levels(df_temp$factor) <- levels_to_replace
  
  ###' (3) Reshape to wide format 
  df_temp <- df_temp %>%
    unite(key, stat, factor) %>%
    spread(key = key, value = value, fill = 0)
  
  ### Add table name & AcademicYear
  df_temp <- df_temp %>%
    mutate(table = table_name, 
           AcademicYear = as.numeric(year))
  
  ### Reorder variables by factor levels
  N_vars <- paste0("N_", levels_to_replace)
  PCT_vars <- paste0("PCT_", levels_to_replace)
  
  N_vars_select <- N_vars[N_vars %in% names(df_temp)]
  PCT_vars_select <- PCT_vars[PCT_vars %in% names(df_temp)]
  
  df_temp <- df_temp %>%
    dplyr::select(ends_with("Code"),
           table, AcademicYear, subtotal,
           N_vars_select, 
           PCT_vars_select)
  
  return(df_temp)
}


###'######################################################################
###'
###' school_composition_sum()
###' 
###' => Generate table for summarizing school-level compositions
###' => Based on "SUMMATION" of variable
###'
###'

school_composition_sum <- function(df, 
                                   var_to_sum, 
                                   factor, 
                                   levels_to_replace = NULL, 
                                   table_name = "Staff Type", 
                                   year = NULL){
  
  ### Enquote variables
  var_to_sum <- enquo(var_to_sum)
  factor <- enquo(factor)
  
  ### Calculate Subtotal
  df_subtotal <- df %>%
    group_by(CountyCode, DistrictCode, SchoolCode) %>%
    summarise(subtotal = sum(!!var_to_sum, na.rm = TRUE)) 
  
  ### Calculate Subgroup Counts
  df_by_subgroup <- df %>% 
    group_by(CountyCode, DistrictCode, SchoolCode, 
             !!factor) %>% 
    summarise(N = sum(!!var_to_sum, na.rm = TRUE)) 
  
  ### Calculate Subgroup Percentages
  df_by_subgroup <- df_by_subgroup %>%
    left_join(df_subtotal, 
              by = c("CountyCode", "DistrictCode", "SchoolCode")) %>%
    mutate(PCT = 100*(N/subtotal))
  
  ###' Reshape from long to wide data format 
  ###' (1) Long data format to spread multiple variables
  df_temp <- df_by_subgroup %>%
    gather(stat, value, N, PCT)  
  
  ###' (2) Recode factor levels to brief names
  df_temp <- df_temp %>%
    rename(factor = !!factor)
  
  levels(df_temp$factor) <- levels_to_replace
  
  ###' (3) Reshape to wide format 
  df_temp <- df_temp %>%
    unite(key, stat, factor) %>%
    spread(key = key, value = value, fill = 0)
  
  ### Add table name & AcademicYear
  df_temp <- df_temp %>%
    mutate(table = table_name, 
           AcademicYear = as.numeric(year))
  
  ### Reorder variables by factor levels
  N_vars <- paste0("N_", levels_to_replace)
  PCT_vars <- paste0("PCT_", levels_to_replace)
  
  N_vars_select <- N_vars[N_vars %in% names(df_temp)]
  PCT_vars_select <- PCT_vars[PCT_vars %in% names(df_temp)]
  
  df_temp <- df_temp %>%
    dplyr::select(ends_with("Code"),
           table, AcademicYear, subtotal,
           N_vars_select, 
           PCT_vars_select)
  
  return(df_temp)
}




###'######################################################################
###'
###' school_summarize()
###' 
###' => Generate descriptive statistics table 
###'    for summarizing distribution of counituous variable
###'
###'

school_summarize <- function(df, 
                             var_to_summarize, 
                             table_name, 
                             year = NULL){
  
  ### Enquote variables
  var_to_summarize <- enquo(var_to_summarize)
  
  ### Generate summary table
  df_sum <- df %>%
    group_by(CountyCode, DistrictCode, SchoolCode) %>%
    summarise(N = sum(!is.na(!!var_to_summarize)),  
              mean = mean(!!var_to_summarize, na.rm = TRUE), 
              median = median(!!var_to_summarize, na.rm = TRUE), 
              sd = sd(!!var_to_summarize, na.rm = TRUE), 
              min = min(!!var_to_summarize, na.rm = TRUE), 
              max = max(!!var_to_summarize, na.rm = TRUE))
  
  ### Add table name and year & reorder variables
  df_sum <- df_sum %>%
    mutate(table = table_name, 
           AcademicYear = as.numeric(year)) %>%
    dplyr::select(CountyCode:SchoolCode, table, AcademicYear, 
           everything())

  return(df_sum)
}



###'######################################################################
###'
###' dualsave()
###' 
###' => save dataframe in both .rda and .dta formats
###'
###'

dualsave <- function(df_to_save, 
                     filename = "temp_name"){
  
  library(foreign)
  save(df_to_save, file = paste0(filename, ".rda"))
  write.dta(df_to_save, file = paste0(filename, ".dta"))
  
}



###'######################################################################
###'
###' full_join_track()
###' 
###' => add "_merge" variables like Stata or Python
###' 
###' 


full_join_track <- function(x, y, by = NULL, suffix = c(".x", ".y"),
                            .merge = FALSE, ...){
  
  # Checking to make sure used variable names are not already in use
  if(".x_tracker" %in% names(x)){
    message("Warning: variable .x_tracker in left data was dropped")
  }
  if(".y_tracker" %in% names(y)){
    message("Warning: variable .y_tracker in right data was dropped")
  }
  if(.merge & (".merge" %in% names(x) | ".merge" %in% names(y))){
    stop("Variable .merge already exists; change name before proceeding")
  }
  
  # Adding simple merge tracker variables to data frames
  x[, ".x_tracker"] <- 1
  y[, ".y_tracker"] <- 1
  
  # Doing full join
  joined <- full_join(x, y, by = by, suffix = suffix,  ...)
  
  # Calculating merge diagnoses 
  matched <- joined %>%
    filter(!is.na(.x_tracker) & !is.na(.y_tracker)) %>%
    NROW()
  unmatched_x <- joined %>%
    filter(!is.na(.x_tracker) & is.na(.y_tracker)) %>%
    NROW()
  unmatched_y <- joined %>%
    filter(is.na(.x_tracker) & !is.na(.y_tracker)) %>%
    NROW()
  
  # Print merge diagnoses
  message(
    unmatched_x, " Rows ONLY from left data frame", "\n",
    unmatched_y, " Rows ONLY from right data frame", "\n",
    matched, " Rows matched"
  )
  
  # Create .merge variable if specified
  if(.merge){
    joined <- joined %>%
      mutate(.merge = 
               case_when(
                 !is.na(.$.x_tracker) & is.na(.$.y_tracker) ~ "left_only",
                 is.na(.$.x_tracker) & !is.na(.$.y_tracker) ~ "right_only",
                 TRUE ~ "matched"
               )
      )
  }
  
  # Dropping tracker variables and returning data frame
  joined <- joined %>%
    dplyr::select(-.x_tracker, -.y_tracker)
  return(joined)
}



###'######################################################################
###'
###' tag_ID_miss_dups()
###' 
###' => Tag missing and duplicated IDs
###'

tag_ID_miss_dups <- function(df, 
                             match_key = NULL){
  
  ### Tag cases with missing IDs
  idx_complete <- df %>%
    dplyr::select(match_key) %>%
    complete.cases()
  
  idx_miss <- !idx_complete
  
  
  ### Tag cases with duplicated IDs
  idx_duplicated <- df %>%
    dplyr::select(match_key) %>%
    duplicated()
  
  
  ### Get a dataframe of duplicated IDs
  tbl_duplicated <- df[idx_duplicated, ] %>%
    dplyr::select(match_key) %>%
    distinct() %>%
    mutate(tag = "ID_duplicated")
  
  
  ### Merge with the original dataset to filter all cases with duplicated IDs
  df_miss_dups <- df %>%
    left_join(tbl_duplicated, by = match_key)
  
  
  ### Tag missing IDs
  df_miss_dups$tag[idx_miss] <- "ID_missing"
  
  
  ### Tag unduplicated IDs
  df_miss_dups$tag[is.na(df_miss_dups$tag)] <- "ID_valid"
  
  
  ### Print diagnoses
  n_tag <- sum(!is.na(df_miss_dups$tag))
  n_miss <- sum(idx_miss)
  n_dups <- sum(df_miss_dups$tag == "ID_duplicated")
  n_valid <- sum(df_miss_dups$tag == "ID_valid")

  
  message(
    paste0(n_miss, " (", round(100*n_miss/n_tag, 1), "%)"), 
    " Rows have MISSING IDs", "\n",
    paste0(n_dups, " (", round(100*n_dups/n_tag, 1), "%)"), 
    " Rows have DUPLICATED IDs", "\n", 
    paste0(n_valid, " (", round(100*n_valid/n_tag, 1), "%)"), 
    " Rows have VALID IDs to use in matching"
  )
  
  return(df_miss_dups)
  
}

