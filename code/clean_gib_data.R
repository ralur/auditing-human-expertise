library(tidyverse)
library(glue)

###
# This script takes the input patient characteristics and simply
  # (1) selects the columns relevant to the GBS calculation and
  # (2) normalizes each column to lie in [0, 1]

  # Note that only the 'gib_score_only.csv' output file is available in the public Github repository
  # The remaining files, which contain a richer set of confidential patient data, were made available to reviewers
  # We will do our best to accommodate individual requests for this data (each request is subject to IRB approval)
###

input_filenames <- c('gib_raw_values', 'gib_gbs_clean_inputs')

for (file in input_filenames) {
  
  d <- read_csv((glue('../data/{file}.csv')))
  if (!('total_score' %in% colnames(d))) {
    d <- d %>% mutate(total_score = 0) # raw_values does not contain total_score column, include 0s to ensure table schemas match
  }
  d %>% filter(ed_disposition == 'Other') %>% nrow # 10 patients for whom we have no admit/discharge decision; we'll remove these
  
  d <- d %>%
    filter(ed_disposition != 'Other') %>% 
    mutate(y_hat = map_dbl(ed_disposition, function(e) e == 'Admit'), y = as.numeric(outcome_composite | return_admission))
  
  d2 <- d %>% select(c(
    BUN,
    HGB,
    SBP,
    PULSE,
    cardiac_failure,
    hepatic_disease,
    melena,
    syncope,
    SEX_C,
    total_score,
    y_hat,
    y
  ))
  
  d2 <- d2 %>%
    mutate(SEX_C = SEX_C - 1) %>%
    mutate(across(-c(total_score, y, y_hat), function(col) col / max(col)))
  
  
  d2 %>% summarize(across(everything(), max))
  d2 %>% summarize(across(everything(), min))
  
  d2 %>% write_csv(glue('../data/{file}_normalized.csv'))
  
}


d2 <- read_csv('../data/gib_gbs_clean_inputs.csv')
d %>% select(total_score, y_hat, y) %>% write_csv('../data/gib_score_only.csv')
