library(data.table)
library(glue)
library(parallel)
library(doParallel)
library(Rfast) # provides a (very!) efficient implementation of the NxN distance calculation
library(xtable)

###
# This script produces Table 2, 3 and 4, which test for physician expertise in various feature spaces
# Which table is produced depends on the 'input_file' chosen in the SETUP and HYPERPARAMS section below
# Note that only the 'gib_score_only.csv' file is available in the public Github repository
# The other two files, which contain a richer set of confidential patient data, were made available to reviewers
# We will do our best to accommodate individual requests for this data (each request is subject to IRB approval)
# It also produces a plot of the pairwise Euclidian distance between each pair of patients chosen for the test
###

source('./helpers.R')
source('./logic.R')
source('./gib_experiment_helpers.R')

### SETUP and HYPERPARAMS ###

use.saved.data <- TRUE
save.tables <- TRUE
include.distance <- TRUE

#input_file <- 'gib_raw_values_normalized.csv'
input_file <- 'gib_gbs_clean_inputs_normalized.csv'
#input_file <- 'gib_score_only.csv'

suffix = ''
title = ''
feature_space_desc = ''

if (input_file == 'gib_raw_values_normalized.csv') {
  suffix = '_raw'
  title <- 'Pairwise Distances in Original Feature Space'
  feature_space_desc = 'nine patient characteristics. The Glasgow-Blatchford score is computed from these nine characteristics, after a pre-processing step which discretizes each feature.'
} else if (input_file == 'gib_gbs_clean_inputs_normalized.csv') {
  suffix = '_clean'
  title <- 'Pairwise Distances in Discretized Feature Space'
  feature_space_desc = 'the nine (discrete) inputs to the Glasgow-Blatchford score.'
} else if (input_file == 'gib_score_only.csv') {
  title <- 'Pairwise Distances in GBS Space'
  suffix = '_score_only'
  feature_space_desc = 'their Glasgow-Blatchford scores.'
} else {
  stop('unrecognized input file!')
}

seed <- 1
set.seed(seed)

### LOGIC ###
if (!use.saved.data) {
  if (!file.exists(input_file)) {
    stop(glue('File: {input_file} does not exist! Only the gib_score_only.csv is available in the public repository (see README)'))
  }
  d2 <- read_csv(glue('../data/{input_file}'))
  
  if (input_file != 'gib_score_only.csv') {
    d2 <- d2 %>% select(-total_score)  
  }
  max.L <- floor(nrow(d2) / 2)
}

if (!file.exists(glue('../data/res.100{suffix}.RData')) || !use.saved.data) {
  res.100 <- generate_p_value(d2, L = 100, K = 1000, loss_type = '01', return.diagnostics = T)
  save(res.100, file = glue('../data/res.100{suffix}.RData'))
}

if (!file.exists(glue('../data/res.250{suffix}.RData')) || !use.saved.data) {
  res.250 <- generate_p_value(d2, L = 250, K = 1000, loss_type = '01', return.diagnostics = T)
  save(res.250, file = glue('../data/res.250{suffix}.RData'))
}

if (!file.exists(glue('../data/res.500{suffix}.RData')) || !use.saved.data) {
  res.500 <- generate_p_value(d2, L = 500, K = 1000, loss_type = '01', return.diagnostics = T)
  save(res.500, file = glue('../data/res.500{suffix}.RData'))
}

if (!file.exists(glue('../data/res.1000{suffix}.RData')) || !use.saved.data) {
  res.1000 <- generate_p_value(d2, L = 1000, K = 1000, loss_type = '01', return.diagnostics = T)
  save(res.1000, file = glue('../data/res.1000{suffix}.RData'))
}

if (!file.exists(glue('../data/res.max{suffix}.RData')) || !use.saved.data) {
  res.max <- generate_p_value(d2, L = max.L, K = 1000, loss_type = '01', return.diagnostics = T)
  save(res.max, file = glue('../data/res.max{suffix}.RData'))
}

load(glue('../data/res.100{suffix}.RData'))
load(glue('../data/res.250{suffix}.RData'))
load(glue('../data/res.500{suffix}.RData'))
load(glue('../data/res.1000{suffix}.RData'))
load(glue('../data/res.max{suffix}.RData'))

max.L <- res.max$paired %>% nrow

table <- print_table(res.100, 100) %>%
  rbind(print_table(res.250, 250)) %>%
  rbind(print_table(res.500, 500)) %>%
  rbind(print_table(res.1000, 1000)) %>%
  rbind(print_table(res.max, max.L)) %>%
  mutate(tau = map_chr(tau, function(t) if (t < .001) glue('<.001') else round(t,3) %>% as.character))

colnames(table) <- c('L',
                     'mismatched pairs',
                     'swaps that increase loss',
                     'swaps that decrease loss', 
                     '$\\tau$')

hline <- c(-1,0,nrow(table))
htype <- c("\\toprule ", "\\midrule ","\\bottomrule ")

label <- str_c('tab:testing for expertise ', substr(suffix, 2, nchar(suffix)))

if (save.tables) {
  caption <- str_c(
    'The results of running ExpertTest, where each pair of patients is chosen to be as similar as possible with respect to ', 
    feature_space_desc)
  
  caption <- str_c(caption,
                   ' $L$ indicates the number of pairs selected for the test, of which `mismatched pairs\' are not identical to each other.
                   Swaps that decrease (respectively, increase) loss indicates how many of the $L$ pairs result in a decrease (respectively, increase) in the 0/1 loss
                   when their corresponding hospitalization decisions are exchanged with each other. $\\tau$ is the p-value obtained
                   from running ExpertTest.')
  
  print(xtable(table, 
               type = "latex",
               caption = caption,
               label = label,
               auto = TRUE),
        include.rownames = FALSE, table.placement = '!htbp', caption.placement = 'bottom', sanitize.text.function=function(x){x},
        add.to.row = list(pos = as.list(hline),
                          command = htype),
        hline.after = NULL,
        file = glue('../tables/testing_expertise{suffix}.tex'))
}


fig <- collect_distances(res.100, 100) %>%
  rbind(collect_distances(res.250, 250)) %>%
  rbind(collect_distances(res.500, 500)) %>%
  rbind(collect_distances(res.1000, 1000)) %>%
  rbind(collect_distances(res.max, max.L))

p1 <- fig %>%
  mutate(L = as.factor(L)) %>%
  ggplot(aes(x = distances, y = L, fill = L)) +
  geom_boxplot(aes(group = L)) +
  geom_vline(xintercept = sqrt(9), color = 'red') +
  theme_bw() +
  xlab('pairwise euclidian distance') +theme(title = element_text(size = 25),
                                             axis.text=element_text(size=25),
                                             axis.title=element_text(size=14,face="bold"),
                                             legend.title = element_text(size = 25),
                                             legend.text = element_text(size = 15)) + ggtitle(title)

ggsave(glue('../figures/pairwise_distances{suffix}.png'), p1, width = 12, height = 6)


