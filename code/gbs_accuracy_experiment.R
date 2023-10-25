library(data.table)
library(glue)
library(parallel)
library(doParallel)
library(Rfast) # provides a (very!) efficient implementation of the NxN distance calculation
library(xtable)

###
# This script produces Table 1, while compares the accuracy of GBS thresholding to physicians' discretionary decisions #
###

source('./helpers.R')
source('./logic.R')
source('./gib_experiment_helpers.R')

use.saved.data <- TRUE

seed <- 1
set.seed(seed)

d2 <- read_csv('../data/gib_score_only.csv')
mean(d2$y)

acc.human <- printable_confusion_matrix(d2) %>% mutate(decision.rule = 'Physician Discretion')

raw.acc.gbs <- lapply(0:22, function(i) {
  calc_confusion_matrix(d2 %>% mutate(y_hat = total_score > i)) %>% mutate(threshold = i)
}) %>% rbindlist %>% as_tibble

acc.gbs <- lapply(0:22, function(i) {
  printable_confusion_matrix(d2 %>% mutate(y_hat = total_score > i)) %>% mutate(threshold = i)
}) %>% rbindlist %>% as_tibble


acc.gbs.filtered <- (acc.gbs %>% filter(threshold <= 2)) %>% rbind(acc.gbs[which.max(raw.acc.gbs$acc), ]) %>%
  mutate(decision.rule = map_chr(threshold, function(t) glue('Admit GBS > {t}'))) %>%
  select(-threshold)

table1 <- rbind(acc.human, acc.gbs.filtered)
table1 <- table1[, c('decision.rule', 'admitted.frac', 'acc', 'tpr', 'tnr')]

colnames(table1) <- c('Decision Rule', 'Fraction Hospitalized', 'Accuracy', 'Sensitivity', 'Specificity')

hline <- c(-1,0,nrow(table1))
htype <- c("\\toprule ", "\\midrule ","\\bottomrule ")

caption <- '
Comparing the accuracy of physician hospitalization decisions (`Physician Discretion\') to those made by thresholding the GBS.
For example, `Admit GBS $> 1$\' hospitalizes patients with a GBS strictly larger than $1$.  
`Fraction Hospitalized\' indicates the fraction of patients hospitalized by each rule.
`Accuracy\' indicates the 0/1 accuracy of each rule, where a decision is correct if it hospitalizes a patient who suffers an adverse outcome (as defined above) or discharges a patient who does not.
`Sensitivity\' indicates the fraction of patients who suffer an adverse outcome that are correctly hospitalized,
and `Specificity\' indicates the fraction of patients who do not suffer an adverse outcome that are correctly discharged.
Results are reported to $\\pm 2$ standard errors. 
'

print(xtable(table1, 
             type = "latex",
             caption = caption,
             label = 'tab:physician and gbs perf',
             auto = TRUE,
             digits = 2),
      include.rownames = FALSE, table.placement = '!htbp', caption.placement = 'bottom', sanitize.text.function=function(x){x},
      add.to.row = list(pos = as.list(hline),
                        command = htype),  
      hline.after = NULL,
      file = '../tables/human_gbs_accuracy.tex')


