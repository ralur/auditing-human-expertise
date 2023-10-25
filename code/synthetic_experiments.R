library(tidyverse)
library(data.table)
library(ggplot2)
library(parallel)
library(doParallel)
library(glue)

source('logic.R')
source('helpers.R')

set.seed(1)
use.saved.data <- TRUE

# TOY EXAMPLE
res <- generate_loss_comparison(1000, 100)

table3 <- data.frame(
  algorithm_mse = format_mean_sd(res$algorithm_mse),
  human_mse = format_mean_sd(res$human_mse),
  rescaled_mse = format_mean_sd(res$human_mse_rescaled)
)

colnames(table3) <- c('Algorithm MSE', 'Human MSE', 'Rescaled Human MSE')

hline <- c(-1,0,nrow(table3))
htype <- c("\\toprule ", "\\midrule ","\\bottomrule ")

print(xtable(table3, 
             type = "latex",
             caption = "Expert vs Algorithm Performance",
             label = 'tab: human v algo performance',
             auto = TRUE),
      include.rownames = FALSE, table.placement = '!htbp', caption.placement = 'top', sanitize.text.function=function(x){x},
      add.to.row = list(pos = as.list(hline),
                        command = htype),
      hline.after = NULL,
      file = '../tables/synthetic_experiment.tex')


p.value.cdf <- function(p.dist, critical_value = .05) {
  discovery_rate <- mean(p.dist <= critical_value)
  
  x <- seq(0, 1, by = .01)
  df <- data.frame(x, y = map_dbl(x, function(x) mean(p.dist <= x)))
  
  df %>%
    ggplot(aes(x, y)) +
    geom_line() +
    geom_abline(slope = 1, linetype = "dashed") +
    geom_vline(xintercept = .05, color = 'red') +
    xlab('\u03b1') +
    ylab('P(\u03c4 <= \u03b1)') +
    theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), text = element_text(size=20)) +
    theme_bw()
}

if (!file.exists('../data/synth_distribution.RData') || !use.saved.data) {
  res <- generate_distribution_synth(reps = 100, n = 1000, L = 100, K = 100, run.serial = TRUE)  
  save(res, file = '../data/synth_distribution.RData')
}
load('../data/synth_distribution.RData')

f2 <- p.value.cdf(res$p.value.h0)
ggsave('../figures/p_value_h0.png', f2 , dpi = 100, width = 8, height = 5)

f3 <- p.value.cdf(res$p.value.h1)
ggsave('../figures/p_value_h1.png', f3 , dpi = 100, width = 8, height = 5)



## POWER ## 
if (!file.exists('../data/power_by_n.csv') || !use.saved.data) {
  r1 <- generate_distribution_power(100, 200, 25, 100)
  r2 <- generate_distribution_power(100, 600, 75, 100)
  r3 <- generate_distribution_power(100, 1200, 150, 100)
  
  rbind(r1, r2, r3) %>% write_csv('../data/power_by_n.csv')
}
if (!file.exists('../data/power_by_L.csv') || !use.saved.data) {
  r4 <- lapply(1:10, function(i) generate_distribution_power(500, 600, i * 20, 100, eps = c(.2))) %>%
    rbindlist %>% as_tibble
  r4 %>% write_csv('../data/power_by_L.csv')
}

f4 <- read_csv('../data/power_by_n.csv') %>%
  mutate(n = factor(n,
                    levels = c(200, 600, 1200),
                    labels = c('n=200, L=25', 'n=600, L=75', 'n=1200, L=150')
                  )) %>%
  rename(delta = epsilon) %>%
  group_by(n, delta) %>% summarize(discovery.rate = mean(p.value <= .05)) %>%
  ggplot(aes(delta, discovery.rate, color = n)) +
    geom_line() +
    geom_hline(yintercept = .80, linetype = 'dashed') +
    theme_bw() +
    ylab('discovery rate') + xlab('\u03b4')

ggsave('../figures/power_by_n.png', f4 , dpi = 100, width = 8, height = 5)


f5 <- read_csv('../data/power_by_L.csv') %>%
  group_by(L) %>%
  summarize(p.value = mean(p.value <= .05)) %>% 
  ggplot() +
    geom_line(aes(L, p.value)) +
    theme_bw() + ylim(0, 1) +
    geom_hline(yintercept = .8, linetype = 'dashed') +
    xlab('L') + ylab('discovery rate')

ggsave('../figures/power_by_L.png', f5 , dpi = 100, width = 8, height = 5)

## EXCESS TYPE I ERROR ##

if (!file.exists('../data/error_by_L.csv') || !use.saved.data) {
  L <- map_dbl(1:10, function(i) i * 25)
  r5 <- generate_distribution_excess_error(50, 500, L, 50)
  r5 %>% write_csv('../data/error_by_L.csv')
}

f6 <- read_csv('../data/error_by_L.csv') %>%
  group_by(L) %>%
  summarize(p.value = mean(p.value <= .05)) %>% 
  ggplot() +
  geom_line(aes(L, p.value)) +
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = .05, linetype = 'dashed') +
  xlab('L') + ylab('false discovery rate')

ggsave('../figures/type_1_error_by_L.png', f6 , dpi = 100, width = 8, height = 5)


