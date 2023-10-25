print_table <- function(res, L) {
  nimp <- res$n_imperfect
  chg <- res$possible_changes_01 %>% nrow
  impr <- res$possible_improvements_01 %>% nrow
  decr <- chg - impr
  tau <- res$p_value
  
  data.frame(L = L,
             `mismatched pairs` = nimp,
             `swaps that increase loss` = decr,
             `swaps that decrease loss` = impr,
             tau = tau)
}

collect_distances <- function(res, L) {
  distances <- res$paired %>% select(-c(y, y_hat, y.2, y_hat.2)) %>%
    apply(1, function(x) sqrt(sum(((x[1:(length(x)/2)] - x[((length(x)/2)+1):length(x)])^2))))
  
  data.frame(
    L = rep(L, length(distances)),
    distances = distances
  )
}

print_distances <- function(res, L) {
  distances <- res$paired %>% select(-c(y, y_hat, y.2, y_hat.2)) %>%
    apply(1, function(x) sqrt(sum(((x[1:(length(x)/2)] - x[((length(x)/2)+1):length(x)])^2)))) 
  
  data.frame(
    L = L,
    min = min(distances),
    tenth = quantile(distances, .10) %>% as.numeric,
    median = quantile(distances, .5)%>% as.numeric,
    ninetieth = quantile(distances, .90 )%>% as.numeric,
    max = max(distances)
  ) %>% mutate(across(-L, function(x) round(x, 3)))
  
}
