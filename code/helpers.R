euclidean_distance <- function(row1, row2) {
  return(sqrt(sum((row1 - row2)^2)))
}

distance_matrix <- function(df, dims = NULL) {
  filtered <- df
  if (!is.null(dims)) {
    filtered <- df[,dims]
  }
  distance_mat <- Dist(filtered, method = "euclidean", square = FALSE, p = 0, vector = FALSE)
  diag(distance_mat) <- Inf
  distance_mat
}

permute_preds <- function(paired, ratio = 1) {
  r <- ratio / (ratio + 1)
  randomness <- runif(nrow(paired))
  to_swap <- randomness > r
  
  out <- paired
  tmp <- out[to_swap,'y_hat']
  
  out[to_swap,'y_hat'] <- out[to_swap,'y_hat.2']
  out[to_swap,'y_hat.2'] <- tmp
  
  return(out)
}

eval_loss_l2 <- function(paired) {
  r1 <- paired$y - paired$y_hat
  r2 <- paired$y.2 - paired$y_hat.2
  
  r <- c(r1, r2)
  
  sum(r^2) / (nrow(paired) * 2)
}

eval_loss_01 <- function(paired) {
  r1 <- paired$y != paired$y_hat
  r2 <- paired$y.2 != paired$y_hat.2
  
  r <- c(r1, r2)
  
  mean(r)
}

eval_loss_FN <- function(paired) {
  r1 <- (paired$y != paired$y_hat) & paired$y == 1
  r2 <- (paired$y.2 != paired$y_hat.2) & paired$y == 1
  
  r <- c(r1, r2)
  
  mean(r)
}

eval_loss <- function(paired, loss_type) {
  if (loss_type == 'l2') {
    eval_loss_l2(paired)
  } else if (loss_type == '01') {
    eval_loss_01(paired)
  } else if (loss_type == 'FN') {
    eval_loss_FN(paired) 
  }else {
    stop(glue('unrecognized loss type: {loss_type}'))
  }
}

generate_forecast_synth <- function(n) {
  
  x1 <- runif(n, min = -2, max = 2)
  x2 <- runif(n, min = -1, max = 1)
  
  y <- x1 + x2 + rnorm(n, mean = 0, sd = 1)
  y_hat <- sign(x1) + sign(x2) + rnorm(n, mean = 0, sd = 1)
  y_model <- x1
  
  mean((y - y_hat)^2)
  mean((y - y_model)^2)
  
  rescale.1 <- lm(y ~ y_hat)
  pred.1 <- predict(rescale.1, data.frame(y_hat = y_hat))
  mean((y - pred.1)^2)
  
  rescale.2 <- lm(y ~ y_model + y_hat)
  pred.2 <- predict(rescale.2, data.frame(y_model = y_model, y_hat = y_hat))
  mean((y - pred.2)^2)
  
  data.frame(x1 = x1, x2 = x2, y = y, y_hat = y_hat)
}

generate_loss_comparison <- function(n, reps) {
  lapply(1:reps, function(i) {
    set.seed(i)
    data <- generate_forecast_synth(n)
    
    y <- data$y
    y_hat <- data$y_hat
    x1 <- data$x1
    
    rescale.1 <- lm(y ~ y_hat)
    pred.1 <- predict(rescale.1, data.frame(y_hat = y_hat))
    
    rescale.2 <- lm(y ~ x1 + y_hat)
    pred.2 <- predict(rescale.2, data.frame(x1 = x1, y_hat = y_hat))
    
    human_mse <- mean((y - y_hat)^2)
    human_mse_rescaled <- mean((y - pred.1)^2)
    algorithm_mse <- mean((x1 - y)^2)
    combined_mse <- mean((y - pred.2)^2)
    
    data.frame(human_mse,
               human_mse_rescaled,
               algorithm_mse,
               combined_mse,
               seed = i)
        
  }) %>% rbindlist %>% as_tibble
  
}

run_experiment <- function(seed, n, L, K) {
  set.seed(seed)
  X <- generate_forecast_synth(n)
  
  x1 <- X$x1
  x2 <- X$x2
  
  y <- X$y
  y_hat <- X$y_hat
  
  rescale.1 <- lm(y ~ y_hat)
  pred.1 <- predict(rescale.1, data.frame(y_hat = y_hat))
  
  rescale.2 <- lm(y ~ x1 + y_hat)
  pred.2 <- predict(rescale.2, data.frame(x1 = x1, y_hat = y_hat))
  
  p.value.h0 <- generate_p_value(X, L, K, loss_type = 'l2', seed = seed)
  p.value.h1 <- generate_p_value(X %>% select(-x2), L, K, loss_type = 'l2', seed = seed)
  
  data.frame(
    n = n,
    L = L,
    K = K,
    p.value.h0 = p.value.h0,
    p.value.h1 = p.value.h1,
    seed = seed,
    l2.y_hat = mean((y - y_hat)^2),
    l2.y_hat.rescale = mean((y - pred.1)^2),
    l2.optimal = mean((y - X$x1)^2),
    l2.comb = mean((y - pred.2)^2)
  )
}

run_power_test <- function(seed, n, L, K, eps) {
  if (n %% 2 != 0) {
    stop('n must be even for simplicity')
  }
  
  set.seed(seed)
  
  lapply(eps, function(e) {
    x <- runif(n / 2)
    y <- rep(1, n/2)
    x <- c(x, x)
    y <- c(y, 1 - y)
    y_hat <- runif(n/2) <= .5 + e %>% as.numeric
    y_hat <- c(y_hat, 1 - y_hat)
    
    d <- data.frame(x = x, y = y, y_hat = y_hat)
    
    p.value <- generate_p_value(d, L, K, loss_type = '01', seed = seed)
    
    data.frame(
      n = n,
      L = L,
      K = K,
      p.value = p.value,
      epsilon = e,
      seed = seed
    )
  }) %>% rbindlist %>% as_tibble
  
}

generate_distribution_synth <- function(reps, n, L, K, run.serial = FALSE) {
  if (reps <= 10 || run.serial) {
    lapply(1:reps, function(r) {
      run_experiment(r, n, L, K)
    }) %>% rbindlist %>% as_tibble
  }
  else {
    cl <- parallel::makeCluster(detectCores())
    doParallel::registerDoParallel(cl)
    
    foreach(r = 1:reps, .combine = 'rbind') %dopar% {
      # reload required libraries in each thread
      library(tidyverse)
      library(data.table)
      library(glue)
      library(parallel)
      library(doParallel)
      library(Rfast) # provides a (very!) efficient implementation of the NxN distance calculation
      
      source('./helpers.R')
      source('./logic.R')
      
      run_experiment(r, n, L, K)
    }
    
    stopCluster(cl)
  }
}

generate_distribution_power <- function(reps, n, L, K, eps = NULL) {
  if (is.null(eps)) {
    eps <- 1:10 / 20 - .05
  }
  
  lapply(1:reps, function(r) {
    run_power_test(r, n, L, K, eps)
  }) %>% rbindlist %>% as_tibble
}

generate_distribution_excess_error <- function(reps, n, L, K) {
  lapply(L, function(l) {
    lapply(1:reps, function(r) {
      seed <- r
      set.seed(seed)
      
      x1 <- runif(n, max = 10)
      x2 <- runif(n, max = 10)
      x3 <- runif(n, max = 10)
      
      eps1 <- rnorm(n, sd = 1)
      eps2 <- rnorm(n, sd = 1)
      
      y <- x1 + x2 + x3 + eps1
      y_hat <- x1 + x2 + x3 + eps2
      
      d <- data.frame(x1 = x1, x2 = x2, x3 = x3, y = y, y_hat = y_hat)
      
      p.value <- generate_p_value(d, l, K, loss_type = 'l2', seed = seed)
      
      data.frame(
        n = n,
        L = l,
        K = K,
        p.value = p.value,
        seed = seed
      )
      
    }) %>% rbindlist %>% as_tibble
  }) %>% rbindlist %>% as_tibble
}



calc_confusion_matrix <- function(df) {
  m.y <- df$y %in% c(0, 1)
  m.y_hat <- df$y_hat %in% c(0, 1)
  if (any(!m.y) || any(!m.y_hat)) {
    stop('confusion matrix only defined for binary outcomes')
  }
  
  data.frame(
    acc = mean(df$y == df$y_hat),
    admitted.frac = mean(df$y_hat == 1),
    tpr = (sum(df$y == 1 & df$y_hat == 1) / sum(df$y == 1)),
    tnr = (sum(df$y == 0 & df$y_hat == 0) / sum(df$y == 0)),
    fpr = (sum(df$y == 0 & df$y_hat == 1) / sum(df$y == 0)),
    fnr = (sum(df$y == 1 & df$y_hat == 0) / sum(df$y == 1))
  ) 
}

printable_confusion_matrix <- function(df) {
  m.y <- df$y %in% c(0, 1)
  m.y_hat <- df$y_hat %in% c(0, 1)
  if (any(!m.y) || any(!m.y_hat)) {
    stop('confusion matrix only defined for binary outcomes')
  }
  
  data.frame(
    acc = format_mean_se(df$y == df$y_hat),
    admitted.frac = format_mean_se(df$y_hat == 1),
    tpr = format_mean_se((df$y_hat == 1)[df$y == 1]),
    tnr = format_mean_se((df$y_hat == 0)[df$y == 0]),
    fpr = format_mean_se((df$y_hat == 1)[df$y == 0]),
    fnr = format_mean_se((df$y_hat == 0)[df$y == 1])
  ) 
}

format_mean_se <- function(data, digits = 2) {
  m <- mean(data) %>% round(digits)
  se <- (sd(data) / sqrt(length(data))) %>% round(digits)
  
  glue('{format(m, nsmall=digits)} \u00B1 {format(2*se, nsmall = digits)}')
}


format_mean_sd <- function(data, digits = 2) {
  m <- mean(data) %>% round(digits)
  s <- sd(data) %>% round(digits)
  
  glue('{format(m, nsmall=digits)} \u00B1 {format(2*s, nsmall = digits)}')
}
