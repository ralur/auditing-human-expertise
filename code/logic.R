
generate_synth_data <- function(n) {
  ### CREATE SYNTHETIC DATA ###
  x1 <- runif(n, -1, 1)
  x2 <- runif(n, -1, 1)
  x3 <- runif(n, -2, 2)
  x4 <- x3 + runif(n, -1, 1)
  
  X <- data.frame(x1, x2, x3, x4)  
  
  eps1 <- rnorm(n, 0, 1)
  
  y <- .5 * x1 + 2 * x2 + 1.75 * x3 - 3*x4 + (x1 * x2) + eps1
  
  model <- lm(y ~ x1 + x2 + x3 + x4)
  
  eps2 <- rnorm(n, 0, 2)
  y_hat <- predict(model, X) + eps2
  
  outcomes <- X %>% mutate(y = y, y_hat = y_hat) %>% as_tibble
  
  outcomes
}

generate_p_value <- function(data, L, K, seed = 1, return.diagnostics = FALSE, loss_type = 'l2') {
  set.seed(seed)
  
  X <- data %>% select(-c(y, y_hat))
  
  if (L > floor(nrow(X) / 2)) {
    stop('L must be at most n/2')
  }
  
  ### PAIR NEAREST NEIGHBORS ###
  
  distance_mat <- distance_matrix(X)
  
  neighbors <- data.frame(i = rep(0, L), j = rep(0, L))
  distances <- rep(0, L)
  
  for (k in 1:L) {
    m <- which(distance_mat == min(distance_mat), arr.ind = TRUE)
    fst <- m[1,1]
    snd <- m[1,2]
    distances[k] <- distance_mat[fst, snd]
    distance_mat[fst,] <- Inf
    distance_mat[snd,] <- Inf
    distance_mat[,fst] <- Inf
    distance_mat[,snd] <- Inf
    
    neighbors[k,] <- c(fst, snd)
  }
  
  
  paired <- apply(neighbors, 1, function(row) {
    
    nearest <- data[row[2],]
    colnames(nearest) <- map_chr(colnames(nearest), function(c) str_c(c, '.2'))
    c(data[row[1],], nearest)
    
  }) %>% rbindlist %>% as_tibble
  
  ### GENERATE P-VALUE ###
  
  actual_loss <- eval_loss(paired, loss_type)
  
  synthetic_losses <- sapply(1:(K+1), function(k) {
    permuted <- permute_preds(paired)
    eval_loss(permuted, loss_type)
  })
  
  lt <- sum(synthetic_losses < actual_loss)
  # break any ties uniformly at random
  r <- runif(length(synthetic_losses)) <= .5
  eq <- sum((synthetic_losses == actual_loss) & r)
  
  p_value <- (1 + lt + eq) / (K+1)
  
  if (return.diagnostics) {
    
    median <- quantile(synthetic_losses, .5)
    critical_value <- quantile(synthetic_losses, .05)
    
    n_imperfect <- paired %>%
      select(-c(y, y_hat, y.2, y_hat.2)) %>%
      apply(1, function(x) any(x[1:(length(x)/2)] != x[((length(x)/2)+1):length(x)])) %>%
      sum
    
    possible_changes_01 <- paired %>%
      filter(y != y.2 & y_hat != y_hat.2)
    
    possible_improvements_01 <- possible_changes_01 %>%
      filter(y_hat.2 == y)
    
    # plot
    p <- data.frame(losses = synthetic_losses) %>%
      ggplot(aes(x=losses)) + geom_histogram() +
      geom_vline(xintercept = actual_loss, color = 'blue') +
      geom_vline(xintercept = critical_value, color = 'red') +
      geom_vline(xintercept = median) +
      annotate("text", x=actual_loss, y = 10, label="actual loss", angle=90) +
      annotate("text", x=median, y = 10, label="median loss", angle=90) +
      annotate("text", x=critical_value, y = 10, label="critical value", angle=90)
    
    list(
      plot = p,
      actual_loss = actual_loss,
      synthetic_losses = synthetic_losses,
      median = median,
      critical_value = critical_value,
      p_value = p_value,
      paired = paired,
      distances = distances,
      n_imperfect = n_imperfect,
      possible_changes_01 = possible_changes_01,
      possible_improvements_01 = possible_improvements_01,
      neighbors = neighbors
    )
  } else {
    p_value
  }
}


