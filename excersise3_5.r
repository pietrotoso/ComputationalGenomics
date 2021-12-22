# salvare la variabile

group1 <- dati_norm[,c(15:19)]
group2 <- dati_norm[,c(1:14)]

mean_by_row <- function(data){
  n_row = dim(data)[1]
  media <- c()
  for (i in 1:n_row){
    row <- as.matrix(data[i,])
    media[i] <- mean(row)
  }
  return (media)
}

t_test <- function(data){
  group1 <- data[,c(15:19)]
  group2 <- data[,c(1:14)]
  mean1 <- mean_by_row(group1)
  mean2 <- mean_by_row(group2)
  mean_diff <- mean1 - mean2
  alpha <- 0.05
  n_rows <- dim(data)[1]
  t_test_result <- c()
  for (i in (1:n_rows)){
    t_test_result[i] <- t.test(group1[i,], group2[i,], alternative = "two.sided", mu = mean_diff[i], conf.level = alpha)
  }
  return (t_test_result)
}

