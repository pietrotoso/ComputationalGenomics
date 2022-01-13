t_test <- function(data){
  group1 <- data[,c(15:19)]
  group2 <- data[,c(1:14)]
  alpha <- 0.05
  n_rows <- length(data[,1])
  t_test_result <- c()
  for (i in (1:n_rows)){
    t_test_result[i] <- t.test(group1[i,], group2[i,], alternative = "two.sided", 0, conf.level = alpha)$p.value
  }
  tabella<-matrix(t_test_result)
  colnames(tabella)<-'P_value'
  rownames(tabella)<-rownames(data)
  return (tabella)
}

wilcoxon_test <- function(data){
  group1 <- data[, c(15:19)]
  group2 <- data[, c(1:14)]
  alpha <- 0.05
  n_rows <- length(data[,1])
  wilcoxon_test_result <- c()
  for (i in (1:n_rows)){
    wilcoxon_test_result[i] <- wilcox.test(as.numeric(group1[i,]), as.numeric(group2[i,]), alternative = "two.sided", mu = 0, paired = FALSE)$p.value
    
  }
  tabella<-matrix(wilcoxon_test_result)
  colnames(tabella)<-'P_value'
  rownames(tabella)<-rownames(data)
  return (tabella)
}