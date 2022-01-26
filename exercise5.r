

# install.packages('ppcor')
# library(ppcor)

# function to compute partial correlations
partial_corr <- function(data){
  # input: matrix of data (genes/clusters on the rows, subjects on columns)
  # output: all informations of partial correlation
  partial_matrix <- pcor(data, method = 'pearson')
  return(partial_matrix)
}

# install.packages('infotheo')
# library(infotheo)

# function to compute mutual information
mutual <- function(data, n_levels){
  # input: data: matrix of objects (rows = variables, columns = samples)
  #        n_levels: number of levels for entropy computation
  # output: matrix of mutual information
  # discretize dataset (we need variables on columns and samples on rows)
  N <- dim(data)[1]
  M <- dim(data)[2]
  discretized <- data
  discretized<-t(discretized)
  discretized<-discretize(discretized,nbins=7,disc='equalwidth')
  rownames(discretized)<-colnames(data)
  mutualm<-matrix(rep(0,N^2),nrow=N)
  rownames(mutualm)<-colnames(discretized)
  colnames(mutualm)<-colnames(discretized)
  mutualm<-mutinformation(discretized)
  return(mutualm)
}

exactmutual<-function(mutual,normali){
  discretized <- normali
  discretized<-t(discretized)
  discretized<-discretize(discretized,nbins=7,disc='equalwidth')
  rownames(discretized)<-colnames(normali)
  names_df <- colnames(discretized)
  values_df <- rep(0, N)
  entropy_df <- data.frame(entropia = values_df, row.names = names_df)
  # computing entropy
  for (i in (1:N)){
    row_curr <- discretized[,i]
    entropy_df[i,] <- entropy(row_curr)
  }
  
  relazioni<-matrix(c(0,0),nrow=1)
  for (i in (1:N)){
    for (j in (1:N)){
      if (i!=j){
        if (entropy_df[i,]==mutual[i,j]){
          relazioni<-rbind(relazioni,c((rownames(discretized))[i],(rownames(discretized))[j]))
        }
      }
    }
  }
  relazioni<-relazioni[-1,]
  colnames(relazioni)<-c('causa','effetto')
  return(relazioni)
}

# function to perform hierarchical clustering
clustering <- function(data){
  # input: data frame of objects to cluster
  # output: object of class hclust (with relative attributes)
  # computing distance matrix
  dist_matrix <- dist(data, method = 'euclidean')
  # average linkage is used to update distance matrix
  hierarchical <- hclust(dist_matrix, method = 'average')
  groups<-cutree(hierarchical,4000)
  gruppi<-list()
  for(i in 1:4000){
    gruppi<-append(gruppi,list(names(groups[groups==i])))
  }
  return(gruppi)
}


