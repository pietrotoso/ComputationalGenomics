<<<<<<< HEAD
# function to perform hierarchical clustering
clustering <- function(data){
  # input: data frame of objects to cluster
  # output: object of class hclust (with relative attributes)
  # computing distance matrix
  dist_matrix <- dist(data, method = 'euclidean')
  # average linkage is used to update distance matrix
  hierarchical <- hclust(dist_matrix, method = 'average')
  return(hierarchical)
}

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
  # output:
  # discretize dataset (we need variables on columns and samples on rows)
  N <- dim(data)[1]
  M <- dim(data)[2]
  discretized <- data
  #discretized <- matrix(0, nrow = N, ncol = M)
  for (i in (1:N)){
    # current row
    row_curr <- data[i,]
    minimum <- min(row_curr)
    maximum <- max(row_curr)
    levels <- seq(minimum , maximum , (maximum - minimum)/n_levels)
    for (j in (1:M)){
      # analising current row
      # set parameters for while loop
      k <- 1
      boolean <- TRUE
      while (boolean == TRUE){
        # find level
        if (k < (n_levels + 1)){
          if ((levels[k] <= row_curr[j]) && (row_curr[j] <= levels[k + 1])){
            boolean <- FALSE
            # update matrix
            discretized[i,j] <- k
            k <- 1
            break
          }else{
            k <- k + 1
          }
        }else{
          discretized[i,j] <- k
          boolean <- FALSE
        }
        #k <- 1
      }
    }
    
  }
  # compute entropy
  # intialize a matrix where each row is a gene and the value is entropy
  names_df <- rownames(discretized)
  values_df <- rep(0, N)
  entropy_df <- data.frame(entropia = values_df, row.names = names_df)
  # computing entropy
  for (i in (1:N)){
    row_curr <- discretized[i,]
    entropy_df[i,] <- entropy(t(row_curr))
  }
  
  #entropy_mat <- entropy(discretized)
  # compute mutual information
  #matrix_discr <- as.matrix(discretized)
  mutual<-matrix(rep(0,N^2),nrow=N)
  rownames(mutual)<-rownames(discretized)
  colnames(mutual)<-rownames(discretized)
  for (i in (1:(N-1))){
    for (j in ((i+1):N)){
      mutual[i,j]<-mutinformation(t(discretized[i,]),t(discretized[j,]))
    }
  }
  mutual<-mutual+t(mutual)
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

# function to compute entropy
entropia <- function(discretized, n_levels){
  # input: discretized: data frame of discretized data (row = genes, columns = samples)
  #        n_levels: nmumber of levels to compute entropy
  # output: data frame with row = genes, column = entropy value
  # number of rows and columns
  N <- dim(discretized)[1]
  M <- dim(discretized)[2]
  # count of many sample for each level
  counts_lev <- rep(0, n_levels)
  for (i in (1:N)){
    for (j in (1:M)){
      current <- discretized[i, j] # value and position in counts_lev
      counts_lev[current] <- counts_lev[current] + 1
    }
  }
  # computation of entropy
  nomi <- rownames(discretized)
  entropia_val <- data.frame(entropia = rep(0, N), row.names = nomi)
  for (i in (1:N)){
    for (j in (1:n_levels)){
      current_count <- counts_lev[j]
      entropia_val[i,] <- entropia_val[i,] - (current_count/M * log(current_count/M,2))
    }
  }
  return(entropia_val)
}



data <- dati_norm[c(1:1000), c(1:5)]
cluster <- clustering(data)
# matrix that describes the merging of clusters at step i (if negative is an
# agglomeration of singletons)
merging <- cluster$merge
# distance between clusters at each merging step
altezza <- cluster$height

# fissare la soglia di altezza: combinaione di massima distanza voluta e numero
# di cluster?














# # function that computes the matrix of euclidean distances between all pairs
# # of genes
# euclidean <- function(data){
#   # input: data frame where columns are subjects and rows are genes
#   # output: matrix of euclidean distances betwwen pairs
#   #size of pairwise distance matrix
#   dimension <- dim(data)[1]
#   # size of objects
#   M <- dim(data)[2]
#   distance <- 0
#   pair_distance <- matrix(0, nrow = dimension, ncol = dimension)
#   # filling matrix
#   for (i in (1:dimension)){
#     for (j in (i:dimension)){
#       gene1 <- data.matrix(data[i,], rownames.force = NA)
#       # elimitating the 'key'
#       gene1 <- gene1[2:M]
#       gene2 <- data.matrix(data[j,], rownames.force = NA)
#       gene2 <- gene2[2:M]
#       # computing euclidean distance
#       distance <- sqrt(sum((gene1 - gene2)^2))
#       pair_distance[i,j] <- distance
#     }
#   }
#   return(pair_distance)
# }
# 
# 
# 
# # function for agglomerative clustering but with a maximum number of iterations
# 
# # function to compute which clusters will be merged
# which_merge <- function(distance){
#   # input: matrix of distances between clusters
#   # output: which clusters are merged next
#   # matrix with boolean values on positive entries
#   pos_dist <- (distance > 0)
#   # indices of positive elements
#   indices <- which(pos_dist == TRUE)
#   # corresponding values
#   values <- distance[indices]
#   minimun = min(values)
#   # position of the minimum, which coordinates corresponds to the clusters to merge
#   position <- which(distance == minimum, arr.ind = T)
#   return(position)
# }
# 
# 
# # function for agglomerative clustering with a maximum number of iterations
# clustering <- function(data, distance, threshold){
#   # input:
#   # output:
#   # clusters to merge
#   clusters <- which_merge(distance)
#   cluster1 <- clusters[1]
#   cluster2 <- clusters[2]
#   distance_curr <- distance[cluster1][cluster2]
#   while (distance_curr < threshold){
#     
#   }
# }
#   
# 
# 
# 
# 
# 
# 
# 
# # function that determines the number of iterations for agglomerative clustering
# iterations <- function(distances, threshold){
#   # input: distances: matrix of pair distances
#   #        threshold: threshold for the maximum distance tolerate (we stop merging
#   #                   when it is reached)
#   # output: number of iterations for agglomerative clustering
#   lower_than <- (distances)
# }
# 




=======
# function to perform hierarchical clustering
clustering <- function(data){
  # input: data frame of objects to cluster
  # output: object of class hclust (with relative attributes)
  # computing distance matrix
  dist_matrix <- dist(data, method = 'euclidean')
  # average linkage is used to update distance matrix
  hierarchical <- hclust(dist_matrix, method = 'average')
  return(hierarchical)
}

data <- dati_norm[c(1:1000), c(1:5)]
cluster <- clustering(data)
# matrix that describes the merging of clusters at step i (if negative is an
# agglomeration of singletons)
merging <- cluster$merge
# distance between clusters at each merging step
altezza <- cluster$height

# fissare la soglia di altezza: combinaione di massima distanza voluta e numero
# di cluster?














# # function that computes the matrix of euclidean distances between all pairs
# # of genes
# euclidean <- function(data){
#   # input: data frame where columns are subjects and rows are genes
#   # output: matrix of euclidean distances betwwen pairs
#   #size of pairwise distance matrix
#   dimension <- dim(data)[1]
#   # size of objects
#   M <- dim(data)[2]
#   distance <- 0
#   pair_distance <- matrix(0, nrow = dimension, ncol = dimension)
#   # filling matrix
#   for (i in (1:dimension)){
#     for (j in (i:dimension)){
#       gene1 <- data.matrix(data[i,], rownames.force = NA)
#       # elimitating the 'key'
#       gene1 <- gene1[2:M]
#       gene2 <- data.matrix(data[j,], rownames.force = NA)
#       gene2 <- gene2[2:M]
#       # computing euclidean distance
#       distance <- sqrt(sum((gene1 - gene2)^2))
#       pair_distance[i,j] <- distance
#     }
#   }
#   return(pair_distance)
# }
# 
# 
# 
# # function for agglomerative clustering but with a maximum number of iterations
# 
# # function to compute which clusters will be merged
# which_merge <- function(distance){
#   # input: matrix of distances between clusters
#   # output: which clusters are merged next
#   # matrix with boolean values on positive entries
#   pos_dist <- (distance > 0)
#   # indices of positive elements
#   indices <- which(pos_dist == TRUE)
#   # corresponding values
#   values <- distance[indices]
#   minimun = min(values)
#   # position of the minimum, which coordinates corresponds to the clusters to merge
#   position <- which(distance == minimum, arr.ind = T)
#   return(position)
# }
# 
# 
# # function for agglomerative clustering with a maximum number of iterations
# clustering <- function(data, distance, threshold){
#   # input:
#   # output:
#   # clusters to merge
#   clusters <- which_merge(distance)
#   cluster1 <- clusters[1]
#   cluster2 <- clusters[2]
#   distance_curr <- distance[cluster1][cluster2]
#   while (distance_curr < threshold){
#     
#   }
# }
#   
# 
# 
# 
# 
# 
# 
# 
# # function that determines the number of iterations for agglomerative clustering
# iterations <- function(distances, threshold){
#   # input: distances: matrix of pair distances
#   #        threshold: threshold for the maximum distance tolerate (we stop merging
#   #                   when it is reached)
#   # output: number of iterations for agglomerative clustering
#   lower_than <- (distances)
# }
# 




>>>>>>> 10580a9ad4ad398ac0807735707fa8c085a5dbf1
