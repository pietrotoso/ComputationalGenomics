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




