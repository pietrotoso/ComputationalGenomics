data <- dati_norm[c(1:1000), c(1:5)]
clusteredGenes <- clustering(data)#usata la function già implementata
clustered_genes <- cutree(clusteredGenes,k = 100)#prendo il livello in cui ho 100 cluster
#quanti geni sono nei vari cluster ?Per esempio per il soggetto Group2_S8 in percentuale ...
#tapply(matrix_b$Group2_S8, clusteredGenes , mean)

#number_cluster <- clustered_genes[14]#dice il numero del cluster a cui appartiene

#per vedere tutti i geni appartenenti a quel cluster specifico
#cluster <-subset(clusteredGenes,clustered_genes == 1) 

#-----------------------------------------------------------------------------------------------------------------
#ma quanti clusters ???
#library('cluster')
#x <- data

#pam1 <- function(x,k) list(cluster = pam(x,k>=2, cluster.only=TRUE))
#n_clusters <- clusGap(matrix_b, FUN=pam1 ,K.max=100)


#library(clusterGenomics)
#n_clusters <- gap(data,Kmax=1000,B=10)
#library(cluster)
#pam() =Partitioning (clustering) of the data into k clusters "around medoids", a more robust version of K-means.
#n_clu <- silhouette(data,dist = distance)


#---------------------------------------------------------------------------------------------------------------------------
#Bayesian network
library(bnlearn)

#1) Learning structure

#create the network:
x <- clusteredGenes
net <-  hc(x) #hill-climbing (HC) greedy search
#oppure net <- tabu(x)
plot(net)
#2) Training 

#learn the conditional probability tables (CPT) at each node (ogni nodo è un soggetto)
#run the Expectation-Maximation (EM) algorithm to learn CPT for different node
prob_cond <- bn.fit(net, data = x)
print(prob_cond$Group1_S1)

prob <- bn.net(prob_cond,debug=FALSE)

print(prob$arcs)


#3) Inference


#cpquery(prob_cond,event =(Group1_S1 ==''),evidence =  ( == ''))


#--------------------------------------------------------------------------------

#revengc è il package per fare reverse engineering che ho trovato

