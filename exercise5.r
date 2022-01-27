

# install.packages('ppcor')
# library(ppcor)

# function to compute partial correlations
partial_corr <- function(dati){
  # input: matrix of data (genes/clusters on the rows, subjects on columns)
  # output: all informations of partial correlation
  dati<-as.matrix(dati)
  partial_matrix <- pcor(t(dati), method = 'pearson')$estimate
  nulla<-dati
  N<-dim(nulla)[1]
  M<-dim(nulla)[2]
  for (i in 1:N){
    nulla[i,]<-sample(dati[i,],M)
  }
  parnull<- pcor(t(nulla), method = 'pearson')$estimate
  valori<-as.numeric(parnull)
  valori<-valori[!(valori==1)]
  vettorepar<-as.numeric(partial_matrix)
  vettorepar<-vettorepar[!(vettorepar==1)]
  totali<-length(valori)
  totcorr<-length(vettorepar)
  soglia<-0.78
  selezionati<-length(vettorecorr[abs(vettorecorr)>=soglia])
  expFP<-(length(valorinulli[abs(valorinulli)>=soglia])/totali)*totcorr
  FDR<-expFP/selezionati
  bool<-abs(partial_matrix)>soglia
  relazioni<-which(bool==TRUE,arr.ind = TRUE)
  relazioni<-as.data.frame(relazioni)
  relazioni<-relazioni[!(relazioni[,1]==relazioni[,2]),]
  return(relazioni)
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

# mutualrelations<-function(mutmatrix, normali){
#   discretized <- normali
#   discretized<-t(discretized)
#   discretized<-discretize(discretized,nbins=7,disc='equalwidth')
#   rownames(discretized)<-colnames(normali)
#   names_df <- colnames(discretized)
#   values_df <- rep(0, N)
#   entropy_df <- data.frame(entropia = values_df, row.names = names_df)
#   # computing entropy
#   for (i in (1:N)){
#     row_curr <- discretized[,i]
#     entropy_df[i,] <- entropy(row_curr)
#   }
#   
#   
# }

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


tabfinale<-function(data,gruppi){
  tabella<-matrix(rep(0,104000),nrow=4000)
  for (i in 1:4000){
    if (length(gruppi[[i]])==1){
      tabella[i,]<-data[gruppi[[i]],]
    }
    else{
      tabella[i,]<-colMeans(data[gruppi[[i]],])
    }
  }
  colnames(tabella)<-colnames(data)
  rownames(tabella)<-1:4000
  return(tabella)
}

dividi<-function(normali){
  lunghezze<-read.delim('raw_trascr_count_annot.txt')
  lunghezze<-lunghezze[,c('Symbol','Length')]
  nomitr<-rownames(normali)
  lunghezze<-lunghezze[lunghezze$Symbol%in%nomitr,]
  lunghezze<-as.data.frame(lunghezze)
  normali<-normali[lunghezze$Symbol,]
  for (i in 1:length(normali[,1])){
    normali[i,]<-normali[i,]/(lunghezze[lunghezze$Symbol==rownames(normali)[i],])$Length
  }
  return(normali*100)
}

coeffvariation<-function(data){
  deviazioni<-apply(data,FUN=sd,1)
  medie<-apply(data,FUN=mean,1)
  coeff<-deviazioni/medie
  return(coeff)
  
}

filtraescala<-function(data,coeff){
  filtrati<-data[coeff>0.3,]
  trasp<-t(filtrati)
  trasp<-scale(trasp)
  filtratiscalati<-t(trasp)
  return(filtratiscalati)
}
pvcorr<-function(dati){
  correlazioni<-cor(t(dati))#genes on the rows
  valorinulli<-nullcorr(dati,4)
  soglie<-sort(abs(valorinulli),decreasing=TRUE)[1000:1050]
  fdrlist<-c()
  vettorecorr<-as.vector(correlazioni)
  vettorecorr<-vettorecorr[!(vettorecorr==1)]
  totali<-length(valorinulli)
  totcorr<-length(vettorecorr)
  soglia<-0.78
  selezionati<-length(vettorecorr[abs(vettorecorr)>=soglia])
  expFP<-(length(valorinulli[abs(valorinulli)>=soglia])/totali)*totcorr
  FDR<-expFP/selezionati
  boole<-abs(correlazioni)>=soglia
  relazioni<-which(boole==TRUE,arr.ind = TRUE)
  relazioni<-as.data.frame(relazioni)
  relazioni<-relazioni[!(relazioni[,1]==relazioni[,2]),]
  M=dim(relazioni)[1]
  segni<-c()
  for (i in 1:M){
    segni<-append(segni,sign(correlazioni[relazioni[i,1],relazioni[i,2]]))
  }
  for (j in (1 :length(segni))){
    if (segni[j]==1){
      segni[j]<-'+'
    }
    else {
      segni[j]<-'-'
    }
  }
  relazioni<-cbind(relazioni,segni)
  rownames(relazioni)<-1:M
  colnames(relazioni)<-c('cluster1','cluster2','sign')
  
  
  return(relazioni)
}
nullcorr<-function(dati,k){
  valori<-c()
  for ( i in 1:k){
    nulla<-dati
    N<-dim(nulla)[1]
    M<-dim(nulla)[2]
    for (i in 1:N){
      nulla[i,]<-sample(dati[i,],M)
    }
    cornull<-cor(t(nulla))
    valcor<-as.vector(cornull)
    valori<-c(valori,valcor)
    
  }
  valori<-valori[!(valori==1)]
  return(valori)
}