

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

entropia<-function(dati){
  N=dim(dati)[1]
  dati<-t(dati)
  dati<-discretize(dati,nbins=7,disc='equalwidth')
  names_df <- colnames(dati)
  values_df <- rep(0, N)
  entropy_df <- data.frame(entropia = values_df, row.names = names_df)
  for (i in (1:N)){
    row_curr <- dati[,i]
    entropy_df[i,] <- entropy(row_curr)
  }
  return(entropy_df)
}



mutualsimilraity<-function(scalnul,scal){#(mutue,dati,nulmut){
  # nulla<-dati
  # N<-dim(nulla)[1]
  # M<-dim(nulla)[2]
  # for (i in 1:N){
  #   nulla[i,]<-sample(dati[i,],M)
  # }
  # normentropy<-entropia(dati)
  # nullentropy<-entropia(nulla)
  # nulla<-t(nulla)
  # nulla<-discretize(nulla,nbins=7,disc='equalwidth')
  # rownames(nulla)<-colnames(dati)
  # nullmutual<-matrix(rep(0,N^2),nrow=N)
  # nullmutual<-mutinformation(nulla)
  # scalnul<-scalamutua(nulmut,nullentropy)
  # scal<-scalamutua(mutue,normentropy)
  soglia=0.63
  selezionati<-length(which(scal>=soglia&scal!=1))
  valori<-length(which(scalnul>=soglia&scalnul!=1))
  alfa<-valori/length(which(scalnul!=1))
  expFP<-alfa*length(which(scal!=1))
  FDR<-expFP/selezionati
  maggiori<-which(scal>soglia &scal!=1, arr.ind = TRUE)
  colnames(maggiori)<-c('cluster1','cluster2')
  rownames(maggiori)<-1:dim(maggiori)[1]
  maggiori<-checkaracne(maggiori)
  return(maggiori)
}

checkaracne<-function(relazioni){
  unici<-unique((relazioni[,1]))
  finali<-c()
  triplette<-c()
  for (i in unici){
    opposite<-relazioni[relazioni[,1]==i,2]
    for (j in opposite){
      opposite2<-relazioni[relazioni[,1]==j,2]
      for (k in opposite2){
        if (k %in% opposite){
          triplette<-append(triplette,c(i,j,k))
          
        }
      }
    }
  }
  if (length(triplette)==0){
    return(relazioni)
  }
}
scalamutua<-function(mutua,entropia){
  N<-dim(mutua)[1]
  for (i in 1:N){
    for (j in 1:N){
      mutua[i,j]<-mutua[i,j]/max(entropia[i,],entropia[j,])
    }
  }
  return(mutua)
  
}
intersezione<-function(subset1,subset2){
  subset1<-subset1[c(1,2),]
  subset2<-subset2[c(1,2),]
  rownames(subset1)<-c('cluster1','cluster2')
  rownames(subset2)<-c('cluster1','cluster2')
  intersez<-intersect(subset1,subset2)
  return(intersez)
}

entrambi<-function() {
  dati<-read.delim('GeneticData.txt',row.names=1)
  annotazione<-read.delim('ensembl_annot.txt')
  cromosomi<-row.names(dati)
  simboli<-annotazione$Symbol
  simboli<-annotazione$Variation.ID
  associazione<-matrix(nrow=length(cromosomi))
  colnames(associazione)<-'Gene'
  geni<-c()
  for (i in cromosomi){
    geni<-append(geni,annotazione$Symbol[grep(unlist(strsplit(i,'_'))[2],simboli)[1]])
  }
  associazione[,1]<-geni
  unici<-unique(geni)
  raw<-read.delim('raw_trascr_count.txt',row.names=1)
  totali<-rownames(raw)
  comuni<-unici[unici%in%totali]
  dati<-cbind(dati,geni)
  comgenetic<-dati[dati$geni%in%comuni,]
  raw<-raw[rownames(raw)%in%comuni,]
  raw<-raw[,1:19]
  raw<-cbind(raw,rownames(raw))
  colnames(raw)[length(colnames(raw))]<-'geni'
  soggettigen<-colnames(comgenetic)
  soggettiraw<-colnames(raw)
  comgenetic<-comgenetic[,c(soggettigen%in%soggettiraw)]
  ultima<-comgenetic[,length(comgenetic[1,])]
  comgenetic<-booleana(comgenetic)
  comgenetic[,length(comgenetic[1,])]<-ultima
  comgenetic<-as.data.frame(comgenetic)
  tabcor<-matrix(rep(0,length(rownames(comgenetic))),nrow=length(rownames(comgenetic)))
  rownames(tabcor)<-rownames(comgenetic)
  comgenetic<-comgenetic[,order(colnames(comgenetic))]
  raw<-raw[,order(colnames(raw))]
  for (i in 1:length(rownames(comgenetic))){
    gene<-comgenetic$geni[i]
    tabcor[i]<-cor(as.numeric(comgenetic[i,2:length(comgenetic[1,])])
                   ,as.numeric(raw[gene,2:length(raw[1,])]))
  }
  posizioni<-which(abs(tabcor)>=0.7)
  tenute<-tabcor[posizioni,]
  pdf('snpimportanti.pdf')
  par(mfrow=c(2,2))
  for (j in 1:length(tenute)){
    snp<-names(tenute)[j]
    gene<-comgenetic[snp,]$geni
    plot(1:19,comgenetic[snp,2:20], ylim=c(-4,4),xlab='subjects',ylab='expression')
    vettore<-raw[gene,2:20]
    vettore<-t(vettore)
    vettore<-scale(vettore)
    vettore<-t(vettore)
    points(1:19,vettore,col='red')
    legend(x='bottomright',c(snp,gene), fill=c('black','red'),bty='n')
  }
  dev.off()
  return(tabcor)
}


booleana <- function(data, file=TRUE){
  # input: data: file name
  if( file==TRUE){
    genetici <- read.delim(data, row.names = 1)
  }
  else{
    genetici<-data
  }
  N <- dim(genetici)[1]
  M <- dim(genetici)[2]
  row_names <- rownames(genetici)
  col_names <- colnames(genetici)
  booleana_tab <- matrix(rep(0, N*M), nrow = N)
  rownames(booleana_tab) <- row_names
  colnames(booleana_tab) <- col_names
  for (i in (1:N)){
    for (j in (1:M)){
      if (genetici[i,j] == '---'){
        booleana_tab[i,j] <- 0
      }
      else{
        booleana_tab[i,j] <- 1
      }
    }
  }
  return(booleana_tab)
}


# elimina righe di tutti 0 e di tutti 1
processing <- function(data){
  N <- dim(data)[1]
  M <- dim(data)[2]
  lista_indici <- c()
  for (i in (1:N)){
    somma <- sum(data[i,])
    if ((somma < 1) | (somma == M)){
      lista_indici <- append(lista_indici, i)
    }
  }
  data <- data[-lista_indici,]
  return(data)
}


#elimina righe uguali
cluster_bool <- function(data){
  dist_matrix <- dist(data, method = 'euclidean')
  cluster_boo <- hclust(dist_matrix, method = 'average')
  altezze <- cluster_boo$height
  rimanenti <- length(altezze[altezze > 0])
  rimanenti <- rimanenti + 1
  clustering <- cutree(cluster_boo, rimanenti)
  gruppi<-list()
  for(i in 1:rimanenti){
    gruppi<-append(gruppi,list(names(clustering[clustering==i])))
  }
  return(gruppi)
}


tabfinale_boo<-function(data,gruppi){
  tabella<-matrix(rep(0,14994),nrow=714)
  for (i in 1:714){
    # righe <- data[gruppi[[i]],]
    # riga <- righe[1,]
    # tabella[i,] <- riga
    
    if (length(gruppi[[i]])==1){
      tabella[i,]<-data[gruppi[[i]],]
    }
    else{
      righe <- data[gruppi[[i]],]
      riga <- righe[1,]
      tabella[i,]<-riga
    }
  }
  colnames(tabella)<-colnames(data)
  rownames(tabella)<-1:714
  return(tabella)
}

rete_booleana <- function(dati){
  N=dim(dati)[1]
  dati<-t(dati)
  names_df <- colnames(dati)
  values_df <- rep(0, N)
  entropy_df <- data.frame(entropia = values_df, row.names = names_df)
  for (i in (1:N)){
    row_curr <- dati[,i]
    entropy_df[i,] <- entropy(row_curr)
  }
  discreti <- discretize(dati, nbins = 2, disc = 'equalwidth')
  mutual_bool<-mutinformation(discreti)
  
  relazioni<-matrix(c(0,0),nrow=1)
  for (i in (1:N)){
    for (j in (1:N)){
      if (i!=j){
        if (entropy_df[i,]==mutual_bool[i,j]){
          relazioni<-rbind(relazioni,c((rownames(mutual_bool))[i],(rownames(mutual_bool))[j]))
        }
      }
    }
  }
  relazioni<-relazioni[-1,]
  colnames(relazioni)<-c('causa','effetto')
  return(relazioni)
}


# TOGLIERE CHI HA I REGOLATORI
# number è il numero di geni senza regolatore
# rete_booleana_doppia <- function(dati){
#   N=dim(dati)[1]
#   dati<-t(dati)
#   names_df <- colnames(dati)
#   values_df <- rep(0, N)
#   entropy_df <- data.frame(entropia = values_df, row.names = names_df)
#   for (i in (1:N)){
#     row_curr <- dati[,i]
#     entropy_df[i,] <- entropy(row_curr)
#   }
#   discreti <- discretize(dati, nbins = 2, disc = 'equalwidth')
#   
#   # inizializzo la matrice
#   relazioni <- matrix(c(0,0), nrow = 1)
#   relazioni <- as.data.frame(relazioni)
#   #rownames(relazioni) <- rownames(dati)
#   coppie <- t(combn((1:N),2))
#   
#   #regolazioni <- matrix(c(0,0), nrow = 1)
# 
#   # mutual information con due regolatori
#   for (k in (1:N)){
#     # colonna corrente
#     corrente <- discreti[,k]
#     nome_corr <- k
#     entropia_corr <- entropy_df[k,]
#     # ora per la colonna cerco tutte le coppie
#     for (i in (1:length(coppie[,1]))){
#       elem1 <- coppie[i,1]
#       elem2 <- coppie[i,2]
#       if ((elem1 != nome_corr) & (elem2 != nome_corr)){
#         matrice <- discreti[,c(elem1,elem2)]
#         mutua_bool <- mutinformation(corrente, matrice)
#         if (mutua_bool == entropia_corr){
#           stringa <- paste(coppie[i,1], coppie[i,2], sep = ',')
#           relazioni<-rbind(relazioni, c(k, stringa))
#         }
#       }
#     }
#   }
#   
# 
#   relazioni<-relazioni[-1,]
#   colnames(relazioni)<-c('causa','effetto')
#   rownames(relazioni) <- NULL
#   return(relazioni)
# }


# BAYESIAN
#rete <- rsmax2(dati)
