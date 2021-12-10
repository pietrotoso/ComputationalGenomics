Normalization<-function(file){
  #normalizes data and plots MvAs and boxplots for both raw and normalized data
  #input: filename containing original data
  #output: normalized logaritmic data
  dati<-read.delim(file,row.names = 1)
  seq_depth<-depth(dati)
  log_counts<-logaritmici(dati)
  plotMvA(log_counts,'MvAs.pdf')
  normalized<-scaling(log_counts)
  plotMvA(normalized,'MvAs_norm.pdf')
  #normalized<-2^(normalized)
  return(normalized)
}

depth<-function(dati){
  #input: raw count data
  #output: vector sequencing depth of every sample
  seq_depth<-c()
  n_samples<-length(dati[1,])
  for (i in 1:n_samples){
    seq_depth<-append(seq_depth,sum(dati[,i]))
  }
  return(seq_depth)
}

logaritmici<-function(dati){
  #input: raw count data
  #output: logaritmic data
  n_samples<-length(dati[1,])
  for (i in 1:n_samples){
    dati[,i]<-dati[,i]+1
    dati[,i]<-log(dati[,i],2)
    }
  return(dati)
}

plotMvA<-function(dati,filename){
  #input: dati<-data of which I want MvAs plots
  #       filename<-name of the pdf in which plots will be saved
  # plots MvAs and boxplots
  n_samples<-length(dati[1,])
  first<-dati[,1]
  genes<-length(dati[,1])
  zeros<-rep(0,genes)
  pdf(filename)
  par(mfrow=c(2,2))
  stringa1<-'M=log2(r*1)-log2(r*'
  stringa2<-'A=(log2(r*1)+log2(r*'
  for (i in 2:n_samples){
    M<-first-dati[,i]
    A<-(first+dati[,i])/2
    plot(A,M,pch='.', main='MvA plot',xlab=paste(stringa2,i,'))/2',sep=''), ylab=paste(stringa1,i,')',sep=''))
    lines(A,zeros,col='red')
    boxplot(first, dati[,i],main='Boxplot', names=c('sample1',paste('sample',i,sep='')), ylab='log(r*i)')
  }
  dev.off()
}

scaling<-function(dati){
  #input: raw count data
  #output: normalized data
  n_samples<-length(dati[1,])
  first<-dati[,1]
  for (i in 2:n_samples){
    M<-first-dati[,i]
    scale_factor<-mean(M,trim=0.1)
    dati[,i]<-dati[,i]+scale_factor
  }
  return(dati)
}