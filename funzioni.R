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
  colnames(relazioni)<-c('gene1','gene2','sign')
  
  
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