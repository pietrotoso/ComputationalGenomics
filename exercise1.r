
# Per  i geni che avevano più di una variante ho inserito nella tabella i conteggi
# di ogni genotipo presente dato che il chisquare test non ha un numero di colonne fissato
# come in questo caso
#         G A REF ACA ACA A A
# group1   3   2       0   0
# group2   2  12       1   1
# Un'alternativa più generale di cui avevo parlato con la prof potrebbe essre racchiudere in un'unica
# categoria Aa tutte le varianti con alleli diversi di quel gene e in aa tutte le varianti con
# alleli uguali, ma credo sarebbe meno informativa




VARIANTanalysis<-function(filename){
  mydata<-read.delim(filename, row.names = 1)
  MAT_results<-ptab(mydata)
  return(MAT_results)}
  

conti<-function(riga){
  #input: list corresponding to a row (gene) in the data.frame
  #output: named list with genotypes and corresponding counts
  
  genotypes<-c()
  conteggio<-list()
  #create genotypes vector and counts list
  for (data in riga){
    if (!(data %in% genotypes)) {
      genotypes<-append(genotypes,data)
      num<- riga[riga==data]
      conteggio<-append(conteggio,length(num))}
  }
  #assigning names to the list
  names(conteggio)<-genotypes
  conteggio<-nomina(conteggio)
  return(conteggio)}


nomina<-function(lista){
  #input: named list with 'raw' genotypes (as written in the file) and counts
  #output: named list with genotypes and counts (REF if genotype is the refernce one,
  #'allele1 allele2' otherwise)
  
  for (i in (1:length(lista))){
    if (names(lista)[i]=='---'){
      names(lista)[i]<-'REF'}
    #obtaining genotypes from raw genotypes and updating them
    else {
      dividi<-strsplit(names(lista)[i],'_')
      dividi<-unlist(dividi)
      geno<-dividi[3]
      varianti<-c(dividi[1],unlist(strsplit(dividi[2],',')))
      alleli<-strsplit(geno,'/')
      alleli<-unlist(alleli)
      names(lista)[i]<-paste(varianti[as.integer(alleli[1])+1],varianti[as.integer(alleli[2])+1])
    }
  }
  return (lista)}

tabestesa<-function(conteggio1,conteggio2){
  #input:conteggio1<-named list with genotypes and relative counts of group 1
  #      conteggio2<-named list with genotypes and relative counts of group 2
  #output: contingency table of genotypes counts for the two groups
  
  #craeting the names of the columns for the table
  columns<-c()
  for (i in (1:length(conteggio1))){
    if (!(names(conteggio1)[i] %in% columns)){
      columns<-append(columns,names(conteggio1)[i])
    }
  }
  for (i in (1:length(conteggio2))){
    if (!(names(conteggio2)[i] %in% columns)){
      columns<-append(columns,names(conteggio2)[i])
    }
  }
  lunghezza<-length(columns)
  #creating an empty table
  group1<-rep(0,lunghezza)
  group2<-rep(0,lunghezza)
  tabella<-rbind(group1,group2)
  colnames(tabella)<-columns
  #filling in the table with genotypes counts
  for (j in (1:length(columns))){
    if (colnames(tabella)[j] %in% names(conteggio1)){
      index<-which(names(conteggio1)==colnames(tabella)[j])
      number<-as.integer(conteggio1[index])
      tabella[1,j]<-sum(number)
    }
    if (colnames(tabella)[j] %in% names(conteggio2)){
      index<-which(names(conteggio2)==colnames(tabella)[j])
      number<-as.integer(conteggio2[index])
      tabella[2,j]<-sum(number)
    }
  }
  
  return(tabella)
}
tabellariga<-function(num,dati){
  #input: num<-number of the row
  #       dati<-data.frame with genetic data
  #output: contingency table of genotypes counts for the two groups
  
  primo<-conti(dati[num,1:5])
  secondo<-conti(dati[num,6:21])
  tabella<-tabestesa(primo,secondo)
  return(tabella)}

ptab<-function(dati){
  #input: data.frame with genetic data
  #output: matrix n*1 with genes and relative p-values of the chisquare test
  
  #create the vector with pvalues and a vector 'indici' of genes indices for which
  #we have just one column, since for these the chisq.test doesn't work
  indici<-c()
  lung<-length(dati[,1])
  valori<-c()
  for (i in (1:lung)){
    pvalue<-as.double(chisq.test(tabellariga(i,dati),simulate.p.value = TRUE)[3])
    valori<-append(valori,pvalue)
    if (length(tabellariga(i,dati)[1,])==1){
      indici<-append(indici,i)
    }
  }
  pvalues<-matrix(valori)
  colnames(pvalues)<-'p-value'
  rownames(pvalues)<-rownames(dati)
  #setting to 1 the p-value for the genes corresponding to the indici vector, since
  #these genes are not informative(all genotypes are equal)
  pvalues[indici,1]<-1
  return(pvalues)
}