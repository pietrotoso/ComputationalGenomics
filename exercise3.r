pvaledger<-function(filename){
  #Questo è un test più semplice di quelli della funzione dopo, è quello che avevo
  #usato inizialmente, lo lascio qui, ma molto probabilmente è meglio edgernuovo
  #
  #input: file containing raw counts
  #output: matrix containing the pvalues for DE analysis of genes
  dati<-read.delim('raw_trascr_count.txt', row.names = 1)
  dati<-dati[,1:19]
  groups<-rep(c('Group2','Group1'),c(14,5))
  y<-DGEList(counts=dati,group=groups)
  keep<-filterByExpr(y)
  y<-y[keep, ,keep.lib.sizes=FALSE]
  y<-calcNormFactors(y)
  y<-estimateDisp(y)
  et<-exactTest(y,pair=c('Group1','Group2'))
  pvalues<-et$table$PValue
  names(pvalues)<-rownames(et$table)
  return(pvalues)
}
edgernuovo<-function(filename){
  #A seconda delle righe che si commentano 43/44 e 46/47 si fa Quasi-likelihood
  #generalized linear model o un normale generalized linear model
  #Il primo dovrebbe essere preciso e dà un grafico migliore, ma seleziona solo 7 geni
  #Il secondo ha un grafico peggiore, ma seleziona più geni
  #
  #input: file containing raw counts
  #output: matrix containing the pvalues for DE analysis of genes
  dati<-read.delim('raw_trascr_count.txt', row.names = 1)
  dati<-dati[,1:19]
  group1<-grep("Group1",colnames(dati))
  group2<-grep("Group2",colnames(dati))
  subjects<-c()
  subjects[group1]<-1
  subjects[group2]<-2
  Group<-factor(subjects)
  design <- model.matrix(~0+Group) 
  rownames(design)<-colnames(dati[,1:19]) 
  y<-DGEList(dati)
  keep<-filterByExpr(y,design)
  y<-y[keep, ,keep.lib.sizes=FALSE]
  y<-calcNormFactors(y)
  y<-estimateDisp(y,design)
  fit<-glmQLFit(y,design)
  #fit<-glmFit(y,design)
  contrasts <- makeContrasts(DE = Group1-Group2,levels=design)
  #RES<-glmLRT(fit,contrast=contrasts[,"DE"])
  RES<-glmQLFTest(fit,contrast=contrasts[,"DE"])
  pvalues<-RES$table$PValue
  names(pvalues)<-rownames(RES$table)
  return(pvalues)
}

stime<-function(alfa,pvalues,G0){
  #input: alfa<-the significance level we are considering
  #       pvalues<- vector containing pvalues of genes
  #       G0<-the value of G0 which we're considering(our estimate)
  #output: list containing expectations of TN, FP, tp, FN
  totgenes<-length(pvalues)
  selected<-length(pvalues[pvalues<=alfa])
  expFP<-min(G0*alfa, selected)
  expTP<-max(0,selected-expFP)
  expTN<-min(G0*(1-alfa),totgenes-selected)
  expFN<-totgenes-selected-expTN
  valori<-c(expTN,expFP,expTP,expFN)
  names(valori)<-c('E[TN]','E[FP]','E[TP]','E[FN]')

  return(valori)
}

G0estimate<-function(pvalues, start=0.02,end=0.99){
  #input: pvalues<-vector containing pvalues of the genes
  #       start<-dove cominciare la ricerca dela retta di stima
  #       end<-dove finire la ricerca della retta di stima
  #start e end sono utili, una volta visto il grafico di G0 e decisa la zona in cui si
  #stabilizza, per trovare il valore di stima di G0
  #output: estimated value of G0
  totgenes<-length(pvalues)
  #setting the lambda grid
  lambdas<-seq(start,end,by=0.01)
  G0s<-c()
  for (i in lambdas){
    selected<-length(pvalues[pvalues<=i])
    predict<-(totgenes-selected)/(1-i)
    G0s<-append(G0s,predict)
  }
  plot(lambdas,G0s,ylim=c(15000,17000))
  #create a vector of possible G0s(from min value obtained to total number of genes)
  stimeG0<-seq(round(min(G0s)),length(pvalues))
  somme<-c()
  #find the value of G0 that minimizes average quadratic error
  for (stima in stimeG0){
    somma<-0
    for (i in G0s){
      somma<-somma+(stima-i)^2
    }
    somma<-somma/length(G0s)
    somme<-append(somme,somma)
  }
  indice<-which(somme==min(somme))
  valore<-stimeG0[indice]
  #plotting on the graph the line corresponding to our estimate
  linea<-rep(valore,length(lambdas))
  lines(lambdas,linea)
  
  return(stimeG0[indice])
}
select_alfa_FDR<-function(pvalues,G0, FDR){
  #input: pvalues<-vector containing pvalues of genes
  #       G0<-the value of G0 which we're considering(our estimate)
  #       FDR<-FDR that we want to have
  #output: alfa selected that gives the closest FDR to the one we want
  ordinati<-sort(pvalues)
  ordinati<-unique(ordinati)
  distanza<-1
  #determine the minimum distance between two pvalues
  for ( i in 1:(length(ordinati)-1)){
    if (ordinati[i+1]-ordinati[i]<distanza){
      distanza<-ordinati[i+1]-ordinati[i]
    }
  }
  #set epsilon as half the minimum distance
  epsilon<-distanza/2
  #create grid of pvalues
  alfatest<-c()
  alfatest<-ordinati+epsilon
  
  #find the alfa which gives the closest FDR to the one requested
  FDRlist<-c()
  for (alfa in alfatest){
    selected<-length(pvalues[pvalues<=alfa])
    expFP<-G0*alfa#round?
    FDRlist<-append(FDRlist,expFP/selected)
  }
  diff<-FDRlist-FDR
  indice<-which(abs(diff)==min(abs(diff)))
  selectedalfa<-alfatest[indice]
  return(c(selectedalfa,FDRlist[indice]))
}

selectgenes<-function(pvalues,alfa){
  selezionati<-pvalues[pvalues<=alfa]
  selezionati<-names(selezionati)
  return(selezionati)
}

t_test <- function(data){
  #input: data counts(normalized ones)
  #output: matrix with pvalue of each gene for the t-student test
  group1 <- data[,c(15:19)]
  group2 <- data[,c(1:14)]
  alpha <- 0.05
  n_rows <- length(data[,1])
  t_test_result <- c()
  for (i in (1:n_rows)){
    t_test_result[i] <- t.test(group1[i,], group2[i,], alternative = "two.sided", 0, conf.level = alpha)$p.value
  }
  names(t_test_result)<-rownames(data)
  return (t_test_result)
}

wilcoxon_test <- function(data){
  #input: data counts(normalized ones)
  #output: matrix with pvalue of each gene for the wilcoxon test
  group1 <- data[, c(15:19)]
  group2 <- data[, c(1:14)]
  alpha <- 0.05
  n_rows <- length(data[,1])
  wilcoxon_test_result <- c()
  for (i in (1:n_rows)){
    wilcoxon_test_result[i] <- wilcox.test(as.numeric(group1[i,]), as.numeric(group2[i,]), alternative = "two.sided", mu = 0, paired = FALSE)$p.value
    
  }
  names(wilcoxon_test_result)<-rownames(data)
  return (wilcoxon_test_result)
}