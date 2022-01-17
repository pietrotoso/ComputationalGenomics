#GO Enrichment analysis FUNCTIONS#

call_db<-function(selezionati,nomitotali){
  #input: selezionati<- names of selected genes from previous exercises
  #       nomitotali<-names of all the genes
  #output: matrix with all genes and relative GOs with a column that tells if the gene is
  #        one of the selected ones
  
  #selection of information relative to our genes
  gene_GO<-select(org.Hs.eg.db,keys=nomitotali,columns=c('ENTREZID','GOALL'),keytype='SYMBOL')
  #eliminate evidenceall column
  gene_GO<-gene_GO[,-4]
  #select lines that are repeated
  doppi<-duplicated(gene_GO[,2:3])
  #eliminate repeated lines
  gene_GO<-gene_GO[doppi==FALSE,]
  
  #add a column to track if a gene is part of the selected
  SELECTED<-""
  gene_GO<-cbind(gene_GO,SELECTED)
  gene_GO$SELECTED<-gene_GO$SYMBOL%in%selezionati
  #eliminate NAs
  gene_GO_NA<-na.omit(gene_GO)
  return(gene_GO_NA)
}

#function that assign the total num of genes and selected genes to each GO term
assigned<-function(value,go,sel)
{
  index<-which(go==value)
  num_go<-length(index)
  num_sel<-length(which(sel[index]==TRUE))
  mat<-rbind(num_go,num_sel)
  return(mat)
}
#Elimination of single gene-GO terms and creation of unique GO term list
best<-function(gene_GO_NA){
  #input: matrix with genes and relative GOs with a column that tells if the gene is
  #       one of the selected ones
  #output: matrix that tells for each GO the total occurrences of the term in the genes and
  #        the occurences in the selected subset
  
  #eliminate duplicates
  GO_present<-duplicated(gene_GO_NA$GOALL)
  GO<-gene_GO_NA$GOALL[GO_present==FALSE] 
  
  #create the output matrix
  count<-sapply(GO,assigned,gene_GO_NA$GOALL,gene_GO_NA$SELECTED)
  count<-t(count)
  colnames(count)<-c("N_GENES","N_SELECTED")
  
  
  
  #remove GOs that have just one count
  index_one<-which(count[,1]==1)
  FINAL<-count[-index_one,]
  return(FINAL)
}

#Total amount of genes and selected genes needed in Fisher Test
total<-function(gene_GO_NA,FINAL){
  #input: gene_go_NA<-matrix with genes and relative GOs with a column that tells if the gene is
  #                   one of the selected ones
  #       FINAL<- matrix that tells for each GO the total occurrences of the term in the genes and
  #               the occurences in the selected subset
  #output:counts of total genes remained and total selected remained
  
  #total genes remained
  gene_GO_NA<-gene_GO_NA[gene_GO_NA$GOALL%in%rownames(FINAL),]
  gene_tot<-duplicated(gene_GO_NA$SYMBOL)
  ind_gene<-which(gene_tot==FALSE)
  genes<-length(ind_gene)
  
  #numero selected totali 
  selected_tot<-length(which(gene_GO_NA$SELECTED[ind_gene]==TRUE))
  return(c(genes,selected_tot))
}

#function to compute Fisher Exact test
fisher_test<-function(data,selected_tot,genes){
  #input:data<-matrix that tells for each GO the total occurrences of the term in the genes and
  #            the occurences in the selected subset
  #      selected_tot<- number of selected genes remained
  #      genes<-number of total genes remained
  #output: p-values of the GOs fisher test for enrichment analysis
  
  #initialize variables for the test
  tot_ac<- data[1]
  a<- data[2]
  c<- tot_ac-a
  b<- selected_tot-a
  d<- genes-c-selected_tot 
  #create matrix
  mat<-matrix(c(a,c,b,d),2,2)
  #perform test
  res<-fisher.test (mat, alternative="greater")  
  return(res$p.value)
}
pvfisher<-function(go_mat){
  #questa funzione racchiude le altre scritte in precedenza, fatta per poter riapplicare alle diverse ontology classes
  #input: matrix with genes and relative GOs with a column that tells if the gene is
  #       one of the selected ones
  #output: p-values of the GOs fisher test for enrichment analysis
  
  #Unique GO terms matrix and associated genes
  FINAL_GO<-best(go_mat)
  FINAL_GO<-FINAL_GO[-which(FINAL_GO[,2]<1),] # elimination of a<1 selected genes
  nGO<-nrow(FINAL_GO)
  #Final matrix to be used in Fisher test
  tot<-total(go_mat,FINAL_GO)
  tot_genes<-tot[1]
  tot_selected<-tot[2]
  # Fisher test on each row of the GO matrix
  p_value<-apply(FINAL_GO,1,fisher_test,tot_selected,tot_genes)
  return(p_value)
}
G0estimate<-function(pvalues, start=0.02,end=0.99,ylim1=100,ylim2=500, disegna=TRUE){
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
  plot(lambdas,G0s,ylim=c(ylim1,ylim2))
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
  if (disegna==TRUE){
    linea<-rep(valore,length(lambdas))
    lines(lambdas,linea)
  }
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
enrichedGO<-function(pvalues,start,end,y1,y2){
  #funzione che dovrà selezionare i GO enriched, da scrivere quando si decide come selezionare
  par(mfrow=c(1,2))
  prima<-G0estimate(pvalues,ylim1 = y1,ylim2=y2, disegna=FALSE)
  stima<-G0estimate(pvalues,start,end,y1,y2)
  alfa<-select_alfa_FDR(pvalues,stima,0.05)[1]
  selezionati<-selectgenes(pvalues,alfa)
  return(selezionati)
}
####################################################################

#Main #
main<-function(selected,totali){
  library(org.Hs.eg.db)#Homo sapiens database call
  library(GO.db)
  
  #Creation of GO database matrix
  go_mat<-call_db(selected,totali)
  #creation of matrices for different ontology classes
  matriceMF<-go_mat[go_mat$ONTOLOGYALL=='MF',]
  matriceCC<-go_mat[go_mat$ONTOLOGYALL=='CC',]
  matriceBP<-go_mat[go_mat$ONTOLOGYALL=='BP',]
  #computation of the pvalues
  pvaluesMF<-pvfisher(matriceMF)
  pvaluesCC<-pvfisher(matriceCC)
  pvaluesBP<-pvfisher(matriceBP)
  #extraction of enriched GOs
  enrichedMF<-enrichedGO(pvaluesMF)
  enrichedCC<-enrichedGO(pvaluesCC,0.4,0.6,100,400)
  enrichedBP<-enrichedGO(pvaluesBP, 0.4,0.6,100,400)
  #extraction of the corresponding GO term
  descriptionMF<-select(GO.db,keys=names(enrichedMF),columns='TERM',keytype='GOID')
  descriptionCC<-select(GO.db,keys=names(enrichedCC),columns='TERM',keytype='GOID')
  descriptionBP<-select(GO.db,keys=names(enrichedBP),columns='TERM',keytype='GOID')
}

