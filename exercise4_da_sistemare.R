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
  doppi<-duplicated(gene_GO[2:3])
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
enrichedGO<-function(Gopvalues){
  #funzione che dovrà selezionare i GO enriched, da scrivere quando si decide come selezionare
  return(Gopvalues[Gopvalues<0.3])
}
####################################################################

#Main #

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
enrichedCC<-enrichedGO(pvaluesCC)
enrichedBP<-enrichedGO(pvaluesBP)
#extraction of the corresponding GO term
descriptionMF<-select(GO.db,keys=names(enrichedMF),columns='TERM',keytype='GOID')
descriptionCC<-select(GO.db,keys=names(enrichedCC),columns='TERM',keytype='GOID')
descriptionBP<-select(GO.db,keys=names(enrichedBP),columns='TERM',keytype='GOID')

