#GO Enrichment analysis FUNCTIONS#

#function that assign the total num of genes and selected genes to each GO term
assigned<-function(value,go,sel)
{
  index<-which(go==value)
  num_go<-length(index)
  num_sel<-length(which(sel[index]==TRUE))
  mat<-rbind(num_go,num_sel)
  return(mat)
}

#Database call to get the GO annotations
call_db<-function(DATA_log,p_value,alfa)
{
  #org.Hs.eg
  gene_symbol<-row.names(DATA_log)
  gene_keys<-mapIds(org.Hs.eg.db,gene_symbol,'ENTREZID','SYMBOL')
  gene_GO<-select(org.Hs.eg.db,gene_keys,columns = c("SYMBOL","GO"))
  
  selected_symbol<-gene_symbol[which(p_value<=alfa)]
  SELECTED<-""
  gene_GO<-cbind(gene_GO,SELECTED)
  gene_GO$SELECTED<-gene_GO$SYMBOL%in%selected_symbol
  gene_GO_NA<-na.omit(gene_GO)
  return(gene_GO_NA)
}

#Elimination of single gene-GO terms and creation of unique GO term list
best<-function(gene_GO_NA)
{
  GO_present<-duplicated(gene_GO_NA$GO)
  GO<-gene_GO_NA$GO[GO_present==FALSE] 
  
  count<-sapply(GO,assigned,gene_GO_NA$GO,gene_GO_NA$SELECTED)
  count<-t(count)
  colnames(count)<-c("N_GENES","N_SELECTED")
  
  index_one<-which(count[,1]==1)
  FINAL<-count[-index_one,]
  return(FINAL)
}

#Total amount of genes and selected genes needed in Fisher Test
total<-function(gene_GO_NA,FINAL)
{
  #numero totale geni
  GO_names<-gene_GO_NA$GO%in%rownames(FINAL)
  ind<-which(GO_names==FALSE)
  gene_GO_NA<-gene_GO_NA[-ind,]
  
  gene_tot<-duplicated(gene_GO_NA$SYMBOL)
  ind_gene<-which(gene_tot==FALSE)
  genes<-length(ind_gene)
  
  #numero selected totali 
  selected_tot<-length(which(gene_GO_NA$SELECTED[ind_gene]==TRUE))
  return(c(genes,selected_tot))
}

#function to compute Fisher Exact test
fisher_test<-function(data,selected_tot,genes)
{
  #initialize variables for the test
  tot_ac<- data[1]
  a<- data[2]
  c<- tot_ac-a
  b<- selected_tot-a
  d<- genes-c-selected_tot 
  
  mat<-matrix(c(a,c,b,d),2,2)
  res<-fisher.test (mat, alternative="greater")  
  return(res$p.value)
}
####################################################################

# Main #

library(org.Hs.eg.db) #Homo sapiens database call

#Creation of GO database matrix
go_mat<-call_db(DATA_log,pvalues,alpha_final)

#Unique GO terms matrix and associated genes
FINAL_GO<-best(go_mat)
FINAL_GO<-FINAL_GO[-which(FINAL_GO[,2]<2),] # elimination of a<2 selected genes
nGO<-nrow(FINAL_GO)

#Final matrix to be used in Fisher test
tot<-total(go_mat,FINAL_GO)
tot_genes<-tot[1]
tot_selected<-tot[2]

# Fisher test on each row of the GO matrix
p_value<-matrix(0,lnGO,1)
p_value<-apply(FINAL_GO,1,fisher_test,tot_selected,tot_genes)

