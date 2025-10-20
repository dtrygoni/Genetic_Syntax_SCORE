library(SNPassoc)
load('SNP dataset')
data<-dataset
snp_data<-data[,48:ncol(data)]
response<-data$Binary.SS
snps<-setupSNP(data=snp_data,colSNPs=1:ncol(snp_data),sep="")
snps$target<-response
levels_dominant<-list()
levels_recessive<-list()
levels_additive<-list()

for(col_name in colnames(snps)[1:228]){
  print(col_name)
  r<-SNPassoc::association(as.formula(paste0('target ~',col_name)),data=snps)
  
  ##Dominant level
  if(length(row.names(r))>10){
    
  dom_levels<-rownames(r[6:7,])
  levels_dominant[[col_name]]<- dom_levels
  
  
  ## Recessive level
  rec_levels<-rownames(r[9:10,])
  levels_recessive[[col_name]]<-rec_levels
  
  ## Additive level
  additiv_levels<-rownames(r[2:4,])
  levels_additive[[col_name]]<-additiv_levels
  }
}




                              
