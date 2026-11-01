######################## SNP filtered per step per fold ####################
##Libraries
library(SNPassoc)

load('./Whole_dataset.Rda')
dataset<-dataset[!colnames(dataset) %in% c("MPV","PDW","BMI")]
load('./Test_inds.Rda')
quality_controlSNP<-function(snp_data_train,snp_data_test,cols,genetic_modelling='Codominant',HWE_threshold=0.05,MAF_threshold=1){
  clinical_cols<-setdiff(colnames(snp_data_train),colnames(snp_data_train)[cols])
  snps<-setupSNP(data=snp_data_train,colSNPs=cols,sep="")
  snps_test<-setupSNP(data=snp_data_test,colSNPs=cols,sep="")
  ########### 1. HWE 0.05 p value
  hwe<-tableHWE(snps,Binary.SS)
  #Multiple testing correction: FDR
  p_vals <- hwe[,2]
  fdr_pvals <- p.adjust(p_vals, method = "fdr")
  threshold_p<-HWE_threshold
  snp_filtered<-rownames(hwe)[fdr_pvals < threshold_p | is.na(fdr_pvals)]
  
  wanted_list<-list("HWE"=snp_filtered)
  
  snp_list<-setdiff(rownames(hwe),snp_filtered)
  snp_df<-snps[colnames(snps) %in% c(clinical_cols,snp_list)]
  
  ## 2. MAF < 1% (https://doi.org/10.1038/s41598-024-53310-x)
  summary_ss<-summary(snp_df,print=FALSE)
  MAF_thresh<-100-MAF_threshold
  snps_after_MAF<-rownames(summary_ss)[summary_ss[,2] < MAF_thresh]
  
  
  snps_filtered<-rownames(summary_ss)[summary_ss[,2] > MAF_thresh]
  
  wanted_list[["MAF"]]<-snps_filtered
  
  #no snp filtered with MAF filter (they were filtered before because those with p-values in HWE test, were eliminated)
  snp_df<-snps[colnames(snps) %in% c(clinical_cols,snps_after_MAF)]
  
  
  ## 3. LD R^2>0.5 (10.1002/gepi.20309)
  res<-SNPassoc::LD(snp_df[cols[1]:ncol(snp_df)],which=c("r"))
  r_matrix <- data.frame(res[["r"]]^2)  # Square to get r²
  snp_ld<-list()
  for(snp in colnames(r_matrix)){
    snp_matrix_temp<-r_matrix[snp]
    if(sum(snp_matrix_temp>0.5,na.rm=TRUE)>1){
      
      snp_names<-row.names(snp_matrix_temp)
      
      snps_in_ld_all<-row.names(snp_matrix_temp)[snp_matrix_temp>0.5]
      
      snps_in_ld<-snps_in_ld_all[!is.na(snps_in_ld_all)]
      snp_ld[[snp]]<-snps_in_ld
    } 
  }
  
  groups_snp <- list()
  i <- 0
  for (snp in names(snp_ld)) {
    group_temp <- c(snp, snp_ld[[snp]])
    matched_groups <- c()
    
    for (group in names(groups_snp)) {
      if (any(group_temp %in% groups_snp[[group]])) {
        matched_groups <- c(matched_groups, group)
      }
    }
    
    if (length(matched_groups) == 0) {
      i <- i + 1
      groups_snp[[as.character(i)]] <- group_temp
    } else {
      # Merge all matching groups into one
      merged_snps <- unique(c(group_temp, unlist(groups_snp[matched_groups])))
      
      # Remove old groups
      groups_snp[matched_groups] <- NULL
      
      # Add new merged group
      i <- i + 1
      groups_snp[[as.character(i)]] <- merged_snps
    }
  }
  clean_snps<-setdiff(colnames(snp_df)[cols[1]:ncol(snp_df)],unlist(groups_snp))
  snps_to_check<-groups_snp
  snp_without_ld<-snps[colnames(snps) %in% c(clinical_cols,clean_snps)]
  #606 samples/ 168 SNP columns / 41 clinical/ 3 target
  snp_ld_df<-snps[colnames(snps) %in% c(clinical_cols,unique(unlist(snps_to_check)))]
  #864 samples/ 42 SNP columns in high-LD (>0.5) / 3 target
  ### selected SNPs with high LD from the groups
  selected_snps<-c()
  if(genetic_modelling=='Codominant'){
    ind<-2
  }else if(genetic_modelling=='Dominant'){
    ind<-6
  }else if(genetic_modelling=='Recessive'){
    ind<-9
  }else if(genetic_modelling=='Optimal'){
    ind<-NA
  }
  
  
  for(group in names(snps_to_check)){
    group_snp<-snps_to_check[[group]]
    snp_temp<-snp_ld_df[colnames(snp_ld_df) %in% c('Binary.SS',group_snp)]
    aic_vals<-c()
    for(snp in group_snp){
      r<-SNPassoc::association(as.formula(paste0('Binary.SS ~ ',snp)),data=snp_ld_df)
      if(!is.na(ind)){
        aic_vals<-c(aic_vals,r[ind,"AIC"])
      }else{
        val<-min(r[,9],na.rm=TRUE)
        aic_vals<-c(aic_vals,val)
      }
    }
    snp_select<-group_snp[which.min(aic_vals)]
    selected_snps<-c(selected_snps,snp_select)
  }
  all_snps<-c(selected_snps,clean_snps)
  
  snps_maybe<-unlist(groups_snp)
  eliminated<-setdiff(snps_maybe,selected_snps)
  wanted_list[['LD']]<-eliminated
  ###### Train/test snps set: 693 train / 171 test  (features depend on the seed)
  train_snps<-snps[colnames(snps) %in% c(clinical_cols,all_snps)]
  test_snps<-snps_test[colnames(snps_test) %in% c(clinical_cols,all_snps)]
  
  
  train_unconverted<-snp_data_train[colnames(snp_data_train) %in% c(clinical_cols,all_snps)]
  test_unconverted<-snp_data_test[colnames(snp_data_test) %in% c(clinical_cols,all_snps)]
  
  
  
  
  ######## Modelling - levels to numerical
  snp_cols<-colnames(train_snps)[cols[1]:ncol(train_snps)]
  levels_modelling<-list()
  snps<-c()
  for(snp in snp_cols){
    form<-as.formula(paste('Binary.SS ~',snp))
    res<-association(form,data=train_snps)
    if(length(res[,1])==15){
      snps<-c(snps,snp)
      levels_co<-names(res[,1])[2:4]
      levels_do<-names(res[,1])[6:7]
      levels_re<-names(res[,1])[9:10]
      levels_over<-names(res[,1])[12:13]
      
      if(genetic_modelling=='Codominant'){
        levels_modelling[[snp]]<-levels_co
      }else if(genetic_modelling=='Dominant'){
        levels_modelling[[snp]]<-levels_do
      }else if(genetic_modelling=='Recessive'){
        levels_modelling[[Snp]]<-levels_re
      }else if(genetic_modelling=='Optimal'){
        index<-which.min(res[1:12,9])
        if(index==2){
          levels_opt<-levels_co
        }else if(index==6){
          levels_opt<-levels_do
        }else if(index==9){
          levels_opt<-levels_re
        }else if(index==12){
          levels_opt<-levels_over
        }
        levels_modelling[[snp]]<-levels_opt
      }
      
      
    }
  }
  
  train_converted<-snp_numeric_convert(train_snps,snp_cols,levels_modelling)
  test_converted<-snp_numeric_convert(test_snps,snp_cols,levels_modelling)
  
  
  return(list(`Train_set`=train_converted,`Test_set`=test_converted,`Levels`=levels_modelling,`Train_unconverted`=train_unconverted,
              `Test_uncovered`=test_unconverted,Step.Filtering=wanted_list))
}
snp_numeric_convert <-function(dataset,snp_columns,levels){
  dataset_no_snp<-dataset[! colnames(dataset) %in% snp_columns]
  dataset_snp<-dataset[colnames(dataset) %in% snp_columns]
  for(snp in snp_columns){
    levels_snp<-levels[[snp]]
    if(length(levels_snp)==2){
      level_0<-strsplit(levels_snp[1],split="-")[[1]]
      level_1<-strsplit(levels_snp[2],split="-")[[1]]
      dataset_snp[[snp]]<-as.factor(ifelse(dataset_snp[[snp]] %in% level_0,0,1))
    }else{
      level_0<-levels_snp[1]
      level_1<-levels_snp[2]
      dataset_snp[[snp]]<-as.factor(ifelse(dataset_snp[[snp]] %in% level_0,0,ifelse(dataset_snp[[snp]] %in% level_1,1,2)))
    }
    
  }
  
  df<-cbind(dataset_no_snp,dataset_snp)
  return(df)
}

genetic_modelling<-'Optimal'#Codominant, Optimal
HWE_threshold<-0.05
MAF_threshold<-0.05
clinical_cols<-1:47
snp_cols<-48:ncol(dataset)


step_filtered<-list()
for(jj in 1:10){
  print(jj)
  train_set<-dataset[-list_indexes_test[[as.character(jj)]],]
  test_set<-dataset[list_indexes_test[[as.character(jj)]],]
  res<-quality_controlSNP(train_set,test_set,snp_cols,genetic_modelling=genetic_modelling,HWE_threshold=HWE_threshold,MAF_threshold=MAF_threshold)
  step_filtered[[jj]]<-res$Step.Filtering 
  
  
}
