###################################### Cross validation for Count_part Classifier ############################
##Libraries
library(missForestPredict)
library(Boruta)
library(caret)
library(caret)
library(ggplot2)
library(dplyr)
library(tidyr)
library(randomForest)
library(xgboost)
library(PRROC)
library(precrec)
library(mltools)
library(precrec)
library(pROC)
library(SNPassoc)
library(tibble)
library(recipes)
library(SNPassoc)
set.seed(123)
xgb.set.config(verbosity = 0)

model_arr<-c()
type_arr<-c()
acc_arr<-c()
prec_arr<-c()
rec_arr<-c()
f1_arr<-c()
aucpr_arr<-c()
aucroc_arr<-c()
fold_arr<-c()
final_models<-list()


train_clinical_ls<-list()
test_clinical_ls<-list()
train_snp_ls<-list()
test_snp_ls<-list()
predictors_clinical_ls<-list()
predictors_snp_ls<-list()


###### load
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Whole_dataset.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Test_inds.Rda')


##### Functions
### Helper functions
summary_fnc <- function(data, lev = NULL, model = NULL) {
  # Ensure prediction and true labels are factors
  obs <- factor(data$obs, levels = lev)
  pred <- factor(data$pred, levels = lev)
  
  cm <- caret::confusionMatrix(pred, obs)
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  
  bal_acc <- mean(c(sens, spec), na.rm = TRUE)
  acc<-MLmetrics::Accuracy(pred,obs)
  f1<-FBeta_Score(obs,pred,positive="1",beta=1)
  out <- c(balAcc = bal_acc,F1Score=f1,Spec=spec,Sens=sens,Accuracy=acc)
  return(out)
}
preprocess_baked_structure<-function(data){
  df_clean <- as.data.frame(lapply(data, function(x) {
    if (is.factor(x)) as.numeric(as.character(x))  # convert factors to numeric (via character to keep level values)
    else if (is.matrix(x)) as.numeric(x)           # flatten any matrix columns
    else x
  }))
  
  # Remove all attributes (like "assign", "contrasts")
  df_clean <- as.data.frame(lapply(df_clean, function(x) {
    attributes(x) <- NULL
    return(x)
  }))
  return(df_clean)
  
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
  snp_list<-setdiff(rownames(hwe),snp_filtered)
  snp_df<-snps[colnames(snps) %in% c(clinical_cols,snp_list)]
  ## 2. MAF < 1% (https://doi.org/10.1038/s41598-024-53310-x)
  summary_ss<-summary(snp_df,print=FALSE)
  MAF_thresh<-100-MAF_threshold
  snps_after_MAF<-rownames(summary_ss)[summary_ss[,2] < MAF_thresh]
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
  ###### Train/test snps set: 693 train / 171 test  (features depend on the seed)
  train_snps<-snps[colnames(snps) %in% c(clinical_cols,all_snps)]
  test_snps<-snps_test[colnames(snps_test) %in% c(clinical_cols,all_snps)]
  
  
  
  
  
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
  
  
  return(list(`Train_set`=train_converted,`Test_set`=test_converted,`Levels`=levels_modelling))
}
quality_controlSNP_reg<-function(snp_data_train,snp_data_test,cols,genetic_modelling='Codominant',HWE_threshold=0.05,MAF_threshold=1){
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
  snp_list<-setdiff(rownames(hwe),snp_filtered)
  snp_df<-snps[colnames(snps) %in% c(clinical_cols,snp_list)]
  ## 2. MAF < 1% (https://doi.org/10.1038/s41598-024-53310-x)
  summary_ss<-summary(snp_df,print=FALSE)
  MAF_thresh<-100-MAF_threshold
  snps_after_MAF<-rownames(summary_ss)[summary_ss[,2] < MAF_thresh]
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
        val<-min(r[,8],na.rm=TRUE)
        aic_vals<-c(aic_vals,val)
      }
    }
    snp_select<-group_snp[which.min(aic_vals)]
    selected_snps<-c(selected_snps,snp_select)
  }
  all_snps<-c(selected_snps,clean_snps)
  ###### Train/test snps set: 693 train / 171 test  (features depend on the seed)
  train_snps<-snp_df[colnames(snp_df) %in% c(clinical_cols,all_snps)]
  test_snps<-snps_test[colnames(snps_test) %in% c(clinical_cols,all_snps)]
  
  
  
  
  
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
        index<-which.min(res[1:12,8])
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
  
  
  return(list(`Train_set`=train_converted,`Test_set`=test_converted,`Levels`=levels_modelling))
}
Univariate_LR_filtering<-function(train,test,predictors_clinical,predictors_snp,target_col,LR_threshold=0.10){
  clinical_predictors<-c()
  snp_predictors<-c()

  for(predictor in predictors_clinical){
    form <- as.formula(paste(target_col, '~', paste0('`', predictor, '`')))
    res<-glm(form,family='binomial',data=train)
    pval <- summary(res)$coefficients[2, 4]
    if(pval<0.10){
      clinical_predictors<-c(clinical_predictors,predictor)
    }
  }
  
  ### For all modelling
  for(predictor in predictors_snp){
    form <- as.formula(paste(target_col, '~', paste0('`', predictor, '`')))
    res<-glm(form,family='binomial',data=train)
    pval<-summary(res)$coefficients[2,4]
    if(pval<LR_threshold){
      snp_predictors<-c(snp_predictors,predictor)
    }
  }
  
  
  #### Clinical sets
  train_clinical<-train[colnames(train) %in% c(clinical_predictors,'Binary.SS','SYNTAX.SCORE')]
  train_clinical<-train_clinical[complete.cases(train_clinical),]
  test_clinical<-test[colnames(test) %in% c(clinical_predictors,'Binary.SS','SYNTAX.SCORE')]
  test_clinical<-test_clinical[complete.cases(test_clinical),]
  #### SNP sets
  train_snp<-train[colnames(train) %in% c(clinical_predictors,snp_predictors,'Binary.SS','SYNTAX.SCORE')]
  train_snp<-train_snp[complete.cases(train_snp),]
  test_snp<-test[colnames(test) %in% c(clinical_predictors,snp_predictors,'Binary.SS','SYNTAX.SCORE')]
  test_snp<-test_snp[complete.cases(test_snp),]
  
  return(list(`Train_clinical`=train_clinical,`Train_SNP`=train_snp,`Test_clinical`=test_clinical,`Test_SNP`=test_snp,`Clinical_predictors`=clinical_predictors,
              `SNP_predictors`=snp_predictors))
    
}
Univariate_filtering_reg<-function(train,test,predictors_clinical,predictors_snp,target_col,LR_threshold=0.10){
  clinical_predictors<-c()
  snp_predictors<-c()
  
  for(predictor in predictors_clinical){
    if(is.numeric(train[[predictor]])){
      r<-cor.test(train[[predictor]], train[[target_col]], method = "spearman")
      pval<-r$p.value
    }else if(length(levels(train[[predictor]]))<=2){
      form <- as.formula(paste(target_col, '~', paste0('`', predictor, '`')))
      r<-wilcox.test(form, data = train)
      pval<-r$p.value
    }else{
      form <- as.formula(paste(target_col, '~', paste0('`', predictor, '`')))
      r<-kruskal.test(form,data = train)
      pval<-r$p.value
    }
    
    if(pval<LR_threshold){
      clinical_predictors<-c(clinical_predictors,predictor)
    }
  }
  
  ### For all modelling
  for(predictor in predictors_snp){
    form <- as.formula(paste(target_col, '~', paste0('`', predictor, '`')))
    if(length(levels(train[[predictor]]))>1){
      if(is.numeric(train[[predictor]])){
        r<-cor.test(train[[predictor]], train[[target_col]], method = "spearman")
        pval<-r$p.value
      }else if(length(levels(train[[predictor]]))<=2){
        form <- as.formula(paste(target_col, '~', paste0('`', predictor, '`')))
        r<-wilcox.test(form, data = train)
        pval<-r$p.value
      }else{
        form <- as.formula(paste(target_col, '~', paste0('`', predictor, '`')))
        r<-kruskal.test(form,data = train)
        pval<-r$p.value
      }
      if(pval<LR_threshold){
        snp_predictors<-c(snp_predictors,predictor)
      }
    }
  }
  
  
  #### Clinical sets
  train_clinical<-train[colnames(train) %in% c(clinical_predictors,'Binary.SS','SYNTAX.SCORE')]
  train_clinical<-train_clinical[complete.cases(train_clinical),]
  test_clinical<-test[colnames(test) %in% c(clinical_predictors,'Binary.SS','SYNTAX.SCORE')]
  test_clinical<-test_clinical[complete.cases(test_clinical),]
  #### SNP sets
  train_snp<-train[colnames(train) %in% c(clinical_predictors,snp_predictors,'Binary.SS','SYNTAX.SCORE')]
  train_snp<-train_snp[complete.cases(train_snp),]
  test_snp<-test[colnames(test) %in% c(clinical_predictors,snp_predictors,'Binary.SS','SYNTAX.SCORE')]
  test_snp<-test_snp[complete.cases(test_snp),]
  
  return(list(`Train_clinical`=train_clinical,`Train_SNP`=train_snp,`Test_clinical`=test_clinical,`Test_SNP`=test_snp,`Clinical_predictors`=clinical_predictors,
              `SNP_predictors`=snp_predictors))
  
}
Boruta_filtering<-function(train,test,clinical_predictors,snp_predictors,target_col,strategy='Seperate',Boruta_threshold=0.05){
  ####Clinical
  train_clinical<-train[,c(clinical_predictors,target_col)]
  boruta_res<-Boruta(as.formula(paste(target_col,' ~ .')),data=train_clinical,pValue=Boruta_threshold,maxRuns=1000)
  clinical_cols_zero <-colnames(train_clinical)[boruta_res$finalDecision=='Confirmed']
  
  if(strategy=='Seperate'){
    print('Here')
    snp_boruta<-train[,c(snp_predictors,target_col)]
    boruta_res<-Boruta(Binary.SS ~ .,data=snp_boruta,pValue=Boruta_threshold,maxRuns=1000)
    snp_boruta_cols_zero<-colnames(snp_boruta)[boruta_res$finalDecision=='Confirmed']
    
    
    clinical_train<-train[colnames(train) %in% c(setdiff(clinical_cols_zero,'Binary.SS'),'Binary.SS')]
    clinical_test<-test[colnames(test) %in% c(setdiff(clinical_cols_zero,'Binary.SS'),'Binary.SS')]
    
    snp_train<-train[colnames(train) %in% c(setdiff(clinical_cols_zero,'Binary.SS'),setdiff(snp_boruta_cols_zero,'Binary.SS'),'Binary.SS')]
    snp_test<-test[colnames(test) %in% c(setdiff(clinical_cols_zero,'Binary.SS'),setdiff(snp_boruta_cols_zero,'Binary.SS'),'Binary.SS')]
    
    
    clinical_train$target <-as.numeric(ifelse(clinical_train$Binary.SS==0,0,1))
    clinical_test$target <-as.numeric(ifelse(clinical_test$Binary.SS==0,0,1))
    snp_train$target <-as.numeric(ifelse(snp_train$Binary.SS==0,0,1))
    snp_test$target <-as.numeric(ifelse(snp_test$Binary.SS==0,0,1))
    
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('Binary.SS')]
    clinical_test<-clinical_test[!colnames(clinical_test) %in% c('Binary.SS')]
    snp_train<-snp_train[!colnames(snp_train) %in% c('Binary.SS')]
    snp_test<-snp_test[!colnames(snp_test) %in% c('Binary.SS')]
    
  }else if(strategy=='All'){
    snp_boruta<-train[,c(clinical_predictors,snp_predictors,target_col)]
    boruta_res<-Boruta(Binary.SS ~ .,data=snp_boruta,pValue=Boruta_threshold,maxRuns=1000)
    snp_boruta_cols_zero<-colnames(snp_boruta)[boruta_res$finalDecision=='Confirmed']
    
    
    clinical_train<-train[colnames(train) %in% c(setdiff(clinical_cols_zero,'Binary.SS'),'Binary.SS')]
    clinical_test<-test[colnames(test) %in% c(setdiff(clinical_cols_zero,'Binary.SS'),'Binary.SS')]
    
    snp_train<-train[colnames(train) %in% c(setdiff(snp_boruta_cols_zero,'Binary.SS'),'Binary.SS')]
    snp_test<-test[colnames(test) %in% c(setdiff(snp_boruta_cols_zero,'Binary.SS'),'Binary.SS')]
    
    
    clinical_train$target <-as.numeric(ifelse(clinical_train$Binary.SS==0,0,1))
    clinical_test$target <-as.numeric(ifelse(clinical_test$Binary.SS==0,0,1))
    snp_train$target <-as.numeric(ifelse(snp_train$Binary.SS==0,0,1))
    snp_test$target <-as.numeric(ifelse(snp_test$Binary.SS==0,0,1))
    
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('Binary.SS')]
    clinical_test<-clinical_test[!colnames(clinical_test) %in% c('Binary.SS')]
    snp_train<-snp_train[!colnames(snp_train) %in% c('Binary.SS')]
    snp_test<-snp_test[!colnames(snp_test) %in% c('Binary.SS')]
  } 
 
  
  return(list(`Train_clinical`=clinical_train,`Train_SNP`=snp_train,`Test_clinical`=clinical_test,`Test_SNP`=snp_test,
              `Clinical_predictors`=clinical_cols_zero,`SNP_predictors`=setdiff(snp_boruta_cols_zero,clinical_cols_zero)))
}
Boruta_filtering_reg<-function(train,test,clinical_predictors,snp_predictors,target_col,strategy='Seperate',Boruta_threshold=0.05){
  ####Clinical
  train_clinical<-train[,c(clinical_predictors,target_col)]
  boruta_res<-Boruta(as.formula(paste(target_col,' ~ .')),data=train_clinical,pValue=Boruta_threshold,maxRuns=1000)
  clinical_cols_zero <-c(colnames(train_clinical)[boruta_res$finalDecision=='Confirmed'],'SYNTAX.SCORE')
  
  if(strategy=='Seperate'){
    print('Here')
    snp_boruta<-train[,c(snp_predictors,target_col)]
    boruta_res<-Boruta(SYNTAX.SCORE ~ .,data=snp_boruta,pValue=Boruta_threshold,maxRuns=1000)
    snp_boruta_cols_zero<-colnames(snp_boruta)[boruta_res$finalDecision=='Confirmed']
    
    
    clinical_train<-train[colnames(train) %in% c(setdiff(clinical_cols_zero,'Binary.SS'),'Binary.SS')]
    clinical_test<-test[colnames(test) %in% c(setdiff(clinical_cols_zero,'Binary.SS'),'Binary.SS')]
    
    snp_train<-train[colnames(train) %in% c(setdiff(clinical_cols_zero,'Binary.SS'),setdiff(snp_boruta_cols_zero,'Binary.SS'),'Binary.SS')]
    snp_test<-test[colnames(test) %in% c(setdiff(clinical_cols_zero,'Binary.SS'),setdiff(snp_boruta_cols_zero,'Binary.SS'),'Binary.SS')]
    
    
    clinical_train$target <-as.numeric(ifelse(clinical_train$Binary.SS==0,0,1))
    clinical_test$target <-as.numeric(ifelse(clinical_test$Binary.SS==0,0,1))
    snp_train$target <-as.numeric(ifelse(snp_train$Binary.SS==0,0,1))
    snp_test$target <-as.numeric(ifelse(snp_test$Binary.SS==0,0,1))
    
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('Binary.SS')]
    clinical_test<-clinical_test[!colnames(clinical_test) %in% c('Binary.SS')]
    snp_train<-snp_train[!colnames(snp_train) %in% c('Binary.SS')]
    snp_test<-snp_test[!colnames(snp_test) %in% c('Binary.SS')]
    
  }else if(strategy=='All'){
    snp_boruta<-train[,c(clinical_predictors,snp_predictors,target_col)]
    boruta_res<-Boruta(as.formula(paste(target_col,'~ .')),data=snp_boruta,pValue=Boruta_threshold,maxRuns=1000)
    snp_boruta_cols_zero<-c(colnames(snp_boruta)[boruta_res$finalDecision=='Confirmed'],'SYNTAX.SCORE')
    
    
    clinical_train<-train[colnames(train) %in% c(setdiff(clinical_cols_zero,c('SYNTAX.SCORE',target_col)),'SYNTAX.SCORE',target_col)]
    clinical_test<-test[colnames(test) %in% c(setdiff(clinical_cols_zero,c('SYNTAX.SCORE',target_col)),'SYNTAX.SCORE',target_col)]
    
    snp_train<-train[colnames(train) %in% c(setdiff(snp_boruta_cols_zero,c('SYNTAX.SCORE',target_col)),'SYNTAX.SCORE',target_col)]
    snp_test<-test[colnames(test) %in% c(setdiff(snp_boruta_cols_zero,c('SYNTAX.SCORE',target_col)),'SYNTAX.SCORE',target_col)]
    
    
    # clinical_train$target <-as.numeric(clinical_train$SYNTAX.SCORE)
    # clinical_test$target <-as.numeric(clinical_test$SYNTAX.SCORE)
    # snp_train$target <-as.numeric(snp_train$SYNTAX.SCORE)
    # snp_test$target <-as.numeric(snp_test$SYNTAX.SCORE)
    # 
    # clinical_train<-clinical_train[!colnames(clinical_train) %in% c('SYNTAX.SCORE')]
    # clinical_test<-clinical_test[!colnames(clinical_test) %in% c('SYNTAX.SCORE')]
    # snp_train<-snp_train[!colnames(snp_train) %in% c('SYNTAX.SCORE')]
    # snp_test<-snp_test[!colnames(snp_test) %in% c('SYNTAX.SCORE')]
  } 
  
  
  return(list(`Train_clinical`=clinical_train,`Train_SNP`=snp_train,`Test_clinical`=clinical_test,`Test_SNP`=snp_test,
              `Clinical_predictors`=clinical_cols_zero,`SNP_predictors`=setdiff(snp_boruta_cols_zero,clinical_cols_zero)))
}


### Parameters
genetic_modelling<-'Optimal'#Codominant, Optimal
HWE_threshold<-0.05
MAF_threshold<-0.05
clinical_cols<-1:47
LR_threshold<-0.10
Boruta_threshold<-0.05
strategy_boruta<-'All'# 'Seperate','All'
target_col<-'Binary.SS' #It is the binary.risk
snp_cols<-48:ncol(dataset)
factor_col_opt='One-hot' #One-hot, L-1_encoding
numerical_col_opt='MinMax' #None, MinMax
optimazation_metric="balAcc" #balAcc,F1Score,Spec.Specificity, Sens.Sensitivity  Accuracy
model_selection='RF'
number_iterations=50
#### select split or all test together with for

for(jj in 1:10){
  print(jj)
  train_set<-dataset[-list_indexes_test[[as.character(jj)]],]
  test_set<-dataset[list_indexes_test[[as.character(jj)]],]

  train_set<-train_set[!is.na(train_set$SYNTAX.SCORE),]
  test_set<-test_set[!is.na(test_set$SYNTAX.SCORE),]
  
  train_set<-train_set[train_set$SYNTAX.SCORE>0,]
  test_set<-test_set[test_set$SYNTAX.SCORE>0,]
  train_set$Binary.SS<-as.factor(ifelse(train_set$SYNTAX.SCORE>=23,1,0))
  test_set$Binary.SS<-as.factor(ifelse(test_set$SYNTAX.SCORE>=23,1,0))
  
  ############### SNP quality control ################
  res_control<-quality_controlSNP_reg(train_set,test_set,snp_cols,genetic_modelling=genetic_modelling,HWE_threshold=HWE_threshold,MAF_threshold=MAF_threshold)
  train_set_all<-res_control$Train_set
  test_set_all<-res_control$Test_set
  ############### First step - Logistic Regression -> selected important features ########
  predictors_clinical<-setdiff(colnames(train_set_all)[clinical_cols],c('SYNTAX.SCORE','SYNTAX.Tertiles','Binary.SS'))
  predictors_snp<-colnames(train_set_all)[snp_cols[1]:ncol(train_set_all)]
  
  res_lr<-Univariate_LR_filtering(train_set_all,test_set_all,predictors_clinical=predictors_clinical,predictors_snp=predictors_snp,
                                  target_col = target_col)
  clinical_predictors<-res_lr$Clinical_predictors
  snp_predictors<-res_lr$SNP_predictors
  
  train_clinical<-res_lr$Train_clinical
  test_clinical<-res_lr$Test_clinical
  
  train_snp<-res_lr$Train_SNP
  test_snp<-res_lr$Test_SNP
  ############### Boruta ################
  print("Boruta")
  res_boruta<-Boruta_filtering_reg(train_snp,test_snp,clinical_predictors = clinical_predictors,snp_predictors = snp_predictors,
                               strategy=strategy_boruta,target_col = target_col)
  
  clinical_train<-res_boruta$Train_clinical
  clinical_test<-res_boruta$Test_clinical
  snp_train<-res_boruta$Train_SNP
  snp_test<-res_boruta$Test_SNP
  cl_preds<-res_boruta$Clinical_predictors
  snp_preds<-res_boruta$SNP_predictors
  
  train_clinical_ls[[jj]]<-clinical_train
  test_clinical_ls[[jj]]<-clinical_test
  train_snp_ls[[jj]]<-snp_train
  test_snp_ls[[jj]]<-snp_test
  predictors_clinical_ls[[jj]]<-cl_preds
  predictors_snp_ls[[jj]]<-snp_preds

  
}


save(train_clinical_ls,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Train_clinical_CV_BA_DM.Rda')
save(test_clinical_ls,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Test_clinical_CV_BA_DM.Rda')
save(train_snp_ls,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Train_snp_CV_BA_DM.Rda')
save(test_snp_ls,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Test_snp_CV_BA_DM.Rda')



######################################### Clinical + SNP predictors Initial --- Count Part ###################

clinical_predictors_all<-c()
snp_predictors_all<-c()
for(jj in 1:10){
  print(jj)
  train_set<-dataset[-list_indexes_test[[as.character(jj)]],]
  test_set<-dataset[list_indexes_test[[as.character(jj)]],]
  
  train_set<-train_set[!is.na(train_set$SYNTAX.SCORE),]
  test_set<-test_set[!is.na(test_set$SYNTAX.SCORE),]
  
  train_set<-train_set[train_set$SYNTAX.SCORE>0,]
  test_set<-test_set[test_set$SYNTAX.SCORE>0,]
  train_set$Binary.SS<-as.factor(ifelse(train_set$SYNTAX.SCORE>=23,1,0))
  test_set$Binary.SS<-as.factor(ifelse(test_set$SYNTAX.SCORE>=23,1,0))
  
  ############### SNP quality control ################
  res_control<-quality_controlSNP_reg(train_set,test_set,snp_cols,genetic_modelling=genetic_modelling,HWE_threshold=HWE_threshold,MAF_threshold=MAF_threshold)
  train_set_all<-res_control$Train_set
  test_set_all<-res_control$Test_set
  ############### First step - Logistic Regression -> selected important features ########
  predictors_clinical<-setdiff(colnames(train_set_all)[clinical_cols],c('SYNTAX.SCORE','SYNTAX.Tertiles','Binary.SS'))
  predictors_snp<-colnames(train_set_all)[snp_cols[1]:ncol(train_set_all)]
  
  res_lr<-Univariate_LR_filtering(train_set_all,test_set_all,predictors_clinical=predictors_clinical,predictors_snp=predictors_snp,
                                  target_col = target_col)
  clinical_predictors<-res_lr$Clinical_predictors
  snp_predictors<-res_lr$SNP_predictors
  
  clinical_predictors_all<-c(clinical_predictors_all,clinical_predictors)
  snp_predictors_all<-c(snp_predictors_all,snp_predictors)
  
}




























