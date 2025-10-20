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

####### Functions #########
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
Preprocessing_fnc<-function(clinical_train,clinical_test,snp_train,snp_test,target_col='target',factor_col_opt='One-hot',numerical_col_opt='MinMax'){
  
  ########### Factor preprocess
  if(factor_col_opt=='One-hot'){
    dmy <- dummyVars(" ~ .", data = clinical_train)
    trsf_train <- data.frame(predict(dmy, newdata = clinical_train))
    trsf_test<- data.frame(predict(dmy, newdata = clinical_test))
    
    
    dmy_snp <- dummyVars(" ~ .", data = snp_train)
    trsf_trainsnp <- data.frame(predict(dmy_snp, newdata = snp_train))
    trsf_testsnp<- data.frame(predict(dmy_snp, newdata = snp_test))
    
  }else if(factor_col_opt=="L-1_encoding"){
    rec <- recipe(~ ., data = clinical_train) %>%
      step_dummy(all_nominal(), one_hot = FALSE)  # L-1 encoding
    
    prep_rec <- prep(rec)
    trsf_train <- data.frame(bake(prep_rec, new_data = clinical_train))
    trsf_train<-preprocess_baked_structure(trsf_train)
    trsf_test <- data.frame(bake(prep_rec, new_data = clinical_test))
    trsf_test<-preprocess_baked_structure(trsf_test)
    
    
    
    rec <- recipe(~ ., data = snp_train) %>%
      step_dummy(all_nominal(), one_hot = FALSE)  # L-1 encoding
    
    prep_rec <- prep(rec)
    trsf_trainsnp <- data.frame(bake(prep_rec, new_data = snp_train))
    trsf_trainsnp<-preprocess_baked_structure(trsf_trainsnp)
    
    trsf_testsnp <- data.frame(bake(prep_rec, new_data = snp_test))
    trsf_testsnp<-preprocess_baked_structure(trsf_testsnp)
    
    
  }
  
  ########## Numerical preprocess / Scale
  if(numerical_col_opt=='None'){
    clinical_train<-trsf_train
    clinical_test<-trsf_test
    
    snp_train<-trsf_trainsnp
    snp_test<-trsf_testsnp
    scalerClinical_obj<-NULL
    scalerSNP_obj<-NULL
    
  }else if(numerical_col_opt=='MinMax'){
    pre_proc <- preProcess(trsf_train[!colnames(trsf_train) %in% c(target_col)], method = c("range"))
    # Transform the dataset
    clinical_train <- predict(pre_proc, newdata = trsf_train)
    clinical_test <- predict(pre_proc, newdata = trsf_test)
    
    scalerClinical_obj<-pre_proc
    
    pre_procsnp <- preProcess(trsf_trainsnp[!colnames(trsf_trainsnp) %in% c(target_col)], method = c("range"))
    # Transform the dataset
    snp_train <- predict(pre_procsnp, newdata = trsf_trainsnp)
    snp_test <- predict(pre_procsnp, newdata = trsf_testsnp)
    
    scalerSNP_obj<-pre_procsnp
    
  }
  
  
  return(list(`ClinicalTrain`=clinical_train,`ClinicalTest`=clinical_test,`SNPtrain`=snp_train,`SNPtest`=snp_test,`ClinicalScaler`=scalerClinical_obj,
              `ClinicalSNPScaler`=scalerSNP_obj))
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
Preprocessing_fnc_count<-function(clinical_train,clinical_test,snp_train,snp_test,target_col='target',factor_col_opt='One-hot',numerical_col_opt='None'){
  ########### Factor preprocess
  if(factor_col_opt=='One-hot'){
    dmy <- dummyVars(" ~ .", data = clinical_train[!colnames(clinical_train) %in% c(target_col)])
    trsf_train <- data.frame(predict(dmy, newdata = clinical_train))
    trsf_test<- data.frame(predict(dmy, newdata = clinical_test))
    
    trsf_train$SYNTAX.SCORE<-clinical_train$SYNTAX.SCORE
    trsf_train$Binary.SS<-clinical_train$Binary.SS
    
    trsf_test$SYNTAX.SCORE<-clinical_test$SYNTAX.SCORE
    trsf_test$Binary.SS<-clinical_test$Binary.SS
    
    dmy_snp <- dummyVars(" ~ .", data = snp_train[!colnames(snp_train) %in% c(target_col)])
    trsf_trainsnp <- data.frame(predict(dmy_snp, newdata = snp_train))
    trsf_testsnp<- data.frame(predict(dmy_snp, newdata = snp_test))
    
    
    trsf_trainsnp$SYNTAX.SCORE<-snp_train$SYNTAX.SCORE
    trsf_trainsnp$Binary.SS<-snp_train$Binary.SS
    
    trsf_testsnp$SYNTAX.SCORE<-snp_test$SYNTAX.SCORE
    trsf_testsnp$Binary.SS<-snp_test$Binary.SS
    
    
  }else if(factor_col_opt=="L-1_encoding"){
    rec <- recipe(~ ., data = clinical_train) %>%
      step_dummy(all_nominal(), one_hot = FALSE)  # L-1 encoding
    
    prep_rec <- prep(rec)
    trsf_train <- data.frame(bake(prep_rec, new_data = clinical_train))
    trsf_train<-preprocess_baked_structure(trsf_train)
    trsf_test <- data.frame(bake(prep_rec, new_data = clinical_test))
    trsf_test<-preprocess_baked_structure(trsf_test)
    
    
    
    rec <- recipe(~ ., data = snp_train) %>%
      step_dummy(all_nominal(), one_hot = FALSE)  # L-1 encoding
    
    prep_rec <- prep(rec)
    trsf_trainsnp <- data.frame(bake(prep_rec, new_data = snp_train))
    trsf_trainsnp<-preprocess_baked_structure(trsf_trainsnp)
    
    trsf_testsnp <- data.frame(bake(prep_rec, new_data = snp_test))
    trsf_testsnp<-preprocess_baked_structure(trsf_testsnp)
    
    
  }
  
  ########## Numerical preprocess / Scale
  if(numerical_col_opt=='None'){
    clinical_train<-trsf_train
    clinical_test<-trsf_test
    
    snp_train<-trsf_trainsnp
    snp_test<-trsf_testsnp
    scalerClinical_obj<-NULL
    scalerSNP_obj<-NULL
    
  }else if(numerical_col_opt=='MinMax'){
    pre_proc <- preProcess(trsf_train[!colnames(trsf_train) %in% c(target_col)], method = c("range"))
    # Transform the dataset
    clinical_train <- predict(pre_proc, newdata = trsf_train)
    clinical_test <- predict(pre_proc, newdata = trsf_test)
    
    scalerClinical_obj<-pre_proc
    
    pre_procsnp <- preProcess(trsf_trainsnp[!colnames(trsf_trainsnp) %in% c(target_col)], method = c("range"))
    # Transform the dataset
    snp_train <- predict(pre_procsnp, newdata = trsf_trainsnp)
    snp_test <- predict(pre_procsnp, newdata = trsf_testsnp)
    
    scalerSNP_obj<-pre_procsnp
    
  }
  
  
  return(list(`ClinicalTrain`=clinical_train,`ClinicalTest`=clinical_test,`SNPtrain`=snp_train,`SNPtest`=snp_test,`ClinicalScaler`=scalerClinical_obj,
              `ClinicalSNPScaler`=scalerSNP_obj))
}



########## Construction of full test sets ########
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Whole_dataset.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Test_inds.Rda')
test_cases_ls<-list()
for(jj in 1:10){
  dt<-dataset[list_indexes_test[[jj]],]
  test_cases_ls[[jj]]<-dt
}


###### Selected Zero part features ########
### Train sets for Zero part
test_zero_clinical_ls<-list()
test_zero_snp_ls<-list()
### Clinical
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Test_clinical_CV_BA_DM.Rda')
### SNP
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Test_snp_CV_BA_DM.Rda')
for(jj in 1:10){
  orig<-test_cases_ls[[jj]]
  
  zero_cl<-test_clinical_ls[[jj]]
  zero_snp<-test_snp_ls[[jj]]
  
  zero_cl_orig<-orig[colnames(orig) %in% c(colnames(zero_cl),'SYNTAX.SCORE')]
  zero_snp_orig<-orig[colnames(orig) %in% c(colnames(zero_snp),'SYNTAX.SCORE')]
  
  test_zero_clinical_ls[[jj]]<-zero_cl_orig
  test_zero_snp_ls[[jj]]<-zero_snp_orig
}


######### Selected features for Count #######
### Train sets for Zero part
test_count_clinical_ls<-list()
test_count_snp_ls<-list()
### Clinical
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Test_clinical_CV_BA_DM.Rda')
### SNP
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Test_snp_CV_BA_DM.Rda')
for(jj in 1:10){
  orig<-test_cases_ls[[jj]]
  
  count_cl<-test_clinical_ls[[jj]]
  count_snp<-test_snp_ls[[jj]]
  
  count_cl_orig<-orig[colnames(orig) %in% c(colnames(count_cl),'SYNTAX.SCORE')]
  count_snp_orig<-orig[colnames(orig) %in% c(colnames(count_snp),'SYNTAX.SCORE')]
  
  test_count_clinical_ls[[jj]]<-count_cl_orig
  test_count_snp_ls[[jj]]<-count_snp_orig
}
















################## Modelling selection for the rs columns (SNP only) ####################
new_test_snp_ls<-list()
new_train_snp_ls<-list()
for(jj in 1:10){
  load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Whole_dataset.Rda')
  train_jj<-dataset[-list_indexes_test[[jj]],]
  train_zero<-train_jj[ colnames(train_jj)%in%colnames(test_zero_snp_ls[[jj]])]
  
  #Find levels -- encoding
  syntax_index<-which(colnames(train_zero)=='SYNTAX.SCORE')
  snp_cols<-syntax_index:ncol(train_zero)
  zero_subset<-train_zero[snp_cols]
  zero_subset$SYNTAX.SCORE<-as.factor(ifelse(zero_subset$SYNTAX.SCORE==0,0,1))
  snp_dataset<-setupSNP(zero_subset,2:ncol(zero_subset),sep="")
  levels_modelling<-list()
  for(i in 2:ncol(zero_subset)){
    snp<-colnames(zero_subset)[i]
    res<-SNPassoc::association(as.formula(paste0('SYNTAX.SCORE ~ ',snp)),data=snp_dataset)
    levels_co<-names(res[,1])[2:4]
    levels_do<-names(res[,1])[6:7]
    levels_re<-names(res[,1])[9:10]
    levels_over<-names(res[,1])[12:13]
    
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
  ind_snp<-syntax_index+1
  test_set<-test_zero_snp_ls[[jj]]
  snp_test<-setupSNP(test_set,ind_snp:ncol(test_set),sep="")
  test_zero_snp<-snp_numeric_convert(snp_test,colnames(snp_test)[ind_snp:ncol(snp_test)],levels_modelling)
  
  train_zero<-setupSNP(train_zero,ind_snp:ncol(snp_test),sep="")
  train_zero_snp<-snp_numeric_convert(train_zero,colnames(snp_test)[ind_snp:ncol(snp_test)],levels_modelling)
  new_train_snp_ls[[jj]]<-train_zero_snp
  
  new_test_zero<-test_zero_snp_ls[[jj]]
  for(col in colnames(snp_dataset)[2:ncol(snp_dataset)]){
    new_test_zero[[col]]<-test_zero_snp[[col]]
    new_test_snp_ls[[jj]]<-new_test_zero
  }
}
test_zero_snp_ls<-new_test_snp_ls
train_zero_snp_ls<-new_train_snp_ls

new_test_snp_ls<-list()
new_train_snp_ls<-list()
for(jj in 1:10){
  load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Whole_dataset.Rda')
  train_jj<-dataset[-list_indexes_test[[jj]],]
  train_count<-train_jj[ colnames(train_jj)%in%colnames(test_count_snp_ls[[jj]])]
  train_count<-train_count[train_count$SYNTAX.SCORE>0,]
  train_count$Binary.SS<-as.factor(ifelse(train_count$SYNTAX.SCORE<23,0,1))
  train_count<-train_count[complete.cases(train_count),]
  #Find levels -- encoding
  syntax_index<-which(colnames(train_count)=='SYNTAX.SCORE')
  snp_cols<-syntax_index:ncol(train_count)
  if(length(snp_cols)>=3){
  count_subset<-train_count[snp_cols]
  snp_dataset<-setupSNP(count_subset,3:ncol(count_subset),sep="")
  levels_modelling<-list()
  for(i in 3:ncol(count_subset)){
    snp<-colnames(count_subset)[i]
    res<-SNPassoc::association(as.formula(paste0('Binary.SS ~ ',snp)),data=snp_dataset)
    levels_co<-names(res[,1])[2:4]
    levels_do<-names(res[,1])[6:7]
    levels_re<-names(res[,1])[9:10]
    levels_over<-names(res[,1])[12:13]
    
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
  ind_snp<-syntax_index+2
  test_set<-test_count_snp_ls[[jj]]
  snp_test<-setupSNP(test_set,ind_snp:ncol(test_set),sep="")
  test_count_snp<-snp_numeric_convert(snp_test,colnames(snp_test)[ind_snp:ncol(snp_test)],levels_modelling)
  train_count<-setupSNP(train_count,ind_snp:ncol(snp_test),sep="")
  
  train_count_snp<-snp_numeric_convert(train_count,colnames(snp_test)[ind_snp:ncol(snp_test)],levels_modelling)
  new_train_snp_ls[[jj]]<-train_count_snp
  
  new_test_count<-test_count_snp_ls[[jj]]
  for(col in colnames(snp_dataset)[2:ncol(snp_dataset)]){
    new_test_count[[col]]<-test_count_snp[[col]]
    new_test_snp_ls[[jj]]<-new_test_count
  }
  }
  else{
    new_test_snp_ls[[jj]]<- test_count_snp_ls[[jj]]
    new_train_snp_ls[[jj]]<-train_count
  }
}
test_count_snp_ls<-new_test_snp_ls
train_count_snp_ls<-new_train_snp_ls













######################### Preprocessing Numerical and Factor variables ####################
####Zero part + Count part
final_test_clinical_zero_ls<-list()
final_test_snp_zero_ls<-list()

final_test_clinical_count_ls<-list()
final_test_snp_count_ls<-list()
interaction_analysis <-function(train,test){
  cols_snp<-which(grepl("rs",colnames(train)))
  ## Interaction
  if(length(cols_snp)>=2){
    
    snp_cols<-colnames(train)[cols_snp]
    
    for(snp1 in snp_cols){
      rest_cols<-snp_cols[!snp_cols %in% snp1]
      
      for(snp2 in rest_cols){
        interaction_string<-paste0(snp1,'_',snp2)
        
        interaction_col<-interaction(train[,snp1],train[,snp2])
        interaction_col_test<-interaction(test[,snp1],test[,snp2])
        
        train[[interaction_string]]<-interaction_col
        test[[interaction_string]]<-interaction_col_test
        
        
        
        snp_cols<-snp_cols[!snp_cols %in% snp1]
        
      }
    }
    
  }else{
    return(list(`Train`=train,`Test`=test))
  }
  
  
  
  return(list(`Train`=train,`Test`=test))
  
}
for(jj in 1:10){
  load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Whole_dataset.Rda')
  clinical_train_zero<-dataset[-list_indexes_test[[jj]],]
  clinical_test_zero<-test_zero_clinical_ls[[jj]]
  clinical_train_zero<-clinical_train_zero[colnames(clinical_train_zero) %in% colnames(clinical_test_zero)]
  
  
  snp_train_zero<-train_zero_snp_ls[[jj]]
  snp_test_zero<-test_zero_snp_ls[[jj]]
  
  
  #Include interactions#
  res_interaction<-interaction_analysis(snp_train_zero,snp_test_zero)
  snp_train_zero<-res_interaction$Train
  snp_test_zero<-res_interaction$Test
  ########### Interaction include ##########
  
  
  res_preprocessing<-Preprocessing_fnc(clinical_train_zero,clinical_test_zero,snp_train_zero,snp_test_zero,target_col='SYNTAX.SCORE',factor_col_opt='One-hot',numerical_col_opt='MinMax')
  clinical_train<-res_preprocessing$ClinicalTrain
  clinical_test<-res_preprocessing$ClinicalTest
  snp_train<-res_preprocessing$SNPtrain
  snp_test<-res_preprocessing$SNPtest
  
  final_test_clinical_zero_ls[[jj]]<-clinical_test
  final_test_snp_zero_ls[[jj]]<-snp_test
  
  
  
  
  
  load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Whole_dataset.Rda')
  clinical_train_count<-dataset[-list_indexes_test[[jj]],]
  clinical_test_count<-test_count_clinical_ls[[jj]]
  clinical_train_count<-clinical_train_count[colnames(clinical_train_count) %in% colnames(clinical_test_count)]
  
  
  snp_train_count<-train_count_snp_ls[[jj]]
  snp_test_count<-test_count_snp_ls[[jj]]
  
  #Include interactions#
  res_interaction<-interaction_analysis(snp_train_count,snp_test_count)
  snp_train_count<-res_interaction$Train
  snp_test_count<-res_interaction$Test
  ########### Interaction include ##########
  
  res_preprocessing<-Preprocessing_fnc_count(clinical_train_count,clinical_test_count,snp_train_count,snp_test_count,target_col=c('SYNTAX.SCORE','Binary.SS'),factor_col_opt='One-hot',numerical_col_opt='MinMax')
  clinical_train<-res_preprocessing$ClinicalTrain
  clinical_test<-res_preprocessing$ClinicalTest
  snp_train<-res_preprocessing$SNPtrain
  snp_test<-res_preprocessing$SNPtest
  
  final_test_clinical_count_ls[[jj]]<-clinical_test
  final_test_snp_count_ls[[jj]]<-snp_test
  }




save(final_test_clinical_zero_ls,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Whole_model_experiments_classifiers/Zero_part/Clinical_test.Rda')
save(final_test_snp_zero_ls,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Whole_model_experiments_classifiers/Zero_part/SNP_test.Rda')

save(final_test_clinical_count_ls,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Whole_model_experiments_classifiers/Count_part/Clinical_test.Rda')
save(final_test_snp_count_ls,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Whole_model_experiments_classifiers/Count_part/SNP_test.Rda')





######## Model check
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Extra/Models_RF_AUC_None.Rda')
model<-final_models$SNP_RF_2
test<-final_test_snp_zero_ls[[2]]
test<-test[complete.cases(test),]
predict(model,test)








load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Whole_model_experiments_classifiers/Predictors/Predictors_count_snp.Rda')
table(predictors_count_snp)


load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Whole_model_experiments_classifiers/Predictors/Predictors_zero_snp.Rda')
table(predictors_zero_snp)
