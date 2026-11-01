rm(list=ls())
######################### Model results ###################
library(caret)
library(PRROC)
library(pROC)
load('./Sensitivity_Analysis_Imputation/Positive_score_component/HighRisk_regressor/Regression_Clinical_High_Models.Rda')
load('./Sensitivity_Analysis_Imputation/Positive_score_component/HighRisk_regressor/Regression_SNP_High_Models.Rda')
#names: regression_cl_model_high, regression_snp_model_high



load('./Sensitivity_Analysis_Imputation/Positive_score_component/LowRisk_regressor/Regression_Clinical_Low_Models.Rda')
load('./Sensitivity_Analysis_Imputation/Positive_score_component/LowRisk_regressor/Regression_SNP_Low_Models.Rda')
#names: regression_cl_model_low, regression_snp_model_low


load('./Sensitivity_Analysis_Imputation/Positive_score_component/Risk_classifier/Classifier_Models.Rda')
#name: classifier_models


load('./Sensitivity_Analysis_Imputation/Positive_score_component/Train_clinical_count.Rda')
load('./Sensitivity_Analysis_Imputation/Positive_score_component/Train_snp_count.Rda')
load('./Sensitivity_Analysis_Imputation/Positive_score_component/Test_clinical_count.Rda')
load('./Sensitivity_Analysis_Imputation/Positive_score_component/Test_snp_count.Rda')






###################### Preprocess test sets #################
######################### Model results ###################
library(caret)
library(PRROC)
library(pROC)




#load('C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_Imputation/Imputed_test_clinical.Rda')
#load('C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_Imputation/Imputed_test_snp.Rda')
load('./Sensitivity_Analysis_Imputation/Positive_score_component/Test_clinical_count.Rda')
load('./Sensitivity_Analysis_Imputation/Positive_score_component/Test_snp_count.Rda')
load('./Sensitivity_Analysis_Imputation/Positive_score_component/Train_clinical_count.Rda')
load('./Sensitivity_Analysis_Imputation/Positive_score_component/Train_snp_count.Rda')
### Prepare clinical sets
Preprocessing_fnc<-function(clinical_train,clinical_test,snp_train,snp_test,target_col='target',factor_col_opt='One-hot',numerical_col_opt='None'){
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
    pre_proc <- preProcess(trsf_train, method = c("range"))
    # Transform the dataset
    clinical_train <- predict(pre_proc, newdata = trsf_train)
    clinical_test <- predict(pre_proc, newdata = trsf_test)
    
    scalerClinical_obj<-pre_proc
    
    pre_procsnp <- preProcess(trsf_trainsnp, method = c("range"))
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




factor_col_opt='One-hot' #One-hot, L-1_encoding
numerical_col_opt='MinMax' #None, MinMax
train_clinical_ls_processed<-list()
test_clinical_ls_processed<-list()
for(i in 1:10){
  res_preprocessing<-Preprocessing_fnc(train_clinical_ls[[i]],test_clinical_ls[[i]],train_snp_ls[[i]],test_snp_ls[[i]],factor_col_opt=factor_col_opt,numerical_col_opt=numerical_col_opt)
  train_set<-res_preprocessing$ClinicalTrain
  test_set<-res_preprocessing$ClinicalTest
  
  train_clinical_ls_processed[[i]]<-train_set
  test_clinical_ls_processed[[i]]<-test_set
}




### Prepare snps sets

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
train_snp_ls_processed<-list()
test_snp_ls_processed<-list()
for(i in 1:10){
  res_interaction<-interaction_analysis(train_snp_ls[[i]],test_snp_ls[[i]])
  snp_train<-res_interaction$Train
  snp_test<-res_interaction$Test
  res_preprocessing<-Preprocessing_fnc(train_clinical_ls[[i]],test_clinical_ls[[i]],snp_train,snp_test,factor_col_opt=factor_col_opt,numerical_col_opt=numerical_col_opt)
  
  train_snp_ls_processed[[i]]<-res_preprocessing$SNPtrain
  test_snp_ls_processed[[i]]<-res_preprocessing$SNPtest
}




target_arr<-c()
fold_arr<-c()
predicted_snp<-c()
predicted_clinical<-c()
library(pROC)
for(fold in 1:10){
  ############ SNP ##################
  classifier_snp<-classifier_models[paste0("SNP_AUC_None_",fold)]
  low_regression_snp<-regression_snp_model_low[paste0("SNP_AUC_None_",fold)]
  high_regression_snp<-regression_snp_model_high[paste0("SNP_AUC_None_",fold)]
  
  test_snp<-test_snp_ls_processed[[fold]]
  #fold_1<-rep(fold,nrow(test_clinical))
  #fold_arr<-c(fold_arr,fold_1)
  
  #target_fold<-imputated_train_clinical[[fold]]$SYNTAX.SCORE
  #target_arr<-c(target_arr,target_fold)
  
  #### Find optimal threshold ##
  target_train<-train_snp_ls[[fold]]$SYNTAX.SCORE
  target_binary<-ifelse(target_train>23,1,0)
  
  train_snp<-train_snp_ls_processed[[fold]]
  preds_probs_snp<-predict(classifier_snp,newdata=train_snp,type="prob")[[1]]
  
  roc_snp<-roc(target_binary,preds_probs_snp[,"pos"])
  best <- pROC::coords(
    roc_snp,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  thresh<-best$threshold
  
  pred_probs<-predict(classifier_snp,newdata=test_snp,type="prob")[[1]]
  preds_class<-ifelse(pred_probs[,2]>thresh,1,0)  
  
  regression_vals_low<-as.numeric(predict(low_regression_snp,newdata=test_snp)[[1]])
  regression_vals_high<-as.numeric(predict(high_regression_snp,newdata=test_snp)[[1]])
  
  
  predicted_val_overall<-round(ifelse(preds_class==0,regression_vals_low,regression_vals_high),2)
  predicted_snp<-c(predicted_snp,predicted_val_overall)
  
  ############ Clinical ##################
  classifier_clinical<-classifier_models[paste0("Clinical_AUC_None_",fold)]
  low_regression_clinical<-regression_cl_model_low[paste0("Clinical_AUC_None_",fold)]
  high_regression_clinical<-regression_cl_model_high[paste0("Clinical_AUC_None_",fold)]

  test_clinical<-test_clinical_ls_processed[[fold]]
  fold_1<-rep(fold,nrow(test_clinical))
  fold_arr<-c(fold_arr,fold_1)
  
  target_fold<-test_clinical_ls[[fold]]$SYNTAX.SCORE
  target_arr<-c(target_arr,target_fold)
  
  #### Find optimal threshold ##
  target_train<-train_clinical_ls[[fold]]$SYNTAX.SCORE
  target_binary<-ifelse(target_train>23,1,0)
  
  train_clinical<-train_clinical_ls_processed[[fold]]
  preds_probs_clinical<-predict(classifier_clinical,newdata=train_clinical,type="prob")[[1]]
  
  roc_clinical<-roc(target_binary,preds_probs_clinical[,"pos"])
  best <- pROC::coords(
    roc_clinical,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  thresh<-0.03
  
  pred_probs<-predict(classifier_clinical,newdata=test_clinical,type="prob")[[1]]
  preds_class<-ifelse(pred_probs[,2]>thresh,1,0)  
  
  regression_vals_low<-as.numeric(predict(low_regression_clinical,newdata=test_clinical)[[1]])
  regression_vals_high<-as.numeric(predict(high_regression_clinical,newdata=test_clinical)[[1]])
  
  
  predicted_val_overall<-round(ifelse(preds_class==0,regression_vals_low,regression_vals_high),2)
  predicted_clinical<-c(predicted_clinical,predicted_val_overall)
  
  
  
  
}


res_count<-data.frame(Target=target_arr,Fold=fold_arr,Predicted.SNP=predicted_snp,
                      Predicted.Clinical=predicted_clinical)
save(res,file='./Sensitivity_Analysis_Imputation/Count_part.Rda')
