rm(list=ls())
outer_fold<-1
######################### Model results ###################
library(caret)
library(PRROC)
library(pROC)


median_clinical<-c()
median_snp<-c()
for(outer_fold in 1:10){
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/LowRisk_regressor/Regression_Clinical_Low_Models.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/LowRisk_regressor/Regression_SNP_Low_Models.Rda'))
  regression_models<-list()
  for(i in names(regression_cl_model_low)){
    regression_models[i]<-regression_cl_model_low[i]
  }
  for(i in names(regression_snp_model_low)){
    regression_models[i]<-regression_snp_model_low[i]
  }
  
  
  ###### load
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/Train_clinical_CV_BA_DM.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/Test_clinical_CV_BA_DM.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/Train_snp_CV_BA_DM.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/Test_snp_CV_BA_DM.Rda'))
  
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
  
  
  test_all<-c()
  folds_arr<-c()
  snp_preds<-c()
  clinical_preds<-c()
  for(jj in 1:10){
    #Clinical
    test_clinical <-test_clinical_ls_processed[[jj]]
    
    test_clinical$target<-test_clinical_ls[[jj]]$SYNTAX.SCORE
    test_clinical<-test_clinical[test_clinical$target<23,]
    model_clinical<-regression_models[[paste0('Clinical_AUC_None_',jj)]]
    
    preds_cl<-predict(model_clinical,test_clinical)
    clinical_preds<-c(clinical_preds,preds_cl)
    
    #SNP
    test_snp <-test_snp_ls_processed[[jj]] 
    test_snp<-test_snp[test_clinical_ls[[jj]]$SYNTAX.SCORE<23,]
    model_snp<-regression_models[[paste0('SNP_AUC_None_',jj)]]
    
    
    
    preds_snp<-predict(model_snp,test_snp)
    
    snp_preds<-c(snp_preds,preds_snp)
    
    
    test_all<-c(test_all,test_clinical$target[test_clinical$target<23])
    folds_arr<-c(folds_arr,rep(jj,length(test_clinical$target)))
  }
  
  
  med_snp<- median(abs(snp_preds-test_all))
  med_clinical<- median(abs(clinical_preds-test_all))
  
  
  median_clinical<-c(median_clinical,med_clinical)
  median_snp<-c(median_snp,med_snp)
}




######### Fold-level performance ###########
folds<-c()
MAE_clinical<-c()
MAE_snp<-c()
for(fold in 1:10){
  ae_clinical<-abs(test_all[which(folds_arr==fold)]-clinical_preds[which(folds_arr==fold)])
  mae_clinical<-median(ae_clinical)
  
  MAE_clinical<-c(MAE_clinical,mae_clinical)
  
  
  ae_snp<-abs(test_all[which(folds_arr==fold)]-snp_preds[which(folds_arr==fold)])
  mae_snp<-median(ae_snp)
  
  MAE_snp<-c(MAE_snp,mae_snp)
  
  folds<-c(folds,fold)
}
res<-data.frame(Fold=folds,MAE.clinical=MAE_clinical,MAE.SNP=MAE_snp)
quantile(res$MAE.clinical) 


########### NEXT #######################




y_true<-test_all
pred_model_A<-clinical_preds
pred_model_B<-snp_preds

n_boot <- 5000
boot_medae_A <- numeric(n_boot)
boot_medae_B <- numeric(n_boot)
boot_delta_medae <- numeric(n_boot)

complete_index <- which(
  complete.cases(y_true, pred_model_A, pred_model_B)
)

for (i in seq_len(n_boot)) {
  
  boot_index <- sample(
    complete_index,
    length(complete_index),
    replace = TRUE
  )
  
  boot_medae_A[i] <- median(
    abs(y_true[boot_index] - pred_model_A[boot_index])
  )
  
  boot_medae_B[i] <- median(
    abs(y_true[boot_index] - pred_model_B[boot_index])
  )
  
  boot_delta_medae[i] <-
    boot_medae_B[i] - boot_medae_A[i]
}

ci_medae_A <- quantile(
  boot_medae_A,
  c(0.025, 0.975),
  na.rm = TRUE
)

ci_medae_B <- quantile(
  boot_medae_B,
  c(0.025, 0.975),
  na.rm = TRUE
)

ci_delta_medae <- quantile(
  boot_delta_medae,
  c(0.025, 0.975),
  na.rm = TRUE
)




