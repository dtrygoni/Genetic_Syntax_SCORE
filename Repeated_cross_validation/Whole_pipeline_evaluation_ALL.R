rm(list=ls())
################################## Evaluation of complete SYNTAX prediction pipeline across repeated cross validationr results ###############

library(caret)
library(PRROC)
library(pROC)
library(SNPassoc)


median_clinical<-c()
median_snp<-c()

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
optimal_encoding_categorical<-function(train,snp_cols,target_col){
  snp_cols_wanted<-which(colnames(train) %in% snp_cols)
  trainSNP<-setupSNP(data=train,colSNPs=snp_cols_wanted,sep="")
  
  levels_modelling<-list()
  snps<-c()
  genetic_modelling<-"Optimal"
  
  
  ind_categorical<-8
  ind_numerical<-7
  
  wanted<-ind_categorical
  
  
  
  for(snp in snp_cols_wanted){
    form<-as.formula(paste(target_col,' ~',colnames(train)[snp]))
    res<-association(form,data=trainSNP)
    if(length(res[,1])==15){
      snps<-c(snps,snp)
      levels_co<-names(res[,1])[2:4]
      levels_do<-names(res[,1])[6:7]
      levels_re<-names(res[,1])[9:10]
      levels_over<-names(res[,1])[12:13]
      
      if(genetic_modelling=='Codominant'){
        levels_modelling[[colnames(train)[snp]]]<-levels_co
      }else if(genetic_modelling=='Dominant'){
        levels_modelling[[colnames(train)[snp]]]<-levels_do
      }else if(genetic_modelling=='Recessive'){
        levels_modelling[[colnames(train)[snp]]]<-levels_re
      }else if(genetic_modelling=='Optimal'){
        index<-which.min(res[1:12,wanted])
        if(index==3){
          levels_opt<-levels_co
        }else if(index==6){
          levels_opt<-levels_do
        }else if(index==9){
          levels_opt<-levels_re
        }else if(index==12){
          levels_opt<-levels_over
        }
        levels_modelling[[colnames(train)[snp]]]<-levels_opt
      }
      
      
    }
  }
  
  return(levels_modelling)
}
optimal_encoding_numeric<-function(train,snp_cols,target_col){
  snp_cols_wanted<-which(colnames(train) %in% snp_cols)
  trainSNP<-setupSNP(data=train,colSNPs=snp_cols_wanted,sep="")
  
  levels_modelling<-list()
  snps<-c()
  genetic_modelling<-"Optimal"
  
  
  ind_categorical<-8
  ind_numerical<-7
  
  wanted<-ind_numerical
  
  
  
  for(snp in snp_cols_wanted){
    form<-as.formula(paste(target_col,' ~',colnames(train)[snp]))
    res<-association(form,data=trainSNP)
    if(length(res[,1])==15){
      snps<-c(snps,snp)
      levels_co<-names(res[,1])[2:4]
      levels_do<-names(res[,1])[6:7]
      levels_re<-names(res[,1])[9:10]
      levels_over<-names(res[,1])[12:13]
      
      if(genetic_modelling=='Codominant'){
        levels_modelling[[colnames(train)[snp]]]<-levels_co
      }else if(genetic_modelling=='Dominant'){
        levels_modelling[[colnames(train)[snp]]]<-levels_do
      }else if(genetic_modelling=='Recessive'){
        levels_modelling[[colnames(train)[snp]]]<-levels_re
      }else if(genetic_modelling=='Optimal'){
        index<-which.min(res[1:12,wanted])
        if(index==3){
          levels_opt<-levels_co
        }else if(index==6){
          levels_opt<-levels_do
        }else if(index==9){
          levels_opt<-levels_re
        }else if(index==12){
          levels_opt<-levels_over
        }
        levels_modelling[[colnames(train)[snp]]]<-levels_opt
      }
      
      
    }
  }
  
  return(levels_modelling)
}
snp_numeric_convert <-function(dataset,snp_columns,levels){
  snp_cols_wanted<-which(colnames(dataset) %in% snp_columns)
  trainSNP<-setupSNP(data=dataset,colSNPs=snp_cols_wanted,sep="")
  
  dataset_no_snp<-trainSNP[! colnames(trainSNP) %in% snp_columns]
  dataset_snp<-trainSNP[colnames(trainSNP) %in% snp_columns]
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

for(outer_fold in 1:10){
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Zero_part_classification/Train_clinical_CV_BA_DM.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Zero_part_classification/Test_clinical_CV_BA_DM.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Zero_part_classification/Train_snp_CV_BA_DM.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Zero_part_classification/Test_snp_CV_BA_DM.Rda'))
  
  
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
  
  train_clinical_ls_processed_zero<-train_clinical_ls_processed
  test_clinical_ls_processed_zero<-test_clinical_ls_processed
  train_snp_ls_processed_zero<-train_snp_ls_processed
  test_snp_ls_processed_zero<-test_snp_ls_processed
  
  ################### Find threshold for zero-part ##################
  thresholds_clinical_zero<-c()
  thresholds_snp_zero<-c()
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Zero_part_classification/Models/Models_RF_AUC_None.Rda'))
  for(jj in 1:10){
    #Clinical
    test_clinical <-train_clinical_ls_processed_zero[[jj]]
    
    model_clinical<-final_models[[paste0('Clinical_RF_',jj)]]
    
    train_clinical<-train_clinical_ls_processed_zero[[jj]]
    probs<-predict(model_clinical,train_clinical,type="prob")[,2]
    roc_curve<-pROC::roc(train_clinical$target,probs)
    best_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity", "ppv", "npv","youden"), best.method = "youden")
    threshold<-as.numeric(best_coords$threshold)
    thresholds_clinical_zero<-c(thresholds_clinical_zero,threshold)
    
    
    #SNP
    test_snp <-test_snp_ls_processed_zero[[jj]] 
    model_snp<-final_models[[paste0('SNP_RF_',jj)]]
    
    train_snp<-train_snp_ls_processed_zero[[jj]]
    probs<-predict(model_snp,train_snp,type="prob")[,2]
    roc_curve<-pROC::roc(train_snp$target,probs)
    best_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity", "ppv", "npv","youden"), best.method = "youden")
    threshold<-as.numeric(best_coords$threshold)
    thresholds_snp_zero<-c(thresholds_snp_zero,threshold)
    
    
  }
  
  
  ################## Find thresholds for count-part classifier #############
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/Risk_classifier/Classifier_Models.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/Train_clinical_CV_BA_DM.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/Test_clinical_CV_BA_DM.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/Train_snp_CV_BA_DM.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/Test_snp_CV_BA_DM.Rda'))
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
  thresholds_clinical_count<-c()
  thresholds_snp_count<-c()
  for(jj in 1:10){
    #Clinical
    test_clinical <-test_clinical_ls_processed[[jj]]
    
    test_clinical$target<-ifelse(test_clinical_ls[[jj]]$SYNTAX.SCORE<23,1,0)
    model_clinical<-classifier_models[[paste0('Clinical_AUC_None_',jj)]]
    
    train_clinical<-train_clinical_ls_processed[[jj]]
    probs<-predict(model_clinical,train_clinical,type="prob")[,2]
    
    train_clinical$target<-ifelse(train_clinical_ls[[jj]]$SYNTAX.SCORE<23,1,0)
    
    roc_curve<-pROC::roc(train_clinical$target,probs)
    best_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity", "ppv", "npv","youden"), best.method = "youden")
    threshold<-as.numeric(best_coords$threshold)
    thresholds_clinical_count<-c(thresholds_clinical_count,threshold)
    
    #SNP
    test_snp <-test_snp_ls_processed[[jj]] 
    model_snp<-classifier_models[[paste0('SNP_AUC_None_',jj)]]
    
    train_snp<-train_snp_ls_processed[[jj]]
    probs<-predict(model_snp,train_snp,type="prob")[,2]
    
    train_snp$target<-train_clinical$target
    
    roc_curve<-pROC::roc(train_snp$target,probs)
    
    best_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity", "ppv", "npv","youden"), best.method = "youden")
    
    threshold<-as.numeric(best_coords$threshold)
    thresholds_snp_count<-c(thresholds_snp_count,threshold)
  }
  
  
  
  
  
  
  
  
  
  
  #################### Prepare test sets ##############################
  ### 1. Need the construction of the same set twice (one for zero predictors and one for count predictors) ###
  load("C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Whole_dataset.Rda")
  load(paste0("C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold",outer_fold,"/Test_inds.Rda"))
  
  
  ### 2. Find the participants with non-NA values in the columns needed for both zero and count part
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Zero_part_classification/Train_clinical_CV_BA_DM.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Zero_part_classification/Test_clinical_CV_BA_DM.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Zero_part_classification/Train_snp_CV_BA_DM.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Zero_part_classification/Test_snp_CV_BA_DM.Rda'))
  train_clinical_ls_zero<-train_clinical_ls
  test_clinical_ls_zero<-test_clinical_ls
  train_snp_ls_zero<-train_snp_ls
  test_snp_ls_zero<-test_snp_ls
  
  
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/Risk_classifier/Classifier_Models.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/Train_clinical_CV_BA_DM.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/Test_clinical_CV_BA_DM.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/Train_snp_CV_BA_DM.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/Test_snp_CV_BA_DM.Rda'))
  
  

  test_sets_clinical_count<-list()
  test_sets_snp_count<-list()
  test_sets_clinical_zero<-list()
  test_sets_snp_zero<-list()
  
  
  train_sets_clinical_count<-list()
  train_sets_clinical_zero<-list()
  train_sets_snp_count<-list()
  train_sets_snp_zero<-list()
  for(jj in 1:10){
    all_predictors<-unique(c(colnames(train_snp_ls_zero[[jj]]),colnames(train_snp_ls[[jj]]),colnames(train_clinical_ls_zero[[jj]]),
                                                                                                     colnames(train_clinical_ls[[jj]])))
    all_predictors<-all_predictors[!all_predictors %in% c('target')]
    train_all_participants<-dataset[-test_indices[[jj]],all_predictors]
    test_all_participants<-dataset[test_indices[[jj]],all_predictors]
    
    train_all_noNA<-train_all_participants[complete.cases(train_all_participants),]
    test_all_noNA<-test_all_participants[complete.cases(test_all_participants),]
    test_all_noNA$target<-ifelse(test_all_noNA$SYNTAX.SCORE>0,1,0)
    train_all_noNA$target<-ifelse(train_all_noNA$SYNTAX.SCORE>0,1,0)
    
    
    test_zero_snp<-test_all_noNA[colnames(train_snp_ls_zero[[jj]])]
    test_zero_clinical<-test_all_noNA[colnames(train_clinical_ls_zero[[jj]])]
    test_sets_clinical_zero[[jj]]<-test_zero_clinical
    test_sets_snp_zero[[jj]]<-test_zero_snp
    
    test_count_snp<-test_all_noNA[colnames(train_snp_ls[[jj]])]
    test_count_clinical<-test_all_noNA[colnames(train_clinical_ls[[jj]])]
    test_sets_clinical_count[[jj]]<-test_count_clinical
    test_sets_snp_count[[jj]]<-test_count_snp
    
    train_zero_snp<-train_all_noNA[colnames(train_snp_ls_zero[[jj]])]
    train_zero_clinical<-train_all_noNA[colnames(train_clinical_ls_zero[[jj]])]
    train_sets_clinical_zero[[jj]]<-train_zero_clinical
    train_sets_snp_zero[[jj]]<-train_zero_snp
    
    train_count_snp<-train_all_noNA[colnames(train_snp_ls[[jj]])]
    train_count_clinical<-train_all_noNA[colnames(train_clinical_ls[[jj]])]
    train_sets_clinical_count[[jj]]<-train_count_clinical
    train_sets_snp_count[[jj]]<-train_count_snp
  }
  
  
  ######## zero sets pre-processing (Train set based) ######
  #### 1. Find SNP levels and transform train/test #####
  for(i in 1:10){
    train_set<-train_sets_snp_zero[[i]]
    test_set<-test_sets_snp_zero[[i]]
    
    rs_wanted<-colnames(train_set)[(grepl("rs",colnames(train_set)))]
    res<-optimal_encoding_categorical(train_set,rs_wanted,target_col='target')
    
    res_convert_train<-snp_numeric_convert(train_set,rs_wanted,res)
    res_convert_test<-snp_numeric_convert(test_set,rs_wanted,res)
    
    train_sets_snp_zero[[i]]<-res_convert_train
    test_sets_snp_zero[[i]]<-res_convert_test
  }
  
  #### 2. Min-Max and One-hot enconding
  factor_col_opt='One-hot' #One-hot, L-1_encoding
  numerical_col_opt='MinMax' #None, MinMax
  train_clinical_ls_processed<-list() #for dataset_clinical_train
  test_clinical_ls_processed<-list() #for dataset_clinical_test
  for(i in 1:10){
    res_preprocessing<-Preprocessing_fnc(train_sets_clinical_zero[[i]],test_sets_clinical_zero[[i]],
                                         train_sets_snp_zero[[i]],test_sets_snp_zero[[i]],factor_col_opt=factor_col_opt,numerical_col_opt=numerical_col_opt)
    train_set<-res_preprocessing$ClinicalTrain
    test_set<-res_preprocessing$ClinicalTest
    
    train_clinical_ls_processed[[i]]<-train_set
    test_clinical_ls_processed[[i]]<-test_set
  }
  
  
  ### Prepare snps sets
  train_snp_ls_processed<-list()
  test_snp_ls_processed<-list()
  for(i in 1:10){
    res_interaction<-interaction_analysis(train_sets_snp_zero[[i]],test_sets_snp_zero[[i]])
    snp_train<-res_interaction$Train
    snp_test<-res_interaction$Test
    res_preprocessing<-Preprocessing_fnc(train_sets_clinical_zero[[i]],test_sets_clinical_zero[[i]],
                                         snp_train,snp_test,factor_col_opt=factor_col_opt,numerical_col_opt=numerical_col_opt)
    
    train_snp_ls_processed[[i]]<-res_preprocessing$SNPtrain
    test_snp_ls_processed[[i]]<-res_preprocessing$SNPtest
  }
  
  train_snp_ls_processed_zero<-train_snp_ls_processed
  test_snp_ls_processed_zero<-test_snp_ls_processed
  train_clinical_ls_processed_zero<-train_clinical_ls_processed
  test_clinical_ls_processed_zero<-test_clinical_ls_processed



  ######## Count sets pre-processing (Train set based) #####
  #### 1. Find SNP levels and transform train/test #####
  for(i in 1:10){
      train_set<-train_sets_snp_count[[i]]
      test_set<-test_sets_snp_count[[i]]
      
      rs_wanted<-colnames(train_set)[(grepl("rs",colnames(train_set)))]
      res<-optimal_encoding_numeric(train_set,rs_wanted,target_col='SYNTAX.SCORE')
      
      res_convert_train<-snp_numeric_convert(train_set,rs_wanted,res)
      res_convert_test<-snp_numeric_convert(test_set,rs_wanted,res)
      
      train_sets_snp_count[[i]]<-res_convert_train
      test_sets_snp_count[[i]]<-res_convert_test
    }
    
  #### 2. Min-Max and One-hot enconding
  factor_col_opt='One-hot' #One-hot, L-1_encoding
  numerical_col_opt='MinMax' #None, MinMax
  train_clinical_ls_processed<-list() #for dataset_clinical_train
  test_clinical_ls_processed<-list() #for dataset_clinical_test
  for(i in 1:10){
      res_preprocessing<-Preprocessing_fnc(train_sets_clinical_count[[i]],test_sets_clinical_count[[i]],
                                           train_sets_snp_count[[i]],test_sets_snp_count[[i]],factor_col_opt=factor_col_opt,numerical_col_opt=numerical_col_opt)
      train_set<-res_preprocessing$ClinicalTrain
      test_set<-res_preprocessing$ClinicalTest
      
      train_clinical_ls_processed[[i]]<-train_set
      test_clinical_ls_processed[[i]]<-test_set
    }
    
    
    
    
    
   ### Prepare snps sets
   train_snp_ls_processed<-list()
   test_snp_ls_processed<-list()
   for(i in 1:10){
      res_interaction<-interaction_analysis(train_sets_snp_count[[i]],test_sets_snp_count[[i]])
      snp_train<-res_interaction$Train
      snp_test<-res_interaction$Test
      res_preprocessing<-Preprocessing_fnc(train_sets_snp_zero[[i]],test_sets_snp_zero[[i]],
                                           snp_train,snp_test,factor_col_opt=factor_col_opt,numerical_col_opt=numerical_col_opt)
      
      train_snp_ls_processed[[i]]<-res_preprocessing$SNPtrain
      test_snp_ls_processed[[i]]<-res_preprocessing$SNPtest
    }
    

    
  
    
    
    
  ############## Prediction module FINALLY ################
  #1. Zero-part classifier
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Zero_part_classification/Models/Models_RF_AUC_None.Rda'))
  #2. Count-part classifier 
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/Risk_classifier/Classifier_Models.Rda'))
  #3. Regression models 
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/LowRisk_regressor/Regression_Clinical_Low_Models.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/LowRisk_regressor/Regression_SNP_Low_Models.Rda'))
  regression_models<-list()
  for(i in names(regression_cl_model_low)){
    regression_models[i]<-regression_cl_model_low[i]
  }
  for(i in names(regression_snp_model_low)){
    regression_models[i]<-regression_snp_model_low[i]
  }
  
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/HighRisk_regressor/Regression_Clinical_High_Models.Rda'))
  load(paste0('C:/Users/dtryg/Desktop/GESS_experiments_code/Repeated_cross_validation/Fold',outer_fold,'/Count_part_regression/HighRisk_regressor/Regression_SNP_High_Models.Rda'))
  high_regression_models<-list()
  for(i in names(regression_cl_model_high)){
    high_regression_models[i]<-regression_cl_model_high[i]
  }
  for(i in names(regression_snp_model_low)){
    high_regression_models[i]<-regression_snp_model_high[i]
  }
  
  
  clinical_preds<-c()
  test_all<-c()
  snp_preds<-c()
  for(jj in 1:10){
    ######## Clinical prediction ###############################################
    #1. Zero classification
    zero_classifier<-final_models[[paste0("Clinical_RF_",jj)]]
    
    #2. define train/test for zero and count
    zero_clinical_train<-train_clinical_ls_processed_zero[[jj]]
    zero_clinical_test<-test_clinical_ls_processed_zero[[jj]]
    count_clinical_train<-train_clinical_ls_processed[[jj]]
    count_clinical_test<-test_clinical_ls_processed[[jj]]
    
    #3. Predict zero-part class
    predicted_class_prob<-predict(zero_classifier,zero_clinical_test,type="prob")[,2]
    predicted_label_zero<-ifelse(predicted_class_prob>thresholds_clinical_zero[jj],"pos","neg")
    
    #4.Predict count-part risk class
    count_part_classifier<-classifier_models[[paste0("Clinical_AUC_None_",jj)]]
    optimal_threshold<-thresholds_clinical_count[jj]
    
    predicted_class_prob<-predict(count_part_classifier,count_clinical_test,type="prob")[,2]
    predicted_label<-ifelse(predicted_class_prob>optimal_threshold,"pos","neg")
    
    #5. Predict regression scores
    low_regression_clinical<-regression_models[[paste0('Clinical_AUC_None_',jj)]]
    high_regression_clinical<-high_regression_models[[paste0('Clinical_AUC_None_',jj)]]
    
    predicted_value_count<-ifelse(predicted_label=="pos",predict(high_regression_clinical,count_clinical_test),predict(low_regression_clinical,count_clinical_test))
    
    
    
    #6. Final prediction ###
    clinical_predicted<-ifelse(predicted_label_zero=="pos",predicted_value_count,0)
    clinical_preds<-c(clinical_preds,clinical_predicted)
    
    
    
    
    ######################### SNP prediction #######################
    #1. Zero classification
    zero_classifier<-final_models[[paste0("SNP_RF_",jj)]]
    
    #2. define train/test for zero and count
    zero_snp_train<-train_snp_ls_processed_zero[[jj]]
    zero_snp_test<-test_snp_ls_processed_zero[[jj]]
    count_snp_train<-train_snp_ls_processed[[jj]]
    count_snp_test<-test_snp_ls_processed[[jj]]
    
    #3. Predict zero-part class
    predicted_class_prob<-predict(zero_classifier,zero_snp_test,type="prob")[,2]
    predicted_label_zero<-ifelse(predicted_class_prob>thresholds_snp_zero[jj],"pos","neg")
    
    #4.Predict count-part risk class
    count_part_classifier<-classifier_models[[paste0("SNP_AUC_None_",jj)]]
    optimal_threshold<-thresholds_snp_count[jj]
    
    predicted_class_prob<-predict(count_part_classifier,count_snp_test,type="prob")[,2]
    predicted_label<-ifelse(predicted_class_prob>optimal_threshold,"pos","neg")
    
    #5. Predict regression scores
    low_regression_snp<-regression_models[[paste0('SNP_AUC_None_',jj)]]
    high_regression_snp<-high_regression_models[[paste0('SNP_AUC_None_',jj)]]
    
    predicted_value_count<-ifelse(predicted_label=="pos",predict(high_regression_snp,count_snp_test),predict(low_regression_snp,count_snp_test))
    
    
    
    #6. Final prediction ###
    snp_predicted<-ifelse(predicted_label_zero=="pos",predicted_value_count,0)
    snp_preds<-c(snp_preds,snp_predicted)
    
    
    
    test_all<-c(test_all,count_clinical_test$SYNTAX.SCORE)
  }
  
  
  median_clinical<-c(median_clinical,median(abs(test_all-clinical_preds)))
  median_snp<-c(median_snp,median(abs(test_all-snp_preds)))
}