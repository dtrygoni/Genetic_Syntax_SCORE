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
library(MLmetrics)
library(reticulate)
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

###### load yours results


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
  acc <- MLmetrics::Accuracy(pred, obs)
  f1 <- FBeta_Score(obs, pred, positive = "pos", beta = 1)
  f2 <- FBeta_Score(obs, pred, positive = "pos", beta = 2)
  f05 <- FBeta_Score(obs, pred, positive = "pos", beta = 0.5)
  
  # AUC Calculation (assumes the column name for positive class is 'pos')
  if ("pos" %in% colnames(data)) {
    auc <- pROC::roc(response = obs, predictor = data[["pos"]], levels = rev(lev), direction = ">")$auc
  } else {
    auc <- NA
  }
  
  out <- c(
    balAcc = bal_acc,
    F1Score = f1,
    Spec = spec,
    Sens = sens,
    Accuracy = acc,
    F2Score = f2,
    F05Score = f05,
    AUC = auc
  )
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
  hwe<-tableHWE(snps,Binary)
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
  }
  
  for(group in names(snps_to_check)){
    group_snp<-snps_to_check[[group]]
    snp_temp<-snp_ld_df[colnames(snp_ld_df) %in% c('Binary',group_snp)]
    aic_vals<-c()
    for(snp in group_snp){
      r<-SNPassoc::association(as.formula(paste0('Binary ~ ',snp)),data=snp_ld_df)
      aic_vals<-c(aic_vals,r[ind,"AIC"])
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
    form<-as.formula(paste('Binary ~',snp))
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
  train_clinical<-train[colnames(train) %in% c(clinical_predictors,'Binary')]
  train_clinical<-train_clinical[complete.cases(train_clinical),]
  test_clinical<-test[colnames(test) %in% c(clinical_predictors,'Binary')]
  test_clinical<-test_clinical[complete.cases(test_clinical),]
  #### SNP sets
  train_snp<-train[colnames(train) %in% c(clinical_predictors,snp_predictors,'Binary')]
  train_snp<-train_snp[complete.cases(train_snp),]
  test_snp<-test[colnames(test) %in% c(clinical_predictors,snp_predictors,'Binary')]
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
    boruta_res<-Boruta(Binary ~ .,data=snp_boruta,pValue=Boruta_threshold,maxRuns=1000)
    snp_boruta_cols_zero<-colnames(snp_boruta)[boruta_res$finalDecision=='Confirmed']
    
    
    clinical_train<-train[colnames(train) %in% c(setdiff(clinical_cols_zero,'Binary'),'Binary')]
    clinical_test<-test[colnames(test) %in% c(setdiff(clinical_cols_zero,'Binary'),'Binary')]
    
    snp_train<-train[colnames(train) %in% c(setdiff(clinical_cols_zero,'Binary'),setdiff(snp_boruta_cols_zero,'Binary'),'Binary')]
    snp_test<-test[colnames(test) %in% c(setdiff(clinical_cols_zero,'Binary'),setdiff(snp_boruta_cols_zero,'Binary'),'Binary')]
    
    
    clinical_train$target <-as.numeric(ifelse(clinical_train$Binary==0,0,1))
    clinical_test$target <-as.numeric(ifelse(clinical_test$Binary==0,0,1))
    snp_train$target <-as.numeric(ifelse(snp_train$Binary==0,0,1))
    snp_test$target <-as.numeric(ifelse(snp_test$Binary==0,0,1))
    
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('Binary')]
    clinical_test<-clinical_test[!colnames(clinical_test) %in% c('Binary')]
    snp_train<-snp_train[!colnames(snp_train) %in% c('Binary')]
    snp_test<-snp_test[!colnames(snp_test) %in% c('Binary')]
    
  }else if(strategy=='All'){
    snp_boruta<-train[,c(clinical_predictors,snp_predictors,target_col)]
    boruta_res<-Boruta(Binary ~ .,data=snp_boruta,pValue=Boruta_threshold,maxRuns=1000)
    snp_boruta_cols_zero<-colnames(snp_boruta)[boruta_res$finalDecision=='Confirmed']
    
    
    clinical_train<-train[colnames(train) %in% c(setdiff(clinical_cols_zero,'Binary'),'Binary')]
    clinical_test<-test[colnames(test) %in% c(setdiff(clinical_cols_zero,'Binary'),'Binary')]
    
    snp_train<-train[colnames(train) %in% c(setdiff(snp_boruta_cols_zero,'Binary'),'Binary')]
    snp_test<-test[colnames(test) %in% c(setdiff(snp_boruta_cols_zero,'Binary'),'Binary')]
    
    
    clinical_train$target <-as.numeric(ifelse(clinical_train$Binary==0,0,1))
    clinical_test$target <-as.numeric(ifelse(clinical_test$Binary==0,0,1))
    snp_train$target <-as.numeric(ifelse(snp_train$Binary==0,0,1))
    snp_test$target <-as.numeric(ifelse(snp_test$Binary==0,0,1))
    
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('Binary')]
    clinical_test<-clinical_test[!colnames(clinical_test) %in% c('Binary')]
    snp_train<-snp_train[!colnames(snp_train) %in% c('Binary')]
    snp_test<-snp_test[!colnames(snp_test) %in% c('Binary')]
  } 
  
  
  return(list(`Train_clinical`=clinical_train,`Train_SNP`=snp_train,`Test_clinical`=clinical_test,`Test_SNP`=snp_test,
              `Clinical_predictors`=clinical_cols_zero,`SNP_predictors`=setdiff(snp_boruta_cols_zero,clinical_cols_zero)))
}
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
Model_training<-function(clinical_train,snp_train,model_selection='XGB',number_iterations=100,optimization_metric='balAcc'){
  set.seed(123)
  
  input_x_clinical <- as.matrix(dplyr::select(clinical_train, -c("target")))
  input_y_clinical <- clinical_train$target
  input_y_clinical <- factor(input_y_clinical, levels = c(0, 1), labels = c("neg", "pos"))
  
  
  input_x_snp<-as.matrix(dplyr::select(snp_train, -c("target")))
  input_y_snp <- snp_train$target
  input_y_snp <- factor(input_y_snp, levels = c(0, 1), labels = c("neg", "pos"))
  
  
  
  ctrl <- trainControl(
    method = "cv",
    number = 4,
    search = "random",
    summaryFunction = summary_fnc,
    classProbs = TRUE,
    verboseIter = TRUE,
    allowParallel = TRUE
  )
  
  # Control for final model training
  tune_control <- trainControl(
    method = "cv",
    number = 4,
    classProbs = TRUE,
    summaryFunction = summary_fnc,
    verboseIter = FALSE,
    allowParallel = TRUE
  )
  if(model_selection=='XGB'){
    # Fit the model with a random subset of the grid (e.g., 20 random combinations)
    xgb_model_clinical <- train(
      x = input_x_clinical,
      y = as.factor(input_y_clinical),          
      method = "xgbTree",
      trControl = ctrl,
      tuneLength =number_iterations,
      metric=optimization_metric
    )
    
    (final_grid<- expand.grid(
      nrounds = xgb_model_clinical$bestTune$nrounds,
      eta = xgb_model_clinical$bestTune$eta,
      max_depth = xgb_model_clinical$bestTune$max_depth,
      gamma = xgb_model_clinical$bestTune$gamma,
      colsample_bytree = xgb_model_clinical$bestTune$colsample_bytree,
      min_child_weight = xgb_model_clinical$bestTune$min_child_weight,
      subsample = xgb_model_clinical$bestTune$subsample
    ))
    
    (model_Clinical <- caret::train(
      x = input_x_clinical,
      y = as.factor(input_y_clinical),
      trControl = tune_control,
      tuneGrid = final_grid,
      method = "xgbTree",
      verbose = TRUE
    ))
    
    
    xgb_model_snp <- train(
      x = input_x_snp,
      y = as.factor(input_y_snp),          
      method = "xgbTree",
      trControl = ctrl,
      tuneLength =number_iterations,
      metric=optimization_metric
    )
    
    (final_grid_snp <- expand.grid(
      nrounds = xgb_model_snp$bestTune$nrounds,
      eta = xgb_model_snp$bestTune$eta,
      max_depth = xgb_model_snp$bestTune$max_depth,
      gamma = xgb_model_snp$bestTune$gamma,
      colsample_bytree = xgb_model_snp$bestTune$colsample_bytree,
      min_child_weight = xgb_model_snp$bestTune$min_child_weight,
      subsample = xgb_model_snp$bestTune$subsample
    ))
    
    (model_SNP <- caret::train(
      x = input_x_snp,
      y = as.factor(input_y_snp),
      trControl = tune_control,
      tuneGrid = final_grid,
      method = "xgbTree",
      verbose = TRUE
    ))
  }else if(model_selection=='RF'){
    rf_model_clinical <- train(
      x = input_x_clinical,
      y = as.factor(input_y_clinical),          
      method = "rf",
      trControl = ctrl,
      tuneLength =number_iterations,
      metric=optimization_metric
    )
    
    (final_grid<- expand.grid(
      mtry=rf_model_clinical$bestTune$mtry
    ))
    
    (model_Clinical <- caret::train(
      x = input_x_clinical,
      y = as.factor(input_y_clinical),
      trControl = tune_control,
      tuneGrid = final_grid,
      method = "rf",
      verbose = TRUE
    ))
    
    
    rf_model_snp <- train(
      x = input_x_snp,
      y = as.factor(input_y_snp),          
      method = "rf",
      trControl = ctrl,
      tuneLength =number_iterations,
      metric=optimization_metric
    )
    
    (final_grid_snp <- expand.grid(
      mtry=rf_model_snp$bestTune$mtry
    ))
    
    (model_SNP <- caret::train(
      x = input_x_snp,
      y = as.factor(input_y_snp),
      trControl = tune_control,
      tuneGrid = final_grid_snp,
      method = "rf",
      verbose = TRUE
    ))}else if(model_selection=='kNN'){
      knn_model_clinical <- train(
        x = input_x_clinical,
        y = as.factor(input_y_clinical),          
        method = "knn",
        trControl = ctrl,
        tuneLength =number_iterations,
        metric=optimization_metric
      )
      
      (final_grid<- expand.grid(
        k=knn_model_clinical$bestTune$k
      ))
      
      (model_Clinical <- caret::train(
        x = input_x_clinical,
        y = as.factor(input_y_clinical),
        trControl = tune_control,
        tuneGrid = final_grid,
        method = "knn"))
      
      knn_model_snp <- train(
        x = input_x_snp,
        y = as.factor(input_y_snp),          
        method = "knn",
        trControl = ctrl,
        tuneLength =number_iterations,
        metric=optimization_metric
      )
      
      (final_grid_snp <- expand.grid(
        k=knn_model_snp$bestTune$k
      ))
      
      (model_SNP <- caret::train(
        x = input_x_snp,
        y = as.factor(input_y_snp),
        trControl = tune_control,
        tuneGrid = final_grid_snp,
        method = "knn"
      ))
      
      
    }else if(model_selection=='SVM'){
      svm_model_clinical <- train(
        x = input_x_clinical,
        y = as.factor(input_y_clinical),          
        method = "svmRadial",
        trControl = ctrl,
        tuneLength =number_iterations,
        metric=optimazation_metric
      )
      
      (final_grid<- expand.grid(
        sigma=svm_model_clinical$bestTune$sigma,
        C=svm_model_clinical$bestTune$C
      ))
      
      (model_Clinical <- caret::train(
        x = input_x_clinical,
        y = as.factor(input_y_clinical),
        trControl = tune_control,
        tuneGrid = final_grid,
        method = "svmRadial",
        probs=TRUE))
      
      svm_model_snp <- train(
        x = input_x_snp,
        y = as.factor(input_y_snp),          
        method = "svmRadial",
        trControl = ctrl,
        tuneLength =number_iterations,
        metric=optimazation_metric
      )
      
      (final_grid_snp <- expand.grid(
        sigma=svm_model_snp$bestTune$sigma,
        C=svm_model_snp$bestTune$C
      ))
      
      (model_SNP <- caret::train(
        x = input_x_snp,
        y = as.factor(input_y_snp),
        trControl = tune_control,
        tuneGrid = final_grid_snp,
        method = "svmRadial",
        probs=TRUE
      ))
      
    }
  
  
  return(list(`ClinicalModel`=model_Clinical,`SNPmodel`=model_SNP))
  
}
Resample_risk<-function(clinical_train,clinical_test,snp_train,snp_test,sampling='under',resample_proportion=0.5,n_neighbors=5){
  clinical_train$group<-clinical_train$target
  clinical_train$orig_index <- seq_len(nrow(clinical_train))
  
  
  
  if(sampling=='under'){
    min_n <- clinical_train %>%
      count(group) %>%
      summarise(min_n = min(n)) %>%
      pull(min_n)
    
    clinical_train_undersampled <- clinical_train %>%
      group_by(group) %>%
      slice_sample(n = min_n, replace = FALSE) %>%
      ungroup()
    
    resampled_indices<-clinical_train_undersampled$orig_index
    
    clinical_train<-clinical_train_undersampled
    clinical_train$target<-clinical_train$group
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('group','orig_index')]
    snp_train<-snp_train[resampled_indices,]
    snp_train$target<-clinical_train$target
    
    
  }else if(sampling=='over'){
    max_n <- clinical_train %>%
      count(group) %>%
      summarise(max_n = max(n)) %>%
      pull(max_n)
    
    clinical_train_oversampled <- clinical_train %>%
      group_by(group) %>%
      slice_sample(n = max_n, replace = TRUE) %>%
      ungroup()
    # Get resampled data and keep original indices
    #resampled <- clinical_train %>%
    #group_by(group) %>%
    #slice_sample(n = 65) %>%
    #ungroup()
    resampled_indices<-clinical_train_oversampled$orig_index
    
    clinical_train<-clinical_train_oversampled
    clinical_train$target<-clinical_train$target
    
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('group','orig_index')]
    snp_train<-snp_train[resampled_indices,]
    snp_train$target<-clinical_train$target
  }else if(sampling=='None'){
    clinical_train$target<-clinical_train$group
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('group','orig_index')]
    snp_train$target<-clinical_train$target
  }else if(sampling=='under_05'){
    group_counts <- clinical_train %>%
      count(group)
    
    # Get min and max class sizes
    min_n <- min(group_counts$n)
    max_n <- max(group_counts$n)
    
    # Target number per group (between min and max)
    target_n <- round(min_n + resample_proportion * (max_n - min_n))
    
    # Sample target_n from each group (or all if fewer than target_n)
    clinical_train_resampled <- clinical_train %>%
      group_by(group) %>%
      group_modify(~ slice_sample(.x, n = min(nrow(.x), target_n))) %>%
      ungroup()
    
    
    resampled_indices<-clinical_train_resampled$orig_index
    
    clinical_train<-clinical_train_resampled
    clinical_train$target<-clinical_train$group
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('group','orig_index')]
    snp_train<-snp_train[resampled_indices,]
    snp_train$target<-clinical_train$target
  }else if(sampling=='over_05'){
    group_counts <- clinical_train %>%
      count(group)
    
    # Get min and max class sizes
    min_n <- min(group_counts$n)
    max_n <- max(group_counts$n)
    
    # Target number per group (between min and max)
    target_n <- round(max_n - resample_proportion * (max_n - min_n))
    
    # Sample target_n from each group (or all if fewer than target_n)
    clinical_train_resampled <- clinical_train %>%
      group_by(group) %>%
      group_modify(~ {
        n_rows <- nrow(.x)
        if (n_rows < target_n) {
          slice_sample(.x, n = target_n, replace = TRUE)
        } else {
          .x  # keep majority class intact
        }
      }) %>%
      ungroup()
    
    
    resampled_indices<-clinical_train_resampled$orig_index
    
    clinical_train<-clinical_train_resampled
    clinical_train$target<-clinical_train$group
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('group','orig_index')]
    snp_train<-snp_train[resampled_indices,]
    snp_train$target<-clinical_train$target
  }else if(sampling=='ENN'){
    ENN <- import("imblearn.under_sampling", convert = TRUE)$EditedNearestNeighbours
    enn <- ENN(sampling_strategy='majority',n_neighbors = as.integer(n_neighbors))
    X <- clinical_train[, -which(names(clinical_train) %in% c("group"))]
    y <- clinical_train$group
    # Use reticulate's conversion to make X and y usable in Python
    resampled <- enn$fit_resample(X, y)
    df_X<-resampled[[1]]
    y<-resampled[[2]]
    
    resampled_indices<-df_X$orig_index
    
    clinical_train<-df_X
    clinical_train$target_risk<-y
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('group','orig_index')]
    
    snp_train<-snp_train[resampled_indices,]
    snp_train$target_risk<-clinical_train$target_risk
    
    
  }else if(sampling=='RENN'){
    
    snp_train$group<-snp_train$target
    snp_train$orig_index <- seq_len(nrow(snp_train))
    
    ENN <- import("imblearn.under_sampling", convert = TRUE)$RepeatedEditedNearestNeighbours
    enn <- ENN(sampling_strategy='majority',n_neighbors = as.integer(n_neighbors),max_iter = as.integer(100))
    X <- snp_train[, -which(names(snp_train) %in% c("group"))]
    y <- snp_train$group
    # Use reticulate's conversion to make X and y usable in Python
    resampled <- enn$fit_resample(X, y)
    df_X<-resampled[[1]]
    y<-resampled[[2]]
    
    resampled_indices<-df_X$orig_index
    snp_train<-df_X
    snp_train$target<-y
    snp_train<-snp_train[!colnames(snp_train) %in% c('group','orig_index')]  
    
    
    X <- clinical_train[, -which(names(clinical_train) %in% c("group",'orig_index'))]
    y <- clinical_train$group
    resampled <- enn$fit_resample(X, y)
    df_X<-resampled[[1]]
    y<-resampled[[2]]
    clinical_train<-df_X
    clinical_train$target<-y
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('group','orig_index')]  
    
    
    
    #clinical_train<-df_X
    #clinical_train$target_risk<-y
    #clinical_train<-clinical_train[!colnames(clinical_train) %in% c('group','orig_index')]
    
    #snp_train<-snp_train[resampled_indices,]
    #snp_train$target_risk<-clinical_train$target_risk
  }
  # Oversample each group up to max_n
  
  
  return(list('Clinical_Train'=clinical_train,'Clinical_Test'=clinical_test,'SNP_train'=snp_train,'SNP_test'=snp_test))
}
####### New functions ######
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
 ### Parameters
genetic_modelling<-'Codominant'
HWE_threshold<-0.05
MAF_threshold<-0.05
clinical_cols<-1:47
LR_threshold<-0.10
Boruta_threshold<-0.05
strategy_boruta<-'Seperate' ## All second experiments
target_col<-'Binary'
snp_cols<-48:ncol(dataset)
factor_col_opt='One-hot' #One-hot, L-1_encoding
numerical_col_opt='MinMax' #None, MinMax
optimization_metric="balAcc" #balAcc,F1Score,Spec.Specificity, Sens.Sensitivity  Accuracy, F2Score, F05Score,"AUC"
model_selection='RF'
number_iterations=50
metrics<-c('AUC')
samplingss<-c('None','under','over','under_05','over_05','RENN3','RENN5','RENN7','RENN9')
samplingss<-c('None')
delong_p<-list()
res_all<-NULL
for(samp in samplingss){
  print(samp)
  if(samp %in% c('RENN3','RENN5','RENN7','RENN9')){
    res_samp<-strsplit(samp,split='RENN')[[1]]
    sampling_opt<-'RENN'
    ENN_n<-as.numeric(res_samp[[2]])}else{
      sampling_opt<-samp
    }
for(metric in metrics){
  print(metric)
  clinical_probs<-c()
  clinical_target<-c()
  snp_probs<-c()
  snp_target<-c()
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
  train_clinical_list<-list()
  test_clinical_list<-list()
  train_snp_list<-list()
  test_snp_list<-c()
  optimization_metric<-metric
for(jj in 1:10){
  print(jj)
  clinical_train<-train_clinical_ls[[jj]]
  clinical_test<-test_clinical_ls[[jj]]
  snp_train<-train_snp_ls[[jj]]
  snp_test<-test_snp_ls[[jj]]
  
  #Include interactions#
  res_interaction<-interaction_analysis(snp_train,snp_test)
  snp_train<-res_interaction$Train
  snp_test<-res_interaction$Test
  
  ########## Pre-processing for factor/numerical variables #########
  res_preprocessing<-Preprocessing_fnc(clinical_train,clinical_test,snp_train,snp_test,factor_col_opt=factor_col_opt,numerical_col_opt=numerical_col_opt)
  clinical_train<-res_preprocessing$ClinicalTrain
  clinical_test<-res_preprocessing$ClinicalTest
  snp_train<-res_preprocessing$SNPtrain
  snp_test<-res_preprocessing$SNPtest
  
  res_resample<-Resample_risk(clinical_train,clinical_test,snp_train,snp_test,sampling = sampling_opt,n_neighbors = ENN_n)
  clinical_train<-res_resample$Clinical_Train
  clinical_test<-res_resample$Clinical_Test
  snp_train<-res_resample$SNP_train
  snp_test<-res_resample$SNP_test
 
  
  train_clinical_list[[jj]]<-clinical_train
  test_clinical_list[[jj]]<-clinical_test
  train_snp_list[[jj]]<-snp_train
  test_snp_list[[jj]]<-snp_test
  
  ########## Model evaluation -- train/hyperparameter tuning ################
  res_model<-Model_training(clinical_train,snp_train,model_selection=model_selection,number_iterations=number_iterations,optimization_metric = optimization_metric)
  clinical_model<-res_model$ClinicalModel
  snp_model<-res_model$SNPmodel
  #### Clinical evaluation ####
  test_x_clinical<- as.matrix(dplyr::select(clinical_test, -c("target")))
  test_y_clinical <- clinical_test$target
  test_y_clinical <- factor(test_y_clinical, levels = c(0, 1), labels = c("neg", "pos"))
  
  xtab <- table(predict(clinical_model,newdata = test_x_clinical),test_y_clinical)
  Measures.Predicted.rf <- confusionMatrix(xtab, positive="pos")
  Accuracy.Predicted.rf <- enframe(round(Measures.Predicted.rf$overall,4),name="Measure", value="Value")
  Measures.Predicted.rf <- enframe(Measures.Predicted.rf$byClass, name="Measure", value="Value")
  Measures.Predicted.rf.df <- rbind(Accuracy.Predicted.rf,Measures.Predicted.rf)
  Measures.Predicted.rf.df <- Measures.Predicted.rf.df[c(1,18,12,13,14),]
  
  fold_arr<-c(fold_arr,jj)
  model_arr<-c(model_arr,model_selection)
  type_arr<-c(type_arr,'Clinical')
  acc_arr<-c(acc_arr,round(Measures.Predicted.rf.df$Value[1],3))
  prec_arr<-c(prec_arr,round(Measures.Predicted.rf.df$Value[3],3))
  rec_arr<-c(rec_arr,round(Measures.Predicted.rf.df$Value[4],3))
  f1_arr<-c(f1_arr,round(Measures.Predicted.rf.df$Value[5],3))
  
  
  res_curves<-precrec::evalmod(mmdata(predict(clinical_model,newdata=test_x_clinical,type="prob")[,2],as.factor(test_y_clinical)))
  aucs_df<-precrec::auc(res_curves)
  roc_auc<-aucs_df$aucs[1]
  pr_auc<-aucs_df$aucs[2]
  
  aucpr_arr<-c(aucpr_arr,round(pr_auc,3))
  aucroc_arr<-c(aucroc_arr,round(roc_auc,3))
  
  
  #Save model
  string<-paste0('Clinical_',model_selection,'_',jj)
  final_models[[string]]<-clinical_model
  clinical_probs<-c(clinical_probs,predict(clinical_model,newdata=test_x_clinical,type="prob")[,2])
  clinical_target<-c(clinical_target,as.factor(test_y_clinical))
  
  ###### SNP model #####
  test_x_snp<- as.matrix(dplyr::select(snp_test, -c("target")))
  test_y_snp <- snp_test$target
  test_y_snp <- factor(test_y_snp, levels = c(0, 1), labels = c("neg", "pos"))
  xtab <- table(predict(snp_model,newdata = test_x_snp),test_y_snp)
  
  Measures.Predicted.rf <- confusionMatrix(xtab, positive="pos")
  Accuracy.Predicted.rf <- enframe(round(Measures.Predicted.rf$overall,4),name="Measure", value="Value")
  Measures.Predicted.rf <- enframe(Measures.Predicted.rf$byClass, name="Measure", value="Value")
  Measures.Predicted.rf.df <- rbind(Accuracy.Predicted.rf,Measures.Predicted.rf)
  Measures.Predicted.rf.df <- Measures.Predicted.rf.df[c(1,18,12,13,14),]
  
  fold_arr<-c(fold_arr,jj)
  model_arr<-c(model_arr,model_selection)
  type_arr<-c(type_arr,'ClinicalSNP')
  acc_arr<-c(acc_arr,round(Measures.Predicted.rf.df$Value[1],3))
  prec_arr<-c(prec_arr,round(Measures.Predicted.rf.df$Value[3],3))
  rec_arr<-c(rec_arr,round(Measures.Predicted.rf.df$Value[4],3))
  f1_arr<-c(f1_arr,round(Measures.Predicted.rf.df$Value[5],3))
  
  
  res_curves<-precrec::evalmod(mmdata(predict(snp_model,newdata=test_x_snp,type="prob")[,2],as.factor(test_y_snp)))
  aucs_df<-precrec::auc(res_curves)
  roc_auc<-aucs_df$aucs[1]
  pr_auc<-aucs_df$aucs[2]
  
  aucpr_arr<-c(aucpr_arr,round(pr_auc,3))
  aucroc_arr<-c(aucroc_arr,round(roc_auc,3))
  
  
  #Save model
  string<-paste0('SNP_',model_selection,'_',jj)
  final_models[[string]]<-snp_model
  snp_probs<-c(snp_probs,predict(snp_model,newdata=test_x_snp,type="prob")[,2])
  snp_target<-c(snp_target,as.factor(test_y_snp))
  
}
  
  
  roc1 <- roc(clinical_target,clinical_probs)
  roc2 <- roc(snp_target,snp_probs)
  ######## Delong
  res<-roc.test(roc1,roc2,method="delong")
  delong_p[[paste0(metric,'_',samp)]]<-res
  
  
  sampling_arr<-rep(samp,length(model_arr))
  res<-data.frame(Model=model_arr,Data_Type=type_arr,Fold=fold_arr,Sampling=sampling_arr,Accuracy=acc_arr,Precision=prec_arr,
                  Recall=rec_arr,F1=f1_arr,AUCROC=aucroc_arr,AUCPR=aucpr_arr)
  res_all<-rbind(res_all,res)
  
  
  str_models<-paste0('./Extra/Models_RF_',metric,'_',samp,'.Rda')
  save(final_models,file=str_models)
  
  str_train_clinical<-paste0('./Extra/Train_clinical_',metric,'_',samp,'.Rda')
  save(train_clinical_list,file=str_train_clinical)
  
  str_test_clinical<-paste0('./Extra/Test_clinical_',metric,'_',samp,'.Rda')
  save(test_clinical_list,file=str_test_clinical)
  
  str_train_snp<-paste0('./Extra/Train_SNP_',metric,'_',samp,'.Rda')
  save(train_snp_list,file=str_train_snp)
  
  str_test_snp<-paste0('./Extra/Test_SNP_',metric,'_',samp,'.Rda')
  save(test_snp_list,file=str_test_snp)
}
}



AUC_clinical<-c()
AUC_SNP<-c()
statistic<-c()
CI_low<-c()
CI_high<-c()
p_val<-c()
for(name in names(delong_p)){
  r<-delong_p[[name]]
  AUC_clinical<-c(AUC_clinical,as.numeric(r$estimate[1]))
  AUC_SNP<-c(AUC_SNP,as.numeric(r$estimate[2]))
  statistic<-c(statistic,r$statistic)
  CI_low<-c(CI_low,r$conf.int[1])
  CI_high<-c(CI_high,r$conf.int[2])
  p_val<-c(p_val,r$p.value)
}

df<-data.frame(Model=names(delong_p),AUC.Clinical=AUC_clinical,AUC.SNP=AUC_SNP,Statistic=statistic,CI.Low=CI_low,CI.High=CI_high,P.value=p_val)
save(df,file='./Results/Delong_results_RF.Rda')
save(res_all,file='./Results/Classifier_results_RF.Rda')






