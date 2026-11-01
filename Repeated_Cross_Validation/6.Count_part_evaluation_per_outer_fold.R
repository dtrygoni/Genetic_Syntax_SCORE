outer_fold<-1 #Change the outer fold to get the values you want
######################### Model results ###################
library(caret)
library(PRROC)
library(pROC)


load(paste0('./Repeated_cross_validation/Fold',outer_fold,'/Positive_score_component/Risk_classifier/Classifier_Models.Rda'))
load(paste0('./Repeated_cross_validation/Fold',outer_fold,'/Positive_score_component/LowRisk_regressor/Regression_Clinical_Low_Models.Rda'))
load(paste0('./Repeated_cross_validation/Fold',outer_fold,'/Positive_score_component/LowRisk_regressor/Regression_SNP_Low_Models.Rda'))
load(paste0('./Repeated_cross_validation/Fold',outer_fold,'/Positive_score_component/HighRisk_regressor/Regression_Clinical_High_Models.Rda'))
load(paste0('./Repeated_cross_validation/Fold',outer_fold,'/Positive_score_component/HighRisk_regressor/Regression_SNP_High_Models.Rda'))



count_classifiers<-classifier_models
count_regressor_low_cl<-regression_cl_model_low
count_regressor_high_cl<-regression_cl_model_high
count_regressor_low_snp<-regression_snp_model_low
count_regressor_high_snp<-regression_snp_model_high


###### load
load(paste0('./Repeated_cross_validation/Fold',outer_fold,'/Positive_score_component/Train_clinical_count.Rda'))
load(paste0('./Repeated_cross_validation/Fold',outer_fold,'/Positive_score_component/Test_clinical_count.Rda'))
load(paste0('./Repeated_cross_validation/Fold',outer_fold,'/Positive_score_component/Train_snp_count.Rda'))
load(paste0('./Repeated_cross_validation/Fold',outer_fold,'/Positive_score_component/Test_snp_count.Rda'))

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



Count_part_prediction<-function(test,classifier,regressor_low,regressor_high,thresh){
  preds_classifiers<-predict(classifier,newdata=test,type="prob")[,2]
  preds_classifiers<-ifelse(preds_classifiers>=thresh,"pos","neg")
  preds_regressor_low<-predict(regressor_low,newdata=test)
  preds_regressor_high<-predict(regressor_high,newdata = test)
  
  preds_f<-round(ifelse(preds_classifiers=='pos',preds_regressor_high,preds_regressor_low),2)
  return(preds_f)
}
Optimal_threshold <- function(
    train, target, classifier, type = "count", seed = 123
) {
  # Accept either a caret model or your one-element model list.
  if (!inherits(classifier, "train")) {
    if (is.list(classifier) &&
        length(classifier) == 1L &&
        inherits(classifier[[1]], "train")) {
      classifier <- classifier[[1]]
    } else {
      stop("classifier must contain a caret::train model.")
    }
  }
  
  if (classifier$method != "rf") {
    stop("This function expects a caret Random Forest classifier.")
  }
  
  # Use exactly the predictors used to fit the original classifier.
  # This excludes outcome columns present in 'train'.
  predictor_names <- classifier$coefnames
  
  if (is.null(predictor_names) ||
      !all(predictor_names %in% names(train))) {
    stop("train must contain all original processed predictor columns.")
  }
  
  x <- as.data.frame(train[, predictor_names, drop = FALSE])
  
  if (nrow(x) != length(target)) {
    stop("train and target must have the same number and order of rows.")
  }
  
  target <- as.character(target)
  
  if (anyNA(target) || !all(target %in% c("0", "1"))) {
    stop("target must contain only 0 and 1, without missing values.")
  }
  
  y <- factor(target, levels = c("0", "1"),
              labels = c("neg", "pos"))
  
  if (anyNA(x)) {
    stop("train must contain complete eligible outer-training rows.")
  }
  
  if (any(table(y) < 4L)) {
    stop("At least four observations per class are needed.")
  }
  
  set.seed(seed)
  
  # Each participant is held out once.
  folds <- caret::createFolds(y, k = 4, returnTrain = FALSE)
  cv_probabilities <- rep(NA_real_, nrow(x))
  
  for (i in seq_along(folds)) {
    validation_rows <- folds[[i]]
    training_rows <- setdiff(seq_len(nrow(x)), validation_rows)
    
    # Refit using the selected mtry and original RF settings.
    # No further hyperparameter search is performed.
    fold_model <- randomForest::randomForest(
      x = x[training_rows, , drop = FALSE],
      y = y[training_rows],
      mtry = classifier$bestTune$mtry,
      ntree = classifier$finalModel$ntree
    )
    
    cv_probabilities[validation_rows] <- predict(
      fold_model,
      newdata = x[validation_rows, , drop = FALSE],
      type = "prob"
    )[, "pos"]
  }
  
  if (any(!is.finite(cv_probabilities))) {
    stop("Some inner-CV probabilities are missing or non-finite.")
  }
  
  roc_obj <- pROC::roc(
    response = y,
    predictor = cv_probabilities,
    levels = c("neg", "pos"),
    direction = "<",
    quiet = TRUE
  )
  
  best <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = "threshold",
    transpose = FALSE
  )
  
  thresholds <- best$threshold
  thresholds <- thresholds[is.finite(thresholds)]
  
  if (!length(thresholds)) {
    stop("No finite optimal threshold was found.")
  }
  
  # Explicit tie rule: larger threshold favors specificity.
  return(max(thresholds))
}
opt_arr<-c("AUC")
resamp_arr<-c("None")
predictions_cl_all<-list()
target_cl_all<-list()
probs_binary_cl_all<-list()

for(opt in opt_arr){
  for(resamp in resamp_arr){
    predictions_clinical<-c()
    target_clinical<-c()
    probs_binary<-c()
    ae_arr<-c()
    fold_arr<-c()
    for(jj in 1:10){
      #snp_zero<-final_test_snp_zero_ls[[jj]]
      snp_count<-test_snp_ls_processed[[jj]]
      
      #complete_snp<-complete.cases(snp_zero)
      complete_snp_count<-complete.cases(snp_count)
      
      #complete_snp_all<-complete_snp & complete_snp_count
      
      
      #clinical_zero<-final_test_clinical_zero_ls[[jj]]
      clinical_count<-test_clinical_ls_processed[[jj]]
      
      #complete_cl<-complete.cases(clinical_zero)
      complete_count<-complete.cases(clinical_count)
      
      #complete_cl<-complete_cl & complete_count
      
      
      complete<-complete_snp_count & complete_count
      
      #clinical_zero<-clinical_zero[complete,]
      clinical_count<-clinical_count[complete,]
      
      
      #zero_jj<-zero_models[paste0('Clinical_RF_',jj)]
      #preds<-predict(zero_jj,clinical_zero)[[1]]
      #preds_probs<-predict(zero_jj,clinical_zero,type='prob')[[1]]
      #probs_binary<-c(probs_binary,preds_probs[,2])
      
      #Clinical
      name<-paste0('Clinical_',opt,'_',resamp,'_',jj)
      classifier<-count_classifiers[[name]]
      thresh<-Optimal_threshold(train_clinical_ls_processed[[jj]],train_clinical_ls[[jj]]$Binary.SS,classifier)
      
      name_regr<-paste0('Clinical_AUC_',resamp,'_',jj)
      regressor_low<-count_regressor_low_cl[[name_regr]]
      regressor_high<-count_regressor_high_cl[[name_regr]]
      preds_count<-Count_part_prediction(clinical_count,classifier,regressor_low,regressor_high,thresh)
      
      
      predictions<-preds_count#ifelse(preds=='neg',0,preds_count)
      predictions_clinical<-c(predictions_clinical,predictions)
      
      
      target<-test_clinical_ls[[jj]]$SYNTAX.SCORE
      target_clinical<-c(target_clinical,target)
      
      ae<-mean(abs(target-predictions))
      ae_arr<-c(ae_arr,ae)
      fold_arr<-c(fold_arr,jj)
    }
    name<-paste0('Clinical_',opt,'_',resamp)
    predictions_cl_all[[name]]<-predictions_clinical
    target_cl_all[[name]]<-target_clinical
    probs_binary_cl_all[[name]]<-probs_binary
  }
}



round(mean(abs(predictions_clinical-target_clinical)),2)
round(median(abs(predictions_clinical-target_clinical)),2)











# Evaluate zero if positive then count model
opt_arr<-c("AUC")
resamp_arr<-c("None")
predictions_snp_all<-list()
target_snp_all<-list()
probs_binary_snp_all<-list()
for(opt in opt_arr){
  for(resamp in resamp_arr){
    predictions_snp<-c()
    target_snp<-c()
    probs_binary<-c()
    
    fold_arr<-c()
    ae_snp<-c()
    for(jj in 1:10){
      #snp_zero<-final_test_snp_zero_ls[[jj]]
      snp_count<-test_snp_ls_processed[[jj]]
      
      #complete_snp<-complete.cases(snp_zero)
      complete_snp_count<-complete.cases(snp_count)
      
      #complete_snp_all<-complete_snp & complete_snp_count
      
      
      #clinical_zero<-final_test_clinical_zero_ls[[jj]]
      #clinical_count<-final_test_clinical_count_ls[[jj]]
      
      #complete_cl<-complete.cases(clinical_zero)
      #complete_count<-complete.cases(clinical_count)
      
      #complete_cl<-complete_cl & complete_count
      
      
      #complete<-complete_cl & complete_snp_all
      
      
      
      #snp_zero<-snp_zero[complete,]
      snp_count<-snp_count[complete_snp_count,]
      
      
      #zero_jj<-zero_models[paste0('SNP_RF_',jj)]
      #preds<-predict(zero_jj,snp_zero)[[1]]
      
      #preds_probs<-predict(zero_jj,snp_zero,type='prob')[[1]]
      #probs_binary<-c(probs_binary,preds_probs[,2])
      #snp
      
      
      
      
      name<-paste0('SNP_',opt,'_',resamp,'_',jj)
      classifier<-count_classifiers[[name]]
      thresh<-Optimal_threshold(train_snp_ls_processed[[jj]],train_snp_ls[[jj]]$Binary.SS,classifier)
      
      name_regr<-paste0('SNP_AUC_',resamp,'_',jj)
      regressor_low<-count_regressor_low_snp[[name_regr]]
      regressor_high<-count_regressor_high_snp[[name_regr]]
      preds_count<-Count_part_prediction(snp_count,classifier,regressor_low,regressor_high,thresh)
      
      
      predictions<-preds_count#ifelse(preds=='neg',0,preds_count)
      predictions_snp<-c(predictions_snp,predictions)
      
      
      target<-test_snp_ls[[jj]]$SYNTAX.SCORE
      target_snp<-c(target_snp,target)
      
      ae<-median(abs(target-predictions))
      
      ae_snp<-c(ae_snp,ae)
    }
    name<-paste0('SNP_',opt,'_',resamp)
    predictions_snp_all[[name]]<-predictions_snp
    target_snp_all[[name]]<-target_snp
    #probs_binary_snp_all[[name]]<-probs_binary
  }
}
round(mean(abs(predictions_snp-target_snp)),2)
round(median(abs(predictions_snp-target_snp)),2)


