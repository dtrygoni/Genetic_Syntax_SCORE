############################################## Model training -- Count Part Random Forest Classifier / Regressors ########################
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
set.seed(123)
xgb.set.config(verbosity = 0)

model_arr<-c()
type_arr<-c()
r2_arr<-c()
mae_arr<-c()
mre_arr<-c()
rmse_arr<-c()
fold_arr<-c()
final_models<-list()
train_list<-list()
test_list<-list()
###### load
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Train_clinical_CV_BA_DM.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Test_clinical_CV_BA_DM.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Train_snp_CV_BA_DM.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Test_snp_CV_BA_DM.Rda')

path_resampling<-'C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Prev/Cross_validation_experiments_Count_Classifier/Resampling_script.py'
library(reticulate)
reticulate::source_python(path_resampling)
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
    bal.Acc = bal_acc,
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
Preprocessing_fnc<-function(clinical_train,clinical_test,snp_train,snp_test,target_col='target',factor_col_opt='One-hot',numerical_col_opt='None'){
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
Model_training_count<-function(clinical_train,snp_train,model_selection='RF',number_iterations=100,optimization_metric='MAE',target_col=c('SYNTAX.SCORE','Binary.SS')){
  set.seed(123)
  input_x_clinical <- as.matrix(dplyr::select(clinical_train, -target_col))
  input_y_clinical <- clinical_train$SYNTAX.SCORE
  #input_y_clinical <- factor(input_y_clinical, levels = c(0, 1), labels = c("neg", "pos"))
  
  
  input_x_snp<-as.matrix(dplyr::select(snp_train, -target_col))
  input_y_snp <- snp_train$SYNTAX.SCORE
  #input_y_snp <- factor(input_y_snp, levels = c(0, 1), labels = c("neg", "pos"))
  
  
  
  ctrl <- trainControl(
    method = "cv",
    number = 4,
    search = "random",
    #summaryFunction = summary_fnc,
    #classProbs = TRUE,
    verboseIter = TRUE,
    allowParallel = TRUE
  )
  
  # Control for final model training
  tune_control <- trainControl(
    method = "cv",
    number = 4,
    #classProbs = TRUE,
    #summaryFunction = summary_fnc,
    verboseIter = FALSE,
    allowParallel = TRUE
  )
  
  
  
  if(model_selection=='XGB'){
    # Fit the model with a random subset of the grid (e.g., 20 random combinations)
    xgb_model_clinical <- train(
      x = input_x_clinical,
      y = input_y_clinical,          
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
      y = input_y_clinical,
      trControl = tune_control,
      tuneGrid = final_grid,
      method = "xgbTree",
      verbose = TRUE
    ))
    
    
    xgb_model_snp <- train(
      x = input_x_snp,
      y = input_y_snp,          
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
      y = input_y_snp,
      trControl = tune_control,
      tuneGrid = final_grid,
      method = "xgbTree",
      verbose = TRUE
    ))
  }else if(model_selection=='RF'){
    rf_model_clinical <- train(
      x = input_x_clinical,
      y = input_y_clinical,          
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
      y =input_y_clinical,
      trControl = tune_control,
      tuneGrid = final_grid,
      method = "rf",
      verbose = TRUE
    ))
    
    
    rf_model_snp <- train(
      x = input_x_snp,
      y = input_y_snp,          
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
      y = input_y_snp,
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
Model_training_classifier<-function(clinical_train,snp_train,model_selection='RF',number_iterations=100,optimization_metric='bal.Acc',target_col=c('SYNTAX.SCORE','Binary.SS')){
  set.seed(123)
  input_x_clinical <- as.matrix(dplyr::select(clinical_train, -target_col))
  input_y_clinical <- clinical_train$Binary.SS
  input_y_clinical <- factor(input_y_clinical, levels = c(0, 1), labels = c("neg", "pos"))
  
  
  input_x_snp<-as.matrix(dplyr::select(snp_train, -target_col))
  input_y_snp <- snp_train$Binary.SS
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
  
  
  
  if(model_selection=='RF'){
    rf_model_clinical <- train(
      x = input_x_clinical,
      y = input_y_clinical,          
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
      y =input_y_clinical,
      trControl = tune_control,
      tuneGrid = final_grid,
      method = "rf",
      verbose = TRUE
    ))
    
    
    rf_model_snp <- train(
      x = input_x_snp,
      y = input_y_snp,          
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
      y = input_y_snp,
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
Preprocess_targetcol<-function(clinical_train,clinical_test,snp_train,snp_test,target_col='target',preprocessor_opt='MinMax'){
  
  if(preprocessor_opt=='MinMax'){
    trsf_train<-clinical_train
    trsf_test<-clinical_test
    trsf_trainsnp<-snp_train
    trsf_testsnp<-snp_test
    
    pre_proc <- preProcess(trsf_train[colnames(trsf_train) %in% c(target_col)], method = c("range"))
    # Transform the dataset
    clinical_train <- predict(pre_proc, newdata = trsf_train)
    clinical_test <- predict(pre_proc, newdata = trsf_test)
    
    scalerClinical_obj<-pre_proc
    
    pre_procsnp <- preProcess(trsf_trainsnp[colnames(trsf_trainsnp) %in% c(target_col)], method = c("range"))
    # Transform the dataset
    snp_train <- predict(pre_procsnp, newdata = trsf_trainsnp)
    snp_test <- predict(pre_procsnp, newdata = trsf_testsnp)
    
    scalerSNP_obj<-pre_procsnp
  }else if(preprocessor_opt=='log2'){
    trsf_train<-clinical_train
    trsf_test<-clinical_test
    trsf_trainsnp<-snp_train
    trsf_testsnp<-snp_test
    
    clinical_train$target<-log(clinical_train$target,base=2)
    clinical_test$target<-log(clinical_test$target,base=2)
    snp_train$target<-log(snp_train$target,base=2)
    snp_test$target<-log(snp_test$target,base=2)
    scalerClinical_obj<-'Log2'
  }else if(preprocessor_opt=='log'){
    trsf_train<-clinical_train
    trsf_test<-clinical_test
    trsf_trainsnp<-snp_train
    trsf_testsnp<-snp_test
    
    clinical_train$target<-log(clinical_train$target)
    clinical_test$target<-log(clinical_test$target)
    snp_train$target<-log(snp_train$target)
    snp_test$target<-log(snp_test$target)
    scalerClinical_obj<-'Log'
  }else if(preprocessor_opt=='BoxCox'){
    trsf_train<-clinical_train
    trsf_test<-clinical_test
    trsf_trainsnp<-snp_train
    trsf_testsnp<-snp_test
    pre_proc <- preProcess(trsf_train[colnames(trsf_train) %in% c(target_col)], method = c("BoxCox"))
    clinical_train <- predict(pre_proc, newdata = trsf_train)
    clinical_test <- predict(pre_proc, newdata = trsf_test)
    
    
    scalerClinical_obj<-pre_proc
    
    pre_procsnp <- preProcess(trsf_trainsnp[colnames(trsf_trainsnp) %in% c(target_col)], method = c("BoxCox"))
    # Transform the dataset
    snp_train <- predict(pre_procsnp, newdata = trsf_trainsnp)
    snp_test <- predict(pre_procsnp, newdata = trsf_testsnp)
  }else if(preprocessor_opt=='None'){
 
    
    scalerClinical_obj<-NA
    
  }
  return(list('Clinical_Train'=clinical_train,'Clinical_Test'=clinical_test,'SNP_train'=snp_train,'SNP_test'=snp_test,'Scaler'=scalerClinical_obj))
}
Resample_risk<-function(clinical_train,clinical_test,snp_train,snp_test,sampling='under',resample_proportion=0.5,n_neighbors=5){
  clinical_train$group<-clinical_train$Binary.SS
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
    clinical_train<-clinical_train[resampled_indices,]
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('group','orig_index')]
    snp_train<-snp_train[resampled_indices,]
    

    
  }else if(sampling=='over'){
    max_n <- clinical_train %>%
      count(group) %>%
      summarise(max_n = max(n)) %>%
      pull(max_n)
    
    clinical_train_oversampled <- clinical_train %>%
      group_by(group) %>%
      slice_sample(n = max_n, replace = TRUE) %>%
      ungroup()
    
    resampled_indices<-clinical_train_oversampled$orig_index
    clinical_train<-clinical_train[resampled_indices,]
    snp_train<-snp_train[resampled_indices,]
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('group','orig_index')]
   
  }else if(sampling=='None'){
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('group','orig_index')]
    
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
    clinical_train$target_risk<-clinical_train$group
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('group','orig_index')]
    snp_train<-snp_train[resampled_indices,]
    snp_train$target_risk<-clinical_train$target_risk
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
    clinical_train$target_risk<-clinical_train$group
    clinical_train<-clinical_train[!colnames(clinical_train) %in% c('group','orig_index')]
    snp_train<-snp_train[resampled_indices,]
    snp_train$target_risk<-clinical_train$target_risk
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
    
    snp_train$group<-snp_train$Binary.SS
    snp_train$orig_index <- seq_len(nrow(snp_train))
    
    ENN <- import("imblearn.under_sampling", convert = TRUE)$RepeatedEditedNearestNeighbours
    enn <- ENN(sampling_strategy='majority',n_neighbors = as.integer(n_neighbors),max_iter = as.integer(10))
    X <- snp_train[, -which(names(snp_train) %in% c("Binary.SS","SYNTAX.SCORE","group","orig_index"))]
    y <- snp_train$group
    # Use reticulate's conversion to make X and y usable in Python
    resampled <- enn$fit_resample(X, y)
    df_X<-resampled[[1]]
    y<-resampled[[2]]
    X_resampled_df <- as.data.frame(df_X)
    # Find logical vector indicating which rows in X are present in X_resampled
    rows_kept_logical <- apply(X, 1, function(row) {
      any(apply(X_resampled_df, 1, function(res_row) all(row == res_row)))
    })
    snp_train<-snp_train[rows_kept_logical,]
    snp_train<-snp_train[!colnames(snp_train) %in% c('group','orig_index')]  
    
    
    X <- clinical_train[, -which(names(clinical_train) %in% c("Binary.SS","SYNTAX.SCORE","group","orig_index"))]
    y <- clinical_train$group
    resampled <- enn$fit_resample(X, y)
    df_X<-resampled[[1]]
    y<-resampled[[2]]
    X_resampled_df <- as.data.frame(df_X)
    # Find logical vector indicating which rows in X are present in X_resampled
    rows_kept_logical <- apply(X, 1, function(row) {
      any(apply(X_resampled_df, 1, function(res_row) all(row == res_row)))
    })
    
    clinical_train<-clinical_train[rows_kept_logical,]
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

### Parameters
target_col<-'SYNTAX.SCORE'
factor_col_opt='One-hot'#One-hot, L-1_encoding
numerical_col_opt='MinMax' #None, MinMax
optimization_metric="balAcc" #balAcc,F1Score,Spec.Specificity, Sens.Sensitivity  Accuracy, F2Score, F05Score,"AUC"
model_selection='RF'
number_iterations=50
metrics<-c('MAE','RMSE','RSquared') ####### Other metrics
preprocessor_target<-'None'#MinMax, log2,log,BoxCox
wilcoxon_test_ls<-list()
preproc_options<-c('None','log')
metrics<-c('AUC')
############## Classifiers model train ###############

classifier_models<-list()
train_list<-list()
test_list<-list()
resampling<-c('None')#,'under','over','RENN3','RENN5','RENN7','RENN9')


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
for(resamp in resampling){
  if(resamp %in% c('RENN3','RENN5','RENN7','RENN9','RENN11')){
    res_samp<-strsplit(resamp,split='RENN')[[1]]
    sampling_opt<-'RENN'
    ENN_n<-as.numeric(res_samp[[2]])}else{
      sampling_opt<-resamp
      ENN_n<-0
    }
  
  
  for(metric in metrics){
  
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
    res_preprocessing<-Preprocessing_fnc(clinical_train,clinical_test,snp_train,snp_test,target_col=c('SYNTAX.SCORE','Binary.SS'),factor_col_opt=factor_col_opt,numerical_col_opt=numerical_col_opt)
    clinical_train<-res_preprocessing$ClinicalTrain
    clinical_test<-res_preprocessing$ClinicalTest
    snp_train<-res_preprocessing$SNPtrain
    snp_test<-res_preprocessing$SNPtest
    
    
    ########## Scale target variable #####################
    #res_scale<-Preprocess_targetcol(clinical_train,clinical_test,snp_train,snp_test,target_col='target',preprocessor_opt = preprocessor_target)
    #scaler<-res_scale$Scaler
    #clinical_train<-res_scale$Clinical_Train
    #clinical_test<-res_scale$Clinical_Test
    #snp_train<-res_scale$SNP_train
    #snp_test<-res_scale$SNP_test
    ############# Resampling train set ##################
    res_resample<-Resample_risk(clinical_train,clinical_test,snp_train,snp_test,sampling = sampling_opt,n_neighbors = ENN_n)
    clinical_train<-res_resample$Clinical_Train
    clinical_test<-res_resample$Clinical_Test
    snp_train<-res_resample$SNP_train
    snp_test<-res_resample$SNP_test
    
    
    
    ########## Model evaluation -- train/hyperparameter tuning ################
    res_model<-Model_training_classifier(clinical_train,snp_train,model_selection=model_selection,number_iterations=number_iterations,optimization_metric = optimization_metric)
    clinical_model<-res_model$ClinicalModel
    snp_model<-res_model$SNPmodel
    
    classifier_models[[paste0('Clinical_',metric,'_',resamp,'_',jj)]]<-clinical_model
    classifier_models[[paste0('SNP_',metric,'_',resamp,'_',jj)]]<-snp_model
    
    train_list[[paste0('Clinical_',metric,'_',resamp,'_',jj)]]<-clinical_train
    train_list[[paste0('SNP_',metric,'_',resamp,'_',jj)]]<-snp_train
    
    test_list[[paste0('Clinical_',metric,'_',resamp,'_',jj)]]<-clinical_test
    test_list[[paste0('SNP_',metric,'_',resamp,'_',jj)]]<-snp_test
    }

  }
}



save(classifier_models,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Classifier_Models.Rda')
save(train_list,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Train_sets.Rda')
save(test_list,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Test_sets.Rda')


clinical_names<-c()
snp_names<-c()
for(resamp in resampling){
  for(jj in 1:10){
    string_wanted<-paste0('Clinical_','AUC_',resamp,'_',jj)
    clinical_names<-c(clinical_names,string_wanted)
    
    string_wanted<-paste0('SNP_','AUC_',resamp,'_',jj)
    snp_names<-c(snp_names,string_wanted)
  }
}

train_clinical_list<-train_list[names(train_list) %in% clinical_names]
train_snp_list<-train_list[names(train_list) %in% snp_names]
test_clinical_list<-test_list[names(test_list) %in% clinical_names]
test_snp_list<-test_list[names(test_list) %in% snp_names]


regression_cl_model_low<-list()
regression_snp_model_low<-list()
regression_cl_model_high<-list()
regression_snp_model_high<-list()

for(i in 1:length(clinical_names)){
  clinical_train<-train_clinical_list[[clinical_names[i]]]
  clinical_test<-test_clinical_list[[clinical_names[i]]]
  
  snp_train<-train_snp_list[[snp_names[i]]]
  
  snp_test<-test_snp_list[[snp_names[i]]]
  
  
  ###### Train - test sets for Regression models == 0
  clinical_train_0<-clinical_train[clinical_train$Binary.SS==0,]
  clinical_test_0 <- clinical_test[clinical_test$Binary.SS==0,]
  snp_train_0<-snp_train[snp_train$Binary.SS==0,]
  snp_test_0<-snp_test[snp_test$Binary.SS==0,]
  
  
  ##### Train -test sets for Regression models == 1
  clinical_train_1<-clinical_train[clinical_train$Binary.SS==1,]
  clinical_test_1<-clinical_test[clinical_test$Binary.SS==1,]
  snp_train_1<-snp_train[snp_train$Binary.SS==1,]
  snp_test_1<-snp_test[snp_test$Binary.SS==1,]
  
  
  
  ######## Model train 0 #####
  res_model<-Model_training_count(clinical_train_0,snp_train_0,number_iterations = number_iterations,optimization_metric="MAE")
  clinical_model_0<-res_model$ClinicalModel
  snp_model_0<-res_model$SNPmodel
  
  regression_cl_model_low[[clinical_names[i]]]<-clinical_model_0
  regression_snp_model_low[[snp_names[[i]]]]<-snp_model_0
  
  ######## Model train 1#####
  res_model<-Model_training_count(clinical_train_1,snp_train_1,number_iterations = number_iterations,optimization_metric="RMSE")
  clinical_model_1<-res_model$ClinicalModel
  snp_model_1<-res_model$SNPmodel
  
  regression_cl_model_high[[clinical_names[i]]]<-clinical_model_1
  regression_snp_model_high[[snp_names[[i]]]]<-snp_model_1
  
  }

save(regression_cl_model_low,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Regression_Clinical_Low_Models.Rda')
save(regression_snp_model_low,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Regression_SNP_Low_Models.Rda')
save(regression_cl_model_high,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Regression_Clinical_High_Models.Rda')
save(regression_snp_model_high,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Regression_SNP_High_Models.Rda')
