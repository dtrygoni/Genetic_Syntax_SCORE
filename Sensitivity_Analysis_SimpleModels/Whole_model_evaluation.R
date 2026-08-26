rm(list=ls())
########################## Construct the test set for the complete evaluation ####################################
library(recipes)
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
### 1. Need the construction of the same set twice (one for zero predictors and one for count predictors) ###
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Whole_dataset.Rda")
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Test_inds.Rda")


#### Zero-part set is the same as the whole dataset test indices (doesn't need SNP levels here) ###
####################### Preparation of test-sets for the zero-component #################################
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Zero_part_classification/Test_snp_zero.Rda")
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Zero_part_classification/Train_snp_zero.Rda")
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Zero_part_classification/Test_clinical_zero.Rda")
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Zero_part_classification/Train_clinical_zero.Rda")
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

train_clinical_ls_processed_zero<-train_clinical_ls_processed
test_clinical_ls_processed_zero<-test_clinical_ls_processed
train_snp_ls_processed_zero<-train_snp_ls_processed
test_snp_ls_processed_zero<-test_snp_ls_processed



############################# Preparation of test-sets for the count part component ##########################
#1. Load the train/test to get the predictors
load('C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Count_part_regression/Test_clinical_count.Rda')
load('C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Count_part_regression/Train_clinical_count.Rda')
load('C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Count_part_regression/Test_snp_count.Rda')
load('C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Count_part_regression/Train_snp_count.Rda')


dataset_clinical_train<-list()
dataset_clinical_test<-list()
dataset_snp_train<-list()
dataset_snp_test<-list()

#### Imputation ###
median_mode_impute <- function(train,test = NULL,exclude = NULL,unique_cutoff = 5) {
  
  train_imputed <- train
  test_imputed <- test
  
  # Columns eligible for imputation
  columns <- setdiff(
    names(train),
    exclude
  )
  
  imputation_info <- data.frame(
    Variable = columns,
    Method = NA_character_,
    Unique_observed = NA_integer_,
    Missing_train = NA_integer_,
    Missing_test = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  # Store imputation values without coercing their types
  imputation_values <- vector(
    mode = "list",
    length = length(columns)
  )
  
  names(imputation_values) <- columns
  
  calculate_mode <- function(x) {
    
    observed <- x[!is.na(x)]
    
    if (length(observed) == 0) {
      return(NA)
    }
    
    frequencies <- table(observed)
    
    mode_name <- names(frequencies)[
      which.max(frequencies)
    ]
    
    # Return the mode using the original variable type
    if (is.factor(x)) {
      return(mode_name)
    }
    
    if (is.logical(x)) {
      return(as.logical(mode_name))
    }
    
    if (is.integer(x)) {
      return(as.integer(mode_name))
    }
    
    if (is.numeric(x)) {
      return(as.numeric(mode_name))
    }
    
    mode_name
  }
  
  for (j in seq_along(columns)) {
    
    variable <- columns[j]
    x_train <- train[[variable]]
    
    observed_values <- x_train[!is.na(x_train)]
    n_unique <- length(unique(observed_values))
    
    imputation_info$Unique_observed[j] <- n_unique
    imputation_info$Missing_train[j] <- sum(
      is.na(x_train)
    )
    
    if (!is.null(test) && variable %in% names(test)) {
      imputation_info$Missing_test[j] <- sum(
        is.na(test[[variable]])
      )
    }
    
    # Cannot estimate an imputation value if all training values are missing
    if (length(observed_values) == 0) {
      
      warning(
        paste(
          "Variable", variable,
          "contains only missing values in the training data",
          "and was not imputed."
        )
      )
      
      imputation_info$Method[j] <- "Not imputed"
      imputation_values[[variable]] <- NA
      
      next
    }
    
    # Mode for variables with <= unique_cutoff values
    # and for all non-numeric variables
    if (
      n_unique <= unique_cutoff ||
      !is.numeric(x_train)
    ) {
      
      method <- "Mode"
      fill_value <- calculate_mode(x_train)
      
    } else {
      
      method <- "Median"
      fill_value <- median(
        x_train,
        na.rm = TRUE
      )
    }
    
    imputation_info$Method[j] <- method
    imputation_values[[variable]] <- fill_value
    
    # Impute training data
    missing_train <- is.na(
      train_imputed[[variable]]
    )
    
    if (is.factor(train_imputed[[variable]])) {
      
      fill_value <- as.character(fill_value)
      
      if (
        !fill_value %in%
        levels(train_imputed[[variable]])
      ) {
        levels(train_imputed[[variable]]) <- c(
          levels(train_imputed[[variable]]),
          fill_value
        )
      }
    }
    
    train_imputed[[variable]][missing_train] <-
      fill_value
    
    # Apply the training-derived value to the test data
    if (
      !is.null(test_imputed) &&
      variable %in% names(test_imputed)
    ) {
      
      missing_test <- is.na(
        test_imputed[[variable]]
      )
      
      if (is.factor(test_imputed[[variable]])) {
        
        fill_value_test <- as.character(fill_value)
        
        if (
          !fill_value_test %in%
          levels(test_imputed[[variable]])
        ) {
          levels(test_imputed[[variable]]) <- c(
            levels(test_imputed[[variable]]),
            fill_value_test
          )
        }
        
      } else {
        
        fill_value_test <- fill_value
      }
      
      test_imputed[[variable]][missing_test] <-
        fill_value_test
    }
  }
  
  list(
    train = train_imputed,
    test = test_imputed,
    imputation_values = imputation_values,
    imputation_summary = imputation_info
  )
}
for(i in 1:10){
  train_clinical_count<-train_clinical_ls[[i]]
  dataset_train_i<-dataset[-list_indexes_test[[i]],colnames(train_clinical_count)]
  dataset_test_i<-dataset[list_indexes_test[[i]],colnames(train_clinical_count)]
  res<-median_mode_impute(dataset_train_i,dataset_test_i)
  dataset_clinical_train[[i]]<-res$train
  dataset_clinical_test[[i]]<-res$test
  
  
  train_snp_count<-train_snp_ls[[i]]
  dataset_train_i<-dataset[-list_indexes_test[[i]],colnames(train_snp_count)]
  dataset_test_i<-dataset[list_indexes_test[[i]],colnames(train_snp_count)]
  res<-median_mode_impute(dataset_train_i,dataset_test_i)
  dataset_snp_train[[i]]<-res$train
  dataset_snp_test[[i]]<-res$test
  
}


### Find the levels of the count sets and apply them to train and test sets of the dataset_snp_train/test ###
library(SNPassoc)
optimal_encoding<-function(train,snp_cols,target_col){
  snp_cols_wanted<-which(colnames(train) %in% snp_cols)
  trainSNP<-setupSNP(data=train,colSNPs=snp_cols_wanted,sep="")
  
  levels_modelling<-list()
  snps<-c()
  genetic_modelling<-"Optimal"
  
  
  ind_categorical<-9
  ind_numerical<-7
  
  if(is.numeric(train[[target_col]])){
    wanted<-ind_numerical
  }else{
    wanted<-ind_categorical
  }
  
  
  
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
        if(index==2){
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
for(i in 1:10){
  train_set<-dataset_snp_train[[i]]
  train_set<-train_set[!is.na(train_set$SYNTAX.SCORE),]
  train_set<-train_set[train_set$SYNTAX.SCORE>0,]
  
  rs_wanted<-colnames(train_set)[(grepl("rs",colnames(train_set)))]
  
  res<-optimal_encoding(train_set,rs_wanted,target_col='SYNTAX.SCORE')
  
  res_convert_train<-snp_numeric_convert(dataset_snp_train[[i]],rs_wanted,res)
  res_convert_test<-snp_numeric_convert(dataset_snp_test[[i]],rs_wanted,res)
  
  dataset_snp_train[[i]]<-res_convert_train
  dataset_snp_test[[i]]<-res_convert_test
}



### Prepare the rest ##
factor_col_opt='One-hot' #One-hot, L-1_encoding
numerical_col_opt='MinMax' #None, MinMax
train_clinical_ls_processed<-list() #for dataset_clinical_train
test_clinical_ls_processed<-list() #for dataset_clinical_test
for(i in 1:10){
  res_preprocessing<-Preprocessing_fnc(dataset_clinical_train[[i]],dataset_clinical_test[[i]],dataset_snp_train[[i]],dataset_snp_test[[i]],factor_col_opt=factor_col_opt,numerical_col_opt=numerical_col_opt)
  train_set<-res_preprocessing$ClinicalTrain
  test_set<-res_preprocessing$ClinicalTest
  
  train_clinical_ls_processed[[i]]<-train_set
  test_clinical_ls_processed[[i]]<-test_set
}


### Prepare snps sets
train_snp_ls_processed<-list()
test_snp_ls_processed<-list()
for(i in 1:10){
  res_interaction<-interaction_analysis(dataset_snp_train[[i]],dataset_snp_test[[i]])
  snp_train<-res_interaction$Train
  snp_test<-res_interaction$Test
  res_preprocessing<-Preprocessing_fnc(dataset_clinical_train[[i]],dataset_clinical_test[[i]],snp_train,snp_test,factor_col_opt=factor_col_opt,numerical_col_opt=numerical_col_opt)
  
  train_snp_ls_processed[[i]]<-res_preprocessing$SNPtrain
  test_snp_ls_processed[[i]]<-res_preprocessing$SNPtest
}






####### train_snp_ls_processed_zero -> Zero part / train_snp_ls_processed -> count part 


#### Prediction module: 1. Classification based on zero component
load("C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Zero_part_classification/Models/Models_LR_AUC_None.Rda")
#### 2. Classification based on count classifier
load('C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Count_part_regression/LowRisk_regressor/Regression_Clinical_Low_Models.Rda')
load('C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Count_part_regression/LowRisk_regressor/Regression_SNP_Low_Models.Rda')
regression_models<-list()
for(i in names(regression_cl_model_low)){
  regression_models[i]<-regression_cl_model_low[i]
}
for(i in names(regression_snp_model_low)){
  regression_models[i]<-regression_snp_model_low[i]
}
load('C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Count_part_regression/HighRisk_regressor/Regression_Clinical_High_Models.Rda')
load('C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Count_part_regression/HighRisk_regressor/Regression_SNP_High_Models.Rda')
high_regression_models<-list()
for(i in names(regression_cl_model_high)){
  high_regression_models[i]<-regression_cl_model_high[i]
}
for(i in names(regression_snp_model_low)){
  high_regression_models[i]<-regression_snp_model_high[i]
}
load('C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Count_part_regression/Risk_classifier/Classifier_Models.Rda')

#### 3. Regression
test_all<-c()
folds_arr<-c()
snp_preds<-c()
clinical_preds<-c()
load('C:/Users/dtryg/Desktop/GESS_experiments_code/Sensitivity_Analysis_SimpleModels/Zero_part_classification/Test_clinical_zero.Rda')
inds<-list()
for(i in 1:10){
  library(dplyr)
  full_with_id <- dataset %>% mutate(full_row = row_number()) 
  matched <- test_clinical_ls[[i]] %>% left_join( full_with_id, by = intersect(names( test_clinical_ls[[i]]), names(dataset)) )
  inds[[i]]<-matched$full_row
}
count_thresholds<-c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
for(jj in 1:10){
  # For clinical!
  dd_init<-dataset[-list_indexes_test[[jj]],]
  no_nan_inds<-which(!is.na(dd_init$SYNTAX.SCORE))
  #1. Zero classification
  zero_classifier<-final_models[[paste0("Clinical_RF_",jj)]]

  #2. define train/test for zero and count
  zero_clinical_train<-train_clinical_ls_processed_zero[[jj]]
  zero_clinical_test<-test_clinical_ls_processed_zero[[jj]]
  count_clinical_train<-train_clinical_ls_processed[[jj]]
  count_clinical_test<-test_clinical_ls_processed[[jj]]
  no_nan_inds<-no_nan_inds[1:nrow(zero_clinical_train)]
  
  #3. Find the optimal threshold for the zero-part classifier
  ## Construct ROC curve
  dataset_clinical_train_wanted<-dataset_clinical_train[[jj]][no_nan_inds,]
  roc_obj <- roc(
    response  = ifelse(dataset_clinical_train_wanted$SYNTAX.SCORE>0,1,0),
    predictor = predict(zero_classifier,zero_clinical_train,type="prob")[,2],
    levels    = c(0, 1),  # negative class first, positive class second
    direction = "<"
  )
  
  ## Threshold maximizing Youden's index:
  ## sensitivity + specificity - 1
  optimal <- coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity", "youden"),
    transpose = FALSE
  )
  optimal_threshold <- as.numeric(optimal$threshold)
  
  predicted_class_prob<-predict(zero_classifier,zero_clinical_test,type="prob")[,2]
  predicted_label_zero<-ifelse(predicted_class_prob>optimal_threshold,"pos","neg")
  
  
  # Count-part implementation
  count_part_classifier<-classifier_models[[paste0("Clinical_AUC_None_",jj)]]
  optimal_threshold<-count_thresholds[jj]
  
  predicted_class_prob<-predict(count_part_classifier,count_clinical_test,type="prob")[,2]
  predicted_label<-ifelse(predicted_class_prob>optimal_threshold,"pos","neg")
  
  
  low_regression_clinical<-regression_models[[paste0('Clinical_AUC_None_',jj)]]
  high_regression_clinical<-high_regression_models[[paste0('Clinical_AUC_None_',jj)]]
  
  predicted_value_count<-ifelse(predicted_label=="pos",predict(high_regression_clinical,count_clinical_test),predict(low_regression_clinical,count_clinical_test))
  
  
  
  #### Final prediction ###
  clinical_predicted<-ifelse(predicted_label_zero=="pos",predicted_value_count,0)
  clinical_preds<-c(clinical_preds,clinical_predicted)
  
  
  
  
  
  
  
  
  
  # For SNP now
  dd_init<-dataset[-list_indexes_test[[jj]],]
  no_nan_inds<-which(!is.na(dd_init$SYNTAX.SCORE))
  #1. Zero classification
  zero_classifier<-final_models[[paste0("SNP_RF_",jj)]]
  
  #2. define train/test for zero and count
  zero_snp_train<-train_snp_ls_processed_zero[[jj]]
  zero_snp_test<-test_snp_ls_processed_zero[[jj]]
  count_snp_train<-train_snp_ls_processed[[jj]]
  count_snp_test<-test_snp_ls_processed[[jj]]
  no_nan_inds<-no_nan_inds[1:nrow(zero_snp_train)]
  
  #3. Find the optimal threshold for the zero-part classifier
  ## Construct ROC curve
  dataset_snp_train_wanted<-dataset_snp_train[[jj]][no_nan_inds,]
  roc_obj <- roc(
    response  = ifelse(dataset_snp_train_wanted$SYNTAX.SCORE>0,1,0),
    predictor = predict(zero_classifier,zero_snp_train,type="prob")[,2],
    levels    = c(0, 1),  # negative class first, positive class second
    direction = "<"
  )
  
  ## Threshold maximizing Youden's index:
  ## sensitivity + specificity - 1
  optimal <- coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity", "youden"),
    transpose = FALSE
  )
  optimal_threshold <- as.numeric(optimal$threshold)
  
  predicted_class_prob<-predict(zero_classifier,zero_snp_test,type="prob")[,2]
  predicted_label_zero<-ifelse(predicted_class_prob>optimal_threshold,"pos","neg")
  
  
  # Count-part implementation
  count_part_classifier<-classifier_models[[paste0("SNP_AUC_None_",jj)]]
  optimal_threshold<-count_thresholds[jj]
  
  predicted_class_prob<-predict(count_part_classifier,count_snp_test,type="prob")[,2]
  predicted_label<-ifelse(predicted_class_prob>optimal_threshold,"pos","neg")
  
  
  low_regression_snp<-regression_models[[paste0('SNP_AUC_None_',jj)]]
  high_regression_snp<-high_regression_models[[paste0('SNP_AUC_None_',jj)]]
  
  predicted_value_count<-ifelse(predicted_label=="pos",predict(high_regression_snp,count_snp_test),predict(low_regression_snp,count_snp_test))
  
  
  
  #### Final prediction ###
  snp_predicted<-ifelse(predicted_label_zero=="pos",predicted_value_count,0)
  snp_preds<-c(snp_preds,snp_predicted)
  
  
  test_truth<-dataset[inds[[jj]],]
  test_truth<-test_truth$SYNTAX.SCORE
  
  test_all<-c(test_all,test_truth)
  folds_arr<-c(folds_arr,rep(jj,length(test_truth)))
}








pred_model_A<-clinical_preds
pred_model_B<-snp_preds

n_boot <- 5000
boot_medae_A <- numeric(n_boot)
boot_medae_B <- numeric(n_boot)
boot_delta_medae <- numeric(n_boot)

complete_index <- 1:799

for (i in seq_len(n_boot)) {
  
  boot_index <- sample(
    complete_index,
    length(complete_index),
    replace = TRUE
  )
  
  boot_medae_A[i] <- median(
    abs(y_true[boot_index] - pred_model_A[boot_index]),na.rm=TRUE
  )
  
  boot_medae_B[i] <- median(
    abs(y_true[boot_index] - pred_model_B[boot_index]),na.rm=TRUE
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



