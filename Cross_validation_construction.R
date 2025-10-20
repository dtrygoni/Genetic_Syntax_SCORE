
###################################### Before splitting all-procedure to test multiple splits ############################
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

dataset<-dataset[!is.na(dataset$SYNTAX.SCORE),]
dataset_positive<-dataset[dataset$SYNTAX.SCORE>0,]
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
              `Test_uncovered`=test_unconverted))
}
Univariate_LR_filtering<-function(train,test,predictors_clinical,predictors_snp,target_col,LR_threshold=0.10){
  clinical_predictors<-c()
  snp_predictors<-c()

  # for(predictor in predictors_clinical){
  #   form <- as.formula(paste(target_col, '~', paste0('`', predictor, '`')))
  #   res<-glm(form,family='binomial',data=train)
  #   pval <- summary(res)$coefficients[2, 4]
  #   if(pval<0.10){
  #     clinical_predictors<-c(clinical_predictors,predictor)
  #   }
  # }
  # 
  # ### For all modelling
  # for(predictor in predictors_snp){
  #   form <- as.formula(paste(target_col, '~', paste0('`', predictor, '`')))
  #   res<-glm(form,family='binomial',data=train)
  #   pval<-summary(res)$coefficients[2,4]
  #   if(pval<LR_threshold){
  #     snp_predictors<-c(snp_predictors,predictor)
  #   }
  # }
  for(predictor in predictors_clinical){
    form <- as.formula(paste(target_col, '~', paste0('`', predictor, '`')))
    res<-glm(form,family='binomial',data=train)
    train_temp<-train[!is.na(train[[predictor]]),]
    res_null<-glm(as.formula(paste0(target_col,' ~ 1')),family='binomial',data=train_temp)
    res_anova<-anova(res_null,res,test="Chisq")
    pval <- res_anova$`Pr(>Chi)`[2]
    if(pval<0.10){
      clinical_predictors<-c(clinical_predictors,predictor)
    }
  }

  ### For all modelling
  for(predictor in predictors_snp){
    form <- as.formula(paste(target_col, '~', paste0('`', predictor, '`')))

    if(is.factor(train[,predictor])){
      if(length(levels(train[,predictor]))<2){
        next
      }
    }
    res<-glm(form,family='binomial',data=train)
    train_temp<-train[!is.na(train[[predictor]]),]
    res_null<-glm(as.formula(paste0(target_col,' ~ 1')),family='binomial',data=train_temp)
    res_anova<-anova(res_null,res,test="Chisq")
    pval <- res_anova$`Pr(>Chi)`[2]
    if(pval<LR_threshold){
      snp_predictors<-c(snp_predictors,predictor)
    }
  }
  
  
  #### Clinical sets
  train_clinical<-train[colnames(train) %in% c(clinical_predictors,'Binary.SS')]
  train_clinical<-train_clinical[complete.cases(train_clinical),]
  test_clinical<-test[colnames(test) %in% c(clinical_predictors,'Binary.SS')]
  test_clinical<-test_clinical[complete.cases(test_clinical),]
  #### SNP sets
  train_snp<-train[colnames(train) %in% c(clinical_predictors,snp_predictors,'Binary.SS')]
  train_snp<-train_snp[complete.cases(train_snp),]
  test_snp<-test[colnames(test) %in% c(clinical_predictors,snp_predictors,'Binary.SS')]
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
  
  
  return(list(`ClinicalTrain`=clinical_train,`ClinicalTest`=clinical_test,`SNPtrain`=snp_train,`SNPtest`=snp_train,`ClinicalScaler`=scalerClinical_obj,
              `ClinicalSNPScaler`=scalerSNP_obj))
}
Model_training<-function(clinical_train,snp_train,model_selection='XGB',number_iterations=100,optimazation_metric='balAcc'){
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
    summaryFunction =summary_fnc,
    verboseIter = TRUE
  )
  tune_control <- caret::trainControl(
    method = "cv", # cross-validation
    number = 4, # with n folds 60% train - 20% validation - 20% test
    #index = createFolds(tr_treated$Id_clean), # fix the folds
    classProbs = TRUE,
    verboseIter = FALSE, # no training log
    allowParallel = TRUE # FALSE for reproducible results 
  )

  if(model_selection=='XGB'){
  # Fit the model with a random subset of the grid (e.g., 20 random combinations)
  xgb_model_clinical <- train(
    x = input_x_clinical,
    y = as.factor(input_y_clinical),          
    method = "xgbTree",
    trControl = ctrl,
    tuneLength =number_iterations,
    metric=optimazation_metric
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
    tuneLength =2,
    metric=optimazation_metric
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
      metric=optimazation_metric
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
      metric=optimazation_metric
    )
    
    (final_grid_snp <- expand.grid(
      mtry=rf_model_snp$bestTune$mtry
    ))
    
    (model_SNP <- caret::train(
      x = input_x_snp,
      y = as.factor(input_y_snp),
      trControl = tune_control,
      tuneGrid = final_grid,
      method = "rf",
      verbose = TRUE
    ))}else if(model_selection=='kNN'){
   knn_model_clinical <- train(
      x = input_x_clinical,
      y = as.factor(input_y_clinical),          
      method = "knn",
      trControl = ctrl,
      tuneLength =number_iterations,
      metric=optimazation_metric
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
     metric=optimazation_metric
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
### Parameters
genetic_modelling<-'Optimal'#Codominant, Optimal
HWE_threshold<-0.05
MAF_threshold<-0.05
clinical_cols<-1:47
LR_threshold<-0.10
Boruta_threshold<-0.05
strategy_boruta<-'All'# 'Seperate','All'
target_col<-'Binary.SS'
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
  ############### SNP quality control ################
  res_control<-quality_controlSNP(train_set,test_set,snp_cols,genetic_modelling=genetic_modelling,HWE_threshold=HWE_threshold,MAF_threshold=MAF_threshold)
  train_set_all<-res_control$Train_set
  test_set_all<-res_control$Test_set
  levels_all<-res_control$Levels
  ### For interaction ###
  train_unconverted<-res_control$Train_unconverted
  test_unconverted<-res_control$Test_uncovered
  
  
  
  
  #### rest of the code #####
  
  
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
  res_boruta<-Boruta_filtering(train_snp,test_snp,clinical_predictors = clinical_predictors,snp_predictors = snp_predictors,
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


save(train_clinical_ls,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Train_clinical_CV_BA_DM.Rda')
save(test_clinical_ls,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Test_clinical_CV_BA_DM.Rda')
save(train_snp_ls,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Train_snp_CV_BA_DM.Rda')
save(test_snp_ls,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Test_snp_CV_BA_DM.Rda')












################################## Clinical + SNP predictors for first step --- Zero part #####################
clinical_predictors_all<-c()

snp_predictors_all<-c()

for(jj in 1:10){
  print(jj)
  train_set<-dataset[-list_indexes_test[[as.character(jj)]],]
  test_set<-dataset[list_indexes_test[[as.character(jj)]],]
  ############### SNP quality control ################
  res_control<-quality_controlSNP(train_set,test_set,snp_cols,genetic_modelling=genetic_modelling,HWE_threshold=HWE_threshold,MAF_threshold=MAF_threshold)
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


names_all<-names(table(snp_predictors_all))
table(snp_predictors_all)['rs2306374' == names_all]











