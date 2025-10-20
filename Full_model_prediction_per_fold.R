############## Prediction of two step models per fold ########
library(caret)


names_clinical<-c()
names_snp<-c()
for(name in names(zero_models)){
  assoc <-strsplit(name,split="_")[[1]][1]
  if(assoc=='Clinical'){
    names_clinical<-c(names_clinical,name)
  }else{
    names_snp<-c(names_snp,name)
  }
}
predictors_zero_clinical<-c()
predictors_zero_snp<-c()
for(name in names_clinical){
  model<-zero_models[[name]]
  predictors_arr<-predictors(model)
  predictors_zero_clinical<-c(predictors_zero_clinical,predictors_arr)
  
}
for(name in names_snp){
  model<-zero_models[[name]]
  predictors_arr<-predictors(model)
  predictors_zero_snp<-c(predictors_zero_snp,predictors_arr)
  
}

save(predictors_zero_clinical,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Whole_model_experiments_classifiers/Predictors/Predictors_zero_clinical.Rda')
save(predictors_zero_snp,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Whole_model_experiments_classifiers/Predictors/Predictors_zero_SNP.Rda')

names_clinical<-c()
names_snp<-c()
for(name in names(count_classifiers)){
  assoc <-strsplit(name,split="_")[[1]][1]
  if(assoc=='Clinical'){
    names_clinical<-c(names_clinical,name)
  }else{
    names_snp<-c(names_snp,name)
  }
}
predictors_count_clinical<-c()
predictors_count_snp<-c()
for(name in names_clinical){
  model<-count_classifiers[[name]]
  predictors_arr<-predictors(model)
  predictors_count_clinical<-c(predictors_count_clinical,predictors_arr)
  
}
for(name in names_snp){
  model<-count_classifiers[[name]]
  predictors_arr<-predictors(model)
  predictors_count_snp<-c(predictors_count_snp,predictors_arr)
  
}
save(predictors_count_clinical,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Whole_model_experiments_classifiers/Predictors/Predictors_count_clinical.Rda')
save(predictors_count_snp,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Whole_model_experiments_classifiers/Predictors/Predictors_count_SNP.Rda')

####### Load models 
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Cross_validation_experiments/Extra/Models_RF_AUC_None.Rda')
zero_models<-final_models
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Classifier_Models.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Regression_Clinical_High_Models.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Regression_SNP_High_Models.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Regression_Clinical_Low_Models.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Regression_SNP_Low_Models.Rda')
count_classifiers<-classifier_models
count_regressor_low_cl<-regression_cl_model_low
count_regressor_high_cl<-regression_cl_model_high
count_regressor_low_snp<-regression_snp_model_low
count_regressor_high_snp<-regression_snp_model_high


#count_scale<-'log'
#count_metric<-'MAE'
####### Load tests
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Whole_model_experiments_classifiers/Zero_part/Clinical_test.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Whole_model_experiments_classifiers/Zero_part/SNP_test.Rda')

load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Whole_model_experiments_classifiers/Count_part/Clinical_test.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Whole_model_experiments_classifiers/Count_part/SNP_test.Rda')




Count_part_prediction<-function(test,classifier,regressor_low,regressor_high){
  preds_classifiers<-predict(classifier,newdata=test)
  preds_regressor_low<-predict(regressor_low,newdata=test)
  preds_regressor_high<-predict(regressor_high,newdata = test)
  
  preds_f<-round(ifelse(preds_classifiers=='pos',preds_regressor_high,preds_regressor_low),2)
  return(preds_f)
}


Count_part_prediction_soft<-function(test,classifier,regressor_low,regressor_high){
  preds_classifiers<-predict(classifier,newdata=test,type="prob")
  preds_regressor_low<-predict(regressor_low,newdata=test)
  preds_regressor_high<-predict(regressor_high,newdata = test)
  preds_f<-round(preds_classifiers[,1]*preds_regressor_low+preds_classifiers[,2]*preds_regressor_high,2)
  return(preds_f)
}

############ Clinical evaluation ##############

# Evaluate zero if positive then count model
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
      snp_zero<-final_test_snp_zero_ls[[jj]]
      snp_count<-final_test_snp_count_ls[[jj]]
      
      complete_snp<-complete.cases(snp_zero)
      complete_snp_count<-complete.cases(snp_count)
      
      complete_snp_all<-complete_snp & complete_snp_count
      
      
      clinical_zero<-final_test_clinical_zero_ls[[jj]]
      clinical_count<-final_test_clinical_count_ls[[jj]]
      
      complete_cl<-complete.cases(clinical_zero)
      complete_count<-complete.cases(clinical_count)
      
      complete_cl<-complete_cl & complete_count
      
      
      complete<-complete_cl & complete_snp_all
      
      clinical_zero<-clinical_zero[complete,]
      clinical_count<-clinical_count[complete,]
      
      
      zero_jj<-zero_models[paste0('Clinical_RF_',jj)]
      preds<-predict(zero_jj,clinical_zero)[[1]]
      preds_probs<-predict(zero_jj,clinical_zero,type='prob')[[1]]
      probs_binary<-c(probs_binary,preds_probs[,2])
      
      #Clinical
      name<-paste0('Clinical_',opt,'_',resamp,'_',jj)
      classifier<-count_classifiers[[name]]
      name_regr<-paste0('Clinical_AUC_',resamp,'_',jj)
      regressor_low<-count_regressor_low_cl[[name_regr]]
      regressor_high<-count_regressor_high_cl[[name_regr]]
      preds_count<-Count_part_prediction(clinical_count,classifier,regressor_low,regressor_high)
      
      
      predictions<-ifelse(preds=='neg',0,preds_count)
      predictions_clinical<-c(predictions_clinical,predictions)
      
      
      target<-clinical_count$SYNTAX.SCORE
      target_clinical<-c(target_clinical,target)
      
      ae<-median(abs(target-predictions))
      ae_arr<-c(ae_arr,ae)
      fold_arr<-c(fold_arr,jj)
    }
    name<-paste0('Clinical_',opt,'_',resamp)
    predictions_cl_all[[name]]<-predictions_clinical
    target_cl_all[[name]]<-target_clinical
    probs_binary_cl_all[[name]]<-probs_binary
  }
}

df<-data.frame(Fold=fold_arr,MdAE.Clinical=ae_clinical)




############ SNP evaluation ##############

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
      snp_zero<-final_test_snp_zero_ls[[jj]]
      snp_count<-final_test_snp_count_ls[[jj]]
      
      complete_snp<-complete.cases(snp_zero)
      complete_snp_count<-complete.cases(snp_count)
      
      complete_snp_all<-complete_snp & complete_snp_count
      
      
      clinical_zero<-final_test_clinical_zero_ls[[jj]]
      clinical_count<-final_test_clinical_count_ls[[jj]]
      
      complete_cl<-complete.cases(clinical_zero)
      complete_count<-complete.cases(clinical_count)
      
      complete_cl<-complete_cl & complete_count
      
      
      complete<-complete_cl & complete_snp_all
      
      
      
      snp_zero<-snp_zero[complete,]
      snp_count<-snp_count[complete,]
      
      
      zero_jj<-zero_models[paste0('SNP_RF_',jj)]
      preds<-predict(zero_jj,snp_zero)[[1]]
      
      preds_probs<-predict(zero_jj,snp_zero,type='prob')[[1]]
      probs_binary<-c(probs_binary,preds_probs[,2])
      #snp
      name<-paste0('SNP_',opt,'_',resamp,'_',jj)
      classifier<-count_classifiers[[name]]
      name_regr<-paste0('SNP_AUC_',resamp,'_',jj)
      regressor_low<-count_regressor_low_snp[[name_regr]]
      regressor_high<-count_regressor_high_snp[[name_regr]]
      preds_count<-Count_part_prediction(snp_count,classifier,regressor_low,regressor_high)
      
      
      predictions<-ifelse(preds=='neg',0,preds_count)
      predictions_snp<-c(predictions_snp,predictions)
      
      
      target<-snp_count$SYNTAX.SCORE
      target_snp<-c(target_snp,target)
      
      ae<-median(abs(target-predictions))
      
      ae_snp<-c(ae_snp,ae)
    }
    name<-paste0('SNP_',opt,'_',resamp)
    predictions_snp_all[[name]]<-predictions_snp
    target_snp_all[[name]]<-target_snp
    probs_binary_snp_all[[name]]<-probs_binary
  }
}

df$AE.SNP<-ae_snp

median(df$MdAE.Clinical) #7.255
median(df$AE.SNP) #9.725













################################ SOFT #############################
opt_arr<-c("AUC")
resamp_arr<-c("None")
predictions_cl_all<-list()
target_cl_all<-list()
for(opt in opt_arr){
  for(resamp in resamp_arr){
    predictions_clinical<-c()
    target_clinical<-c()
    fold_arr<-c()
    ae_clinical<-c()
    for(jj in 1:10){
      snp_zero<-final_test_snp_zero_ls[[jj]]
      snp_count<-final_test_snp_count_ls[[jj]]
      
      complete_snp<-complete.cases(snp_zero)
      complete_snp_count<-complete.cases(snp_count)
      
      complete_snp_all<-complete_snp & complete_snp_count
      
      
      clinical_zero<-final_test_clinical_zero_ls[[jj]]
      clinical_count<-final_test_clinical_count_ls[[jj]]
      
      complete_cl<-complete.cases(clinical_zero)
      complete_count<-complete.cases(clinical_count)
      
      complete_cl<-complete_cl & complete_count
      
      
      complete<-complete_cl & complete_snp_all
      
      clinical_zero<-clinical_zero[complete,]
      clinical_count<-clinical_count[complete,]
      
      
      zero_jj<-zero_models[paste0('Clinical_RF_',jj)]
      preds_probs<-predict(zero_jj,clinical_zero,type='prob')[[1]]
      
      
      #Clinical
      name<-paste0('Clinical_',opt,'_',resamp,'_',jj)
      classifier<-count_classifiers[[name]]
      name_regr<-paste0('Clinical_AUC_',resamp,'_',jj)
      regressor_low<-count_regressor_low_cl[[name_regr]]
      regressor_high<-count_regressor_high_cl[[name_regr]]
      preds_count<-Count_part_prediction_soft(clinical_count,classifier,regressor_low,regressor_high)
      
      
      predictions<-round(preds_probs[,2]*preds_count,2)
      predictions_clinical<-c(predictions_clinical,predictions)
      
      
      target<-clinical_count$SYNTAX.SCORE
      target_clinical<-c(target_clinical,target)
      
      fold_arr<-c(fold_arr,jj)
      ae<-median(abs(target-predictions))
      ae_clinical<-c(ae_clinical,ae)
    }
    name<-paste0('Clinical_',opt,'_',resamp)
    predictions_cl_all[[name]]<-predictions_clinical
    target_cl_all[[name]]<-target_clinical
  }
}
df_soft<-data.frame(Fold=fold_arr,AE.Clinical=ae_clinical)


opt_arr<-c("AUC")
resamp_arr<-c("None")
predictions_snp_all<-list()
target_snp_all<-list()
for(opt in opt_arr){
  for(resamp in resamp_arr){
    predictions_snp<-c()
    target_snp<-c()
    ae_snp<-c()
    for(jj in 1:10){
      snp_zero<-final_test_snp_zero_ls[[jj]]
      snp_count<-final_test_snp_count_ls[[jj]]
      
      complete_snp<-complete.cases(snp_zero)
      complete_snp_count<-complete.cases(snp_count)
      
      complete_snp_all<-complete_snp & complete_snp_count
      
      
      clinical_zero<-final_test_clinical_zero_ls[[jj]]
      clinical_count<-final_test_clinical_count_ls[[jj]]
      
      complete_cl<-complete.cases(clinical_zero)
      complete_count<-complete.cases(clinical_count)
      
      complete_cl<-complete_cl & complete_count
      
      
      complete<-complete_cl & complete_snp_all
      
      
      
      snp_zero<-snp_zero[complete,]
      snp_count<-snp_count[complete,]
      
      
      zero_jj<-zero_models[paste0('SNP_RF_',jj)]
      preds_probs<-predict(zero_jj,snp_zero,type="prob")[[1]]
      
      
      #snp
      name<-paste0('SNP_',opt,'_',resamp,'_',jj)
      classifier<-count_classifiers[[name]]
      name_regr<-paste0('SNP_AUC_',resamp,'_',jj)
      regressor_low<-count_regressor_low_snp[[name_regr]]
      regressor_high<-count_regressor_high_snp[[name_regr]]
      preds_count<-Count_part_prediction_soft(snp_count,classifier,regressor_low,regressor_high)
      
      
      predictions<-round(preds_probs[,2]*preds_count,2)
      predictions_snp<-c(predictions_snp,predictions)
      
      
      target<-snp_count$SYNTAX.SCORE
      target_snp<-c(target_snp,target)
      
      ae<-median(abs(target-predictions))
      ae_snp<-c(ae_snp,ae)
    }
    name<-paste0('SNP_',opt,'_',resamp)
    predictions_snp_all[[name]]<-predictions_snp
    target_snp_all[[name]]<-target_snp
  }
}

df_soft$AE.SNP<-ae_snp


median(df_soft$AE.Clinical) #7.69
median(df_soft$AE.SNP) #8.0475
