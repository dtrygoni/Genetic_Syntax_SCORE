############## Prediction of two step models ########
library(caret)

########## Latest change for different threshold ##########
######### thresholds clinical+SNP are in Zero part experiments folder in CI_ROC.R ##############
######## thresholds for count part clinical+SNP are Metrics_evaluation_per_fold.R ###############
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

save(predictors_zero_clinical,file='./Predictors_zero_clinical.Rda')
save(predictors_zero_snp,file='./Predictors_zero_SNP.Rda')

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
save(predictors_count_clinical,file='./Predictors_count_clinical.Rda')
save(predictors_count_snp,file='./Predictors_count_SNP.Rda')

####### Load models 
load('./Models_RF_AUC_None.Rda')
zero_models<-final_models
load('./Classifier_Models.Rda')
load('./Regression_Clinical_High_Models.Rda')
load('./Regression_SNP_High_Models.Rda')
load('./Regression_Clinical_Low_Models.Rda')
load('./Regression_SNP_Low_Models.Rda')
count_classifiers<-classifier_models
count_regressor_low_cl<-regression_cl_model_low
count_regressor_high_cl<-regression_cl_model_high
count_regressor_low_snp<-regression_snp_model_low
count_regressor_high_snp<-regression_snp_model_high


#count_scale<-'log'
#count_metric<-'MAE'
####### Load tests
load('./Zero_part/Clinical_test.Rda')
load('./Zero_part/SNP_test.Rda')

load('./Count_part/Clinical_test.Rda')
load('./Count_part/SNP_test.Rda')




Count_part_prediction<-function(test,classifier,regressor_low,regressor_high){
  preds_classifiers<-predict(classifier,newdata=test)
  preds_regressor_low<-predict(regressor_low,newdata=test)
  preds_regressor_high<-predict(regressor_high,newdata = test)
  
  preds_f<-round(ifelse(preds_classifiers=='pos',preds_regressor_high,preds_regressor_low),2)
  return(preds_f)
}
Count_part_prediction_new<-function(test,classifier,regressor_low,regressor_high,threshold){
  preds_classifiers<-predict(classifier,newdata=test,type='prob')[,2]
  preds_regressor_low<-predict(regressor_low,newdata=test)
  preds_regressor_high<-predict(regressor_high,newdata = test)
  
  preds_f<-round(ifelse(preds_classifiers>=threshold,preds_regressor_high,preds_regressor_low),2)
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
fold_arr<-c()
for(opt in opt_arr){
  for(resamp in resamp_arr){
    predictions_clinical<-c()
    target_clinical<-c()
    probs_binary<-c()
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
      preds_probs<-predict(zero_jj,clinical_zero,type='prob')[[1]][,2]
      probs_binary<-c(probs_binary,preds_probs)
      
      #Clinical
      name<-paste0('Clinical_',opt,'_',resamp,'_',jj)
      classifier<-count_classifiers[[name]]
      name_regr<-paste0('Clinical_AUC_',resamp,'_',jj)
      regressor_low<-count_regressor_low_cl[[name_regr]]
      regressor_high<-count_regressor_high_cl[[name_regr]]
      #preds_count<-Count_part_prediction(clinical_count,classifier,regressor_low,regressor_high)
      preds_count<-Count_part_prediction_new(clinical_count,classifier,regressor_low,regressor_high,threshold=threshold_classifier[jj])
      
      
      
      
      predictions<-ifelse(preds_probs>=threshold_clinical[jj],preds_count,0)
      predictions_clinical<-c(predictions_clinical,predictions)
      
      
      target<-clinical_count$SYNTAX.SCORE
      target_clinical<-c(target_clinical,target)
      fold_arr<-c(fold_arr,rep(jj,length(target)))
    }
    name<-paste0('Clinical_',opt,'_',resamp)
    predictions_cl_all[[name]]<-predictions_clinical
    target_cl_all[[name]]<-target_clinical
    probs_binary_cl_all[[name]]<-probs_binary
  }
}

df<-data.frame(Actual=target_cl_all$Clinical_AUC_None,Fold=fold_arr,Predicted.Whole.Clinical.Hard=predictions_cl_all$Clinical_AUC_None,
               Zero.ProbsPositive.Clinical=probs_binary_cl_all$Clinical_AUC_None)

res_list<-list()
for(name in names(probs_binary_cl_all)){
  df<-data.frame(Actual=target_cl_all[[name]],Predicted.Whole.Clinical.Hard=predictions_cl_all[[name]],
                                  Zero.ProbsPositive.Clinical=probs_binary_cl_all[[name]])
  
  components<-strsplit(name,split="_")[[1]]
  name_res<-paste0(components[2],'_',components[3])
  res_list[[name_res]]<-df
}


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
      
      preds_probs<-predict(zero_jj,snp_zero,type='prob')[[1]][,2]
      probs_binary<-c(probs_binary,preds_probs)
      #snp
      name<-paste0('SNP_',opt,'_',resamp,'_',jj)
      classifier<-count_classifiers[[name]]
      name_regr<-paste0('SNP_AUC_',resamp,'_',jj)
      regressor_low<-count_regressor_low_snp[[name_regr]]
      regressor_high<-count_regressor_high_snp[[name_regr]]
      preds_count<-Count_part_prediction_new(snp_count,classifier,regressor_low,regressor_high,threshold_classifier_snp[jj])
      
      
      predictions<-ifelse(preds_probs>=threshold_clinical[jj],preds_count,0)
      predictions_snp<-c(predictions_snp,predictions)
      
      
      target<-snp_count$SYNTAX.SCORE
      target_snp<-c(target_snp,target)
    }
    name<-paste0('SNP_',opt,'_',resamp)
    predictions_snp_all[[name]]<-predictions_snp
    target_snp_all[[name]]<-target_snp
    probs_binary_snp_all[[name]]<-probs_binary
  }
}



for(name in names(probs_binary_snp_all)){
  predictions_snp<-predictions_snp_all[[name]]
  probs_snp<-probs_binary_snp_all[[name]]
  
  components<-strsplit(name,split="_")[[1]]
  name_list<-paste0(components[2],'_',components[3])
  
  df<-res_list[[name_list]]
  
  df$`Predicted.Whole.SNP.Hard`<-predictions_snp
  df$`Zero.ProbsPositive.SNP`<-probs_snp
  
  res_list[[name_list]]<-df
}


df<-res_list$AUC_None
df$Absolute.Error.Clinical<-abs(df$Actual-df$Predicted.Whole.Clinical.Hard)
df$Relative.Error.Clinical<-abs(df$Actual-df$Predicted.Whole.Clinical.Hard)/df$Actual


df$Absolute.Error.SNP<-abs(df$Actual-df$Predicted.Whole.SNP.Hard)
df$Relative.Error.SNP<-abs(df$Actual-df$Predicted.Whole.SNP.Hard)/df$Actual


median(df$Absolute.Error.Clinical) #10.04
median(df$Absolute.Error.SNP) #9.02

res<-data.frame(Fold=fold_arr,Clinical.Predicted=predictions_cl_all$Clinical_AUC_None,
                SNP.Predicted=df$Predicted.Whole.SNP.Hard,Target=target_cl_all$Clinical_AUC_None)
save(res,file='./Count_part_classifiers/Full_model.Rda')
median(abs(res$Target-res$Clinical.Predicted))
median(abs(res$Target-res$SNP.Predicted))

save(df,file='./Results/Dataframe_Hard.Rda')

############# Save results ################

save(predictions_cl_all,file='./Results/Clinical_predictions_hard.Rda')
save(target_cl_all,file='./Results/Clinical_target.Rda')

save(predictions_snp_all,file='./Results/SNP_predictions_hard.Rda')
save(target_snp_all,file='./Results/SNP_target.Rda')



########################################## SOFT EVALUATION ########################################
############ Clinical evaluation ##############

# Evaluate zero if positive then count model
opt_arr<-c("AUC")
resamp_arr<-c("None")
predictions_cl_all<-list()
target_cl_all<-list()
for(opt in opt_arr){
  for(resamp in resamp_arr){
    predictions_clinical<-c()
    target_clinical<-c()
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
    }
    name<-paste0('Clinical_',opt,'_',resamp)
    predictions_cl_all[[name]]<-predictions_clinical
    target_cl_all[[name]]<-target_clinical
  }
}


for(name in names(predictions_cl_all)){
  predictions_soft<-predictions_cl_all[[name]]
  #probs_snp<-probs_binary_snp_all[[name]]
  
  components<-strsplit(name,split="_")[[1]]
  name_list<-paste0(components[2],'_',components[3])
  
  df<-res_list[[name_list]]
  
  df$`Predicted.Whole.Clinical.Soft`<-predictions_soft
  #df$`Zero.ProbsPositive.SNP`<-probs_snp
  
  res_list[[name_list]]<-df
}

#df$Predicted.Whole.Clinical.Soft <-predictions_cl_all$Clinical_AUC_None







############ SNP evaluation ##############

# Evaluate zero if positive then count model
opt_arr<-c("AUC")
resamp_arr<-c("None")
predictions_snp_all<-list()
target_snp_all<-list()
for(opt in opt_arr){
  for(resamp in resamp_arr){
    predictions_snp<-c()
    target_snp<-c()
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
    }
    name<-paste0('SNP_',opt,'_',resamp)
    predictions_snp_all[[name]]<-predictions_snp
    target_snp_all[[name]]<-target_snp
  }
}


for(name in names(predictions_snp_all)){
  predictions_soft<-predictions_snp_all[[name]]
  #probs_snp<-probs_binary_snp_all[[name]]
  
  components<-strsplit(name,split="_")[[1]]
  name_list<-paste0(components[2],'_',components[3])
  
  df<-res_list[[name_list]]
  
  df$`Predicted.Whole.SNP.Soft`<-predictions_soft
  #df$`Zero.ProbsPositive.SNP`<-probs_snp
  
  res_list[[name_list]]<-df
}

#df$Predicted.Whole.SNP.Soft <-predictions_snp_all$SNP_AUC_None


#save(df,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Whole_model_experiments_classifiers/Results/Results_predictions.Rda')


df<-res_list$AUC_None
df<-df[c("Actual","Predicted.Whole.Clinical.Soft","Predicted.Whole.SNP.Soft","Zero.ProbsPositive.Clinical","Zero.ProbsPositive.SNP")]
df$Absolute.Error.Clinical<-abs(df$Actual-df$Predicted.Whole.Clinical.Soft)
df$Relative.Error.Clinical<-abs(df$Actual-df$Predicted.Whole.Clinical.Soft)/df$Actual


df$Absolute.Error.SNP<-abs(df$Actual-df$Predicted.Whole.SNP.Soft)
df$Relative.Error.SNP<-abs(df$Actual-df$Predicted.Whole.SNP.Soft)/df$Actual
save(df,file='./Results/Dataframe_Soft.Rda')

save(res_list,file='./Results_predictions_allresamplings.Rda')




############# Save results ################

save(predictions_cl_all,file='./Results/Clinical_predictions_soft.Rda')
save(target_cl_all,file='./Results/Clinical_target.Rda')

save(predictions_snp_all,file='./Results/SNP_predictions_soft.Rda')
save(target_snp_all,file='./Results/SNP_target.Rda')
























######################## Wilcoxon results for res_list ##################################

wilcoxon_res<-NULL
for(name in names(res_list)){
  df<-res_list[[name]]
  
  ae_clinical_hard<-abs(df$Predicted.Whole.Clinical.Hard-df$Actual)
  ae_clinical_soft<-abs(df$Predicted.Whole.Clinical.Soft-df$Actual)
  
  ae_snp_hard<-abs(df$Predicted.Whole.SNP.Hard-df$Actual)
  ae_snp_soft<-abs(df$Predicted.Whole.SNP.Soft-df$Actual)
  

  res_hard_ae<-wilcox.test(ae_snp_hard,ae_clinical_hard)
  res_soft_ae<-wilcox.test(ae_snp_soft,ae_clinical_soft)
  
  relesterror_clinical_soft<-abs(df$Predicted.Whole.Clinical.Soft-df$Actual)/df$Predicted.Whole.Clinical.Soft
  relesterror_snp_soft<-abs(df$Predicted.Whole.SNP.Soft-df$Actual)/df$Predicted.Whole.SNP.Soft
  res_relestierr_soft<-wilcox.test(relesterror_clinical_soft,relesterror_snp_soft)
  
  
  components<-strsplit(name,split="_")[[1]]
  opt_metric<-components[1]
  resamp<-components[2]
  
  
  res_df<-data.frame(Optimization.Metric=opt_metric,
                     Resampling=resamp,
                     Med.AE.Clinical.Hard=median(ae_clinical_hard),
                     Med.AE.SNP.Hard=median(ae_snp_hard),
                     Wilcox.p.val.Hard=round(res_hard_ae$p.value,4),
                     Med.AE.Clinical.Soft=median(ae_clinical_soft),
                     Med.AE.SNP.Soft=median(ae_snp_soft),
                     Wilcox.p.val.Soft=round(res_soft_ae$p.value,4),
                     Med.Rel.Estimate.Error.Clinical.Soft=median(relesterror_clinical_soft),
                     Med.Rel.Estimate.Error.SNP.Soft=median(relesterror_snp_soft),
                     Wilcox.p.val.REE=round(res_relestierr_soft$p.value,4))
  wilcoxon_res<-rbind(wilcoxon_res,res_df)
  }
save(wilcoxon_res,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Whole_model_experiments_classifiers/Results/Wilcoxon_Results_allresamplings.Rda')

