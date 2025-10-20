load('./DataFrame_AE_Hard.Rda')

library(caret)

median(ae_hard_df[ae_hard_df$Resampling=='None','Clinical'])
median(ae_hard_df[ae_hard_df$Resampling=='None','SNP'])

## AE hard clinical (6.905 ) SNP (6.9075)


load('./DataFrame_RAE_Hard.Rda')
median(rae_hard_df[rae_hard_df$Resampling=='None','Clinical'])
median(rae_hard_df[rae_hard_df$Resampling=='None','SNP'])
## RAE hard clinical (0.4825 ) SNP (0.492)



load('./DataFrame_AE_Soft.Rda')
median(ae_soft_df[ae_soft_df$Resampling=='None','Clinical'])
median(ae_soft_df[ae_soft_df$Resampling=='None','SNP'])
##AE hard clinical (7.0795 ) SNP (7.2485)



load('./DataFrame_RAE_Soft.Rda')
median(rae_soft_df[rae_soft_df$Resampling=='None','Clinical'])
median(rae_soft_df[rae_soft_df$Resampling=='None','SNP'])
##AE hard clinical (0.4125 ) SNP (0.4455)



















######################## Again calculation #######################
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Classifier_Models.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Regression_Clinical_High_Models.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Regression_Clinical_Low_Models.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Regression_SNP_High_Models.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Regression_SNP_Low_Models.Rda')


load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Test_sets.Rda')
load('C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Extra/Train_sets.Rda')


type<-'hard'
ae_clinical<-c()
ae_snp<-c()

for(fold in 1:10){
  test_clinical<-test_list[[paste0('Clinical_AUC_None_',fold)]]
  
  test_snp<-test_list[[paste0('SNP_AUC_None_',fold)]]
  
  target<-test_clinical$SYNTAX.SCORE
  
  #Clinical_prediction
  classifier<-classifier_models[[paste0('Clinical_AUC_None_',fold)]]
  classif_risk<-predict(classifier,test_clinical)
  regressor_low<-regression_cl_model_low[[paste0('Clinical_AUC_None_',fold)]]
  regressor_high<-regression_cl_model_high[[paste0('Clinical_AUC_None_',fold)]]
  
  preds_low<-predict(regressor_low,test_clinical)
  preds_high<-predict(regressor_high,test_clinical)
  
  if(type=='hard'){
    preds_clinical<-round(ifelse(classif_risk=='neg',preds_low,preds_high),2)
  }else{
    classif_probs<-predict(classifier,test_clinical,type="prob")
    preds_clinical<-round(classif_probs[,1]*preds_low+classif_probs[,2]*preds_high,2)
  }
  
  ae_cl<-median(abs(target-preds_clinical))
  ae_clinical<-c(ae_clinical,ae_cl)
  
  #SNP_prediction
  classifier<-classifier_models[[paste0('SNP_AUC_None_',fold)]]
  classif_risk<-predict(classifier,test_snp)
  regressor_low<-regression_snp_model_low[[paste0('SNP_AUC_None_',fold)]]
  regressor_high<-regression_snp_model_high[[paste0('SNP_AUC_None_',fold)]]
  
  preds_low<-predict(regressor_low,test_snp)
  preds_high<-predict(regressor_high,test_snp)
  
  if(type=='hard'){
    preds_snp<-round(ifelse(classif_risk=='neg',preds_low,preds_high),2)
  }else{
    classif_probs<-predict(classifier,test_clinical,type="prob")
    preds_snp<-round(classif_probs[,1]*preds_low+classif_probs[,2]*preds_high,2)
  }
  
  ae_sn<-median(abs(target-preds_snp))
  ae_snp<-c(ae_clinical,ae_sn)
}


median(ae_clinical) #7.00
median(ae_snp) #7.26


type<-'soft'
ae_clinical<-c()
ae_snp<-c()

for(fold in 1:10){
  test_clinical<-test_list[[paste0('Clinical_AUC_None_',fold)]]
  
  test_snp<-test_list[[paste0('SNP_AUC_None_',fold)]]
  
  target<-test_clinical$SYNTAX.SCORE
  
  #Clinical_prediction
  classifier<-classifier_models[[paste0('Clinical_AUC_None_',fold)]]
  classif_risk<-predict(classifier,test_clinical)
  regressor_low<-regression_cl_model_low[[paste0('Clinical_AUC_None_',fold)]]
  regressor_high<-regression_cl_model_high[[paste0('Clinical_AUC_None_',fold)]]
  
  preds_low<-predict(regressor_low,test_clinical)
  preds_high<-predict(regressor_high,test_clinical)
  
  if(type=='hard'){
    preds_clinical<-round(ifelse(classif_risk=='neg',preds_low,preds_high),2)
  }else{
    classif_probs<-predict(classifier,test_clinical,type="prob")
    preds_clinical<-round(classif_probs[,1]*preds_low+classif_probs[,2]*preds_high,2)
  }
  
  ae_cl<-median(abs(target-preds_clinical))
  ae_clinical<-c(ae_clinical,ae_cl)
  
  #SNP_prediction
  classifier<-classifier_models[[paste0('SNP_AUC_None_',fold)]]
  classif_risk<-predict(classifier,test_snp)
  regressor_low<-regression_snp_model_low[[paste0('SNP_AUC_None_',fold)]]
  regressor_high<-regression_snp_model_high[[paste0('SNP_AUC_None_',fold)]]
  
  preds_low<-predict(regressor_low,test_snp)
  preds_high<-predict(regressor_high,test_snp)
  
  if(type=='hard'){
    preds_snp<-round(ifelse(classif_risk=='neg',preds_low,preds_high),2)
  }else{
    classif_probs<-predict(classifier,test_snp,type="prob")
    preds_snp<-round(classif_probs[,1]*preds_low+classif_probs[,2]*preds_high,2)
  }
  
  ae_sn<-median(abs(target-preds_snp))
  ae_snp<-c(ae_clinical,ae_sn)
}
median(ae_clinical) #7.26
median(ae_snp) #7.30



type<-'hard'
rae_clinical<-c()
rae_snp<-c()

for(fold in 1:10){
  test_clinical<-test_list[[paste0('Clinical_AUC_None_',fold)]]
  
  test_snp<-test_list[[paste0('SNP_AUC_None_',fold)]]
  
  target<-test_clinical$SYNTAX.SCORE
  
  #Clinical_prediction
  classifier<-classifier_models[[paste0('Clinical_AUC_None_',fold)]]
  classif_risk<-predict(classifier,test_clinical)
  regressor_low<-regression_cl_model_low[[paste0('Clinical_AUC_None_',fold)]]
  regressor_high<-regression_cl_model_high[[paste0('Clinical_AUC_None_',fold)]]
  
  preds_low<-predict(regressor_low,test_clinical)
  preds_high<-predict(regressor_high,test_clinical)
  
  if(type=='hard'){
    preds_clinical<-round(ifelse(classif_risk=='neg',preds_low,preds_high),2)
  }else{
    classif_probs<-predict(classifier,test_clinical,type="prob")
    preds_clinical<-round(classif_probs[,1]*preds_low+classif_probs[,2]*preds_high,2)
  }
  
  rae_cl<-median(abs(target-preds_clinical)/target)
  rae_clinical<-c(rae_clinical,rae_cl)
  
  #SNP_prediction
  classifier<-classifier_models[[paste0('SNP_AUC_None_',fold)]]
  classif_risk<-predict(classifier,test_snp)
  regressor_low<-regression_snp_model_low[[paste0('SNP_AUC_None_',fold)]]
  regressor_high<-regression_snp_model_high[[paste0('SNP_AUC_None_',fold)]]
  
  preds_low<-predict(regressor_low,test_snp)
  preds_high<-predict(regressor_high,test_snp)
  
  if(type=='hard'){
    preds_snp<-round(ifelse(classif_risk=='neg',preds_low,preds_high),2)
  }else{
    classif_probs<-predict(classifier,test_clinical,type="prob")
    preds_snp<-round(classif_probs[,1]*preds_low+classif_probs[,2]*preds_high,2)
  }
  
  rae_sn<-median(abs(target-preds_snp)/target)
  rae_snp<-c(rae_clinical,rae_sn)
}


median(rae_clinical) #0.494
median(rae_snp) #0.500


type<-'soft'
rae_clinical<-c()
rae_snp<-c()

for(fold in 1:10){
  test_clinical<-test_list[[paste0('Clinical_AUC_None_',fold)]]
  
  test_snp<-test_list[[paste0('SNP_AUC_None_',fold)]]
  
  target<-test_clinical$SYNTAX.SCORE
  
  #Clinical_prediction
  classifier<-classifier_models[[paste0('Clinical_AUC_None_',fold)]]
  classif_risk<-predict(classifier,test_clinical)
  regressor_low<-regression_cl_model_low[[paste0('Clinical_AUC_None_',fold)]]
  regressor_high<-regression_cl_model_high[[paste0('Clinical_AUC_None_',fold)]]
  
  preds_low<-predict(regressor_low,test_clinical)
  preds_high<-predict(regressor_high,test_clinical)
  
  if(type=='hard'){
    preds_clinical<-round(ifelse(classif_risk=='neg',preds_low,preds_high),2)
  }else{
    classif_probs<-predict(classifier,test_clinical,type="prob")
    preds_clinical<-round(classif_probs[,1]*preds_low+classif_probs[,2]*preds_high,2)
  }
  
  rae_cl<-median(abs(target-preds_clinical)/target)
  rae_clinical<-c(rae_clinical,rae_cl)
  
  #SNP_prediction
  classifier<-classifier_models[[paste0('SNP_AUC_None_',fold)]]
  classif_risk<-predict(classifier,test_snp)
  regressor_low<-regression_snp_model_low[[paste0('SNP_AUC_None_',fold)]]
  regressor_high<-regression_snp_model_high[[paste0('SNP_AUC_None_',fold)]]
  
  preds_low<-predict(regressor_low,test_snp)
  preds_high<-predict(regressor_high,test_snp)
  
  if(type=='hard'){
    preds_snp<-round(ifelse(classif_risk=='neg',preds_low,preds_high),2)
  }else{
    classif_probs<-predict(classifier,test_snp,type="prob")
    preds_snp<-round(classif_probs[,1]*preds_low+classif_probs[,2]*preds_high,2)
  }
  
  rae_sn<-median(abs(target-preds_snp)/target)
  rae_snp<-c(rae_clinical,rae_sn)
}


median(rae_clinical) #0.469
median(rae_snp) #0.458




#################### Pooled cases ################

### threshold  investigation ####
threshold_classifier<-c()
threshold_classifier_snp<-c()
library(pROC)
for(fold in 1:10){
  ### Clinical ###
  train<-train_list[[paste0('Clinical_AUC_None_',fold)]]
  model<-classifier_models[[paste0('Clinical_AUC_None_',fold)]]
  probs<-predict(model,train,type="prob")[,2]
  roc_curve<-pROC::roc(train$Binary.SS,probs)
  best_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity", "ppv", "npv","youden"), best.method = "youden")
  threshold<-as.numeric(best_coords$threshold)
  threshold_classifier<-c(threshold_classifier,threshold)
  
  
  
  #### SNP ###
  train<-train_list[[paste0('SNP_AUC_None_',fold)]]
  model<-classifier_models[[paste0('SNP_AUC_None_',fold)]]
  probs<-predict(model,train,type="prob")[,2]
  roc_curve<-pROC::roc(train$Binary.SS,probs)
  best_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity", "ppv", "npv","youden"), best.method = "youden")
  threshold<-as.numeric(best_coords$threshold)
  threshold_classifier_snp<-c(threshold_classifier_snp,threshold)
}

######## NEW threshold ###########
type<-'hard'
ae_clinical<-c()
ae_snp<-c()
target_all<-c()
clinical_preds<-c()
snp_preds<-c()
fold_arr<-c()
for(fold in 1:10){
  
  test_clinical<-test_list[[paste0('Clinical_AUC_None_',fold)]]
  
  test_snp<-test_list[[paste0('SNP_AUC_None_',fold)]]
  
  target<-test_clinical$SYNTAX.SCORE
  target_all<-c(target_all,target)
  #Clinical_prediction
  classifier<-classifier_models[[paste0('Clinical_AUC_None_',fold)]]
  classif_risk<-predict(classifier,test_clinical,type="prob")[,2]
  regressor_low<-regression_cl_model_low[[paste0('Clinical_AUC_None_',fold)]]
  regressor_high<-regression_cl_model_high[[paste0('Clinical_AUC_None_',fold)]]
  
  preds_low<-predict(regressor_low,test_clinical)
  preds_high<-predict(regressor_high,test_clinical)
  
  
  preds_clinical<-round(ifelse(classif_risk>threshold_classifier[fold],preds_high,preds_low),2)
  
  
  ae_cl<-median(abs(target-preds_clinical))
  ae_clinical<-c(ae_clinical,ae_cl)
  
  clinical_preds<-c(clinical_preds,preds_clinical)
  #SNP_prediction
  classifier<-classifier_models[[paste0('SNP_AUC_None_',fold)]]
  classif_risk<-predict(classifier,test_snp,type='prob')[,2]
  regressor_low<-regression_snp_model_low[[paste0('SNP_AUC_None_',fold)]]
  regressor_high<-regression_snp_model_high[[paste0('SNP_AUC_None_',fold)]]
  
  preds_low<-predict(regressor_low,test_snp)
  preds_high<-predict(regressor_high,test_snp)
  
  
  preds_snp<-round(ifelse(classif_risk>threshold_classifier_snp[fold],preds_high,preds_low),2)

  
  ae_sn<-median(abs(target-preds_snp))
  ae_snp<-c(ae_clinical,ae_sn)
  
  snp_preds<-c(snp_preds,preds_snp)
  
  fold_arr<-c(fold_arr,rep(fold,length(preds_snp)))
}

res<-data.frame(Fold=fold_arr,Clinical.Predicted=clinical_preds,SNP.Predicted=snp_preds,Target=target_all)
save(res,file='C:/Users/30697/Desktop/PhD/Working/Review_ML_SNP/Code/Count_part_classifiers/Count_part.Rda')


median(abs(res$Target-res$Clinical.Predicted))
median(abs(res$Target-res$SNP.Predicted))

median(abs(target_all-clinical_preds)) #9.70
median(abs(target_all-snp_preds)) #9.17

median(abs(target_all-clinical_preds)/target_all) #0.559
median(abs(target_all-snp_preds)/target_all) #0.550


####################################










type<-'hard'
ae_clinical<-c()
ae_snp<-c()
target_all<-c()
clinical_preds<-c()
snp_preds<-c()
for(fold in 1:10){
  test_clinical<-test_list[[paste0('Clinical_AUC_None_',fold)]]
  
  test_snp<-test_list[[paste0('SNP_AUC_None_',fold)]]
  
  target<-test_clinical$SYNTAX.SCORE
  target_all<-c(target_all,target)
  #Clinical_prediction
  classifier<-classifier_models[[paste0('Clinical_AUC_None_',fold)]]
  classif_risk<-predict(classifier,test_clinical)
  regressor_low<-regression_cl_model_low[[paste0('Clinical_AUC_None_',fold)]]
  regressor_high<-regression_cl_model_high[[paste0('Clinical_AUC_None_',fold)]]
  
  preds_low<-predict(regressor_low,test_clinical)
  preds_high<-predict(regressor_high,test_clinical)
  
  if(type=='hard'){
    preds_clinical<-round(ifelse(classif_risk=='neg',preds_low,preds_high),2)
  }else{
    classif_probs<-predict(classifier,test_clinical,type="prob")
    preds_clinical<-round(classif_probs[,1]*preds_low+classif_probs[,2]*preds_high,2)
  }
  
  ae_cl<-median(abs(target-preds_clinical))
  ae_clinical<-c(ae_clinical,ae_cl)
  
  clinical_preds<-c(clinical_preds,preds_clinical)
  #SNP_prediction
  classifier<-classifier_models[[paste0('SNP_AUC_None_',fold)]]
  classif_risk<-predict(classifier,test_snp)
  regressor_low<-regression_snp_model_low[[paste0('SNP_AUC_None_',fold)]]
  regressor_high<-regression_snp_model_high[[paste0('SNP_AUC_None_',fold)]]
  
  preds_low<-predict(regressor_low,test_snp)
  preds_high<-predict(regressor_high,test_snp)
  
  if(type=='hard'){
    preds_snp<-round(ifelse(classif_risk=='neg',preds_low,preds_high),2)
  }else{
    classif_probs<-predict(classifier,test_clinical,type="prob")
    preds_snp<-round(classif_probs[,1]*preds_low+classif_probs[,2]*preds_high,2)
  }
  
  ae_sn<-median(abs(target-preds_snp))
  ae_snp<-c(ae_clinical,ae_sn)
  
  snp_preds<-c(snp_preds,preds_snp)
}

median(abs(target_all-clinical_preds)) #7.27
median(abs(target_all-snp_preds)) #7.13

median(abs(target_all-clinical_preds)/target_all) #0.499
median(abs(target_all-snp_preds)/target_all) #0.495



type<-'soft'
ae_clinical<-c()
ae_snp<-c()
target_all<-c()
clinical_preds<-c()
snp_preds<-c()
for(fold in 1:10){
  test_clinical<-test_list[[paste0('Clinical_AUC_None_',fold)]]
  
  test_snp<-test_list[[paste0('SNP_AUC_None_',fold)]]
  
  target<-test_clinical$SYNTAX.SCORE
  target_all<-c(target_all,target)
  #Clinical_prediction
  classifier<-classifier_models[[paste0('Clinical_AUC_None_',fold)]]
  classif_risk<-predict(classifier,test_clinical)
  regressor_low<-regression_cl_model_low[[paste0('Clinical_AUC_None_',fold)]]
  regressor_high<-regression_cl_model_high[[paste0('Clinical_AUC_None_',fold)]]
  
  preds_low<-predict(regressor_low,test_clinical)
  preds_high<-predict(regressor_high,test_clinical)
  
  if(type=='hard'){
    preds_clinical<-round(ifelse(classif_risk=='neg',preds_low,preds_high),2)
  }else{
    classif_probs<-predict(classifier,test_clinical,type="prob")
    preds_clinical<-round(classif_probs[,1]*preds_low+classif_probs[,2]*preds_high,2)
  }
  
  ae_cl<-median(abs(target-preds_clinical))
  ae_clinical<-c(ae_clinical,ae_cl)
  
  clinical_preds<-c(clinical_preds,preds_clinical)
  #SNP_prediction
  classifier<-classifier_models[[paste0('SNP_AUC_None_',fold)]]
  classif_risk<-predict(classifier,test_snp)
  regressor_low<-regression_snp_model_low[[paste0('SNP_AUC_None_',fold)]]
  regressor_high<-regression_snp_model_high[[paste0('SNP_AUC_None_',fold)]]
  
  preds_low<-predict(regressor_low,test_snp)
  preds_high<-predict(regressor_high,test_snp)
  
  if(type=='hard'){
    preds_snp<-round(ifelse(classif_risk=='neg',preds_low,preds_high),2)
  }else{
    classif_probs<-predict(classifier,test_snp,type="prob")
    preds_snp<-round(classif_probs[,1]*preds_low+classif_probs[,2]*preds_high,2)
  }
  
  ae_sn<-median(abs(target-preds_snp))
  ae_snp<-c(ae_clinical,ae_sn)
  
  snp_preds<-c(snp_preds,preds_snp)
}

median(abs(target_all-clinical_preds)) #7.42
median(abs(target_all-snp_preds)) #7.35

median(abs(target_all-clinical_preds)/target_all) #0.457
median(abs(target_all-snp_preds)/target_all) #0.459

