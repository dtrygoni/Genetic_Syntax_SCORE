########################## Evaluation of different classifiers/ regressors ########################################
library(caret)
library(MLmetrics)


###### load models/tests
load('./Classifier_Models.Rda')
load('./Test_sets.Rda')
load('./Train_sets.Rda')



classifier_results<-NULL
type_arr<-c()
opt_arr<-c()
resamp_arr<-c()
fold_arr<-c()
balAcc_arr<-c()
sens_arr<-c()
spec_arr<-c()
for(i in names(test_list)){
  comp<-strsplit(i,split="_")[[1]]
  type<-comp[1]
  optimization<-comp[2]
  resampling<-comp[3]
  fold<-comp[4]
  
  model<-classifier_models[[i]]
  test<-test_list[[i]]
  
  predicted<-predict(model,newdata = test)
  predicted<-as.factor(ifelse(predicted=='pos',1,0))
  actual<-test$Binary.SS
  
  
  cm<-confusionMatrix(data=predicted,reference = actual,positive="1")
  balAcc<-round(as.numeric(cm$byClass['Balanced Accuracy']),3)
  sens<-round(as.numeric(cm$byClass['Sensitivity']),3)
  spec<-round(as.numeric(cm$byClass['Specificity']),3)
  
  
  
  type_arr<-c(type_arr,type)
  opt_arr<-c(opt_arr,optimization)
  resamp_arr<-c(resamp_arr,resampling)
  fold_arr<-c(fold_arr,fold)
  balAcc_arr<-c(balAcc_arr,balAcc)
  sens_arr<-c(sens_arr,sens)
  spec_arr<-c(spec_arr,spec)
  }
df<-data.frame(Type=type_arr,Optimization=opt_arr,Resampling=resamp_arr,Fold=fold_arr,`Balanced Accuracy`=balAcc_arr,Sensitivity=sens_arr,Specificity=spec_arr)
save(df,file='./Classifier_results_per_fold.Rda')

types<-unique(type_arr)
optimizations<-unique(opt_arr)
resamplings<-unique(resamp_arr)
confusion_matrices<-list()
for(type in types){
  for(opt in optimizations){
    for(resamp in resamplings){
      test_cases<-c()
      predicted_all<-c()
      for(jj in 1:10){
        name<-paste0(type,"_",opt,"_",resamp,"_",jj)
        test<-test_list[[name]]
        model<-classifier_models[[name]]
        
        test_cases<-c(test_cases,test$Binary.SS)
        predicted<-predict(model,newdata=test)
        predicted<-ifelse(predicted=='pos',1,0)
       
        predicted_all<-c(predicted_all,predicted)
        
        
      }
      test_cases<-as.factor(ifelse(test_cases==1,0,1))
      predicted_all<-as.factor(predicted_all)
      cm<-confusionMatrix(data=predicted_all,reference=test_cases,positive="1")
      
      names_f<-paste0(type,'_',opt,'_',resamp)
      confusion_matrices[[names_f]]<-cm
    }
  }
}
save(confusion_matrices,file='./Classifier_results_confusion_matrices_allfolds.Rda')




################################## Regressors
##########Low
load('./Regression_Clinical_Low_Models.Rda')
load('./Regression_SNP_Low_Models.Rda')
clinical_names<-c()
snp_names<-c()
for(resamp in resamplings){
  for(jj in 1:10){
    string_wanted<-paste0('Clinical_','bal.Acc_',resamp,'_',jj)
    clinical_names<-c(clinical_names,string_wanted)
    
    string_wanted<-paste0('SNP_','bal.Acc_',resamp,'_',jj)
    snp_names<-c(snp_names,string_wanted)
  }
}
test_clinical_list<-test_list[names(test_list) %in% clinical_names]
test_snp_list<-test_list[names(test_list) %in% snp_names]


res_low<-NULL
for(resamp in resamplings){
  target_clinical<-c()
  target_snp<-c()
  
  predicted_clinical<-c()
  predicted_snp<-c()
  for(jj in 1:10){
  clinical_name<-paste0('Clinical_bal.Acc_',resamp,'_',jj)
  snp_name<-paste0('SNP_bal.Acc_',resamp,'_',jj)
  
  test_clinical<-test_clinical_list[[clinical_name]]
  test_snp<-test_snp_list[[snp_name]]
  
  
  test_clinical<-test_clinical[test_clinical$Binary.SS==0,]
  test_snp<-test_snp[test_snp$Binary.SS==0,]
  
  model_clinical<-regression_cl_model_low[[clinical_name]]
  model_snp<-regression_snp_model_low[[snp_name]]
  
  target_syntax_cl<-test_clinical$SYNTAX.SCORE
  target_syntax_snp<-test_snp$SYNTAX.SCORE
  
  predicted_cl_syntax<-predict(model_clinical,newdata=test_clinical)
  predicted_snp_syntax<-predict(model_snp,newdata=test_snp)
  
  target_clinical<-c(target_clinical,target_syntax_cl)
  target_snp<-c(target_snp,target_syntax_snp)
  
  predicted_clinical<-c(predicted_clinical,round(predicted_cl_syntax,3))
  predicted_snp<-c(predicted_snp,round(predicted_snp_syntax,3))
  }
  
  ae_cl<-abs(target_clinical-predicted_clinical)
  rae_cl<-abs(target_clinical-predicted_clinical)/target_clinical
  
  ae_snp<-abs(target_snp-predicted_snp)
  rae_snp<-abs(target_snp-predicted_snp)/target_snp
  
  res_ae<-wilcox.test(ae_cl,ae_snp)
  res_rae<-wilcox.test(rae_cl,rae_snp)
  
  df<-data.frame(Resampling=c(resamp,resamp),Metric=c('AE','RAE'),Clinical=c(round(median(ae_cl),3),round(median(rae_cl),3)),
                 SNP=c(round(median(ae_snp),3),round(median(rae_snp),3)),P.Val=c(res_ae$p.value,res_rae$p.value))
  
  res_low<-rbind(res_low,df)
                 
}
save(res_low,file='./Wilcoxon_results_Regression_low.Rda')
######### High
load('./Regression_Clinical_High_Models.Rda')
load('./Regression_SNP_High_Models.Rda')
res_high<-NULL
for(resamp in resamplings){
  target_clinical<-c()
  target_snp<-c()
  
  predicted_clinical<-c()
  predicted_snp<-c()
  for(jj in 1:10){
    clinical_name<-paste0('Clinical_bal.Acc_',resamp,'_',jj)
    snp_name<-paste0('SNP_bal.Acc_',resamp,'_',jj)
    
    test_clinical<-test_clinical_list[[clinical_name]]
    test_snp<-test_snp_list[[snp_name]]
    
    
    test_clinical<-test_clinical[test_clinical$Binary.SS==1,]
    test_snp<-test_snp[test_snp$Binary.SS==1,]
    
    model_clinical<-regression_cl_model_high[[clinical_name]]
    model_snp<-regression_snp_model_high[[snp_name]]
    
    target_syntax_cl<-test_clinical$SYNTAX.SCORE
    target_syntax_snp<-test_snp$SYNTAX.SCORE
    
    predicted_cl_syntax<-predict(model_clinical,newdata=test_clinical)
    predicted_snp_syntax<-predict(model_snp,newdata=test_snp)
    
    target_clinical<-c(target_clinical,target_syntax_cl)
    target_snp<-c(target_snp,target_syntax_snp)
    
    predicted_clinical<-c(predicted_clinical,round(predicted_cl_syntax,3))
    predicted_snp<-c(predicted_snp,round(predicted_snp_syntax,3))
  }
  
  ae_cl<-abs(target_clinical-predicted_clinical)
  rae_cl<-abs(target_clinical-predicted_clinical)/target_clinical
  
  ae_snp<-abs(target_snp-predicted_snp)
  rae_snp<-abs(target_snp-predicted_snp)/target_snp
  
  res_ae<-wilcox.test(ae_cl,ae_snp)
  res_rae<-wilcox.test(rae_cl,rae_snp)
  
  df<-data.frame(Resampling=c(resamp,resamp),Metric=c('AE','RAE'),Clinical=c(round(median(ae_cl),3),round(median(rae_cl),3)),
                 SNP=c(round(median(ae_snp),3),round(median(rae_snp),3)),P.Val=c(res_ae$p.value,res_rae$p.value))
  
  res_high<-rbind(res_high,df)
  
}
save(res_high,file='./Wilcoxon_results_Regression_High.Rda')




####### Together with classifier + regressors
res_all<-NULL
types<-unique(type_arr)
optimizations<-unique(opt_arr)
resamplings<-unique(resamp_arr)
for(opt in optimizations){
    for(resamp in resamplings){
      target_clinical<-c()
      target_snp<-c()
      
      predicted_clinical<-c()
      predicted_snp<-c()
      
      for(jj in 1:10){
        clinical_name<-paste0('Clinical_bal.Acc_',resamp,'_',jj)
        snp_name<-paste0('SNP_bal.Acc_',resamp,'_',jj)
        
        test_clinical<-test_clinical_list[[clinical_name]]
        test_snp<-test_snp_list[[snp_name]]
        
        ########## Clinical ####################
        classifier_cl_name<-paste0('Clinical_',opt,'_',resamp,'_',jj)
        classifier_clinical<-classifier_models[[classifier_cl_name]]
        
        regressor_cl_low<-regression_cl_model_low[[clinical_name]]
        regressor_cl_high<-regression_cl_model_high[[clinical_name]]
        
        labels_predicted<-predict(classifier_clinical,newdata=test_clinical)
        predictions_low<-predict(regressor_cl_low,newdata=test_clinical)
        predictions_high<-predict(regressor_cl_high,newdata=test_clinical)
        
        #predictions_low_prob<-predict(regressor_cl_low,newdata=test_clinical,type="prob")
        
        final_predictions<-round(ifelse(labels_predicted=='pos',predictions_high,predictions_low),3)
        target<-test_clinical$SYNTAX.SCORE
        target_clinical<-c(target_clinical,target)
        predicted_clinical<-c(predicted_clinical,final_predictions)
        
        
        
        ########## SNP ####################
        classifier_snp_name<-paste0('SNP_',opt,'_',resamp,'_',jj)
        classifier_snp<-classifier_models[[classifier_snp_name]]
        
        regressor_snp_low<-regression_snp_model_low[[snp_name]]
        regressor_snp_high<-regression_snp_model_high[[snp_name]]
        
        labels_predicted<-predict(classifier_snp,newdata=test_snp)
        predictions_low<-predict(regressor_snp_low,newdata=test_snp)
        predictions_high<-predict(regressor_snp_high,newdata=test_snp)
        
        final_predictions<-round(ifelse(labels_predicted=='pos',predictions_high,predictions_low),3)
        target<-test_snp$SYNTAX.SCORE
        target_snp<-c(target_snp,target)
        predicted_snp<-c(predicted_snp,final_predictions)
        
      }
      
      ae_cl<-abs(target_clinical-predicted_clinical)
      rae_cl<-abs(target_clinical-predicted_clinical)/target_clinical
      
      ae_snp<-abs(target_snp-predicted_snp)
      rae_snp<-abs(target_snp-predicted_snp)/target_snp
      
      res_ae<-wilcox.test(ae_cl,ae_snp)
      res_rae<-wilcox.test(rae_cl,rae_snp)
      
      df<-data.frame(Resampling=c(resamp,resamp),Optimization=c(opt,opt),Metric=c('AE','RAE'),Clinical=c(round(median(ae_cl),3),round(median(rae_cl),3)),
                     SNP=c(round(median(ae_snp),3),round(median(rae_snp),3)),P.Val=c(res_ae$p.value,res_rae$p.value))
      res_all<-rbind(res_all,df)
    }
  }
save(res_all,file='./Wilcoxon_results_Regression_All.Rda')





