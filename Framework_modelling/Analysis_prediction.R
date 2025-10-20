################## Regression result evaluation --- Wilcoxon, Multiclassification ##############
load('./Clinical_predictions.Rda')
load('./Clinical_target.Rda')

load('./SNP_predictions.Rda')
load('./SNP_target.Rda')


res_all<-NULL
for(name in names(predictions_cl_all)){
  predictions_cl<-predictions_cl_all[[name]]
  target_cl<-target_cl_all[[name]]
  ae_cl<-abs(predictions_cl-target_cl)
  
  
  name_spl<-strsplit(name,split='_')[[1]]
  name_snp<-paste0('SNP_',name_spl[2],'_',name_spl[3])
  predictions_snp<-predictions_snp_all[[name_snp]]
  target_snp<-target_snp_all[[name_snp]]
  ae_snp<-abs(predictions_snp - target_snp)
  
  res<-wilcox.test(ae_cl,ae_snp)
  med_cl<-median(ae_cl)
  med_snp<-median(ae_snp)
  
  df<-data.frame(Optimization=c(name_spl[2]),Resampling=c(name_spl[3]),MedAE.Clinical=med_cl,MedAE.SNP=med_snp,Wilc.p.val=res$p.value)
  res_all<-rbind(res_all,df)
}
save(res_all,file='./Wilcoxon_results_Alltogether.Rda')
conf_matrix_list<-list()
for(name in names(predictions_cl_all)){
  predictions_cl<-predictions_cl_all[[name]]
  target_cl<-target_cl_all[[name]]
  ae_cl<-abs(predictions_cl-target_cl)
  predicted_classes_clinicial<-as.factor(ifelse(predictions_cl==0,0,ifelse(predictions_cl<23,1,2)))#ifelse(predictions_clinical_all<32,2,3))))
  target_classes<-as.factor(ifelse(target_cl==0,0,ifelse(target_cl<23,1,2)))#ifelse(target_clinical_all<32,2,3))))
  cli_cm<-caret::confusionMatrix(predicted_classes_clinicial,target_classes)
  
  
  name_spl<-strsplit(name,split='_')[[1]]
  name_snp<-paste0('SNP_',name_spl[2],'_',name_spl[3])
  predictions_snp<-predictions_snp_all[[name_snp]]
  target_snp<-target_snp_all[[name_snp]]
  ae_snp<-abs(predictions_snp - target_snp)
  predicted_classes_snp<-as.factor(ifelse(predictions_snp==0,0,ifelse(predictions_snp<23,1,2)))#ifelse(predictions_snp_all<32,2,3))))
  target_classes_snp<-as.factor(ifelse(target_snp==0,0,ifelse(target_snp<23,1,2)))#ifelse(target_snp_all<32,2,3))))
  snp_cm<-caret::confusionMatrix(predicted_classes_snp,target_classes_snp)
  
  
  conf_matrix_list[[name]]<-cli_cm
  conf_matrix_list[[name_snp]]<-snp_cm
}
save(conf_matrix_list,file='./ConfusionMatrix_results_Alltogether_3class.Rda')


conf_matrix_list4<-list()
for(name in names(predictions_cl_all)){
  predictions_cl<-predictions_cl_all[[name]]
  target_cl<-target_cl_all[[name]]
  ae_cl<-abs(predictions_cl-target_cl)
  predicted_classes_clinicial<-as.factor(ifelse(predictions_cl==0,0,ifelse(predictions_cl<23,1,ifelse(predictions_cl<32,2,3))))
  target_classes<-as.factor(ifelse(target_cl==0,0,ifelse(target_cl<23,1,ifelse(target_cl<32,2,3))))
  cli_cm<-caret::confusionMatrix(predicted_classes_clinicial,target_classes)
  
  
  name_spl<-strsplit(name,split='_')[[1]]
  name_snp<-paste0('SNP_',name_spl[2],'_',name_spl[3])
  predictions_snp<-predictions_snp_all[[name_snp]]
  target_snp<-target_snp_all[[name_snp]]
  ae_snp<-abs(predictions_snp - target_snp)
  predicted_classes_snp<-as.factor(ifelse(predictions_snp==0,0,ifelse(predictions_snp<23,1,ifelse(predictions_snp<32,2,3))))
  target_classes_snp<-as.factor(ifelse(target_snp==0,0,ifelse(target_snp<23,1,ifelse(target_snp<32,2,3))))
  snp_cm<-caret::confusionMatrix(predicted_classes_snp,target_classes_snp)
  
  
  conf_matrix_list4[[name]]<-cli_cm
  conf_matrix_list4[[name_snp]]<-snp_cm
}
save(conf_matrix_list4,file='./ConfusionMatrix_results_Alltogether_4class.Rda')





